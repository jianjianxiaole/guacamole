package org.hammerlab.guacamole.commands

import org.apache.spark.SparkContext
import org.apache.spark.rdd.RDD
import org.bdgenomics.adam.rdd.ADAMContext
import org.bdgenomics.formats.avro.DatabaseVariantAnnotation
import org.hammerlab.guacamole.alignment.AffineGapPenaltyAlignment
import org.hammerlab.guacamole.assembly.AssemblyUtils
import org.hammerlab.guacamole.distributed.LociPartitionUtils.{LociPartitioning, partitionLociAccordingToArgs}
import org.hammerlab.guacamole.distributed.WindowFlatMapUtils.windowFlatMapWithState
import org.hammerlab.guacamole.logging.{DelayedMessages, LoggingUtils}
import org.hammerlab.guacamole.pileup.Pileup
import org.hammerlab.guacamole.reads.MappedRead
import org.hammerlab.guacamole.readsets.{InputFilters, ReadSets, SomaticCallerArgs}
import org.hammerlab.guacamole.reference.{ContigSequence, ReferenceBroadcast}
import org.hammerlab.guacamole.variants._
import org.hammerlab.guacamole.windowing.SlidingWindow
import org.kohsuke.args4j.{Option => Args4jOption}

object SomaticAssemblyCaller {

  class Arguments extends SomaticCallerArgs with Serializable  {

    @Args4jOption(name = "--kmer-size", usage = "Length of kmer used for DeBruijn Graph assembly")
    var kmerSize: Int = 45

    @Args4jOption(name = "--snv-window-range", usage = "Number of bases before and after to check for additional matches or deletions")
    var snvWindowRange: Int = 20

    @Args4jOption(name = "--min-alignment-quality", usage = "Minimum alignment qualities of the read")
    var minAlignmentQuality: Int = 30

    @Args4jOption(name = "--reference-fasta", required = true, usage = "Local path to a reference FASTA file")
    var referenceFastaPath: String = null

    @Args4jOption(name = "--min-area-vaf", required = false, usage = "Minimum variant allele frequency to investigate area")
    var minAreaVaf: Int = 5

    @Args4jOption(name = "--min-occurrence", required = false, usage = "Minimum occurrences to include a kmer ")
    var minOccurrence: Int = 3

    @Args4jOption(name = "--dbsnp-vcf", required = false, usage = "VCF file to identify DBSNP variants")
    var dbSnpVcf: String = ""

    @Args4jOption(name = "--odds", usage = "Minimum log odds threshold for possible variant candidates")
    var oddsThreshold: Int = 1000

    @Args4jOption(name = "--min-likelihood", usage = "Minimum Phred-scaled likelihood. Default: 0 (off)")
    var minLikelihood: Int = 30

    @Args4jOption(name = "--shortcut-assembly", required = false, usage = "Skip assembly process in inactive regions")
    var shortcutAssembly: Boolean = false

    @Args4jOption(name = "--min-read-depth", usage = "Minimum number of reads for a genotype call")
    var minReadDepth: Int = 8

  }

  object Caller extends SparkCommand[Arguments] {
    override val name = "somatic-assembly"
    override val description = "call somatic variants by assembling the surrounding region of reads"

    override def run(args: Arguments, sc: SparkContext): Unit = {

      VariantUtils.validateArguments(args)
      val loci = args.parseLoci(sc.hadoopConfiguration)
      val filters = InputFilters(
        overlapsLoci = loci,
        nonDuplicate = true,
        passedVendorQualityChecks = true
      )

      val reference = ReferenceBroadcast(args.referenceFastaPath, sc)

      val (tumorReads, normalReads, contigLengths) =
        ReadSets.loadTumorNormalReads(
          args,
          sc,
          filters
        )

      val minAlignmentQuality = args.minAlignmentQuality

      def filterToQualityReads(reads: RDD[MappedRead]) = {
        reads.filter(_.alignmentQuality > minAlignmentQuality)
      }

      val normalQualityReads = filterToQualityReads(normalReads.mappedReads)
      val tumorQualityReads = filterToQualityReads(tumorReads.mappedReads)

      val lociPartitions =
        partitionLociAccordingToArgs(
          args,
          loci.result(contigLengths),
          Vector(tumorReads.mappedReads, normalReads.mappedReads)
        )

      val kmerSize = args.kmerSize
      var potentialGenotypes = discoverSomaticGenotypes(
        normalQualityReads,
        tumorQualityReads,
        kmerSize = kmerSize,
        snvWindowRange = args.snvWindowRange,
        minOccurrence = args.minOccurrence,
        minAreaVaf = args.minAreaVaf / 100.0f,
        reference,
        lociPartitions,
        minDepth = args.minReadDepth,
        minLikelihood = args.minLikelihood,
        oddsThreshold = args.oddsThreshold,
        minAltReads = args.minOccurrence,
        shortcutAssembly = args.shortcutAssembly
      )

      potentialGenotypes.persist()
      LoggingUtils.progress("Computed %,d potential genotypes".format(potentialGenotypes.count))

      if (args.dbSnpVcf != "") {
        val adamContext = new ADAMContext(sc)
        val dbSnpVariants = adamContext.loadVariantAnnotations(args.dbSnpVcf)
        potentialGenotypes = potentialGenotypes
          .keyBy(_.adamVariant)
          .leftOuterJoin(dbSnpVariants.keyBy(_.getVariant))
          .map(_._2).map({
          case (calledAllele: CalledSomaticAllele, dbSnpVariant: Option[DatabaseVariantAnnotation]) =>
            calledAllele.copy(rsID = dbSnpVariant.map(_.getDbSnpId))
        })
      }

      VariantUtils.writeVariantsFromArguments(
        args,
        potentialGenotypes.flatMap(AlleleConversions.calledSomaticAlleleToADAMGenotype)
      )

      DelayedMessages.default.print()
    }

    /**
     *
     * @param normalReads    Mapped reads of the normal sample
     * @param tumorReads     Mapped reads of the tumor sample
     * @param kmerSize       Length of subsequence to use for assembly
     * @param snvWindowRange Number of bases to consider left and right of the locus
     * @param minOccurrence  Minimum appearances of kmers to include in the graph
     * @param minAreaVaf     Minimum variant allele frequency in a region to compute variants
     * @param reference      Broadcasted reference sequences
     * @param lociPartitions
     * @return RDD of called somatic variants
     */
    def discoverSomaticGenotypes(normalReads: RDD[MappedRead],
                                 tumorReads: RDD[MappedRead],
                                 kmerSize: Int,
                                 snvWindowRange: Int,
                                 minOccurrence: Int,
                                 minAreaVaf: Float,
                                 reference: ReferenceBroadcast,
                                 lociPartitions: LociPartitioning,
                                 minLikelihood: Int,
                                 minDepth: Int,
                                 oddsThreshold: Int,
                                 minAltReads: Int = 3,
                                 shortcutAssembly: Boolean = false): RDD[CalledSomaticAllele] = {

      val genotypes: RDD[CalledSomaticAllele] =
        windowFlatMapWithState[MappedRead, CalledSomaticAllele, Option[Long]](
          Vector(normalReads, tumorReads),
          lociPartitions,
          skipEmpty = true,
          halfWindowSize = snvWindowRange,
          initialState = None,
          (lastCalledLocus, windows) => {
            val currentLocus = windows.head.currentLocus
            val normalWindow = windows(0)
            val tumorWindow = windows(1)
            val referenceName = normalWindow.referenceName

            val referenceContig = reference.getContig(referenceName)
            val currentLocusNormalReads =
              normalWindow
                .currentRegions()
                .filter(_.overlapsLocus(currentLocus))

            val currentLocusTumorReads =
              tumorWindow
                .currentRegions()
                .filter(_.overlapsLocus(currentLocus))

            def hasMinDepth(locus: Long): Boolean = {
              normalWindow.currentRegions().count(_.overlapsLocus(locus)) >= minDepth &&
                tumorWindow.currentRegions().count(_.overlapsLocus(locus)) >= minDepth
            }
            val referenceStart = (currentLocus - normalWindow.halfWindowSize).toInt
            val referenceEnd = (currentLocus + normalWindow.halfWindowSize).toInt

            if (
              normalWindow.currentRegions().isEmpty ||
              tumorWindow.currentRegions().isEmpty ||
              !hasMinDepth(referenceStart) || !hasMinDepth(referenceEnd)
            )
              (lastCalledLocus, Iterator.empty)
            else if (
              shortcutAssembly &&
                !AssemblyUtils.isActiveRegion(currentLocusNormalReads, referenceContig, minAreaVaf) &&
                  !AssemblyUtils.isActiveRegion(currentLocusTumorReads, referenceContig, minAreaVaf)
            ) {
              callPileupBasedSomaticVariants(
                Pileup(currentLocusNormalReads, referenceName, currentLocus, referenceContig),
                Pileup(currentLocusTumorReads, referenceName, currentLocus, referenceContig),
                minLikelihood,
                oddsThreshold,
                minAltReads,
                lastCalledLocus
              )
            } else {
              assembleRegionAndCallVariants(
                normalWindow,
                tumorWindow,
                kmerSize,
                minOccurrence,
                reference,
                lastCalledLocus,
                referenceContig)
            }
          }
        )
      genotypes
    }


    /**
     *
     * @param normalWindow
     * @param tumorWindow
     * @param kmerSize
     * @param minOccurrence
     * @param reference
     * @param lastCalledLocus
     * @param referenceContig
     * @return
     */
    private def assembleRegionAndCallVariants(normalWindow: SlidingWindow[MappedRead],
                                              tumorWindow: SlidingWindow[MappedRead],
                                              kmerSize: Int,
                                              minOccurrence: Int,
                                              reference: ReferenceBroadcast,
                                              lastCalledLocus: Option[Long],
                                              referenceContig: ContigSequence): (Option[Long], Iterator[CalledSomaticAllele]) = {

      val currentLocus = normalWindow.currentLocus
      val halfWindowSize = normalWindow.halfWindowSize

      val referenceStart = (currentLocus - halfWindowSize).toInt
      val referenceEnd = (currentLocus + halfWindowSize).toInt
      log.warn(s"Performing variant calling from assembly in ${normalWindow.referenceName}:${referenceStart}-$referenceEnd")

      val topNormalPaths = AssemblyUtils.discoverHaplotypes(
        normalWindow,
        kmerSize,
        reference,
        minOccurrence
      )

      lazy val tumorPaths = AssemblyUtils.discoverHaplotypes(
        tumorWindow,
        kmerSize,
        reference,
        minOccurrence
      )

      if (topNormalPaths.isEmpty || tumorPaths.isEmpty) {
        (lastCalledLocus, Iterator.empty)
      } else {

        def buildVariant(variantLocus: Int,
                         referenceBases: Array[Byte],
                         alternateBases: Array[Byte]): CalledSomaticAllele = {

          // This is just placeholder evidence
          val evidence = AlleleEvidence(
            likelihood = 1,
            readDepth = 1,
            alleleReadDepth = 1,
            forwardDepth = 1,
            alleleForwardDepth = 1,
            meanMappingQuality = 1,
            medianMappingQuality = 1,
            meanBaseQuality = 1,
            medianBaseQuality = 1,
            medianMismatchesPerRead = 0
          )

          val allele = Allele(referenceBases, alternateBases)
          CalledSomaticAllele(
            tumorWindow.currentRegions().head.sampleName,
            normalWindow.referenceName,
            variantLocus,
            allele,
            somaticLogOdds = 10, //placeholder
            normalReferenceEvidence = evidence,
            tumorVariantEvidence = evidence
          )
        }

        val variants =
          (for {
            tumorPath <- tumorPaths
            alignmentsToNormal = normalPaths.map(normalPath => (normalPath, AffineGapPenaltyAlignment.align(tumorPath, normalPath)))
            if !alignmentsToNormal.exists(_._2.nonVariant)
            topAlignmentToNormal = alignmentsToNormal.minBy(_._2.alignmentScore)
            pathVariants = AssemblyUtils.buildVariantsFromPath[CalledSomaticAllele](
              tumorPath,
              referenceStart,
              referenceContig,
              path => topAlignmentToNormal._2,
              buildVariant
            )
            variant <- pathVariants
          } yield {
            variant
          }).toSet
            .filter(v => v.allele.altBases.nonEmpty && v.allele.refBases.nonEmpty)
            .filter(variant => lastCalledLocus.forall(_ < variant.start)) // Filter variants before last called

        val lastCalledVariantLocus = variants.view.map(_.start).reduceOption(_ max _).orElse(lastCalledLocus)

        // Advance both windows
        normalWindow.setCurrentLocus(referenceEnd - kmerSize)
        tumorWindow.setCurrentLocus(referenceEnd - kmerSize)

        (lastCalledVariantLocus, variants.iterator)
      }
    }

    private def callPileupBasedSomaticVariants(normalPileup: Pileup,
                                               tumorPileup: Pileup,
                                               minLikelihood: Int,
                                               oddsThreshold: Int,
                                               minAltReads: Int,
                                               lastCalledLocus: Option[Long]): (Option[Long], Iterator[CalledSomaticAllele]) = {
      val tumorAltReads = tumorPileup.depth - tumorPileup.referenceDepth
      if (tumorAltReads >= minAltReads) {
        val variants: Seq[CalledSomaticAllele] = SomaticStandard.callSomaticPileupVariant(
          normalPileup = normalPileup,
          tumorPileup = tumorPileup,
          oddsThreshold = oddsThreshold
        ).filter(v => v.allele.altBases.nonEmpty && v.allele.refBases.nonEmpty)
          .filter(_.phredScaledSomaticLikelihood > minLikelihood)

        (variants.lastOption.map(_.start).orElse(lastCalledLocus), variants.iterator)
      } else {
        (lastCalledLocus, Iterator.empty)
      }
    }
  }
}
