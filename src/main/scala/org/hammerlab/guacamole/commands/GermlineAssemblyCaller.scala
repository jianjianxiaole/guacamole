package org.hammerlab.guacamole.commands

import breeze.linalg.DenseVector
import breeze.stats.{ mean, median }
import htsjdk.samtools.CigarOperator
import org.apache.spark.SparkContext
import org.apache.spark.rdd.RDD
import org.hammerlab.guacamole.Common.Arguments.GermlineCallerArgs
import org.hammerlab.guacamole._
import org.hammerlab.guacamole.alignment.AffineGapPenaltyAlignment
import org.hammerlab.guacamole.assembly.DeBruijnGraph
import org.hammerlab.guacamole.reads.{ MappedRead, Read }
import org.hammerlab.guacamole.reference.{ ReferenceBroadcast, ReferenceGenome }
import org.hammerlab.guacamole.variants.{ Allele, AlleleConversions, AlleleEvidence, CalledAllele }
import org.hammerlab.guacamole.windowing.SlidingWindow
import org.kohsuke.args4j.{ Option => Args4jOption }

import scala.collection.JavaConversions._

/**
 * Simple assembly based germline variant caller
 *
 * Overview:
 * - Find areas where > `min-area-vaf` %  of the reads show a variant
 * - Re-examine those areas with N base window around it
 * - Place the reads in the `snv-window-range`-base window into a DeBruijn graph
 *   and find paths between the start and end of the reference sequence that covers that region
 * - Align those paths to the reference sequence to discover variants
 *
 */
object GermlineAssemblyCaller {

  class Arguments extends GermlineCallerArgs {

    @Args4jOption(name = "--kmer-size", usage = "Length of kmer used for DeBrujin Graph assembly")
    var kmerSize: Int = 45

    @Args4jOption(name = "--snv-window-range", usage = "Number of bases before and after to check for additional matches or deletions")
    var snvWindowRange: Int = 20

    @Args4jOption(name = "--min-average-base-quality", usage = "Minimum average of base qualities in the read")
    var minAverageBaseQuality: Int = 20

    @Args4jOption(name = "--min-alignment-quality", usage = "Minimum alignment qualities of the read")
    var minAlignmentQuality: Int = 30

    @Args4jOption(name = "--reference-fasta", required = false, usage = "Local path to a reference FASTA file")
    var referenceFastaPath: String = null

    @Args4jOption(name = "--min-area-vaf", required = false, usage = "Minimum variant allele frequency to investigate area")
    var minAreaVaf: Int = 5

    @Args4jOption(name = "--min-occurrence", required = false, usage = "Minimum occurrences to include a kmer ")
    var minOccurrence: Int = 3

  }

  object Caller extends SparkCommand[Arguments] {
    override val name = "germline-assembly"
    override val description = "call germline variants by assembling the surrounding region of reads"

    /**
     *
     * @param graph An existing DeBruijn graph of the reads
     * @param currentWindow Window of reads that overlaps the current loci
     * @param kmerSize Length kmers in the DeBruijn graph
     * @param minOccurrence Minimum times a kmer must appear to be in the DeBruijn graph
     * @param expectedPloidy Expected ploidy, or expected number of valid paths through the graph
     * @param maxPathsToScore Number of paths to align to the reference to score them
     * @return A set of variants in the window
     *         TODO: This currently passes along a graph as state, but rebuilds it each time
     *         This can be updated to use the existing graph and just update it
     */
    def discoverHaplotypes(graph: Option[DeBruijnGraph],
                           currentWindow: SlidingWindow[MappedRead],
                           kmerSize: Int,
                           reference: ReferenceGenome,
                           minOccurrence: Int = 3,
                           expectedPloidy: Int = 2,
                           maxPathsToScore: Int = 6): (Option[DeBruijnGraph], Iterator[CalledAllele]) = {

      val locus = currentWindow.currentLocus
      val reads = currentWindow.currentRegions()

      // TODO: Should update graph instead of rebuilding it
      // Need to keep track of reads removed from last update and reads added
      val currentGraph: DeBruijnGraph = DeBruijnGraph(
        reads.map(_.sequence),
        kmerSize,
        minOccurrence,
        mergeNodes = true
      )

      val referenceStart = (reads.head.start).toInt
      val referenceEnd = (reads.last.end).toInt

      val currentReference: Array[Byte] = reference.getReferenceSequence(
        currentWindow.referenceName,
        referenceStart,
        referenceEnd
      )

      val referenceKmerSource = currentReference.take(kmerSize)
      val referenceKmerSink = currentReference.takeRight(kmerSize)

      // If the current window size doesn't cover the kmer size we won't be able to find the reference start/end
      if (currentReference.size < kmerSize)
        return (graph, Iterator.empty)

      val paths = currentGraph.depthFirstSearch(
        referenceKmerSource,
        referenceKmerSink
      )

      val sampleName = reads.head.sampleName
      val referenceContig = reads.head.referenceContig

      // Score up to the maximum number of paths, by aligning them against the reference
      // Take the best aligning `expectedPloidy` paths
      val pathAndAlignments =
        if (paths.length <= maxPathsToScore) {
          paths.map(path => {
            val mergedPath = DeBruijnGraph.mergeOverlappingSequences(path, kmerSize)
            (mergedPath, AffineGapPenaltyAlignment.align(mergedPath, currentReference))
          })
            .sortBy(_._2.alignmentScore)
            .take(expectedPloidy)
        } else {
          log.warn(s"In window ${referenceContig}:${referenceStart}-$referenceEnd " +
            s"there were ${paths.length} paths found, all variants skipped")
          List.empty
        }

      // Build a variant using the current offset and
      // read evidence
      def buildVariant(referenceOffset: Int,
                       referenceBases: Array[Byte],
                       alternateBases: Array[Byte]) = {
        val allele = Allele(
          referenceBases,
          alternateBases
        )

        val depth = reads.length
        val mappingQualities = DenseVector(reads.map(_.alignmentQuality.toFloat).toArray)
        val baseQualities = DenseVector(reads.flatMap(_.baseQualities).map(_.toFloat).toArray)
        CalledAllele(
          sampleName,
          referenceContig,
          referenceStart + referenceOffset,
          allele,
          AlleleEvidence(
            likelihood = 1,
            readDepth = depth,
            alleleReadDepth = depth,
            forwardDepth = depth,
            alleleForwardDepth = depth,
            meanMappingQuality = mean(mappingQualities),
            medianMappingQuality = median(mappingQualities),
            meanBaseQuality = mean(baseQualities),
            medianBaseQuality = median(baseQualities),
            medianMismatchesPerRead = 0
          )
        )
      }

      val variants =
        pathAndAlignments.flatMap(kv => {
          val path = kv._1
          val alignment = kv._2
          var referenceIndex = 0
          var pathIndex = 0

          // Find the alignment sequences using the CIGAR
          val cigarElements = alignment.toCigar.getCigarElements
          cigarElements.flatMap(cigarElement => {
            val cigarOperator = cigarElement.getOperator
            val referenceLength = CigarUtils.getReferenceLength(cigarElement)
            val pathLength = CigarUtils.getReadLength(cigarElement)

            // Yield a resulting variant when there is a mismatch, insertion or deletion
            val possibleVariant = cigarOperator match {
              case CigarOperator.X =>
                val referenceAllele = currentReference.slice(referenceIndex, referenceIndex + referenceLength)
                val alternateAllele = path.slice(pathIndex, pathIndex + pathLength)
                Some(buildVariant(referenceIndex, referenceAllele, alternateAllele.toArray))
              case (CigarOperator.I | CigarOperator.D) =>
                // For insertions and deletions, report the variant with the last reference base attached
                val referenceAllele = currentReference.slice(referenceIndex - 1, referenceIndex + referenceLength)
                val alternateAllele = path.slice(pathIndex - 1, pathIndex + pathLength)
                Some(buildVariant(referenceIndex - 1, referenceAllele, alternateAllele.toArray))
              case _ => None
            }
            referenceIndex += referenceLength
            pathIndex += pathLength

            possibleVariant

          })
        })

      (graph, variants.iterator)
    }

    override def run(args: Arguments, sc: SparkContext): Unit = {

      val loci = Common.loci(args)
      val readSet = Common.loadReadsFromArguments(
        args,
        sc,
        Read.InputFilters(mapped = true, nonDuplicate = true))

      val lociPartitions = DistributedUtil.partitionLociAccordingToArgs(
        args,
        loci.result(readSet.contigLengths),
        readSet.mappedReads
      )

      val kmerSize = args.kmerSize
      val minAlignmentQuality = args.minAlignmentQuality
      val minAverageBaseQuality = args.minAverageBaseQuality
      val reference = ReferenceBroadcast(args.referenceFastaPath, sc)

      val qualityReads = readSet
        .mappedReads
        .filter(_.alignmentQuality > minAlignmentQuality)

      val minAreaVaf = args.minAreaVaf
      // Find loci where there are variant reads
      val lociOfInterest = DistributedUtil.pileupFlatMap[VariantLocus](
        qualityReads,
        lociPartitions,
        skipEmpty = true,
        function = pileup => VariantLocus(pileup).iterator,
        reference = Some(reference)
      ).filter(_.variantAlleleFrequency > (minAreaVaf / 100.0))

      Common.progress(s"Found ${lociOfInterest.count} loci with non-reference reads")

      // Re-examine loci with variant read
      val lociOfInterestSet = LociSet.union(
        lociOfInterest.map(
          locus => LociSet(locus.contig, locus.locus, locus.locus + 1))
          .collect(): _*
      )
      val lociOfInterestPartitions = DistributedUtil.partitionLociUniformly(
        args.parallelism,
        lociOfInterestSet
      )

      val minOccurrence = args.minOccurrence
      val genotypes: RDD[CalledAllele] =
        DistributedUtil.windowFlatMapWithState[MappedRead, CalledAllele, Option[DeBruijnGraph]](
          Seq(qualityReads),
          lociOfInterestPartitions,
          skipEmpty = true,
          halfWindowSize = args.snvWindowRange,
          initialState = None,
          (graph, windows) => {
            discoverHaplotypes(
              graph,
              windows(0),
              kmerSize,
              reference,
              minOccurrence
            )
          }
        )
      genotypes.persist()

      Common.progress(s"Found ${genotypes.count} variants")

      val outputGenotypes =
        genotypes.flatMap(AlleleConversions.calledAlleleToADAMGenotype)

      Common.writeVariantsFromArguments(args, outputGenotypes)
      DelayedMessages.default.print()
    }

  }
}
