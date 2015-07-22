package org.hammerlab.guacamole.commands

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics
import org.apache.spark.SparkContext
import org.apache.spark.mllib.clustering.{ GaussianMixture, GaussianMixtureModel }
import org.apache.spark.mllib.linalg.Vectors
import org.apache.spark.rdd.RDD
import org.apache.spark.storage.StorageLevel
import org.hammerlab.guacamole._
import org.hammerlab.guacamole.filters.PileupFilter.PileupFilterArguments
import org.hammerlab.guacamole.pileup.Pileup
import org.hammerlab.guacamole.reads.Read.InputFilters
import org.hammerlab.guacamole.reads.{ MappedRead, Read }
import org.kohsuke.args4j.{ Argument, Option => Args4jOption }

/**
 * VariantLocus is locus and the variant allele frequency at that locus
 * @param locus Position of non-reference alleles
 * @param variantAlleleFrequency Frequency of non-reference alleles
 */
case class VariantLocus(locus: Long, variantAlleleFrequency: Float)

object VariantLocus {

  /**
   * Construct VariantLocus from a pileup
   * @param pileup Pileup of reads at a given locus
   * @return VariantLocus at reference position, locus
   */
  def apply(pileup: Pileup): Option[VariantLocus] = {
    if (pileup.referenceDepth != pileup.depth) {
      Some(VariantLocus(pileup.locus, (pileup.depth - pileup.referenceDepth).toFloat / pileup.depth))
    } else {
      None
    }
  }

}

object VAFHistogram {

  protected class Arguments extends DistributedUtil.Arguments with PileupFilterArguments {

    @Args4jOption(name = "--out", required = false,
      usage = "Path to save the variant allele frequency histogram. (Print to screen if not provided)")
    var output: String = ""

    @Args4jOption(name = "--bins", required = false,
      usage = "Number of bins for the variant allele frequency histogram (Default: 20)")
    var bins: Int = 20

    @Args4jOption(name = "--cluster", required = false,
      usage = "Cluster the variant allele frequencies using a Gaussian mixture model")
    var cluster: Boolean = false

    @Args4jOption(name = "--num-clusters", required = false, depends = Array("--cluster"),
      usage = "Number of clusters for the Gaussian mixture model (Default: 3)")
    var numClusters: Int = 3

    @Args4jOption(name = "--samplePercent", usage = "Percent of variant to use for the calculations (Default: 25)")
    var samplePercent: Int = 25

    @Argument(required = true, multiValued = true,
      usage = "BAMs")
    var bams: Array[String] = Array.empty

  }

  object Caller extends SparkCommand[Arguments] {
    override val name = "vaf-histogram"
    override val description = "Compute and cluster the variant allele frequencies"

    override def run(args: Arguments, sc: SparkContext): Unit = {

      val filters = Read.InputFilters(mapped = true, nonDuplicate = true, passedVendorQualityChecks = true)
      val samplePercent = args.samplePercent

      val readSets: Seq[ReadSet] = args.bams.zipWithIndex.map(
        bamFile =>
          ReadSet(
            sc,
            bamFile._1,
            requireMDTagsOnMappedReads = false,
            InputFilters.empty,
            token = bamFile._2,
            contigLengthsFromDictionary = false)
      )

      val loci = Common.loci(args, readSets(0))
      val lociPartitions = DistributedUtil.partitionLociAccordingToArgs(
        args,
        loci,
        readSets(0).mappedReads // Use the first set of reads as a proxy for read depth
      )

      val variantLoci = readSets.map(readSet =>
        variantLociFromReads(
          readSet.mappedReads,
          lociPartitions,
          samplePercent)
      )

      val bins = args.bins
      val variantAlleleHistograms =
        variantLoci.map(variantLoci => generateVAFHistogram(variantLoci, bins)) // Sort by variant allele frequency

      val sampleNames = readSets.map(_.mappedReads.take(1)(0).sampleName)
      if (args.output != "") {
        sampleNames.zip(variantAlleleHistograms).foreach(kv => {
          val sampleName = kv._1
          val histogram = kv._2
          // Parallelize histogram and save on HDFS
          val histogramStrings = histogram.toSeq.sortBy(_._1).map(kv => f"VAF: ${kv._1} -> ${kv._2}").toSeq
          sc.parallelize(histogramStrings).saveAsTextFile(args.output + '-' + sampleName)
        })
      } else {
        // Print histograms to standard out
        variantAlleleHistograms.foreach(histogram =>
          histogram.toSeq.sortBy(_._1).foreach(kv => println(f"VAF: ${kv._1} -> ${kv._2}"))
        )
      }

      if (args.cluster) {
        val numClusters = args.numClusters
        variantLoci.foreach(buildMixtureModel(_, numClusters))
      }

    }

  }

  /**
   * Generates a count of loci in each variant allele frequency bins
   * @param variantAlleleFrequencies RDD of loci with variant allele frequency > 0
   * @param bins Number of bins to group the VAFs into
   * @return Map of rounded variant allele frequency to number of loci with that value
   */
  def generateVAFHistogram(variantAlleleFrequencies: RDD[VariantLocus], bins: Int): Map[Int, Long] = {

    def roundToBin(variantAlleleFrequency: Float) = {
      val variantPercent = (variantAlleleFrequency * 100).toInt
      variantPercent - (variantPercent % (100 / bins))
    }
    variantAlleleFrequencies.keyBy(vaf => roundToBin(vaf.variantAlleleFrequency)).countByKey().toMap
  }

  /**
   * Find all non-reference loci in the sample
   * @param reads RDD of mapped reads for the sample
   * @param lociPartitions Positions which to examine for non-reference loci
   * @param samplePercent Percent of non-reference loci to use for descriptive statistics
   * @return RDD of VariantLocus, which contain the locus and non-zero variant allele frequency
   */
  def variantLociFromReads(reads: RDD[MappedRead],
                           lociPartitions: LociMap[Long],
                           samplePercent: Int = 100): RDD[VariantLocus] = {
    val sampleName = reads.take(1)(0).sampleName
    val variantLoci = DistributedUtil.pileupFlatMap[VariantLocus](
      reads,
      lociPartitions,
      skipEmpty = true,
      pileup => VariantLocus(pileup).iterator
    )
    variantLoci.persist(StorageLevel.MEMORY_ONLY)

    val numVariantLoci = variantLoci.count
    Common.progress("%d non-zero variant loci in sample %s".format(numVariantLoci, sampleName))

    // Sample variant loci to compute descriptive statistics
    val sampledVAFs =
      if (samplePercent < 100)
        variantLoci
          .sample(withReplacement = false, fraction = samplePercent / 100.0)
          .collect()
      else
        variantLoci.collect()

    val stats = new DescriptiveStatistics()
    sampledVAFs.foreach(v => stats.addValue(v.variantAlleleFrequency))

    // Print out descriptive statistics for the variant allele frequency distribution
    Common.progress("Variant loci stats for %s (min: %f, max: %f, median: %f, mean: %f, 25Pct: %f, 75Pct: %f)".format(
      sampleName,
      stats.getMin,
      stats.getMax,
      stats.getPercentile(50),
      stats.getMean,
      stats.getPercentile(25),
      stats.getPercentile(75)
    ))

    variantLoci
  }

  /**
   * Fit a Gaussian mixture model to the distribution of variant allele frequencies
   * @param variantAlleleFrequencies RDD of loci with variant allele frequency > 0
   * @param numClusters Number of Gaussian distributions to fit
   * @param maxIterations Maximum number of iterations to run EM
   * @param convergenceTol Largest change in log-likelihood before convergence
   * @return GaussianMixtureModel
   */
  def buildMixtureModel(variantAlleleFrequencies: RDD[VariantLocus],
                        numClusters: Int,
                        maxIterations: Int = 50,
                        convergenceTol: Double = 1e-2): GaussianMixtureModel = {
    val vafVectors = variantAlleleFrequencies.map(vaf => Vectors.dense(vaf.variantAlleleFrequency))
    val model = new GaussianMixture()
      .setK(numClusters)
      .setConvergenceTol(convergenceTol)
      .setMaxIterations(maxIterations)
      .run(vafVectors)

    for (i <- 0 until model.k) {
      println(s"Cluster $i: mean=${model.gaussians(i).mu(0)}, std. deviation=${model.gaussians(i).sigma}, weight=${model.weights(i)}")
    }

    model
  }

}
