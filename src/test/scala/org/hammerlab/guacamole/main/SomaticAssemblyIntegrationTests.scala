package org.hammerlab.guacamole.main

import org.apache.spark.SparkContext
import org.hammerlab.guacamole.commands.SomaticAssemblyCaller.Arguments
import org.hammerlab.guacamole.commands.jointcaller.{InputCollection, SomaticJoint}
import org.hammerlab.guacamole.commands.{SomaticAssemblyCaller, SparkCommand}
import org.hammerlab.guacamole.data.CancerWGSTestUtil
import org.hammerlab.guacamole.distributed.LociPartitionUtils
import org.hammerlab.guacamole.distributed.LociPartitionUtils._
import org.hammerlab.guacamole.loci.set.LociParser
import org.hammerlab.guacamole.readsets.ReadsRDD
import org.hammerlab.guacamole.reference.ReferenceBroadcast
import org.hammerlab.guacamole.util.TestUtil
import org.hammerlab.guacamole.variants.{AlleleConversions, VariantComparisonTest, VariantUtils}

object SomaticAssemblyIntegrationTests extends SparkCommand[Arguments] with VariantComparisonTest {

  override val name: String = "somatic-assembly-integration-test"
  override val description: String = "output various statistics to stdout"

  val callerArgs = new Arguments()

  // Read loading config
  callerArgs.bamReaderAPI = "hadoopbam"
  callerArgs.parallelism = 1

  // Somatic assembly config
  callerArgs.kmerSize = 37
  callerArgs.assemblyWindowRange = 120
  callerArgs.minAreaVaf = 30
  callerArgs.minOccurrence = 3
  callerArgs.minLikelihood = 50

  callerArgs.shortcutAssembly = false

  def main(args: Array[String]): Unit = run(callerArgs)

  override def run(args: Arguments, sc: SparkContext): Unit = {

    val loci = LociParser(((1).until(22).map(i => "chr%d".format(i)) ++ Seq("chrX", "chrY")).mkString(","))

    val inputs = InputCollection(CancerWGSTestUtil.bams)
    val readSets = SomaticJoint.inputsToReadSets(sc, inputs, loci)
    val normalReadSet = readSets(0)
    val primaryTumorReadSet = readSets(1)
    val recurrenceTumorReadSet = readSets(2)

    val lociPartitions = LociPartitionUtils.partitionLociUniformly(
      numPartitions = args.parallelism,
      loci = loci.result(readSets.contigLengths)
    )

    val partialFasta = TestUtil.testDataPath("hg19.partial.fasta")
    val reference = ReferenceBroadcast(partialFasta, sc, partialFasta = true)

    args.variantOutput = "/tmp/somatic-assembly-cancer-wgs-primary-guacamole-tests.vcf"
    makeCallsAndCompare(sc, args, normalReadSet, primaryTumorReadSet, lociPartitions, reference, "primary")

//    args.variantOutput = "/tmp/somatic-assembly-cancer-wgs-recurrence-guacamole-tests.vcf"
//    makeCallsAndCompare(sc, args, normalReadSet, recurrenceTumorReadSet, lociPartitions, reference, "recurrence")

  }

  def makeCallsAndCompare(sc: SparkContext,
                          args: Arguments,
                          normalReads: ReadsRDD,
                          tumorReads: ReadsRDD,
                          lociPartitions: LociPartitioning,
                          reference: ReferenceBroadcast,
                          comparisonSample: String): Unit = {
    val calls =
      SomaticAssemblyCaller.Caller.discoverSomaticGenotypes(
        normalReads.mappedReads.filter(_.alignmentQuality > args.minAlignmentQuality),
        tumorReads.mappedReads.filter(_.alignmentQuality > args.minAlignmentQuality),
        kmerSize = args.kmerSize,
        assemblyWindowRange = args.assemblyWindowRange,
        minOccurrence = args.minOccurrence,
        minAreaVaf = (args.minAreaVaf / 100.0f),
        reference = reference,
        lociPartitions = lociPartitions,
        minDepth = args.minReadDepth,
        minMeanKmerQuality = args.minMeanKmerQuality,
        minLikelihood = args.minLikelihood,
        oddsThreshold = args.oddsThreshold,
        minAltReads = 3,
        shortcutAssembly = args.shortcutAssembly
      )

    VariantUtils.writeVariantsFromArguments(
      args,
      calls.flatMap(AlleleConversions.calledSomaticAlleleToADAMGenotype)
    )

    println("************* CANCER WGS1 SOMATIC CALLS *************")
    compareToCSV(
      args.variantOutput + "/part-r-00000",
      CancerWGSTestUtil.expectedSomaticCallsCSV,
      CancerWGSTestUtil.referenceBroadcast(sc),
      Set(comparisonSample)
    )
  }
}
