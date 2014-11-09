package org.hammerlab.guacamole.commands

import org.apache.spark.SparkContext._
import org.apache.spark.rdd.RDD
import org.apache.spark.{ Logging, SparkContext }
import org.hammerlab.guacamole.reads.Read
import org.hammerlab.guacamole.{ ReadSet, DistributedUtil, SparkCommand }
import org.kohsuke.args4j.{ Argument, Option }

object ReadDepthHist {

  protected class Arguments extends DistributedUtil.Arguments {

    @Argument(metaVar = "OUTFILE", required = true, usage = "BAM file to compute histogram for.")
    var inputFile: String = null

    @Option(name = "-repartition", metaVar = "N", required = false, usage = "Number of shards to repartition to.")
    var repartition: Int = -1
  }

  object Cmd extends SparkCommand[Arguments] with Serializable with Logging {

    override val name = "read-depth-hist"
    override val description = "Compute a histogram of read depths"

    override def run(args: Arguments, sc: SparkContext): Unit = {

      val readSet = ReadSet(sc, args.inputFile, Read.InputFilters.empty, 1, true)
      val reads = readSet.reads

      val unmappedReads = reads.filter(!_.isMapped)
      val mappedReads =
        if (args.repartition > 0)
          reads.flatMap(_.getMappedReadOpt).coalesce(args.repartition, shuffle = true)
        else
          reads.flatMap(_.getMappedReadOpt)

      val readDepthPerLocus: RDD[((String, Long), Long)] =
        mappedReads.flatMap(read => {
          (0 until read.sequence.size).map(offset =>
            ((read.referenceContig, read.start + offset), 1L)
          )
        }).reduceByKey(_ + _)

      val lociPerReadDepth: collection.Map[Long, Long] =
        readDepthPerLocus.map({
          case (locus, count) => (count, 1L)
        }).reduceByKeyLocally(_ + _)

      println("Loci per read depth:\n\n%s".format(
        lociPerReadDepth.toList.sortBy(_._1).map({
          case (depth, numLoci) => "%8d: %d".format(depth, numLoci)
        }).mkString("\t", "\n\t", "")
      ))
    }
  }
}
