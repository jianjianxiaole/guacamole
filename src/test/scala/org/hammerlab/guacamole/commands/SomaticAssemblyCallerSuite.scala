package org.hammerlab.guacamole.commands

import org.hammerlab.guacamole.commands.SomaticAssemblyCaller.Arguments
import org.hammerlab.guacamole.commands.jointcaller.{InputCollection, SomaticJoint}
import org.hammerlab.guacamole.data.CancerWGSTestUtil
import org.hammerlab.guacamole.distributed.LociPartitionUtils
import org.hammerlab.guacamole.loci.set.LociParser
import org.hammerlab.guacamole.reference.ReferenceBroadcast
import org.hammerlab.guacamole.util.{Bases, GuacFunSuite, TestUtil}
import org.hammerlab.guacamole.variants.CalledSomaticAllele
import org.scalatest.{BeforeAndAfterAll, Matchers}

class SomaticAssemblyCallerSuite extends GuacFunSuite with Matchers with BeforeAndAfterAll {

  val args = new Arguments
  args.parallelism = 1

  var reference: ReferenceBroadcast = _

  val partialFasta = TestUtil.testDataPath("hg19.partial.fasta")

  override def beforeAll() {
    super.beforeAll()
    reference = ReferenceBroadcast(partialFasta, sc, partialFasta = true)
  }

  val inputs = InputCollection(CancerWGSTestUtil.bams)
  def verifyVariantsAtLocus(locus: Int,
                            contig: String = "chr1",
                            kmerSize: Int = 31,
                            snvWindowRange: Int = 45,
                            minOccurrence: Int = 3,
                            minVaf: Float = 0.1f,
                            shortcutAssembly: Boolean = false)(
                             expectedVariants: (String, Int, String, String)*
                           ) = {

    val windowStart = locus - snvWindowRange
    val windowEnd = locus + snvWindowRange

    val lociBuilder = LociParser(s"$contig:$windowStart-$windowEnd")

    val readSets = SomaticJoint.inputsToReadSets(sc, inputs, lociBuilder)
    val normalReadSet = readSets(0)
    val tumorReadSet = readSets(1)

    val lociPartitions = LociPartitionUtils.partitionLociUniformly(
      numPartitions = args.parallelism,
      loci = lociBuilder.result(readSets.contigLengths)
    )

    val variants =
      SomaticAssemblyCaller.Caller.discoverSomaticGenotypes(
        normalReadSet.mappedReads,
        tumorReadSet.mappedReads,
        kmerSize = kmerSize,
        snvWindowRange = snvWindowRange,
        minOccurrence = minOccurrence,
        minAreaVaf = minVaf,
        reference = reference,
        lociPartitions = lociPartitions,
        minDepth = 8,
        minLikelihood = 30,
        oddsThreshold = 20,
        shortcutAssembly = shortcutAssembly
      ).collect().sortBy(_.start)

    val actualVariants =
      for {
        CalledSomaticAllele(_, contig, start, allele, _, _, _, _, _) <- variants
      } yield {
        (contig, start, Bases.basesToString(allele.refBases), Bases.basesToString(allele.altBases))
      }

    actualVariants should be(expectedVariants)

  }


  test("test somatic assembly caller: simple snv 1") {
    verifyVariantsAtLocus(65857040, "chr12")(
      ("chr12", 65857040, "G", "C")
    )
  }

  test("test somatic assembly caller: simple snv 2") {
    verifyVariantsAtLocus(82833487, "chr5")(
      ("chr5", 82833487, "G", "A")
    )
  }

  test("test somatic assembly caller: simple snv 3") {
    verifyVariantsAtLocus(22672428, "chr7")(
      ("chr7", 22672428, "A", "G")
    )
  }

  test("test somatic assembly caller: simple snv 4") {
    verifyVariantsAtLocus(217142481, "chr2")(
      ("chr2", 217142481, "A", "G")
    )
  }

  test("test somatic assembly caller: simple snv 5") {
    verifyVariantsAtLocus(68691575, "chr4")(
      ("chr4", 68691575, "C", "A")
    )
  }

  test("test somatic assembly caller: simple deletion") {
    verifyVariantsAtLocus(82649007, "chr5")(
      ("chr5", 82649007, "TCTTTAGAAA", "T")
    )
  }

  test("test somatic assembly caller: simple deletion 2") {
    verifyVariantsAtLocus(158015819, "chr3")(
      ("chr3", 158015819, "TA", "T")
    )
  }

  test("test somatic assembly caller: exclude germline variants") {

    args.variantOutput = TestUtil.tmpPath(".vcf")
    verifyVariantsAtLocus(38019746, "chr1")(
      // No somatic variants in this region
    )

    verifyVariantsAtLocus(120270638, "chr12")(
      // No somatic variants in this region
    )
  }

  test("test somatic assembly caller: call somatic variant near germline variant") {
    // Germline variant G>A at 4:120371251
    // Somatic variant G>C at 4:120371197
    verifyVariantsAtLocus(120371197, "chr4", snvWindowRange = 60)(
      ("chr4", 120371197, "G", "C")
    )
  }
}
