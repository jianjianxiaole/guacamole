/**
 * Licensed to Big Data Genomics (BDG) under one
 * or more contributor license agreements.  See the NOTICE file
 * distributed with this work for additional information
 * regarding copyright ownership.  The BDG licenses this file
 * to you under the Apache License, Version 2.0 (the
 * "License"); you may not use this file except in compliance
 * with the License.  You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package org.hammerlab.guacamole.reads

import htsjdk.samtools.TextCigarCodec
import org.hammerlab.guacamole.util.{AssertBases, GuacFunSuite, TestUtil}
import org.hammerlab.guacamole.util.TestUtil.Implicits._
import org.hammerlab.guacamole.util.TestUtil.makeRead

class MappedReadSuite extends GuacFunSuite {

  test("mappedread is mapped") {
    val read = MappedRead(
      "read1",
      "TCGACCCTCGA",
      Array[Byte]((10 to 20).map(_.toByte): _*),
      true,
      "some sample name",
      "chr5",
      50,
      325352323,
      TextCigarCodec.decode(""),
      false,
      isPositiveStrand = true,
      isPaired = true
    )

    read.isMapped should be(true)
    read.asInstanceOf[Read].isMapped should be(true)

  }

  test("mixed collections mapped and unmapped reads") {
    val uread = UnmappedRead(
      "read1",
      "TCGACCCTCGA",
      Array[Byte]((10 to 20).map(_.toByte): _*),
      isDuplicate = true,
      "some sample name",
      failedVendorQualityChecks = false,
      isPaired = true
    )

    val mread = MappedRead(
      "read1",
      "TCGACCCTCGA",
      Array[Byte]((10 to 20).map(_.toByte): _*),
      isDuplicate = true,
      "some sample name",
      "chr5",
      50,
      325352323,
      TextCigarCodec.decode(""),
      failedVendorQualityChecks = false,
      isPositiveStrand = true,
      isPaired = true
    )

    val collectionMappedReads: Seq[Read] = Seq(uread, mread)
    collectionMappedReads(0).isMapped should be(false)
    collectionMappedReads(1).isMapped should be(true)
  }

  // This must only be accessed from inside a spark test where SparkContext has been initialized
  def reference = TestUtil.makeReference(sc,
    Seq(
      ("chr1", 8, "GGTCGATCGATCAA")
    ))

  test("slice read matching read") {
    val chr1Contig = reference.getContig("chr1")
    val readLength = 10
    val qualityScores = (0 until readLength).map( _ => 30)

    val simpleRead = makeRead(
      "TCGATCGATC",
      start = 10,
      cigar = "10M",
      qualityScores = Some(qualityScores)
    )

    val sliceAll = simpleRead.slice(10L, 20L, chr1Contig).get
    sliceAll should be(simpleRead)

    val sliceNone = simpleRead.slice(20L, 30L, chr1Contig)
    sliceNone should be (None)

    val sliceFirstFive = simpleRead.slice(10L, 15L, chr1Contig).get
    sliceFirstFive.start should be (10L)
    AssertBases(sliceFirstFive.sequence, "TCGAT")
    sliceFirstFive.cigar.toString should be ("5M")
    sliceFirstFive.end should be (15L)

    val sliceLastFive = simpleRead.slice(15L, 20L, chr1Contig).get
    sliceLastFive.start should be (15L)
    AssertBases(sliceLastFive.sequence, "CGATC")
    sliceLastFive.cigar.toString should be ("5M")
    sliceLastFive.end should be (20L)

  }

  test("slice read with deletion") {
    val chr1Contig = reference.getContig("chr1")
    val readLength = 10
    val qualityScores = (0 until readLength).map( _ => 30)

    val deletionRead = makeRead(
      "GGTCGATCAA",
      start = 8,
      cigar = "6M4D4M",
      qualityScores = Some(qualityScores)
    )

    val sliceBeforeDeletion = deletionRead.slice(11L, 20L, chr1Contig).get
    sliceBeforeDeletion.start should be(11L)
    AssertBases(sliceBeforeDeletion.sequence, "CGATC")
    sliceBeforeDeletion.cigar.toString should be ("3M4D2M")
    sliceBeforeDeletion.end should be (20L)

    val sliceInDeletion = deletionRead.slice(16L, 20L, chr1Contig).get
    sliceInDeletion.start should be(16L)
    AssertBases(sliceInDeletion.sequence, "TC")
    sliceInDeletion.cigar.toString should be ("2D2M")
    sliceInDeletion.end should be (20L)

  }

  test("slice read with insertion") {
    val chr1Contig = reference.getContig("chr1")
    val readLength = 15
    val qualityScores = (0 until readLength).map( _ => 30)

    val insertionRead = makeRead(
      "TCGACCCCCTCGATC",
      start = 10,
      cigar = "4M5I6M",
      qualityScores = Some(qualityScores)
    )

    val sliceBeforeInsertion = insertionRead.slice(12L, 18L, chr1Contig).get
    sliceBeforeInsertion.start should be(12L)
    AssertBases(sliceBeforeInsertion.sequence, "GACCCCCTCGA")
    sliceBeforeInsertion.cigar.toString should be ("2M5I4M")
    sliceBeforeInsertion.end should be (18L)

    val sliceAfterInsertion = insertionRead.slice(14L, 18L, chr1Contig).get
    sliceAfterInsertion.start should be(14L)
    AssertBases(sliceAfterInsertion.sequence, "TCGA")
    sliceAfterInsertion.cigar.toString should be ("4M")
    sliceAfterInsertion.end should be (18L)
  }

}
