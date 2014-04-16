package org.bdgenomics.guacamole

import org.bdgenomics.adam.avro.ADAMRecord
import scala.collection.mutable
import org.bdgenomics.adam.rich.{DecadentRead, RichADAMRecord}


case class SlidingReadWindow(windowSize: Long, rawSortedReads: Iterator[ADAMRecord]) {
  var currentLocus = -1
  private var referenceName: Option[String] = None
  private var mostRecentReadStart: Long = 0
  private val sortedReads: Iterator[DecadentRead] = rawSortedReads.map(read => {
    require(read.getReadMapped, "Reads must be mapped")
    if (referenceName.isEmpty) referenceName = Some(read.getReferenceName.toString)
    require(read.getReferenceName == referenceName.get, "Reads must have the same reference name")
    require(read.getStart >= mostRecentReadStart, "Reads must be sorted by start locus")
    require(read.getCigar.length > 1, "Reads must have a CIGAR string")
    DecadentRead(read)
  })

  val currentReads = {
    // Order reads by end locus, increasing.
    def readOrdering = new Ordering[DecadentRead] {
      def compare(first: DecadentRead, second: DecadentRead) = second.record.end.get.compare(first.record.end.get)
    }
    new mutable.PriorityQueue[DecadentRead]()(readOrdering)
  }

  def setCurrentLocus(locus: Long): Iterator[DecadentRead] = {
    assume(locus >= currentLocus, "Pileup window can only move forward in locus")

    def overlaps(read: DecadentRead) = {
      (read.record.getStart >= locus - windowSize && read.record.getStart <= locus + windowSize) ||
      (read.record.end.get >= locus - windowSize && read.record.end.get <= locus + windowSize)
    }

    // Remove reads that are no longer in the window.
    while (!currentReads.isEmpty && currentReads.head.record.end.get < locus - windowSize) {
      val dropped = currentReads.dequeue()
      assert(!overlaps(dropped))
    }
    // Add new reads that are now in the window.
    val newReads = sortedReads.takeWhile(_.record.getStart <= locus + windowSize).filter(overlaps)
    currentReads.enqueue(newReads.toSeq :_*)
    assert(currentReads.forall(overlaps))  // Correctness check.
    newReads // We return the newly added reads.
  }
}
object SlidingReadWindow {
  

}