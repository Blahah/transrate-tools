#include <stdio.h>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include "api/BamReader.h"
#include "pileup.h"
#include "atomicops.h"
#include "readerwriterqueue.h"

using namespace BamTools;

using namespace std;

using namespace moodycamel;

typedef vector<BamAlignment> Batch;
typedef BlockingReaderWriterQueue<Batch> BatchQueue;

double nullprior = 0.7;

class BamRead
{
public:
  int estimate_fragment_size(std::string file);
  void load_bam(std::string, int n_threads);
  void load_bam_header();
  void process_alignment(BamAlignment alignment);
  int bar;
  int seq_count;
  uint32_t nm_tag;
  int ldist;
  int realistic_distance;
  int refid;
  BamReader reader;
  std::vector<TransratePileup> array;
};

// function that each thread runs
// this function gets passed a queue and just processes it until
// it gets an empty batch of alignments
void process_queue(BamRead &br, BatchQueue &queue, bool done) {

  Batch batch;
  queue.wait_dequeue(batch);
  while (batch.size() > 0) {
    for (auto &alignment : batch) {
      br.process_alignment(alignment);
    }
    // free up the memory for this batch
    batch.clear();
    queue.wait_dequeue(batch);
  }

}

// do final summary processing of each contig
void process_contigs(BamRead &br, int thread_id) {
  for (int i = 0; i < br.seq_count; ++i) {
    if (i % thread_id == 0) {
      br.array[i].calculateUncoveredBases();
      br.array[i].setPNotSegmented();
    }
  }
}
