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

typedef ReaderWriterQueue<BamAlignment> AlnQueue;

double nullprior = 0.7;

class BamRead
{
public:
  BamRead(void);
  BamRead(const BamRead &in) {
    bar = in.bar;
    seq_count = in.seq_count;
    nm_tag = in.nm_tag;
    ldist = in.ldist;
    realistic_distance = in.realistic_distance;
    array = in.array;
  }
  ~BamRead(void);
  int estimate_fragment_size(std::string file);
  void load_bam(std::string, int n_threads);
  void load_bam_header();
  void process_alignment(BamAlignment alignment);
  int bar;
  int seq_count;
  uint32_t nm_tag;
  int ldist;
  int realistic_distance;
  BamReader reader;
  std::vector<TransratePileup> array;
};

// function that each thread runs
// this function gets passed a queue and just processes it until
// it gets an empty batch of alignments
void process_queue(BamRead &br, AlnQueue &queue) {

  BamAlignment alignment;
  while (queue.try_dequeue(alignment)) {
    br.process_alignment(alignment);
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
