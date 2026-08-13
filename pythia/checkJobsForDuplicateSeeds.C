#include <vector>
#include <iostream>
#include <typeinfo>

#include "TCanvas.h"
#include "TError.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH3D.h"
#include "THnSparse.h"
#include "TString.h"
#include "TLegend.h"

struct SeedStruct {
  private:
    int _batch;
    int _job;
    int _seed;
  public:
    int getBatch() { return _batch; }
    int getJob() { return _job; }
    int getSeed() { return _seed; }

    void print() { string s = TString::Format("%d_%d: %d", _batch, _job, _seed).Data(); cout << s << "\n"; }

    SeedStruct(int b, int j, int s) : _batch(b), _job(j), _seed(s) {}
};

bool lineContainsSeed(string line) {
  string prepending = "Using seed ";
  if (line.find(prepending) != std::string::npos) {
    return true;
  }
  return false;
}

int getSeedFromLine(string line) {
  string prepending = "Using seed ";
  string seedString = line;
  seedString.replace(0, prepending.size(), ""); // Isolate the seed
  return std::stoi(seedString);
}

// For each batch, get all jobs with corresponding seeds
vector<SeedStruct> getSeedStructs(vector<int> batches) {
  vector<SeedStruct> seedStructs;
  int nJobsMax = 1e3; // Just to have a max for the loop
  for (int batch : batches) {
    for (int job = 0; job < nJobsMax; job++) {
      string fn = TString::Format("V0JetClustering/Batch%d/%d_%d.out", batch, batch, job).Data();
      ifstream f(fn);
      string line;
      if (f.is_open()) {
        while (getline(f, line)) {
          if (lineContainsSeed(line)) {
            int seed = getSeedFromLine(line);
            SeedStruct s(batch, job, seed);
            seedStructs.push_back(s);
          }
        }
      } else {
        // cout << "Can\'t open file " << fn << "\n";
        break;
      }
    } // for job
  } // for batch
  return seedStructs;
}

// Given a set of jobs and seeds, find the duplicates
// Label them so we know exactly which seeds/jobs are duplicates of which
void printDuplicates(vector<SeedStruct> seedStructs) {
  int nSeeds = seedStructs.size();
  
  // Find the duplicates
  int uniqueSeedId = -1;
  int duplicateSetID = uniqueSeedId + 1;
  vector<int> duplicateIterators(nSeeds, uniqueSeedId);

  // Don't need to double-check, so j can start at i+1
  // i = n - 1 leaves no values for j, so can be skipped
  for (int i = 0; i < nSeeds - 2; i++) {
    if (duplicateIterators[i] != uniqueSeedId) // Already matched with something
      continue;

    SeedStruct ss1 = seedStructs[i];
    for (int j = i + 1; j < nSeeds - 1; j++) {
      if (duplicateIterators[j] != uniqueSeedId) // Already matched with something
        continue;

      SeedStruct ss2 = seedStructs[j];
      if (ss1.getSeed() == ss2.getSeed()) {
        duplicateIterators[i] = duplicateSetID;
        duplicateIterators[j] = duplicateSetID;
      }
    }

    // Change duplicateSetID only if i is a duplicate, to keep them identifiable
    if (duplicateIterators[i] != uniqueSeedId)
      duplicateSetID++;
  }

  if (false) {
    cout << "[";
    for (int dit : duplicateIterators)
      cout << dit << ", ";

    cout << "]\n";
  }

  // Print the duplicates
  duplicateSetID = uniqueSeedId + 1; // Reset to loop over all sets
  for (int did : duplicateIterators) {
    if (did < duplicateSetID) { // Filter out uniques and sets we've passed already
      continue;
    } else {
      vector<int> duplicateIndices;
      // Collect duplicates in a vector
      for (int index = 0; index < nSeeds - 1; index++) {
        if (duplicateIterators[index] == did)
          duplicateIndices.push_back(index);
      }
      cout << "Duplicate seeds found in jobs: (duplicate set " << duplicateSetID <<  ")\n";
      for (int index : duplicateIndices) {
        seedStructs[index].print();
      }
      cout << endl;
      duplicateSetID++;
    }
  }
}

void checkBatchesForDuplicates(vector<int> batches) {
  vector<SeedStruct> seedStructs = getSeedStructs(batches);
  printDuplicates(seedStructs);
}

void checkBatchesEscheme() {
  vector<int> batches = { 4861164, 5190393, 5190859 };
  checkBatchesForDuplicates(batches);
}
void checkBatchesPtscheme() {
  vector<int> batches = { 4861167, 5190394, 5190860 };
  checkBatchesForDuplicates(batches);
}