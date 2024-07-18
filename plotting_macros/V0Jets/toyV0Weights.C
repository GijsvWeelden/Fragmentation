
#include <vector>
#include <iostream>
#include <typeinfo>

#include "TCanvas.h"
#include "TError.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"
#include "TH3D.h"
#include "THnSparse.h"
#include "TString.h"
#include "TLegend.h"

// #include "../histUtils.C"

// Sourcre: https://stackoverflow.com/questions/1700079/howto-create-combinations-of-several-vectors-without-hardcoding-loops-in-c
/*
void printAll(const vector<vector<string> > &allVecs, size_t vecIndex, string strSoFar, vector<string>& result)
{
    if (vecIndex >= allVecs.size()) {
        result.push_back(strSoFar);
        return;
    }
    for (size_t i=0; i<allVecs[vecIndex].size(); i++) {
        printAll(allVecs, vecIndex+1, strSoFar+allVecs[vecIndex][i], result);
    }
}

void printAll(const vector<vector<string> > &allVecs, size_t vecIndex, vector<string> strSoFar, vector<vector<string>>& result)
{
    if (vecIndex >= allVecs.size()) {
        // This is where the result lives
        for (auto s : strSoFar) {
            // cout << s << ", ";
        }
        result.push_back(strSoFar);
        return;
    }
    for (size_t i=0; i<allVecs[vecIndex].size(); i++) {
      // cout << "vecIndex: " << vecIndex << ", i: " << i << ", " << allVecs[vecIndex][i] << endl;
      strSoFar.push_back(allVecs[vecIndex][i]);
      printAll(allVecs, vecIndex+1, strSoFar, result);
    }
}

void printA()
{
  vector<vector<string> > allVecs;
  int n = 3; // number of vectors

  for (int i = 0; i < n; i++) {
    vector<string> vec = {"p"+to_string(i), "q"+to_string(i)};
    allVecs.push_back(vec);
  }

  vector<string> result;
  vector<vector<string> > splitResult;
  printAll(allVecs, 0, "", result);

  cout << "Size of result: " << result.size() << endl;
  cout << "Expected size: " << pow(2, n) << endl;
  cout << "Result: " << endl;
  for (int i = 0; i < result.size(); i++) {
    string s = result[i];

    if (i == result.size()-1) {
      cout << s << endl;
    }
    else {
      cout << s << ", ";
    }

    vector<string> v;
    string x;
    strncpy(x, s.c_str(), 2);
    v.push_back(x);
    strncpy(x, &s[2], 2);
    v.push_back(x);
    strncpy(x, &s[4], 2);
    v.push_back(x);
    splitResult.push_back(v);
  }

  cout << "Size of splitResult: " << splitResult.size() << endl;
  cout << "Expected size: " << 3 * pow(2, n) << endl;
  cout << "Split result: " << endl;
  for (auto v : splitResult) {
    for (auto s : v) {
      cout << s << ", ";
    }
    cout << endl;
  }
}

void printB()
{
  vector<vector<string> > allVecs;
  int n = 3; // number of vectors

  for (int i = 0; i < n; i++) {
    vector<string> vec = {"p"+to_string(i), "q"+to_string(i)};
    allVecs.push_back(vec);
  }
  vector<vector<string> > result;
  vector<string> strSoFar = {"", "", ""};
  printAll(allVecs, 0, strSoFar, result);
  cout << "Final check: " << endl;
  for (auto &s : strSoFar) {
    cout << s << ", ";
  }
  for (auto v : result) {
    cout << endl;
    for (auto s : v) {
      cout << s << ", ";
    }
  }
  cout << endl;
  cout << "Size of result: " << result.size() << endl;
  cout << "Expected size: " << pow(2, n) << endl;
}

// bool IsBkg(int particle)
// {
//   return particle == 1;
// }
// void rescalings()
// {
//   int n = 3; // number of particles
//   double pt = 5.;

//   double rescaleFactor = 0;
//   for (int i = 0; i < n; i++) {
//     if (IsBkg(particle[i])) {
//       double z = 0.;
//       r += z;
//       z = 2;
//     }

//     pt *= (1 - r);
//     for (int i = 0; i < n; i++) {
//       if (particle[i].z() < 1.1) {
//         particle[i].z() /= (1-r);
//       }
//     }

//   }
// }
// */

bool isPowerOfTwo(int x)
{
  return (x & (x - 1)) == 0;
}
vector<int> getState(uint32_t state, int nClasses = 4, int nParticles = 2)
{
  vector<int> particles(nParticles, nClasses);
  int nStates = pow(nClasses, nParticles);
  int nBitsPerParticle = round(log2(nClasses));
  int nBitsPerInt = sizeof(uint32_t) * 8;

  if (!isPowerOfTwo(nClasses)) {
    cout << "Number of classes (" << nClasses << ") must be a power of 2!" << endl;
    return particles;
  }
  if (nStates <= 0) {
    cout << "Illegal number of states (" << nStates << ")! " << ((nStates == 0) ? "" : "Max = 2^31") << endl;
    return particles;
  }
  if (nParticles * nBitsPerParticle > nBitsPerInt) {
    cout << "Number of bits required to parse the state (" << nParticles << "*" << nBitsPerParticle << " = " << nParticles * nBitsPerParticle << ") is too large for " << nBitsPerInt << " bits per int!" << endl;
    return particles;
  }
  if (state < 0 || state >= nStates) {
    cout << "Illegal state! State " << state << ((state < 0) ? " < 0" : (" >= " + to_string(nStates))) << endl;
    return particles;
  }

  for (int ip = 0; ip < nParticles; ip++) {
    int value = 0;
    int startBit = nBitsPerParticle * ip;
    for (int iBit = 0; iBit < nBitsPerParticle; iBit++) {
      int bit = startBit + iBit;
      int bitVal = ((state & (1 << bit)) > 0);
      value += bitVal * pow(2, iBit);
    }
    particles[ip] = value;
  }
  return particles;
}
vector<double> doCorrection(vector<int> particles, vector<vector<double>> weights, vector<double> values)
{
  // values = (z1, z2, ..., zn, ptjet)
  vector<double> output = values;
  double r = 0.;
  int nParticles = particles.size();
  if (values.size() != nParticles + 1) {
    cout << "Number of values (" << values.size() << ") must be equal to the number of particles (" << nParticles << ") + 1!" << endl;
    return values;
  }

  for (int ip = 0; ip < nParticles; ip++) {
    if (particles[ip] == 0) { // background
      r += values[ip];
    }
  }
  for (int ip = 0; ip < nParticles; ip++) {
    if (particles[ip] != 0) {
      output[ip] = values[ip] / (1 - r);
    }
  }
  output[nParticles] *= (1 - r);
  return output;
}
vector<vector<double> > getWeights(int nClasses, int nParticles)
{
  vector<vector<double> > weights;
  for (int ip = 0; ip < nParticles; ip++) {
    vector<double> w(nClasses, 1./(2*nClasses));
    w[0] = 1. - (nClasses-1) * w[0];
    weights.push_back(w);
  }
  return weights;
}
vector<double> getRandomValues(int nParticles, double ptjet)
{
  vector<double> values;
  double sum = 0.; double total = 0.9;
  for (int ip = 0; ip < nParticles; ip++) {
    double v = (double) rand()/RAND_MAX;
    sum += v;
    values.push_back(v);
  }
  if (sum > total) {
    for (int ip = 0; ip < nParticles; ip++) {
      values[ip] *= total/sum;
    }
  }
  values.push_back(ptjet);
  return values;
}
void processState(uint32_t state, int nClasses = 4, int nParticles = 3)
{
  vector<vector<double> > weights = getWeights(nClasses, nParticles);
  vector<int> particles = getState(state, nClasses, nParticles);
  vector<double> values = getRandomValues(nParticles, 25.);
  vector<double> corrected = doCorrection(particles, weights, values);

  cout << "weights: (";
  for (int iw = 0; iw < nClasses; iw++) {
    cout << weights[0][iw];
    if (iw < nClasses-1) {
      cout << ", ";
    }
  }
  cout << ")" << endl;

  cout << "State " << state << ": ";
  for (int ip = 0; ip < nParticles; ip++) {
    cout << "p" << ip << "[" << particles[ip] << "]";
    if (ip < nParticles-1) {
      cout << " * ";
    }
  }
  cout << " = ";
  double ws = 1.;
  for (int ip = 0; ip < nParticles; ip++) {
    int iw = particles[ip];
    if (iw >= nClasses) {
      // particles[i] >= nClasses means something went wrong, set to background class
      iw = 0;
    }
    cout << weights[ip][iw];
    if (ip < nParticles-1) {
      cout << " * ";
    }
    ws *= weights[ip][particles[ip]];
  }
  cout << " = " << ws << endl;

  cout << "values: (";
  for (int iv = 0; iv < values.size(); iv++) {
    cout << values[iv];
    if (iv < values.size()-1) {
      if (particles[iv] == 0) {
        cout << "\'";
      }
      cout << ", ";
    }
  }
  cout << ")" << endl << "-> (";
  for (int ic = 0; ic < corrected.size(); ic++) {
    cout << corrected[ic];
    if (ic < corrected.size()-1) {
      if (particles[ic] == 0) {
        cout << "\'";
      }
      cout << ", ";
    }
  }
  cout << ")" << endl;
}

void checkStateSum(int nClasses = 4, int nParticles = 3)
{
  vector<vector<double> > weights = getWeights(nClasses, nParticles);
  double totalWeight = 0.;
  for (int is = 0; is < pow(nClasses, nParticles); is++) {
    cout << ".";
    vector<int> particles = getState(is, nClasses, nParticles);
    double stateWeight = 1.;
    for (int ip = 0; ip < nParticles; ip++) {
      stateWeight *= weights[ip][particles[ip]];
    }
    totalWeight += stateWeight;
  }
  cout << endl << "Total weight: " << totalWeight << " (" << totalWeight - 1 << ")" << endl;
}

