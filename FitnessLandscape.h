/* -*- C++ -*- */
#ifndef FITNESS_LANDSCAPE_H
#define FITNESS_LANDSCAPE_H

#include "BitString.h"
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>

class FitnessLandscape
{
public:
  // fitness is the main thing to be overloaded
  // fitnesses should always be between 0 and 1
  virtual double fitness(const BitString&x)
    const = 0;
  // sometimes useful
  virtual bool isLocalPeak(const BitString&x)
    const;
  // draw hypercube graph with fitnesses (be careful if graph is large)
  virtual void drawLandscapeGraph(std::ostream &os)
    const;
  // draw the part of the hypercube that is occupied
  //virtual void drawQuasispeciesGraph(std::ostream &os)
  //  const;
protected:
  virtual void writeLabel(const BitString &s, std::ostream &os) const;
};

// Perelson + Macken 1995, correlated fitness landscape implemented using blocks
class BlockFitnessLandscape : public FitnessLandscape
{
public:
  string seed;
  double waterline;
  BlockFitnessLandscape(string s, double water=0);
  double blockFitness(const BitString &x, int blockno) const;
  virtual double fitness(const BitString&x) const;
};

// an abstraction to create correlation among multiple landscapes by 
// combining a shared fitness landscape with a unique private one.
template <class ExtrinsicClass, class IntrinsicClass>
class CompositeFitnessLandscape : public IntrinsicClass
{
protected:
  ExtrinsicClass *extrinsicLandscape;
  double weighting;
public:
  //Constructor for when IntrinsicClass is BlockFitnessLandscape
  // I don't think it would be good to do this using a partial
  // template specialization because it seems really inconvenient.
  CompositeFitnessLandscape(string s, double water = 0) :
    IntrinsicClass(s, water) {}

  CompositeFitnessLandscape<ExtrinsicClass, IntrinsicClass> &
    setExtrinsicLandscape(ExtrinsicClass&els)
  { extrinsicLandscape = &els;
    return *this;
  }
  CompositeFitnessLandscape<ExtrinsicClass, IntrinsicClass> &
    setWeighting(double w)
  { weighting = w;
    return *this;
  }
  virtual double fitness(const BitString &x)
  { double fitness = 0;
    if (weighting < 1)
      fitness += (1-weighting) * IntrinsicClass::fitness(x);
    if (weighting > 0)
      fitness += weighting*extrinsicLandscape->fitness(x);
    return fitness;
  }
};

class expFitnessLandscape : public FitnessLandscape
{
private:
  static BitString &refBitString;
public:
  virtual double fitness(const BitString&x) const;
};

// MurrayFitnessLandscape assumes nBlocks is 1 and wildType is all 0's
class MurrayFitnessLandscape: public FitnessLandscape
{
public:
  virtual double fitness(const BitString&x) const;
};

class SinglePeakFitnessLandscape: public FitnessLandscape
{
public:
  virtual double fitness(const BitString&x) const;
};

#endif //FITNESS_LANDSCAPE_H
