/* -*- C++ -*- */
#ifndef FITNESS_LANDSCAPE_H
#define FITNESS_LANDSCAPE_H

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include "BitString.h"
class Collective;

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
  virtual void drawQuasispeciesGraph(std::ostream &os)
    const;
protected:
  virtual void writeLabel(const BitString &s, std::ostream &os) const;
};

class expFitnessLandscape : public FitnessLandscape
{
private:
  static BitString &refBitString;
public:
  virtual double fitness(const BitString&x) const;
};

class BlockFitnessLandscape : public FitnessLandscape
{
public:
  BlockFitnessLandscape(double water=0);
  double waterline;
  double blockFitness(int *blockp, int blockno) const;
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
