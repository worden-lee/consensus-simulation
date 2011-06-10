/* -*- C++ -*- */
#ifndef FITNESS_LANDSCAPE_H
#define FITNESS_LANDSCAPE_H

#include <fstream.h>
#include <stdlib.h>

class Genotype;
class ParticleCommunity;
class Strain;

class FitnessLandscape
{
public:
  // fitness is the main thing to be overloaded
  // fitnesses should always be between 0 and 1
  virtual double fitness(const Genotype&x, const ParticleCommunity& context)
    const = 0;
  // fecundity is scaled fitness, possibly including space limitation
  //  need not be overloaded as long as fitness is
  virtual double fecundity(const Genotype&x, const ParticleCommunity& context)
    const;
  // mortality need not be overloaded
  //  it includes limitation due to immune response
  virtual double mortality(const Genotype&x, const ParticleCommunity& context)
    const;
  // sometimes useful
  virtual bool isLocalPeak(const Genotype&x, const ParticleCommunity& context)
    const;
  // draw hypercube graph with fitnesses (be careful if graph is large)
  virtual void drawLandscapeGraph(ostream &os,
				  const ParticleCommunity& context)
    const;
  // draw the part of the hypercube that is occupied
  virtual void drawQuasispeciesGraph(ostream &os,
				     const ParticleCommunity& context)
    const;
protected:
  virtual void writeLabel(const Strain &s, const ParticleCommunity& context,
			  ostream &os) const;
};

class expFitnessLandscape : public FitnessLandscape
{
private:
  static Genotype *refGenotype;
public:
  double fitness(const Genotype&x, const ParticleCommunity& context) const;
};

#include <des.h>

class BlockFitnessLandscape : public FitnessLandscape
{
public:
  BlockFitnessLandscape();
  double blockFitness(int *blockp, int blockno) const;
  double fitness(const Genotype&x, const ParticleCommunity& context) const;
private:
  des_key_schedule keyschedule; // internal form of encryption key
};

class BlockFitnessLandscapeWithTreatment : public BlockFitnessLandscape
{
public:
  double fitness(const Genotype&x, const ParticleCommunity& context) const;
protected:
  void writeLabel(const Strain&s, const ParticleCommunity& context,
		  ostream &os) const;
private:
  double treatmentEffect(const Genotype&x) const;
};

// MurrayFitnessLandscape assumes nBlocks is 1 and wildType is all 0's
class MurrayFitnessLandscape: public FitnessLandscape
{
public:
  double fitness(const Genotype&x, const ParticleCommunity& context) const;
};

class SinglePeakFitnessLandscape: public FitnessLandscape
{
public:
  double fitness(const Genotype&x, const ParticleCommunity& context) const;
};

#endif FITNESS_LANDSCAPE_H
