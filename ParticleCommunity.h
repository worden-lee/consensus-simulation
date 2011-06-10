/* -*- C++ -*-
 *
 * immune system quasispecies model - individual-based version
 *
 * 10/2000  Lee Worden
 *
*/
#ifndef LOCAL_COMMUNITY_H
#define LOCAL_COMMUNITY_H

#include "Community.h"
#include "Genotype.h"
#include "util.h"

class Strain : public ListMember
{
public:
  Genotype &genotype;
  int population;
  int newpop;
  //double delta;
  double deriv, avg;
  double movingTime;
  double cacheFecundity;
  bool scratch;

  const static double K_FLAG = -HUGE;

  Strain(Genotype*x);
  virtual ~Strain() { /*delete &genotype;*/ }
  ListMember *addToList(ListMember*);  // in ParticleCommunity.cpp
  virtual bool operator==(const ListMember &other)
  { return genotype == ((Strain&)other).genotype; }
};

class ParticleCommunity : public Community
{
public:
  Strain *strains;
  unsigned totalPop; // for easy axxess

  ParticleCommunity();
  void initialize(void);
  void recalcTotalPop();
  Strain *dominantStrain(void);
  Strain *fittestStrain(void);
  bool isVariableInUse(const VariableIndex &);
  
  friend class FitnessLandscape;
protected:
  void printForMathematica(ostream &o);
};

#endif LOCAL_COMMUNITY_H

