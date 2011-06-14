/* -*- C++ -*-
 *
 * Collective = community of individuals seeking consensus
 *
 * 6/2011  Lee Worden
 *
*/
#ifndef LOCAL_COMMUNITY_H
#define LOCAL_COMMUNITY_H

#include "Community.h"
#include "Individual.h"
#include "util.h"
#include <set>

class Collective : public Community
{
public:
  vector<Individual> individuals;
  BitString currentProposal;
  bool finished;

  Collective();
  ~Collective();
  void checkAllocation(void);
  void initialize(void);
  bool isVariableInUse(const Index &);
  void calcNextState(double, const VectorAccess<double> *, VectorAccess<double>*);
};

#endif //LOCAL_COMMUNITY_H

