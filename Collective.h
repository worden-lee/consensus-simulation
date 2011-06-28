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
#include <map>

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
  map<string, double> outcomeStats();
protected:
  int n_failures;
};

#endif //LOCAL_COMMUNITY_H

