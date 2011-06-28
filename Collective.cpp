/*
 * A group of individuals for consensus dynamics model
 *
 * 6/2011  Lee Worden
*/
#include "Collective.h"
#include "BitString.h"
#include "LParameters.h"
#include "Site.h"
#include "ParticleIntegrator.h"
#include "Communicator.h"
//#include "ValueWriter.h"
#include "rand.h"
#include "util.h"
#include <algorithm>
#include <numeric>

Collective::Collective(void) : n_failures(0)
{
}

Collective::~Collective()
{
}

void Collective::initialize(void)
{ _nSpecies = nX = lparameters.groupSize();
  finished = false;
  checkAllocation();
  fill(_alive.begin(), _alive.end(), true);
  currentProposal = BitString::wildType(
      lparameters.nBlocks(), lparameters.nBits() / lparameters.nBlocks());
}

void Collective::checkAllocation(void)
{ Community::checkAllocation();
  unsigned as = alive.size();
  if (as < individuals.size())
    individuals.resize(as);
  else while (individuals.size() < as)
    // do this strangely to make sure none of them are duplicates
    individuals.resize(individuals.size() + 1);
}

bool Collective::isVariableInUse(const Index &n)
{ return !finished; //true;
}

void Collective::calcNextState(double t, const VectorAccess<double>*x,
                               VectorAccess<double>*nx)
{ SiteOutputController *oc = site->outputcontroller;
  oc->log( "at time %g, seeking consensus on proposal %s\n", 
                               t, currentProposal.hexString() );
  string strat = lparameters.facilitationStrategy();
  //oc->log("facilitation strategy is '%s'\n", strat.c_str());
  for (int i = 0; i < nX; ++i)
    oc->log( "Individual %d values %s at %g\n", i, currentProposal.hexString(),
             individuals[i].evaluate(currentProposal) );
  if (strat == "each proposes one")
  { vector<BitString> amendments(nX);
    vector< set<unsigned> > blocks(nX);
    copyTo(x, nx);
    for (int i = 0; i < nX; ++i)
    { // ask each individual to make proposed improvements
      amendments[i] = individuals[i].makeProposal(currentProposal);
      if (amendments[i] != currentProposal)
      { oc->log( "Individual %d proposes %s\n", i, amendments[i].hexString() );
        oc->log( "Individual %d values %s at %g\n", i, amendments[i].hexString(),
               individuals[i].evaluate(amendments[i]) );
        // is the new proposal acceptable to everyone (or to anyone)?
        for (int j = 0; j < nX; ++j)
          if (j != i)
          { oc->log( "Individual %d values %s at %g\n", j, amendments[i].hexString(),
               individuals[j].evaluate(amendments[i]) );
            if ( ! individuals[j].acceptableReplacement(amendments[i], currentProposal) )
            { oc->log( "Individual %d blocks %s\n", j, amendments[i].hexString() );
              blocks[i].insert(j);
            }
          }
      }
      else
      oc->log( "Individual %d has no proposal\n", i );
    }
    // any proposals get no blocks?
    finished = true;
    for (int i = 0; i < nX; ++i)
      if (blocks[i].size() == 0 && amendments[i] != currentProposal)
      { oc->log( "Proposal %s has no objections\n", amendments[i].hexString() );
        // FIXME: we can't just take the first one each time
        currentProposal = amendments[i];
        finished = false;
        break;
      }
  }
  else if (strat == "one proposal at a time")
  { unsigned ind = rand_index(individuals.size());
    finished = false;
    BitString newProposal = individuals[ind].makeProposal(currentProposal);
    int mnpf = lparameters.maxNumberOfProposalFailures();
    if (mnpf == -1)
      mnpf = 2*currentProposal.nBits();
    if (newProposal == currentProposal)
    { oc->log("Individual %d has no proposal\n", ind);
      ++n_failures;
      if (n_failures > mnpf)
      { finished = true;
      }
    }
    else
    { oc->log("Individual %d proposes %s\n", ind, newProposal.hexString());
      oc->log("Individual %d values %s at %g\n", ind, newProposal.hexString(),
          individuals[ind].evaluate(newProposal));
      bool accept = true;
      for (unsigned judge = 0; judge < individuals.size(); ++judge)
        if (judge != ind)
        { oc->log("Individual %d values %s at %g\n", judge, 
              newProposal.hexString(), individuals[judge].evaluate(newProposal));
          if ( ! individuals[judge].acceptableReplacement(newProposal, currentProposal) )
          { oc->log("Individual %d blocks %s\n", judge, newProposal.hexString());
            accept = false;
            break;
          }
        }
      if (accept)
      { oc->log("Proposal %s is accepted\n", newProposal.hexString());
        currentProposal = newProposal;
      }
      else
      { oc->log("Proposal %s is rejected\n", newProposal.hexString());
        ++n_failures;
        if (n_failures > mnpf)
        { finished = true;
        }
      }
    }
  }
  else
  { oc->log("Unknown strategy!\n");
    finished = true;
  }
}

map<string, double> Collective::outcomeStats(void)
{ map<string, double> stats;
  stats["n.individuals"] = individuals.size();
  vector<double> values(individuals.size());
  int nsatisfied = 0;
  for (unsigned i = 0; i < individuals.size(); ++i)
  { values[i] = individuals[i].evaluate(currentProposal);
    if (individuals[i].acceptable(currentProposal))
      ++nsatisfied;
  }
  stats["min.value"] = *min_element(values.begin(), values.end());
  stats["mean.value"] = 
    accumulate(values.begin(), values.end(), 0.0) / individuals.size();
  stats["max.value"] = *max_element(values.begin(), values.end());
  stats["n.satisfied"] = nsatisfied;
  return stats;
}
