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

Collective::Collective(void)
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
  currentProposal = BitString::wildType(); // a.k.a. 0.
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

void Collective::calcNextState(double t,
                               const VectorAccess<double>*x, VectorAccess<double>*nx)
{ SiteOutputController *oc = site->outputcontroller;
  oc->log( "at time t, seeking consensus on proposal %s\n", 
                               currentProposal.hexString() );
  for (int i = 0; i < nX; ++i)
    oc->log( "Individual %d values %s at %g\n", i, currentProposal.hexString(),
             individuals[i].evaluate(currentProposal) );
  vector<BitString> amendments(nX);
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
          if ( ! individuals[j].isAnImprovement(amendments[i], currentProposal) )
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
