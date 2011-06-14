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
{ _nSpecies = nX = 1;
  finished = false;
  checkAllocation();
  alive[0] = true;
  currentProposal = BitString::wildType(); // a.k.a. 0.
}

void Collective::checkAllocation(void)
{ Community::checkAllocation();
  if (alive.size() != individuals.size())
    individuals.resize(alive.size());
}

bool Collective::isVariableInUse(const Index &n)
{ return !finished; //true;
}

void Collective::calcNextState(double t,
                               const VectorAccess<double>*x, VectorAccess<double>*nx)
{ SiteOutputController *oc = site->outputcontroller;
  oc->log( "at time t, seeking consensus on proposal %s\n", 
                               currentProposal.hexString() );
  vector<BitString> amendments(nX);
  vector< set<unsigned> > blocks(nX);
  copyTo(x, nx);
  for (int i = 0; i < nX; ++i)
  { // ask each individual to make proposed improvements
    amendments[i] = individuals[i].makeProposal(currentProposal);
    if (amendments[i] != currentProposal)
    { oc->log( "Individual %d proposes %s\n", i, amendments[i].hexString() );
      // is the new proposal acceptable to everyone (or to anyone)?
      for (int j = 0; j < nX; ++j)
        if (j != i && individuals[j].isAnImprovement(amendments[i], currentProposal) < 0)
        { oc->log( "Individual %d blocks %s\n", i, amendments[i].hexString() );
          blocks[i].insert(j);
        }
    }
    else
      oc->log( "Individual %d has no proposal\n", i );
  }
  // any proposals get no blocks?
  finished = true;
  for (int i = 0; i < nX; ++i)
    if (blocks[i].size() == 0 && amendments[i] != currentProposal)
    { oc->log( "proposal %s gets no objections\n", amendments[i].hexString() );
      currentProposal = amendments[i];
      finished = false;
    }
}
