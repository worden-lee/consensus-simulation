#include "ConsensusSimulation.h"
#include "LParameters.h"
#include "Collective.h"
#include "Iterator.h"
#include "LOutputController.h"
//#include "LCommunicator.h"

/* have to populate the global simulation object */
ConsensusSimulation simobj;
Simulation &simulation = simobj;

void ConsensusSimulation::createNode(void)
{
  node = new Node;
  node->simulation = this;
  node->community = 0;
  node->communicator = new Communicator;
  node->outputcontroller = new NodeOutputController();
  node->integrator = NULL;
  const int nsites_per =
    lparameters.totalGridSize() / lparameters.nProcessors();
  //  Site *sites = new Site[nsites_per];
  vector<Site*> sites;
  for (int i = 0; i < nsites_per; ++i)
    sites.push_back(new Site);
  node->assignSites(sites); // which calls populateSite
}

void ConsensusSimulation::populateSite(Site *s)
{
  s->integrator = new Iterator;
  s->community = new Collective;
  s->communicator = new Communicator;
  s->outputcontroller = new LOutputController;
}

#if 0
void ConsensusSimulation::doSimulation(double endtime)
{ //SiteOutputController *soc
  //  = (SiteOutputController*)node->sites[0]->outputcontroller;
  //soc->recordCommunity();
  
}
#endif

#if 0
bool ConsensusSimulation::finished(void)
{
  return ((Iterator*)node->sites[0]->integrator)->timeToQuit;
}
#endif
