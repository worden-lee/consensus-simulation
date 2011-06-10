#include "LSimulation.h"
#include "LParameters.h"
#include "ParticleCommunity.h"
#include "ParticleIntegrator.h"
#include "LValueWriter.h"
#include "LCommunicator.h"

/* have to populate the global simulation object */
LSimulation simobj;
Simulation &simulation = simobj;

void LSimulation::createNode(void)
{
  node = new Node;
  node->simulation = this;
  node->community = 0;
  node->communicator = new Communicator;
  node->valuewriter = new NodeValueWriter();
  node->integrator = NULL;
  const int nsites_per = 
    parameters.totalGridSize / parameters.nProcessors;
  Site *sites = new Site[nsites_per];
  node->assignSites(sites, nsites_per); // which calls populateSite
}

void LSimulation::populateSite(Site *s)
{
  s->integrator = new ParticleIntegrator;
  s->community = new ParticleCommunity;
  s->communicator = new LCommunicator;
  s->valuewriter = new LValueWriter;
}

bool LSimulation::finished(void)
{
  return ((ParticleIntegrator*)node->sites[0].integrator)->timeToQuit;
}
