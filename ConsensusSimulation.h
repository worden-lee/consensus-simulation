/* -*- C++ -*- */

#ifndef LOCAL_SIMULATION_H
#define LOCAL_SIMULATION_H

#include "Simulation.h"
#include "Site.h"
#include "Node.h"

class ConsensusSimulation : public Simulation
{
public:
  void createNode(void);
  void populateSite(Site *);
  //void doSimulation(double endtime);
  //bool finished();
};

#endif //LOCAL_SIMULATION_H
