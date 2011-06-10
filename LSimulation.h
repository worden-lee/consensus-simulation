/* -*- C++ -*- */

#ifndef LOCAL_SIMULATION_H
#define LOCAL_SIMULATION_H

#include "Simulation.h"
#include "Site.h"
#include "Node.h"

class LSimulation : public Simulation
{
public:
  void createNode(void);
  void populateSite(Site *);
  bool finished();
};

#endif LOCAL_SIMULATION_H
