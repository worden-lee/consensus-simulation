/* -*- C++ -*- */
#ifndef PARTICLE_INTEGRATOR_H
#define PARTICLE_INTEGRATOR_H
#include "Integrator.h"
class Strain;

class ParticleIntegrator : public Integrator
{
protected:
  Strain *dominant;
  Strain *lastStrain;
  double dominantSince;
  double equilibriumSince;
public:
  bool timeToQuit;
  void integrate(double t0, double t1);
  ParticleIntegrator(void):
    dominant(NULL), lastStrain(NULL),
    dominantSince(0), equilibriumSince(0),
    timeToQuit(false) {}
};

#endif PARTICLE_INTEGRATOR_H
