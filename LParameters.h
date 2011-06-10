/* -*- C++ -*- */
#ifndef LOCAL_PARAMETERS_H
#define LOCAL_PARAMETERS_H

#include "Parameters.h"
#include "FitnessLandscape.h"
#include <math.h>
#include <string.h>
#include <unistd.h>

#ifdef MPI
#define ROW_LENGTH 4
#define SITES_PER_NODE 4
#else MPI
#define ROW_LENGTH 1
#define SITES_PER_NODE 1
#endif MPI

class LParameters;
// global parameters object, with local type
extern LParameters &lparameters;

// model-specific parameter values

// values for Genotype.h
// being sensitive to the fact these might be defined from outside.
#ifndef NBLOCKS
#define NBLOCKS    1
#endif
#ifndef BLOCKSIZE
#define BLOCKSIZE  16
#endif

class LParameters : public Parameters
{
public:
  double suspendAt;  // if >0, program will suspend itself at this time

  FitnessLandscape *fitnessLandscape;
  char *desKey; // used by BlockFitnessLandscape
  double maxFecundity;
  int initialPopulation;
  int n;
  //int n1;
  double k0, k1;
  double gamma;
  //static const int max_distance = 2;//6;
  //int feasible_n[1+max_distance];
  double p_mutation;
  bool densityDependentMortality, spaceLimitation;
  double delta_0, delta_1;
  double totalSpace;

  bool applyDrugs;
  //double startDrugs, stopDrugs;
  bool applyDrugsAtEquilibrium;
  bool stopDrugsAtEquilibrium;
  bool endRunAtEquilibrium;
  double waitingTimeForTreatment, treatmentDuration;
  double absoluteEquilibriumPopThreshold;
  double absoluteEquilibriumDerivThreshold;
  double equilibriumDecayFactor; // for iir smoothing
  double minMovingTime;  // how long out of eq. before it counts
  int significantPopulationThreshold;
  
  double treatmentStarted;
  bool treating;
  Genotype *wildType; // center of treatment curve
  double treatmentEffectAtWT;
  double treatmentEffectSD;
  
  double delta_t; // make this <1, makes time go slower (smaller timesteps)
  double cellSize;  // scales the immune response and space limitation
  bool runDot;  
  
  double applyingDrugs(double t)
    //{ return applyDrugs && (treatmentStarted <= t &&
    //		  t < treatmentStarted + treatmentDuration); }
  { return applyDrugs && treating; }
    
  LParameters()
    {
      randSeed = time(0) + getpid(); // program can get run twice in 1 second
      desKey = (char *)malloc(50);
      snprintf(desKey,50,"There is a rose in Spanish Harlem\n%u\n",
               (unsigned)randSeed);
      //snprintf(desKey,50,"There is a rose in Spanish Harlem 1\n");
      fitnessLandscape = new BlockFitnessLandscapeWithTreatment();
      n = NBLOCKS * BLOCKSIZE;  // # of positions
      gamma = 3.4e-5;      // that's per capita mutation rate per unit time
      k0 = k1 = log(2) / 1.6;

      //fitnessLandscape = new SinglePeakFitnessLandscape();
      //#define N_BITS 1680  // used by MurrayFitnessLandscape 
      //n1 = 10; // # of positions that are viable 1-mutations
/*    feasible_n[0] = 0;
      feasible_n[1] = 9;
      feasible_n[2] = 3;
*/
/*    feasible_n[1] = 1000;
      feasible_n[2] = 100;
      feasible_n[3] = 20;
      feasible_n[4] = 15;
      feasible_n[5] = 10; // these mean combinations in the first 
      feasible_n[6] = 8;  // n positions are feasible
*/
      p_mutation = 1 - pow(1-gamma,n); // prob that a given offspring is mutant

      densityDependentMortality = true;
      spaceLimitation = false;
      delta_1 = log(2) / 12.8;
      delta_0 = 5 * delta_1;
      maxFecundity = k1;
      totalSpace = 1000;
      applyDrugs = true;
      applyDrugsAtEquilibrium = true;
      stopDrugsAtEquilibrium = true;
      endRunAtEquilibrium = true;
      significantPopulationThreshold = 15;
      
      //      startDrugs = 10000;
      //      stopDrugs = 15000;//3000;
      //waitingTimeForTreatment = 5000;
      //treatmentDuration = 5000;
      absoluteEquilibriumPopThreshold = 100;
      absoluteEquilibriumDerivThreshold = 0;//100;
      equilibriumDecayFactor = 0.999;
      //treatmentStarted = HUGE; // will be changed by the integrator
      treating = false; // will be changed by the integrator
				// wildType will be set too
      //wildType = new Genotype(0xc163fe40); // center of treatment curve
      //wildType = new Genotype(0x2d);
      treatmentEffectAtWT = 0.95;
      treatmentEffectSD = 0;

      delta_t = 1;//0.1; // resolution of time stepper
      cellSize = 1e-7; //cellSize = 0.00001; // scaling for population values
      initialPopulation = (int)(1 / cellSize);
	// equilibrium population is ~3 without scaling (?)
      //initialPopulation = (int)(.3 / cellSize);
      runDot = false;//true;
              
      // override values from Parameters class

      runLength = 50000;//-1;//10000;
      suspendAt = -1;//4000;
      rowLength = ROW_LENGTH;
      sitesPerProcessor = SITES_PER_NODE;
      outputTimeStep = 50;//delta_t * 0.99;
      equilibriumThreshold = 0.0005;
      equilibriumTime = 5000;
      maxWaitForEquilibrium = 10000;
      minMovingTime = 10;
      //speciateAtEquilibrium = false;
      //speciationRate = 0;  // we use gamma for mutation rate
      doSpeciation = false;
      doExtinction = false; // keep empty strains around
      outputExtinctions = false;
      outputTiming = false;
      useAuxDataFiles = false;//true;
#ifdef RUN_GNUPLOT
      runGnuplot = RUN_GNUPLOT;
#else
      runGnuplot = true;
#endif      
      finishConstruct();
    }
};

#endif LOCAL_PARAMETERS_H
