/* -*- C++ -*- */
#ifndef LOCAL_PARAMETERS_H
#define LOCAL_PARAMETERS_H

#include "Parameters.h"
#include <math.h>
#include <string.h>
#include <unistd.h>

#ifdef MPI
#define ROW_LENGTH 4
#define SITES_PER_NODE 4
#else //MPI
#define ROW_LENGTH 1
#define SITES_PER_NODE 1
#endif //MPI

class LParameters;
// global parameters object, with local type
extern LParameters &lparameters;

// model-specific parameter values

// values for BitString.h
// being sensitive to the fact these might be defined from outside.
#ifndef NBLOCKS
#define NBLOCKS    1
#endif
#ifndef BLOCKSIZE
#define BLOCKSIZE  8
#endif

class LParameters : public Parameters
{
public:
  // number of Individuals
  DECLARE_PARAM(int, groupSize)

  // params for BlockFitnessLandscape
  DECLARE_PARAM(string, hashKey)
  DECLARE_PARAM(int, nBits)
  DECLARE_PARAM(int, nBlocks)

  DECLARE_PARAM(string, individualProposalStrategy)
  DECLARE_PARAM(string, facilitationStrategy)
  DECLARE_PARAM(string, blockStrategy)

  DECLARE_PARAM(int, maxNumberOfProposalFailures)

  DECLARE_PARAM(bool, runDot)
  
  LParameters() { }
};

#endif //LOCAL_PARAMETERS_H
