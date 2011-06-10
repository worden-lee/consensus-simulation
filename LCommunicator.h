/* -*- C++ -*- */
#ifndef LCOMMUNICATOR_H
#define LCOMMUNICATOR_H

#include "Communicator.h"
#include "ParticleCommunity.h"
#include <strstream.h>
class StrainUID : public UniqueID
{
public:
  inline StrainUID(Strain *s) { *(Strain**)&id = s; }
  virtual const char *filename(void) const
  { /*static*/ char n[FILENAME_MAX];
    /*static*/ ostrstream ns(n, FILENAME_MAX,ios::trunc);
    /*ns.seekp(0,ios::beg);*/
  ns << "X " << ((Strain*)id)->genotype << ".dat" << ends;
    return n;
  }
};

class LCommunicator: public Communicator
{
public:
  const UniqueID &uniqueIDOf(const VariableIndex i)
  { return StrainUID( ((StrainVariableIndex&)i).strain() ); }
};

#endif LCOMMUNICATOR_H
