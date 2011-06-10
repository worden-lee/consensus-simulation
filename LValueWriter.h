/* -*- C++ -*- */
#ifndef LVALUEWRITER_H
#define LVALUEWRITER_H

#include "ValueWriter.h"
#include "Integrator.h"

class Strain;

class LValueWriter: public ValueWriter
{
protected:
  Strain *dominant, *fittest;
  double lastDidGnuplot;
public:
  void recordValues(double t);
  void recordValuesForSure(double t);
  bool doWrite(const VariableIndex&);
  void recordCommunity(void);
  void showTimeSeriesInGnuplot(void);
  virtual char *filename(const UniqueID &uniqueID) // doesn't
    { static char n[20]; sprintf(n,"X%ld.dat",(long)uniqueID); return n; }
  LValueWriter(void): dominant(NULL), lastDidGnuplot(0)
  { startGnuplot(); writeToGnuplot("set log y\n"); }
};

#include "ParticleCommunity.h"
class StrainVariableIndex : public VariableIndex
{
public:
  inline StrainVariableIndex(Strain *s) { *((Strain **)&i) = s; }
  inline Strain* strain() const { return (Strain*)i; }
  inline StrainVariableIndex& operator++()
  { *((Strain**)&i) = (Strain*)(((Strain*)i)->next); return *this; }
};

#endif LVALUEWRITER_H
