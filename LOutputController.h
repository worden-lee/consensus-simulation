/* -*- C++ -*- */
#ifndef LOUTPUTCONTROLLER_H
#define LOUTPUTCONTROLLER_H

#include "OutputController.h"

class LOutputController: public SiteOutputController
{
protected:
public:
  void recordCommunity(void);
  void finish(void);
};

#endif //LOUTPUTCONTROLLER_H
