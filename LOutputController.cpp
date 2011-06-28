#include "LOutputController.h"
#include "Site.h"
#include "Collective.h"
#include "LParameters.h"
#include "CSVDisplay.h"

void LOutputController::recordCommunity(void)
{
  bool yes = false;
  ofstream lsd;

  log( "Proposal on table is %s\n",
       ((Collective *)(site->community))->currentProposal.hexString() );
  if ( lparameters.nBits() <= 6 )
  {
    Collective *cv = (Collective *)site->community;
    for (int i = 0; (unsigned)i < cv->individuals.size(); ++i )
    { yes = true;
      string filename = lparameters.outputDirectory() + DIR_SEP
        + "landscape-" + int_to_string(i) + ".dot";
      lsd.open(filename.c_str());
      if ( !lsd.is_open() )
        perror("Failed to open landscape graph file");
      cv->individuals[i].fitnesslandscape.drawLandscapeGraph(lsd);
      lsd.close();
    }
  }
#if 0
  else if (((Collective*)site->community)->totalPop > 0)
  {
    yes = true;
    sprintf(buf, "%s%c%s%cquasispecies.dot", 
	    parameters.outputDirectory, DIR_SEP, dir, DIR_SEP);  
    lsd.open(buf);
    if ( !lsd.is_open() )
    {
      closeAnOpenFile();
      lsd.open(buf);
      if (!lsd.is_open())
	perror("Failed to open quasispecies graph file");
    }
    lparameters.fitnessLandscape->
      drawQuasispeciesGraph(lsd, *((Collective*)site->community));
  }
#endif
  if ( yes )
  {
    if ( lparameters.runDot() )
    {
      //system("make dot");
      system("for d in out/*.dot; "
	     "do echo dot -Tps -o ${d}.eps $d; "
	     "dot -Tps -o ${d}.eps $d; done");
    }
  }
}

void LOutputController::finish(void)
{ CSVDisplay csvstats(lparameters.outputDirectory() + "/outcome.csv");
  csvstats.writeLine(
      "n.individuals,min.value,mean.value,max.value,n.satisfied\n");
  map<string, double> stats = ((Collective*)site->community)->outcomeStats();
  csvstats << stats["n.individuals"] << stats["min.value"] 
    << stats["mean.value"] << stats["max.value"]
    << stats["n.satisfied"];
  csvstats.newRow();
  SiteOutputController::finish();
}
