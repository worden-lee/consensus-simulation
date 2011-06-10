#include "LValueWriter.h"
#include "Node.h"

// // just record community file, not individual things
void LValueWriter::recordValuesForSure(double t)
{
  //  recordCommunity();
  for (Strain *str = ((ParticleCommunity*)(site->community))->strains;
       str != NULL; str = (Strain*)(str->next))
  {
    recordValue((StrainVariableIndex)str,t,str->population);
    if ( parameters.useAuxDataFiles )
    {
      streamholder *as = findOrMakeStream((StrainVariableIndex)str);
      as->aux << t << '\t' << (str->avg==0 ? 0 : str->deriv/str->avg) << endl;
      as->aux.flush();
    }
  }
  if ( parameters.outputValues )
  {
    totalbiomass << t << " "
		 << ((ParticleCommunity*)(site->community))->totalPop << "\n";
    //speciescount << t << " " << site->community->speciesCount() << "\n";
  }
}

void LValueWriter::recordValues(double t)
{
  possiblyReopenFiles(t);
  if ( t - lastDidRecordValues >= parameters.outputTimeStep )
  {
    recordValuesForSure(t);
    lastDidRecordValues = t;
  }
  if ( t - lastDidGnuplot >= 100 )
  {
    showTimeSeriesInGnuplot();
    lastDidGnuplot = t;
  }
}

void LValueWriter::recordCommunity(void)
{
  char buf[FILENAME_MAX];
  bool yes = false;
  ofstream lsd;

  if ( Genotype::blockSize * Genotype::nBlocks <= 6 )
  {
    yes = true;
    sprintf(buf, "%s%c%s%clandscape.dot", 
	    parameters.outputDirectory, DIR_SEP, dir, DIR_SEP);  
    lsd.open(buf);
    if ( !lsd.is_open() )
    {
      closeAnOpenFile();
      lsd.open(buf);
      if (!lsd.is_open())
	perror("Failed to open landscape graph file");
    }
    lparameters.fitnessLandscape->
      drawLandscapeGraph(lsd, *((ParticleCommunity*)site->community));
  }
  else if (((ParticleCommunity*)site->community)->totalPop > 0)
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
      drawQuasispeciesGraph(lsd, *((ParticleCommunity*)site->community));
  }
  if ( yes )
  {
    lsd.close();
    if ( lparameters.runDot )
    {
      //system("make dot");
      system("for d in out/0/*.dot; "
	     "do echo dot -Tps -o out/${d#out/0/}.ps $d; "
	     "dot -Tps -o out/${d#out/0/}.ps $d; done");
    }
  }
}

void LValueWriter::showTimeSeriesInGnuplot(void)
{
  startGnuplot();
  flush();
  writeToGnuplot("set nokey\n");
  writeToGnuplot("set data style lines\n");
  bool started=false;
  for (Strain *strain = ((ParticleCommunity*)site->community)->strains;
       strain != NULL; strain = (Strain*)strain->next)
    if (strain->population > 10)
    {
      if(started)
	writeToGnuplot(", '");
      else
      {
	writeToGnuplot("plot '");
	started = true;
      }
      writeToGnuplot("out/0/X ");
      writeToGnuplot(strain->genotype.hexString());
      writeToGnuplot(".dat'");
    }
  if (started)
    writeToGnuplot("\n");
  fflush(gnuplot);
}

bool LValueWriter::doWrite(const VariableIndex &i)
{
  if(lparameters.fitnessLandscape->
     fitness(((StrainVariableIndex&)i).strain()->genotype,
	     *((ParticleCommunity*)(site->community)))
     > 0)
    return true;
  else
    return false;
}

    
