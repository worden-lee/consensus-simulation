#include "ParticleIntegrator.h"
#include "Collective.h"
//#include "LValueWriter.h"
#include "LParameters.h"
#include "Node.h"
#include "rand.h"
#include "util.h"

void ParticleIntegrator::integrate(double t0, double t1)
{
  Collective &comm = *(Collective*)(site->community);
    
  for ( t = t0; t1 < 0 || t < t1; t += lparameters.delta_t )
  {
    //if (t > t0)
    //  ((LValueWriter*)(site->valuewriter))->recordValues(t);

    comm.recalcTotalPop();
    if ( comm.totalPop == 0 )
    {
      //site->node->valuewriter->log("No more population.\n");
      //      t = t1;
      //cout << "No more population.\n";
      //site->valuewriter->flush();
      timeToQuit = true;
      return;
    }
    
    {
      Strain *currentDom = comm.dominantStrain();
      static Strain *keepFittest = 0;
      Strain *fittest = comm.fittestStrain();
      if ( currentDom != dominant )
      {
	dominant = currentDom;
	dominantSince = t;
	//if (t > 0)
	  //site->node->valuewriter->log("(t = %g) new dominant strain %s"
		//		       " (%g)\n",
		//		       t, dominant->genotype.hexString(),
		//		       lparameters.fitnessLandscape
		//		       ->fitness(dominant->genotype,
		//				 comm));
	//if (dominant == keepFittest)
	//  site->node->valuewriter->log("(t = %g) possible equilibrium\n", t);
      }
      if ( fittest->population > lparameters.significantPopulationThreshold )
      {
	if (fittest != keepFittest )
	{
	  //site->node->valuewriter->log("(t = %g) new fittest %s (%g)\n",
		//		       t, fittest->genotype.hexString(),
		//		       lparameters.fitnessLandscape->
		//		       fitness(fittest->genotype,
		//			       comm));
	  keepFittest = fittest;

	  if ( fittest != dominant )
	  {
	    //if (t - equilibriumSince > 1)
	    //site->node->valuewriter->log("(t = %g) no equilibrium\n", t);
	  }
	  //else //if (equilibriumSince >= t - lparameters.delta_t)  
	    //site->node->valuewriter->log("(t = %g) possible equilibrium\n", t);
	}
	if ( fittest != dominant )
	  equilibriumSince = t;
      }
    }
    //if (t - equilibriumSince > lparameters.equilibriumTime)
      //site->node->valuewriter->log("(t = %g) found equilibrium\n", t);

    if ( lparameters.applyDrugs )
    {
      static bool started = false;
      static bool stopped = false;
      if ( !started )
      {
	if ((lparameters.applyDrugsAtEquilibrium &&
	     (t - equilibriumSince > lparameters.equilibriumTime)) ||
	    (!lparameters.applyDrugsAtEquilibrium &&
	     (t - dominantSince > lparameters.waitingTimeForTreatment)))
	{
	  //site->node->valuewriter->log("wild type for treatment = %s (%g)\n",
				       //dominant->genotype.hexString(),
				       //lparameters.fitnessLandscape->
				       //fitness(dominant->genotype,
					       //comm));
	  //site->node->valuewriter->log("(t = %g) start treatment ----------\n",
				       //t);
	  lparameters.wildType = &dominant->genotype;
	  //lparameters.treatmentStarted = t;
	  lparameters.treating = true;
	  started = true;
	  equilibriumSince = dominantSince = t - 1;
	  lparameters.treatmentStarted = t;
	}
	else if ((lparameters.applyDrugsAtEquilibrium &&
		  (t - dominantSince > lparameters.maxWaitForEquilibrium)))
	{
	  //site->node->valuewriter->log("(t = %g) waiting too long "
				       //"before treatment\n", t);
	  timeToQuit = true;
	  return;
	}
      }
      else if ( !stopped )
      {
	if (((lparameters.stopDrugsAtEquilibrium &&
		  (t - equilibriumSince > lparameters.equilibriumTime)) ||
		 (!lparameters.stopDrugsAtEquilibrium &&
		  (t > lparameters.treatmentStarted +
		   lparameters.treatmentDuration))) )
	{
	  //site->node->valuewriter->log("new wild type = %s (%g)\n",
				       //dominant->genotype.hexString(),
				       //lparameters.fitnessLandscape->
				       //fitness(dominant->genotype, comm));
	  //site->node->valuewriter->log("(t = %g) stop treatment ----------\n",
				       //t);
	  stopped = true;
	  lparameters.treating = false;
	  equilibriumSince = dominantSince = t - 1;
	}
	else if ((lparameters.stopDrugsAtEquilibrium &&
		  (t - dominantSince > lparameters.maxWaitForEquilibrium)))
	{
	  //site->node->valuewriter->log("(t = %g) waiting too long "
				       //"in treatment\n", t);
	  timeToQuit = true;
	  return;
	}	  
      }
      else if ( stopped )
      {
	if (((lparameters.endRunAtEquilibrium
	      && (t - equilibriumSince > lparameters.equilibriumTime))
	     || (!lparameters.endRunAtEquilibrium
		 && (t - lparameters.treatmentStarted -
		     lparameters.treatmentDuration >
		     lparameters.waitingTimeForTreatment)) ))
	{
	  //site->node->valuewriter->log("final wild type = %s\n",
				       //dominant->genotype.hexString());
	  //site->node->valuewriter->log("hamming distance "
				       //"from wild type = %u\n",
				       //lparameters.wildType->
				       //hammingDistance(dominant->genotype));
	  //site->node->valuewriter->log("(t = %g) end of run.\n", t);
	  timeToQuit = true;
	  return;
	}
	else if ((lparameters.endRunAtEquilibrium &&
		  (t - dominantSince > lparameters.maxWaitForEquilibrium)))
	{
	  //site->node->valuewriter->log("(t = %g) waiting too long "
				       //"after treatment\n", t);
	  timeToQuit = true;
	  return;
	}
      }
    }  // if applyDrugs
    
    /*    for ( Strain *x = comm.strains; x != NULL; x = (Strain*)(x->next) )
    {
      //x->delta *= lparameters.equilibriumDecayFactor;
      x->deriv = 0;
    }
    */
    
      // fecundity first
    for ( Strain *x = comm.strains; x != NULL; x = (Strain*)(x->next) )
    {
      if ( x->population > 0 )
      {
	double k = lparameters.fitnessLandscape->fecundity(x->genotype, comm);
	// k is per capita offspring per unit time
	k *= lparameters.delta_t;
	
	// generate the mutants
	const int nPos = BitString::nBlocks * BitString::blockSize;
	// prob of mutation at given position
	double pgamma = lparameters.p_mutation / nPos;
	for ( int bl = 0; bl < BitString::nBlocks; bl++ )
	  for ( int bi = 0; bi < BitString::blockSize; bi++ )
	  {
	    //int nmut = (int)bnldev(pgamma, nOffspring);
				// # of this mutant generated
	    int nmut = (int)poidev(pgamma * k * x->population);
	    if ( nmut > 0 )
	    {
	      //nOffspring -= nmut;
	      BitString *mutant = x->genotype.mutate(bl,bi);
	      //Strain *ylast = NULL;
	      bool inserted = false;
	      for ( Strain *y = comm.strains; y != NULL;
		    y = (Strain*)(y->next) )
		if ( y->genotype == *mutant )
		{  // is already in the list
		  y->newpop += nmut;
		  //y->delta += nmut;
		  y->deriv += pgamma * k * x->population
		    * (1 - lparameters.equilibriumDecayFactor);
		  delete mutant;
		  //ylast = NULL;
		  inserted = true;
		  break;
		}
	      /*
		else if ( *mutant < y->genotype )
		{  // insert into list
		  Strain *ms = new Strain(mutant);
		  //ms->delta = nmut;
		  ms->newpop = nmut;
		  ms->deriv = pgamma * k * x->population
		    * (1 - lparameters.equilibriumDecayFactor);
		  ms->next = y;
		  ms->prev = y->prev;
		  if ( y == comm.strains )
		    comm.strains = ms;
		  else
		    y->prev->next = ms;
		  y->prev = ms;
		  ylast = NULL;
		  break;
		}
		else
		  ylast = y;
	      */
	      //if ( ylast != NULL )
	      if (!inserted)
	      {   // put mutant at the beginning
		Strain *ms = new Strain(mutant);
		//ms->delta = nmut;
		ms->newpop = nmut;
		ms->deriv = pgamma * k * x->population
		  * (1 - lparameters.equilibriumDecayFactor);
		ms->next = comm.strains;
		ms->next->prev = ms;
		ms->prev = NULL;
		comm.strains = ms;
		//ylast->next = ms;
	      }
	    }   // end of if (nmut > 0)
	  }   // end of for ( bi )
	// generate the non-mutants
	int nTrue =
	  (int)poidev((1 - lparameters.p_mutation) * k * x->population);
	x->newpop += nTrue;  // # of nonmutant offspring
	//x->delta += nTrue;
	x->deriv += (1 - lparameters.p_mutation) * k * x->population
	  * (1 - lparameters.equilibriumDecayFactor);
      }   // end of if population > 0
    }  // end of for (Strain *x)
      
    // includes making newPop old
    comm.recalcTotalPop();
  
      // mortality second
    Strain *nextx;
    Strain *lastLiving = NULL;
    for ( Strain *x = comm.strains; x != NULL; x = nextx )
    {
      nextx = (Strain*)(x->next);// in case x gets moved to the end
      if ( x->population > 0 )
      {
	lastLiving = x;
	double m = lparameters.fitnessLandscape->mortality(x->genotype, comm);
	// lifetimes are exponentially distr. (given survival up to now)
	//  with expected lifetime 1/m
	//  so prob. of death in next delta_t is:
	double prDeath = 1 - exp(-lparameters.delta_t * m);
	int nDead = (int)bnldev(prDeath,x->population);
	if ( nDead > x->population ) cerr << "Too many dead in integrate()!\n";
	x->deriv -= prDeath * x->population*(1 - lparameters.equilibriumDecayFactor);
	x->population -= nDead;
	//x->delta -= nDead;
	if ( x->population > 0 )
	{
	  //float dlt = x->delta * (1 - lparameters.equilibriumDecayFactor);
	  float dlt = x->deriv;
	  if ( dlt < 0 ) dlt = -dlt;
	  if ( x->population > lparameters.absoluteEquilibriumPopThreshold
	       && dlt > lparameters.absoluteEquilibriumDerivThreshold
	                 * lparameters.delta_t
	       && dlt / x->avg > lparameters.equilibriumThreshold )
	    x->movingTime++;
	  else
	    x->movingTime = 0;
	  if ( t - x->movingTime < equilibriumSince )
	    x->movingTime = 0;
	  //if ( x->movingTime > lparameters.minMovingTime )
	  //  equilibriumSince = t;
	}
#define eqx (1 - 0.99*(1-lparameters.equilibriumDecayFactor))
	x->avg *= eqx;
	x->avg += x->population*(1 - eqx);
	x->deriv *= lparameters.equilibriumDecayFactor;

	if ( parameters.doExtinction )
	{
	  if ( x->population <= 0 && x->newpop <= 0 )
	  {
	    comm.strains = (Strain *)x->removeFromList(comm.strains);
	    //site->valuewriter->extinction(t,(StrainVariableIndex)x);
	    delete x;
	  }
	}
      }  // end of if population > 0
				// lastLiving test keeps us from re-shifting
				// those already dead
      if ( x->population == 0 && x->newpop == 0 && x == lastLiving )
      { // move to the end of the list
	if (!lastStrain)
	  lastStrain = (Strain*)comm.strains->endOfList();
	if ( x->next ) // if not, no need to move
	{
	  x->next->prev = x->prev;
	  if ( x->prev )
	    x->prev->next = x->next;
	  else // for first one on the list
	    comm.strains = (Strain *)x->next;
	  lastStrain->next = x;
	  x->prev = lastStrain;
	  x->next = NULL;
	  lastStrain = x;
	}
      }
    }  // end of for ( x )
  }  // end of for (t)
  
  //((LValueWriter*)(site->valuewriter))->recordValues(t);
  //site->valuewriter->flush();
}

