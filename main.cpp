/*
 * individual based hiv quasispecies program
 * adapted from ODE parallel main program
 *
 * 4/28/2000  Lee Worden
 * 10/20/2000 Lee Worden
*/
#include "Collective.h"
//#include "ParticleIntegrator.h"
//#include "LValueWriter.h"
#include "LParameters.h"
#include "LCommunicator.h"
#include "LNode.h"
#include "rand.h"
#include <stdlib.h>
#include <iostream.h>
#include <signal.h>
#include <sys/resource.h>

#ifdef MPI
#include <mpi.h>
#endif

/* ---------- set up site class -------------------- */
//  this may be the only thing to change in this file
//  for a new project
void LSite::finishInitialize(void)
{
  integrator = new Iterator;
  community = new Collective;
  communicator = new LCommunicator;
  valuewriter = new LValueWriter;
  Site::finishInitialize();
}

 // accessed by signal handler
LNode *globalNode;

void handler_terminate(int sig_num); // handler for Ctrl-C - end of this file
void handler_dumpncont(int sig_num); // handler for Ctrl-\ - end of this file
void handler_abort(int sig_num);     // handler for seg violation - end

int main(int argc, char **argv)
{
  int i;

/*--------- first do mpi setup stuff -------------------*/
  double start_clock, end_clock;
  int mpi_rank, mpi_size;                 /* of MPI_COMM_WORLD */

#ifdef MPI
  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank );
  MPI_Comm_size( MPI_COMM_WORLD, &mpi_size );
  start_clock = MPI_Wtime();
#else
  mpi_rank = 0;
  mpi_size = parameters.rowLength;
  start_clock = time(0);
#endif

  {  // enable dumping core file in case of crash
    struct rlimit rl = { RLIM_INFINITY, RLIM_INFINITY };
    if (setrlimit(RLIMIT_CORE,&rl))
      perror("error in setrlimit");
  }

#ifdef __sun__
#define SIGNAL sigset
#else
#define SIGNAL signal
#endif
  SIGNAL(SIGQUIT, handler_dumpncont); // ^\ to dump out and continue
  SIGNAL(SIGINT,  handler_terminate); // ^C to dump out and quit
  SIGNAL(SIGPIPE, handler_terminate); // if pipe broken, dump + quit
  SIGNAL(SIGSEGV, handler_abort);

/*---------- now do the program -----------------------*/

  // seed the random numbers (before init'ing community, etc)
  //  seed should always be odd according to chave
  long rand_seed = parameters.randSeed + 2*mpi_rank;
  sgenrand2(rand_seed);

  //  cout.form("%lu blocks, %lu bits per block\n",
  //	    (unsigned long)Genotype::nBlocks,
  //	    (unsigned long)Genotype::blockSize);

  parameters.rankAndSizeAre(mpi_rank, mpi_size);
  
  globalNode = new LNode;
  globalNode->finishInitialize();

  if (lparameters.outputODEMessages)
  {
    globalNode->valuewriter->log("rand_seed %lu\n", rand_seed);
  }
  globalNode->valuewriter->log("%lu blocks, %lu bits per block\n",
			       (unsigned long)Genotype::nBlocks,
			       (unsigned long)Genotype::blockSize);
  globalNode->valuewriter->flush();

  double t = 0;
  double endtime = parameters.runLength;
  
  //    globalNode->sites[0].valuewriter->logCurrentState();
  //globalNode->sites[0].valuewriter->recordCommunity();
  ((LValueWriter*)(globalNode->sites[0].valuewriter))->recordValues(t);
  /*
  if ( lparameters.applyDrugs && !lparameters.applyDrugsAtEquilibrium )
  {
    globalNode->sites[0].integrator->integrate(t, lparameters.startDrugs);
    //    globalNode->sites[0].valuewriter->logCurrentState();
    t = lparameters.startDrugs;
    lparameters.treating = true;
    cout << "start treatment\n";
    globalNode->sites[0].integrator->integrate(t, lparameters.stopDrugs);
    //    globalNode->sites[0].valuewriter->logCurrentState();
    t = lparameters.stopDrugs;
    lparameters.treating = false;
    cout << "stop treatment\n";
  }
  */
  while ( (endtime < 0 || t < endtime)
	  && !(((ParticleIntegrator*)globalNode->sites[0].integrator)
	       ->timeToQuit) )
  {
    if ( lparameters.suspendAt > 0 && lparameters.suspendAt <= t )
    {
      handler_dumpncont(SIGSTOP);
      cout << "suspend" << endl;
      cout.flush();
      raise(SIGSTOP);
    }
    int pause;
    if ( t  > 5000 )
      pause = 1000;
    else if ( t > 2500 )
      pause = 500;
    else
      pause = 100;
    globalNode->sites[0].integrator->integrate(t, t+pause); 
    globalNode->sites[0].valuewriter->showTimeSeriesInGnuplot();
    double nt = globalNode->sites[0].integrator->currentTime();
    if ( nt <= t )
      break;
    t = nt;
  }
  //  globalNode->sites[0].valuewriter->logCurrentState();

  /*  
  bool alldead = false;
  // now integrate all sites for a long time
  while ( (t < endtime || endtime == -1) && !alldead )
  {
    double t1 = t + parameters.diffusionTimeStep;
    alldead = true;

     // integrate each one 'second' then do synchronous diffusion between sites
    for ( i = 0; i < globalNode->nSites; i++ )
    {
#ifdef MPI
      ((NodeCommunicator*)(globalNode->communicator))->handleRequests();
#endif MPI
      globalNode->sites[i].integrator->integrate(t, t1);
      if ( globalNode->sites[i].community->speciesCount() > 0 )
	alldead = false;
    }
    t = t1;
    
    if (lparameters.outputMPIMessages)
    {
      globalNode->valuewriter->log("t = %g\n", t);
    }
    
    for ( i = 0; i < globalNode->nSites; i++ )
    {
      globalNode->sites[i].community->possiblyReindex();
      globalNode->sites[i].communicator->doDiffusion();
    }
  } // end while (t < endtime...)
  */

#ifdef MPI
  end_clock = MPI_Wtime();
#else
  end_clock = time(0);
#endif

  for ( i = 0; i < globalNode->nSites; i++ )
    globalNode->sites[i].valuewriter->recordCommunity();

  if (lparameters.outputTiming)
    globalNode->valuewriter->log( "running time = %g\n", 
				  end_clock - start_clock );
    
  for ( i = 0; i < globalNode->nSites; i++ )
    globalNode->sites[i].valuewriter->finish();

#ifdef MPI
  MPI_Finalize();
#endif MPI
  return 0;
}

void handler_dumpncont(int sig_num)
{
  ((NodeValueWriter*)globalNode->valuewriter)->
    logWithAllCommunities("interrupted!\n");

  if ( !globalNode )
  {
    fputs("Error - globalNode pointer is null\n", stderr);
    fflush(stderr);
  }
  else
  {
    if ( globalNode->integrator )
    {
      globalNode->valuewriter->writeTag(cout);
      cout << globalNode->integrator->currentTime() << endl;
    }
    //globalNode->valuewriter->recordCommunity();
    globalNode->valuewriter->flush();
    for ( int i = 0; i < globalNode->nSites; i++ )
    {
      if ( globalNode->sites[i].integrator )
      {
	globalNode->sites[i].valuewriter->writeTag(cout);
	cout << globalNode->sites[i].integrator->currentTime() << endl;
      }
      globalNode->sites[i].valuewriter->recordCommunity();
      globalNode->sites[i].valuewriter->printPhylo();
      globalNode->sites[i].valuewriter->flush();
    }
  }

    //cout << "interrupted!\n";
  cout.flush();
}

void handler_terminate(int sig_num)
{
  handler_dumpncont(sig_num);
  if ( !globalNode )
  {
    fputs("Error - globalNode pointer is null\n", stderr);
    fflush(stderr);
  }
  else
    for ( int i = 0; i < globalNode->nSites; i++ )
    {
      globalNode->sites[i].valuewriter->finish();
    }
#ifdef MPI
  MPI_Finalize();
#endif MPI
  exit(-1);
}

void handler_abort(int sig_num)
{
  fputs("Received SIGSEGV, attempting to dump core...\n", stderr);
  fflush(stderr);
  abort();
}

