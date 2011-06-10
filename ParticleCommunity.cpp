/*
 * HIV quasispecies model - individual-based version
 *
 * assumptions here:
 *   ??
 *
 * 10/2000  Lee Worden
*/
#include "ParticleCommunity.h"
#include "Genotype.h"
#include "LParameters.h"
#include "Site.h"
#include "ParticleIntegrator.h"
#include "Communicator.h"
#include "ValueWriter.h"
#include "rand.h"
#include "util.h"

Strain::Strain(Genotype*g) : genotype(*g), population(0), newpop(0)
{
  cacheFecundity = K_FLAG;
  avg = deriv = population = 0;
  movingTime = 0;
}

ListMember *Strain::addToList(ListMember*l)
{
    // insert, keep list sorted
  if ( l == NULL || this->genotype < ((Strain*)l)->genotype )
  { // insert at beginning of list
    this->next = l;
    if ( this->next )
      this->next->prev = this;
    this->prev = NULL;
    l = this;
  }
  else
  {  // find item to put self after
    Strain *y;
    for (y = (Strain*)l; y->next != NULL; y = (Strain*)(y->next))
      if ( genotype < ((Strain*)(y->next))->genotype )
	break;
    this->next = y->next;
    if ( this->next )
      this->next->prev = this;
    this->prev = y;
    y->next = this;
  }
  return l;
}

ParticleCommunity::ParticleCommunity(void)
{
  strains = NULL;
}

void ParticleCommunity::initialize(void)
{
  nX = 1; // nX used?
  checkAllocation(nX);
  alive[0] = true;

#if 1
  //Genotype *it = new Genotype(0xc163fe40);
  //Genotype *it = new Genotype(0xc1677e40);
  Genotype *it = new Genotype;
  it->randomize();
  //Strain *initStrain = new Strain(lparameters.wildType);
  Strain *initStrain = new Strain(it);
  initStrain->population = lparameters.initialPopulation;
  strains = (Strain*)initStrain->addToList(strains);
  cout << "init strain " << initStrain->genotype << " ("
       << lparameters.fitnessLandscape->fitness(initStrain->genotype, *this)
       << ")\n";
#else
  Genotype *wt = Genotype::wildType();  // strain 0x0
  Strain *initStrain = new Strain(wt);
  initStrain->population = lparameters.initialPopulation;
  strains = (Strain*)initStrain->addToList(strains);

  initStrain = new Strain(wt->mutate(0,0));    // strain 0x1
  initStrain->population = lparameters.initialPopulation;
  strains = (Strain*)initStrain->addToList(strains);

//#else // equilibrium values for murray experiment
  struct { BLOCK_WORD g; double v; } inits[] =
    { { 0x00000000, 0.53 },
      { 0x00000001, 0.177 },      
      { 0x00000002, 0.177 },      
      { 0x00000004, 0.157 },
      { 0x00000008, 0.157 },
      { 0x00000010, 0.157 },
      { 0x00000020, 0.157 },
      { 0x00000040, 0.157 },
      { 0x00000080, 0.157 },
      { 0x00000100, 0.157 },
      { 0x00000200, 0.157 },
      { 0x00000400, 0.157 },
      { 0x00000801, 0.0537 },
      { 0x00001002, 0.0537 },
      { 0, -1 }
    };
  for ( int i = 0; inits[i].v != -1; i++ )
  {
    Strain *s = new Strain(new Genotype(inits[i].g));
    s->population = (int)(inits[i].v / lparameters.cellSize);
    strains = (Strain*)s->addToList(strains);
  }
#endif
  //  site->valuewriter->recordCommunity();
}

void ParticleCommunity::recalcTotalPop(void)
{
  totalPop = 0;
  for (Strain *ss = strains; ss != NULL; ss = (Strain*)(ss->next))
  {
    if ( ss->newpop )
    {
      ss->population += ss->newpop;
      ss->newpop = 0;
    }
    totalPop += ss->population;
  }
}

Strain *ParticleCommunity::dominantStrain(void)
{
  Strain *dom = strains;
  for ( Strain *ss = strains; ss != NULL; ss = (Strain*)(ss->next))
  {
    if ( dom->population < ss->population )
      dom = ss;
  }
  return dom;
}

Strain *ParticleCommunity::fittestStrain(void)
{
  Strain *fs = NULL;
  double ff;
  for ( Strain *ss = strains; ss != NULL; ss = (Strain*)(ss->next))
    if (ss->population > 0)
    {
      double sf = lparameters.fitnessLandscape->fitness(ss->genotype,*this);
      if ( !fs || ff < sf )
      {
	fs = ss;
	ff = sf;
      }
    }
  return fs;
}

bool ParticleCommunity::isVariableInUse(const VariableIndex &n)
{
  return true;
}

void ParticleCommunity::printForMathematica(ostream &o)
{
  const int dim = Genotype::blockSize * Genotype::nBlocks;
  int dist[dim] = {0};
  double cf[dim] = {0};
  const Strain *ref = NULL;
  for (const Strain *ss = strains; ss != NULL; ss = (Strain*)(ss->next))
    if (ref == NULL || ss->population > ref->population)
      ref = ss;
  for (const Strain *ss = strains; ss != NULL; ss = (Strain*)(ss->next))
    if (ss->population > 0)
    {
      int hd = ref->genotype.hammingDistance(ss->genotype);
      double fit = lparameters.fitnessLandscape->fitness(ss->genotype,*this);
      dist[hd] += ss->population;
      cf[hd] += ss->population * fit;
    }
#define o cout
  o << "cx =\n{ ";
  if ( site->integrator )
  {
    o << "t -> " << site->integrator->currentTime() << ",\n  ";
  }
  if (ref)
  {
    o << "mode -> " << ref->genotype << ",\n  ";
    o << "mode_fitness -> "
      << lparameters.fitnessLandscape->fitness(ref->genotype,*this) << ",\n  ";
  }
  else
  {
    o << "mode = {},\n  ";
  }
  int last;
  o << "dist_hist -> {";
  {
    for (last = dim-1; last >= 0; last--)
      if ( dist[last] > 0 )
	break;
    bool fi = true;
    for (int i = 0; i <= last; i++)
      o << (fi? ((fi=false),""):", ")
	<< dist[i];
  }
  o << "}";
  o << ",\n  mf -> {";
  {
    bool fi = true;
    for (int i = 0; i <= last; i++)
      o << (fi ? ((fi=false),""): ", ") << (dist[i] ? cf[i] / dist[i] : 0);
  }
  o << "}";
  addToPrintForMathematica(o);
  o << "\n}\n";
}

#if 0
void ParticleCommunity::printForMathematica(ostream &o)
{
  o << "cx =\n{ ";
  if ( site->integrator )
  {
    o << "t -> " << site->integrator->currentTime() << ",\n  ";
  }
  o << "species -> {";
  { bool fi = true;
    for ( Strain *ss = strains; ss; ss = (Strain*)(ss->next) )
      if ( ss->population )
	o << (fi? ((fi=false),""):", ")
	  << ss->genotype;
  }
  o << "}";
  if ( site->integrator )
  {
    o << ",\n  state -> {";
    {
      bool fi = true;
      for ( Strain *ss = strains; ss; ss = (Strain*)(ss->next) )
	if ( ss->population )
	{ o << (fi ? ((fi=false),""): ", ") << ss->population; }
    }
    o << "}";
  }
  addToPrintForMathematica(o);
  o << "\n}\n";
}
#endif 0
