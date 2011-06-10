#include "FitnessLandscape.h"
#include "Genotype.h"
#include "ParticleCommunity.h"
#include "LParameters.h"
#include "Node.h"

double FitnessLandscape::mortality(const Genotype &gen,
				   const ParticleCommunity &context) const
{
  double m;
  if (lparameters.densityDependentMortality)
    m = lparameters.delta_0
          + lparameters.delta_1 * context.totalPop * lparameters.cellSize;
  else
    m = lparameters.delta_0;
  return m;
}

double FitnessLandscape::fecundity(const Genotype &gen,
				   const ParticleCommunity &context) const
{
  double fec = fitness(gen, context) * lparameters.maxFecundity;
  if (lparameters.spaceLimitation)
    fec *= 1 - context.totalPop / (double)lparameters.totalSpace;
  return fec;
}

bool FitnessLandscape::isLocalPeak(const Genotype &x,
				   const ParticleCommunity &context) const
{
  double kx = fecundity(x,context);
  Genotype y;
  for (int j = 0; j < Genotype::nBlocks; j++)
    for (int k = 0; k < Genotype::blockSize; k++)
    {
      x.mutate(&y,j,k);
      if (fecundity(y,context) >= kx)  // what to do if one is == ?
	return false;
    }
  return true;
}

void FitnessLandscape::drawLandscapeGraph(ostream &os,
					  const ParticleCommunity& context)
  const
{
  // you will want to use something like 'dot -Tps -ols-1-6.ps' on this file
  Genotype g = *(Genotype::wildType());
  Genotype g1;
  os << "graph \"ls-" << Genotype::nBlocks
     << 'x' << Genotype::blockSize << "\" {\n"
    //     << "  page=\"8.5,11\";\n  ratio=auto;\n  ordering=out;\n";
     << "  page=\"8.5,11\";\n  rotate=90;\n  size=\"10,7.5\";\n"
     << "  ratio=fill;\n  ordering=out;\n";
  do
  { // first all vertices
    double f = fitness(g,context);
    os << "  \"" << g << "\" [label = \"" << g << ": " << f
	     << "\"];\n";
    ++g;
  }  while ( g != 0 );
  do
  { // second all edges
    for (int j = 0; j < Genotype::nBlocks; j++)
      for (int k = 0; k < Genotype::blockSize; k++)
      {
	g.mutate(&g1,j,k);
	if ( g < g1 )
	  os << "  \"" << g << "\" -- \"" << g1 << "\";\n";
      }
    ++g;
  }  while ( g != 0 );
  os << "}\n";
}

void FitnessLandscape::writeLabel(const Strain &s,
				  const ParticleCommunity& context,
				  ostream &os) const
{
  double f = fitness(s.genotype,context);
  os << "\" " << s.genotype << "\\n " << f << "\\n " << s.population << "\"";
}

void BlockFitnessLandscapeWithTreatment::writeLabel
	(const Strain &s, const ParticleCommunity& context, ostream &os)
  const
{
  double f = fitness(s.genotype,context);
  os << "\" " << s.genotype << "\\n " << f << "\\n ";
  if ( lparameters.applyingDrugs(context.site->integrator->currentTime()) )
    os << '(' << BlockFitnessLandscape::fitness(s.genotype,context) << ")\\n ";
  os << s.population << "\"";
}

void FitnessLandscape::drawQuasispeciesGraph(ostream &os,
					     const ParticleCommunity& context) 
  const
{
  double totalPop = context.totalPop; 
  // you will want to use something like 'dot -Tps -ols-1-6.ps' on this file
  os << "graph \""
     << Genotype::nBlocks << " block" << (Genotype::nBlocks==1?"":"s") << ", "
     << Genotype::blockSize << " bit" << (Genotype::blockSize==1?"":"s")
     << ", t=" << context.site->integrator->currentTime() << "\" {\n"
    //     << "  page=\"8.5,11\";\n  ratio=auto;\n  ordering=out;\n";
     << "  page=\"8.5,11\";\n  rotate=90;\n  size=\"10,7.5\";\n"
     << "  ratio=fill;\n  ordering=out;\n";
  for (Strain *walk = context.strains; walk != NULL;
       walk = (Strain*)walk->next)  
  { // first name all vertices
    if ( walk->population > 0 )
    {
      Genotype &g = walk->genotype;
      os << "  \"" << g << "\" [label = ";
      writeLabel(*walk,context,os);
      bool ilp = isLocalPeak(g,context);
      if ( ilp )
        os << ", style=bold, style=filled, fontcolor=white"
	  /*<< ", fontname=Helvetica"*/;
      if (ilp || (walk->population/totalPop > 0.001))
	os << ", color=\"0 1 "
	   << (1-log(walk->population/totalPop)/log(0.001)) << "\"";
      os << "];\n";
    }
  }
  Genotype g1;
  Strain s1(&g1);
  for (Strain *walk = context.strains; walk != NULL;
       walk = (Strain*)walk->next)  
    if (walk->population > 0)
    { // now draw all edges
      Genotype &g = walk->genotype;
      for (int j = 0; j < Genotype::nBlocks; j++)
	for (int k = 0; k < Genotype::blockSize; k++)
	{
	  g.mutate(&g1,j,k);
	  Strain *fs1 = (Strain*)s1.findOnList(walk->next);
	  if ( (fs1 != NULL) && (fs1->population > 0) )
	  {
	    double bigger = fs1->population;
	    if (walk->population > bigger) bigger = walk->population;
	    os << "  \"" << g << "\" -- \"" << g1 << "\"";
	    if (bigger/totalPop > 0.001)
	      os << " [color=\"0 1 "
		 << (1-log(bigger/totalPop)/log(0.001)) << "\"]";
	    os << ";\n";
	  }
	}
    }
  os << "}\n";
}

Genotype *expFitnessLandscape::refGenotype = Genotype::wildType();

double expFitnessLandscape::fitness(const Genotype &gen,
				      const ParticleCommunity &context) const
{
  const double lambda = 1;
  int dist = refGenotype->hammingDistance(gen);
  return exp(- lambda * dist);
}

/* BlockFitnessLandscape produces uncorrelated fitness function for each
   block by using DES encryption
*/
BlockFitnessLandscape::BlockFitnessLandscape()
{
  des_cblock keyblock;
  des_string_to_key(lparameters.desKey,&keyblock);
  des_check_key = true;
  if ( des_set_key(&keyblock,keyschedule) )
    cerr << "DES Error setting key!!\n";
}

/* returns fitness from 0 to maxFecundity for each block
   thus these values should be averaged, not summed
*/
double BlockFitnessLandscape::blockFitness(int *block, int blockno) const
{
  des_cblock outblock;
  // ivec is a kind of seed for the encryption
  //  include blockno so each block's landscape is different
  des_cblock ivec = {0xfe,0xdc,0xba,0x98,0x76,0x54,0x32,0x10+blockno};
  // compute des checksum from the block value
  des_cbc_cksum((des_cblock*)block, &outblock, (Genotype::blockSize+7)/8,
		(des_ks_struct*)keyschedule, &ivec);
  // reduce checksum into a single int using xor
  // this code will be in trouble if cblock stops being the size of 2 ints
  unsigned int *cksump = (unsigned int*)&outblock;
  unsigned int cksum = cksump[0] ^ cksump[1];

  static double factor = 1 / (1 + (double)UINT_MAX);
  //  if (factor == 0) 
  //    factor = 1 / (1 + (double)UINT_MAX);
  return cksum * factor;
}

double BlockFitnessLandscape::fitness(const Genotype&x,
				      const ParticleCommunity& context) const
{
  double fitness = 0;
  for ( int i = 0; i < Genotype::nBlocks; i++ )
    fitness += blockFitness((int*)&x.genome[i], i);
  fitness /= Genotype::nBlocks;
  return fitness;
}

double BlockFitnessLandscapeWithTreatment::fitness(const Genotype&x,
						   const ParticleCommunity&
						     context) const
{
  static bool debugging = false;
  double fitness = BlockFitnessLandscape::fitness(x, context);
  if (debugging) cout << "raw fitness " << fitness << endl
		       << "distance = "
		       << lparameters.wildType->hammingDistance(x) << endl
		       << "effect " << treatmentEffect(x) << endl;
  if ( lparameters.applyingDrugs(context.site->integrator->currentTime()) )
    fitness *= treatmentEffect(x);
  if ( debugging) cout << "new fitness " << fitness << endl
		       << endl;
  return fitness;
}

double BlockFitnessLandscapeWithTreatment::treatmentEffect(const Genotype&x)
  const
{
  int dist = lparameters.wildType->hammingDistance(x);
  if (lparameters.treatmentEffectSD > 0)
    return 1 - (1 - lparameters.treatmentEffectAtWT) * 
                    exp( - (dist*dist) / (lparameters.treatmentEffectSD*
					  lparameters.treatmentEffectSD) );
  else if (lparameters.treatmentEffectSD == 0)
    return (dist == 0) ? lparameters.treatmentEffectAtWT : 1;
  else // ??
    return 0;
}

// built in:
// viable strains are
//   1-mutations in positions 0 .. 10
//   2-mutations in (0,11), (1,12)
// fitnesses are distorted by drugs
//   1-mutations 2 .. 10 multiplied by 0.95 .. 1.03
//   2-mutation (1,12) by 1.03
#if 1
//assumes genotype is defined as 1 block
double MurrayFitnessLandscape::fitness(const Genotype &x,
				       const ParticleCommunity& context) const
{
  static Genotype *wt = Genotype::wildType();
  int distance = wt->hammingDistance(x);
  switch(distance)
  {
  case 0: // wild type
    return 1;
  case 1: // 1-mutation
      // feasible only if all the set bits are in the first feasible_n bits
    for (int i = 1; i < Genotype::wordsPerBlock; i++)
      if (x.genome[0].words[i])
	return 0;
#if USE_EXTRA_BITS
    if (x.genome[0].extraBits)
      return 0;
#endif
    for (int i = 0; i <= 10; i++)
    {
      BLOCK_WORD_TYPE mask = ((BLOCK_WORD_TYPE)1 << i);
      if (x.genome[0].words[0] == mask)
      {
/*
	double t = context.site->integrator->currentTime();
       	if (lparameters.applyDrugs && lparameters.startDrugs <= t &&
	    t < lparameters.stopDrugs && i >= 2)
	  return 0.95 + (i-2) * (1.03 - 0.95) / (10 - 2);
	else
*/	  return 1;
      }
    }
    // if not one of the 'approved' 1-mutations
    return 0;
  case 2: // 2-mutation
    for (int i = 1; i < Genotype::wordsPerBlock; i++)
      if (x.genome[0].words[i])
	return 0;
#if USE_EXTRA_BITS
    if (x.genome[0].extraBits)
      return 0;
#endif
    switch(x.genome[0].words[0])
    {
    case 0x00000801: // (0,11)
      return 1;
    case 0x00001002: // (1,12)
      {
/*	double t = context.site->integrator->currentTime();
	if(lparameters.applyDrugs && lparameters.startDrugs <= t
	   && t < lparameters.stopDrugs)
	  return 1.03;
	else
*/	  return 1;
      }
    default: // wrong 2-mutation
      return 0;
    }
  default: // > 2 mutations
    return 0;
  }
}
#else 1
double MurrayFitnessLandscape::fitness(const Genotype &x,
				       const ParticleCommunity& context) const
{
  // for now, without drugs
  static Genotype &wt = Genotype::wildType();
  int distance = wt.hammingDistance(x);
  if (distance == 0)
    return 1;
  else if (distance <= lparameters.max_distance)
  {  // feasible only if all the set bits are in the first feasible_n bits
    int feas = lparameters.feasible_n[distance];
    int word = 0; // assume 1 block here

    // skip any words that are all feasible
    while (feas >= Genotype::bitsPerWord)
    {
      feas -= Genotype::bitsPerWord;
      word++;
    }
    // test word that is partially feasible
    if (feas > 0)
    {
      BLOCK_WORD_TYPE mask = ~((((BLOCK_WORD_TYPE)1) << feas) - 1);
#if USE_EXTRA_BITS
      if (word == Genotype::wordsPerBlock)
      { // it's in the extra bits
	if (x.genome[0].extraBits & mask)
	  return 0;
      }
      else
#endif
      { // it's in a normal word
	if (x.genome[0].words[word] & mask)
	  return 0; // if it has bits outside the feasible range
      }
      word++;
    }

    // test any remaining words in entirety
    while (word < Genotype::wordsPerBlock)
    {
      if (x.genome[0].words[word])
	return 0;
      word++;
    }    
    // test extra bits if haven't yet
#if USE_EXTRA_BITS
    if (word == Genotype::wordsPerBlock)
      if (x.genome[0].extraBits)
	return 0;
#endif
    // if no stray bits, no problem
    return 1;
  }
  else // if distance is too great
    return 0;
}
#endif 0

// assumes 1 block of 1 word
double SinglePeakFitnessLandscape::fitness(const Genotype &x,
					   const ParticleCommunity& context)
  const
{
  switch ( x.genome[0].words[0] )
  {
  case 0:
    //  case 1: // for double peak
    return 1;
  default:
    return 0.95;
  }
}
