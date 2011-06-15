#include "FitnessLandscape.h"
#include "BitString.h"
#include "Collective.h"
#include "LParameters.h"
#include "Node.h"

bool FitnessLandscape::isLocalPeak(const BitString &x) const
{
  double kx = fitness(x);
  BitString y;
  for (unsigned j = 0; j < BitString::nBlocks; j++)
    for (unsigned k = 0; k < BitString::blockSize; k++)
    {
      x.mutate(&y,j,k);
      if (fitness(y) >= kx)  // what to do if one is == ?
        return false;
    }
  return true;
}

void FitnessLandscape::drawLandscapeGraph(ostream &os)
  const
{
  // you will want to use something like 'dot -Tps -ols-1-6.ps' on this file
  BitString g = BitString::wildType();
  BitString g1;
  os << "graph \"ls-" << BitString::nBlocks
     << 'x' << BitString::blockSize << "\" {\n"
     << "  page=\"8.5,11\";\n  rotate=90;\n  size=\"10,7.5\";\n"
     << "  ratio=fill;\n  ordering=out;\n";
  do
  { // first all vertices
    double f = this->fitness(g);
    os << "  \"" << g << "\" [label = \"" << g << ": " << f
	     << "\"];\n";
    ++g;
  }  while ( g != 0 );
  do
  { // second all edges
    for (unsigned j = 0; j < BitString::nBlocks; j++)
      for (unsigned k = 0; k < BitString::blockSize; k++)
      {
	g.mutate(&g1,j,k);
	if ( g < g1 )
	  os << "  \"" << g << "\" -- \"" << g1 << "\";\n";
      }
    ++g;
  }  while ( g != 0 );
  os << "}\n";
}

void FitnessLandscape::drawQuasispeciesGraph(ostream &os) 
  const
{
  // you will want to use something like 'dot -Tps -ols-1-6.ps' on this file
  os << "graph \""
     << BitString::nBlocks << " block" << (BitString::nBlocks==1?"":"s") << ", "
     << BitString::blockSize << " bit" << (BitString::blockSize==1?"":"s")
     << "\" {\n"
    //     << "  page=\"8.5,11\";\n  ratio=auto;\n  ordering=out;\n";
     << "  page=\"8.5,11\";\n  rotate=90;\n  size=\"10,7.5\";\n"
     << "  ratio=fill;\n  ordering=out;\n";
#if 0
  for (set<Strain*>::iterator walk = context.strains.begin();
       walk != context.strains.end(); ++walk)
  { // first name all vertices
    BitString &g = (*walk)->genotype;
    os << "  \"" << g << "\" [label = ";
    writeLabel(**walk,context,os);
    bool ilp = isLocalPeak(g,context);
    if ( ilp )
      os << ", style=bold, style=filled, fontcolor=white";
    os << ", color=\"0 1 1\"";
    os << "];\n";
  }
  BitString g1;
  Strain s1(&g1);
  for (set<Strain*>::iterator walk = context.strains.begin();
       walk != context.strains.end(); ++walk)
  { // now draw all edges
    BitString &g = (*walk)->genotype;
    for (int j = 0; j < BitString::nBlocks; j++)
     for (int k = 0; k < BitString::blockSize; k++)
     {
       g.mutate(&g1,j,k);
       if (g1 > (*walk)->genotype)
       { // FIXME this probably doesn't work
         set<Strain*>::iterator fs1 = context.strains.find(&s1);
         if ( (fs1 != context.strains.end()) )
         {
           os << "  \"" << g << "\" -- \"" << g1 << "\"";
           os << " [color=\"0 1 1\"]";
           os << ";\n";
         }
       }
    }
  }
#endif
  os << "}\n";
}
 
void FitnessLandscape::writeLabel(const BitString &s,
				  ostream &os) const
{
  double f = fitness(s);
  os << "\" " << s << "\\n " << f << "\"";
}

BitString &expFitnessLandscape::refBitString = BitString::wildType();

double expFitnessLandscape::fitness(const BitString &gen) const
{
  const double lambda = 1;
  int dist = refBitString.hammingDistance(gen);
  return exp(- lambda * dist);
}

/* BlockFitnessLandscape produces uncorrelated fitness function for each
   block by using DES encryption
*/
BlockFitnessLandscape::BlockFitnessLandscape(string s, double water)
  : seed(s), waterline(water)
{}

/* returns fitness from 0 to maxFecundity for each block
   thus these values should be averaged, not summed
*/
#include <openssl/sha.h>
#include <numeric>
double BlockFitnessLandscape::blockFitness(int *block, int blockno) const
{
#if 0
  des_cblock outblock;
  // ivec is a kind of seed for the encryption
  //  include blockno so each block's landscape is different
  des_cblock ivec = {0xfe,0xdc,0xba,0x98,0x76,0x54,0x32,0x10+blockno};
  // compute des checksum from the block value
  des_cbc_cksum((des_cblock*)block, &outblock, (BitString::blockSize+7)/8,
		(des_ks_struct*)keyschedule, &ivec);
  // reduce checksum into a single int using xor
  // this code will be in trouble if cblock stops being the size of 2 ints
  unsigned int *cksump = (unsigned int*)&outblock;
  unsigned int cksum = cksump[0] ^ cksump[1];
#endif
  string hash_string =
    seed + char('1' + blockno)
      + string((char *)block, (BitString::blockSize+7)/8);
  unsigned char *hash = SHA1((const unsigned char*)hash_string.c_str(), 
      hash_string.length(), NULL);
  unsigned int *hash_ints_begin = (unsigned int *)hash;
  unsigned int n_hash_ints = 
    (SHA_DIGEST_LENGTH*sizeof(unsigned char))/sizeof(unsigned int);
  unsigned int cksum =
    std::accumulate(hash_ints_begin, hash_ints_begin + n_hash_ints, 0);

  static double factor = 1 / (1 + (double)UINT_MAX);
  return cksum * factor;
}

double BlockFitnessLandscape::fitness(const BitString&x) const
{
  double fitness = 0;
  for ( unsigned i = 0; i < BitString::nBlocks; i++ )
    fitness += blockFitness((int*)&x.genome[i], i);
  fitness /= BitString::nBlocks;
  return fitness - waterline;
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
double MurrayFitnessLandscape::fitness(const BitString &x) const
{
  static BitString &wt = BitString::wildType();
  int distance = wt.hammingDistance(x);
  switch(distance)
  {
  case 0: // wild type
    return 1;
  case 1: // 1-mutation
      // feasible only if all the set bits are in the first feasible_n bits
    for (int i = 1; i < BitString::wordsPerBlock; i++)
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
    for (int i = 1; i < BitString::wordsPerBlock; i++)
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
#else //1
double MurrayFitnessLandscape::fitness(const BitString &x,
				       const Collective& context) const
{
  // for now, without drugs
  static BitString &wt = BitString::wildType();
  int distance = wt.hammingDistance(x);
  if (distance == 0)
    return 1;
  else if (distance <= lparameters.max_distance())
  {  // feasible only if all the set bits are in the first feasible_n bits
    unsigned feas = lparameters.feasible_n()[distance];
    int word = 0; // assume 1 block here

    // skip any words that are all feasible
    while (feas >= BitString::bitsPerWord)
    {
      feas -= BitString::bitsPerWord;
      word++;
    }
    // test word that is partially feasible
    if (feas > 0)
    {
      BLOCK_WORD_TYPE mask = ~((((BLOCK_WORD_TYPE)1) << feas) - 1);
#if USE_EXTRA_BITS
      if (word == BitString::wordsPerBlock)
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
    while (word < BitString::wordsPerBlock)
    {
      if (x.genome[0].words[word])
	return 0;
      word++;
    }    
    // test extra bits if haven't yet
#if USE_EXTRA_BITS
    if (word == BitString::wordsPerBlock)
      if (x.genome[0].extraBits)
	return 0;
#endif
    // if no stray bits, no problem
    return 1;
  }
  else // if distance is too great
    return 0;
}
#endif //0

// assumes 1 block of 1 word
double SinglePeakFitnessLandscape::fitness(const BitString &x)
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
