#include "BitString.h"
#include "rand.h"
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>

BitString::BitString(void)
{
}

BitString::BitString(BLOCK_WORD_TYPE b0, ...)
{
  va_list va;
  va_start(va,b0);
  bool rb0 = true;
  for(int b = 0; b < nBlocks; b++)
  {
    for(int w = 0; w < wordsPerBlock; w++)
      genome[b].words[w] = (rb0 ? ((rb0=false),b0) :
			    va_arg(va,BLOCK_WORD_TYPE));
#if USE_EXTRA_BITS
    genome[b].extraBits = (rb0 ? ((rb0=false),b0) :
			   va_arg(va,BLOCK_WORD_TYPE));
#endif
  }
  va_end(va);
}

BitString &BitString::wildType(void)
{
  static BitString swt;
  static bool didInit = false;
  if ( !didInit )
  {
    for ( int i = 0; i < nBlocks; i++ )
    {
      for ( int j = 0; j < wordsPerBlock; j++ )
	swt.genome[i].words[j] = 0;
#if USE_EXTRA_BITS
      swt.genome[i].extraBits = 0;
#endif
    }
    didInit = true;
  }
  return swt;
}

unsigned int BitString::hammingDistance(const BitString&other) const
{
  int distance = 0;
  for ( int i = 0; i < nBlocks; i++ )
  {
    for ( int j = (USE_EXTRA_BITS?-1:0); j < wordsPerBlock; j++ )
    {
      unsigned int x, y;
#if USE_EXTRA_BITS
      if ( j < 0 )
      { x = genome[i].extraBits; y = other.genome[i].extraBits; }
      else
#endif
      {	x = genome[i].words[j]; y = other.genome[i].words[j]; }
      while ( (x != 0) || (y != 0) )
      {
	if ( (x & 1) != (y & 1) )
	  distance++;
	x >>= 1;
	y >>= 1;
      }
    }
  }
  return distance;
}

BitString *BitString::mutate(void) const
{
  BitString *mutant = new BitString();
  mutate(mutant);
  return mutant;
}

BitString *BitString::mutate(int bl, int bi) const
{
  BitString *mutant = new BitString();
  mutate(mutant, bl, bi);
  return mutant;
}

void BitString::mutate(BitString *mutant) const
{
  int mbl = rand_index(nBlocks);
  int mbi = rand_index(blockSize);
  mutate(mutant, mbl, mbi);
}

void BitString::mutate(BitString *mutant, int mbl, int mbi) const
{
  //mutant->genome = this->genome;
  std::copy(this->genome, this->genome + NBLOCKS, mutant->genome);
  int mi = mbi / BITS_PER_WORD;
  mbi -= (mi * BITS_PER_WORD);
#if USE_EXTRA_BITS
  if ( mi == wordsPerBlock )
    mutant->genome[mbl].extraBits ^= 1U << mbi;
  else
#endif
    mutant->genome[mbl].words[mi]  ^= 1U << mbi;
}

void BitString::randomize(void)
{
    // initialize random genotype
  for (int j = 0; j < BitString::nBlocks; j++)
  {
    for ( int k = 0; k < BitString::wordsPerBlock; k++ )
      genome[j].words[k] = genrand2i();
#if USE_EXTRA_BITS
    // careful with this 1U if sizeof extraBits gets big
    const static BLOCK_WORD_TYPE mask = (1U << BitString::nExtraBits) - 1;
    genome[j].extraBits = genrand2i() & mask;
#endif
  }
}

BitString &BitString::operator++(void)
{
  for(int i = 0; i < nBlocks; i++ )
  {
    for ( int j = (USE_EXTRA_BITS?-1:0); j < wordsPerBlock; j++ )
    {  // keep incrementing things as long as they turn over
#if USE_EXTRA_BITS
      if ( j < 0 )
      { if (++genome[i].extraBits != 0) return *this; }
      else
#endif
      {	if (++genome[i].words[j] != 0) return *this; }
    }
  }
  return *this;  // if we get here we've turned over to 0
}
      
bool BitString::operator==(const BitString& other) const
{
  for ( int i = 0; i < nBlocks; i++ )
  {
#if USE_EXTRA_BITS
    if (genome[i].extraBits != other.genome[i].extraBits)
      return false;
#endif
    for ( int j = 0; j < wordsPerBlock; j++ )
      if (genome[i].words[j] != other.genome[i].words[j])
	return false;
  }
  return true;
}

bool BitString::operator<(const BitString& other) const
{
  for ( int i = 0; i < nBlocks; i++ )
  {
    for ( int j = 0; j < wordsPerBlock; j++ )
      if (genome[i].words[j] < other.genome[i].words[j])
	return true;
      else if (genome[i].words[j] > other.genome[i].words[j])
	return false;
#if USE_EXTRA_BITS
    if (genome[i].extraBits < other.genome[i].extraBits)
      return true;
    else if (genome[i].extraBits > other.genome[i].extraBits)
      return false;
#endif
  }
  return false;
}

// assumes 32 bit BLOCK_WORD_TYPE
const char *BitString::hexString(void) const
{
  const int length = BitString::nBlocks * (
#if USE_EXTRA_BITS
					  (1+(N_EXTRA_BITS+3)/4) + 
#endif
					  9*BitString::wordsPerBlock);
  static char name[length];
  name[0] = '\0';
  bool fi = true;
  for ( int i = 0; i < BitString::nBlocks; i++ )
  {
    if ( fi ) fi = false;
    else strcat(name," ");
    
    bool fj = true;
    for ( int j = 0; j < BitString::wordsPerBlock; j++ )
      sprintf(name+strlen(name),
	      fj?((fj=false),"%.8x"):" %.8x",genome[i].words[j]);
#if USE_EXTRA_BITS
    {
      static char xb_form[32] = "";
      if ( !xb_form[0] )
	snprintf(xb_form, sizeof xb_form, " %%.%dx", (N_EXTRA_BITS+3)/4);
      sprintf(name+strlen(name),
	      fj?xb_form+1:xb_form,(unsigned int)genome[i].extraBits);
    }
#endif
  }
  return name;
}

ostream& operator<< (ostream &o, const BitString &gen)
{
  return o << gen.hexString();
}
