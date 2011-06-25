#include "BitString.h"
#include "rand.h"
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>

BitString::BitString() :
  nBlocks(0), blockSize(0), blocks(0)
{}
       
BitString::BitString(unsigned int nbl, unsigned int bls) :
  nBlocks(nbl), blockSize(bls)
{ blocks = new BLOCK_WORD_TYPE*[nBlocks];
  for (unsigned int bl = 0; bl < nBlocks; ++bl)
  { blocks[bl] = new BLOCK_WORD_TYPE[totalWordsPerBlock()];
    for (unsigned int w = 0; w < totalWordsPerBlock(); ++w)
      blocks[bl][w] = 0;
  }
}

BitString::BitString(const BitString&other) :
  nBlocks(other.nBlocks), blockSize(other.blockSize)
{ blocks = new BLOCK_WORD_TYPE*[nBlocks];
  for (unsigned int bl = 0; bl < nBlocks; ++bl)
  { blocks[bl] = new BLOCK_WORD_TYPE[totalWordsPerBlock()];
    for (unsigned int w = 0; w < totalWordsPerBlock(); ++w)
      blocks[bl][w] = other.blocks[bl][w];
  }
}

BitString &BitString::operator=(const BitString &other)
{ if (blocks != 0 && (other.nBlocks != nBlocks || other.blockSize != blockSize))
  { for (int i = 0; i < nBlocks; ++i)
      delete blocks[i];
    delete blocks;
  }
  this->nBlocks = other.nBlocks;
  this->blockSize = other.blockSize;
  blocks = new BLOCK_WORD_TYPE*[nBlocks];
  for (unsigned int bl = 0; bl < nBlocks; ++bl)
  { blocks[bl] = new BLOCK_WORD_TYPE[totalWordsPerBlock()];
    for (unsigned int w = 0; w < totalWordsPerBlock(); ++w)
      blocks[bl][w] = other.blocks[bl][w];
  }
  return *this;
}


BitString::~BitString()
{ for (unsigned int bl = 0; bl < nBlocks; ++bl)
    delete blocks[bl];
  delete blocks;
}

// wild type is 0, the default
// note this is carelessly written, if there are multiple values of
// nbl or bls hell will break loose
const BitString &BitString::wildType(unsigned int nbl, unsigned int bls)
{ static const BitString *wt = 0;
  if (wt == 0)
    wt = new BitString(nbl,bls);
  return *wt;
}

static unsigned int hammingDistance(unsigned int x, unsigned int y)
{ unsigned int distance = 0;
  while ( (x != 0) || (y != 0) )
  { if ( (x & 1) != (y & 1) )
      distance++;
    x >>= 1;
    y >>= 1;
  }
  return distance;
}

unsigned int BitString::hammingDistance(const BitString&other) const
{
  int distance = 0;
  for ( unsigned i = 0; i < nBlocks; ++i )
  { for ( unsigned j = 0; j < wordsPerBlock(); ++j )
      distance += ::hammingDistance(blocks[i][j], other.blocks[i][j]);
    if (useExtraBits())
      distance += ::hammingDistance( blocks[i][wordsPerBlock()],
             other.blocks[i][wordsPerBlock()] );
  }
  return distance;
}

BitString *BitString::mutate(void) const
{ BitString *mutant = new BitString();
  mutate(mutant);
  return mutant;
}

BitString *BitString::mutate(unsigned bl, unsigned bi) const
{ BitString *mutant = new BitString();
  mutate(mutant, bl, bi);
  return mutant;
}

void BitString::mutate(BitString *mutant) const
{ unsigned mbl = rand_index(nBlocks);
  unsigned mbi = rand_index(blockSize);
  mutate(mutant, mbl, mbi);
}

void BitString::mutate(BitString *mutant, unsigned mbl, unsigned mbi) const
{ *mutant = *this;
  unsigned mi = mbi / BITS_PER_WORD;
  mbi -= (mi * BITS_PER_WORD);
  mutant->blocks[mbl][mi]  ^= 1U << mbi;
}

void BitString::randomize(void)
{ // initialize random genotype
  for (unsigned j = 0; j < nBlocks; j++)
  { for ( unsigned k = 0; k < wordsPerBlock(); k++ )
      blocks[j][k] = genrand2i();
    if (useExtraBits())
    { // careful with this 1U if sizeof extraBits gets big
      const static BLOCK_WORD_TYPE mask = (1U << nExtraBits()) - 1;
      blocks[j][wordsPerBlock()] = genrand2i() & mask;
    }
  }
}

BitString &BitString::operator++(void)
{ for(unsigned i = 0; i < nBlocks; i++ )
  {  // keep incrementing things as long as they turn over
    if (useExtraBits())
    { if (++blocks[i][wordsPerBlock()] != (1U << nExtraBits()))
        return *this;
      else
        blocks[i][wordsPerBlock()] = 0;
    }
    for ( unsigned j = 0; j < wordsPerBlock(); j++ )
      if (++blocks[i][j] != 0) 
        return *this;
  }
  return *this;  // if we get here we've turned over to 0
}
      
bool BitString::operator==(const BitString& other) const
{ if (other.nBlocks != nBlocks || other.blockSize != blockSize)
    return false;
  for ( unsigned i = 0; i < nBlocks; i++ )
  { for ( unsigned j = 0; j < wordsPerBlock(); j++ )
      if (blocks[i][j] != other.blocks[i][j])
        return false;
    if (useExtraBits())
    { if (blocks[i][wordsPerBlock()] != other.blocks[i][wordsPerBlock()])
        return false;
    }
  }
  return true;
}

bool BitString::operator<(const BitString& other) const
{
  for ( unsigned i = 0; i < nBlocks; i++ )
  {
    for ( unsigned j = 0; j < wordsPerBlock(); j++ )
      if (blocks[i][j] < other.blocks[i][j])
        return true;
      else if (blocks[i][j] > other.blocks[i][j])
        return false;
    if (useExtraBits())
    { if (blocks[i][wordsPerBlock()] < other.blocks[i][wordsPerBlock()])
        return true;
      else if (blocks[i][wordsPerBlock()] > other.blocks[i][wordsPerBlock()])
        return false;
    }
  }
  return false;
}

// assumes 32 bit BLOCK_WORD_TYPE
const char *BitString::hexString(void) const
{
  const unsigned length = 1 + 
    nBlocks * ((wordsPerBlock()*(1+BITS_PER_WORD/4)) + (1+(nExtraBits()+3)/4));
  static char *name = 0;
  static unsigned name_len = 0;
  if (name_len < length)
  { if(name) delete name;
    name = new char[name_len = length];
  }
  name[0] = '\0';
  bool fi = true;
  for ( unsigned i = 0; i < nBlocks; i++ )
  { if ( fi ) fi = false;
    else strcat(name," ");
    
    bool fj = true;
    for ( unsigned j = 0; j < wordsPerBlock(); j++ )
    { static char b_form[32] = "";
      if ( !b_form[0] )
        snprintf(b_form, sizeof b_form, " %%.%dx", BITS_PER_WORD/4);
      sprintf(name+strlen(name),
	      fj?((fj=false),"%.8x"):" %.8x",blocks[i][j]);
    }
    if (useExtraBits())
    { static char xb_form[32] = "";
      if ( !xb_form[0] )
        snprintf(xb_form, sizeof xb_form, " %%.%dx", (nExtraBits()+3)/4);
      sprintf(name+strlen(name),
	      fj?xb_form+1:xb_form,(unsigned int)blocks[i][wordsPerBlock()]);
    }
  }
  return name;
}

ostream& operator<< (ostream &o, const BitString &gen)
{ return o << gen.hexString();
}

