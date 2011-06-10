#ifndef GENOTYPE_H
#define GENOTYPE_H

#include "LParameters.h"
#include <limits.h>
#include <fstream.h>

/*
 * Genotype class embodies a list of equal-size blocks of bits
 *
 * The BlockFitnessLandscape class is well prepared to evaluate fitness
 *  (aka fecundity/infectivity) of these objects
 *
 * WARNING:  Make sure there aren't any stray bits in the padding
 *  section of the block structs, as they may get into the fitness
 *  function, which uses DES encryption with a granularity of (I think)
 *  bytes
 */

// -----------------------------------------------------------
// NBLOCKS and BLOCKSIZE should be the only things to change,
// to alter the landscape structure.
//  multiple blocks == correlation.
//  BLOCKSIZE is number of bits per block (not bytes). Can be large.

// -----------------------------------------------------------

#if UINT_MAX == 4294967295U
#  define BITS_PER_WORD 32
#  define BLOCK_WORD_TYPE unsigned int
#else
#  error Genotype.h has to be rewritten for different size int!
#endif

#define WORDS_PER_BLOCK (BLOCKSIZE / BITS_PER_WORD)
#define N_EXTRA_BITS (BLOCKSIZE - WORDS_PER_BLOCK * BITS_PER_WORD)

#if N_EXTRA_BITS > 0
#  define USE_EXTRA_BITS 1
#else
#  define USE_EXTRA_BITS 0
#endif

class Genotype
{
public:
  // class stuff
  const static int nBlocks = NBLOCKS, blockSize = BLOCKSIZE;
  const static int wordsPerBlock = WORDS_PER_BLOCK;
  const static int bitsPerWord = BITS_PER_WORD;
  const static int nExtraBits = N_EXTRA_BITS;
  typedef struct block_str
  { BLOCK_WORD_TYPE words[WORDS_PER_BLOCK]; 
#if USE_EXTRA_BITS
    BLOCK_WORD_TYPE extraBits:N_EXTRA_BITS;
#endif
  } block;

  // instance stuff
  block genome[NBLOCKS];

  // functions
  Genotype(void);
    // arguments must be ((wordsPerBlock + (useExtraBits?1:0)) * nBlocks)
    //  number of BLOCK_WORD_TYPEs
  Genotype(BLOCK_WORD_TYPE b0, ...);

  unsigned int hammingDistance(const Genotype&) const;

  // create a mutation of *this at random
  Genotype *mutate(void) const;
  // create a mutation at block 'bl', bit 'bi'
  Genotype *mutate(int bl, int bi) const;
  // store a mutation of *this in *destination
  void mutate(Genotype *destination) const;
  void mutate(Genotype *destination, int bl, int bi) const;

  void randomize(void);

  static Genotype *wildType(void);

  Genotype &operator++(void);

  bool operator == (const Genotype&other) const;
  bool operator != (const Genotype&other) const { return !(*this==other); }
  bool operator < (const Genotype&other) const; // for sorting

  const char *hexString(void) const;
  friend ostream& operator<< (ostream &o, const Genotype &comm);
 protected:
};

#endif GENOTYPE_H
