#ifndef BITSTRING_H
#define BITSTRING_H

#include "LParameters.h"
#include <fstream>
#include <limits.h>

/*
 * BitString class embodies a list of equal-size blocks of bits
 *
 * The BlockFitnessLandscape class is well prepared to evaluate fitness
 *  (aka fecundity/infectivity) of these objects
 *
 * WARNING:  Make sure there aren't any stray bits in the padding
 *  section of the block structs, as they may get into the fitness
 *  function, which uses a SHA1 hash with a granularity of (I think)
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
#  error BitString.h has to be rewritten for different size int!
#endif

class BitString
{
public:
  unsigned int nBlocks, blockSize;
  BLOCK_WORD_TYPE **blocks;

  // functions
    // default constructor initializes to zero bits
  BitString();
    // nontrivial constructor initializes all bits to zero
  BitString(unsigned int nbl, unsigned int bls);
    // copy constructor copies
  BitString(const BitString&other);
    // arguments must be ((wordsPerBlock + (useExtraBits?1:0)) * nBlocks)
    //  number of BLOCK_WORD_TYPEs
  //BitString(BLOCK_WORD_TYPE b0, ...);
  ~BitString();

  BitString &operator=(const BitString &other);

  inline unsigned int wordsPerBlock(void) const
  { return blockSize / BITS_PER_WORD; }
  inline unsigned int nExtraBits(void) const
  { return blockSize - (wordsPerBlock() * BITS_PER_WORD); }
  inline bool useExtraBits(void) const
  { return nExtraBits() > 0; }
  inline unsigned int totalWordsPerBlock(void) const
  { return wordsPerBlock() + (useExtraBits() ? 1 : 0); }

  unsigned int hammingDistance(const BitString&) const;

  // create a mutation of *this at random
  BitString *mutate(void) const;
  // create a mutation at block 'bl', bit 'bi'
  BitString *mutate(unsigned bl, unsigned bi) const;
  // store a mutation of *this in *destination
  void mutate(BitString *destination) const;
  void mutate(BitString *destination, unsigned bl, unsigned bi) const;

  void randomize(void);

  static const BitString &wildType(unsigned int nbl, unsigned int bls);

  BitString &operator++(void);

  bool operator == (const BitString&other) const;
  bool operator != (const BitString&other) const { return !(*this==other); }
  bool operator < (const BitString&other) const; // for sorting
  bool operator > (const BitString&other) const { return (other < *this); }

  const char *hexString(void) const;
  friend ostream& operator<< (ostream &o, const BitString &comm);
 protected:
};

#endif //BITSTRING_H
