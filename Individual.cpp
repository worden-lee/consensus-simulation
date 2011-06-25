/*
 * Individual represents a person in the consensus dynamics.
 *
 * Each individual has a fitness function, can evaluate a proposal,
 * and can propose alternatives to any proposal.
 */
#include "Individual.h"
#include "rand.h"

Individual::Individual()
  : fitnesslandscape(lparameters.hashKey() + ' ' + individual_seed(), 0.5)
{ cout << "New individual, seed = " << fitnesslandscape.seed << endl;
}

const char *Individual::individual_seed(void)
{ //static const unsigned int nchars = sizeof(unsigned int) / sizeof(char);
  static char buf[3];
  buf[0] = (unsigned char)(' ' + rand_index(80));
  buf[1] = (unsigned char)(' ' + rand_index(80));
  buf[2] = '\0';
  return &buf[0];
}

double Individual::evaluate(const BitString &proposal)
{ return fitnesslandscape.fitness(proposal);
}

BitString Individual::makeProposal(const BitString &proposal)
{ BitString curr = proposal, nxt = proposal;
  do
  { curr = nxt;
    double nfitness = evaluate(nxt);
    for (unsigned i = 0; i < proposal.nBlocks; ++i)
      for (unsigned j = 0; j < proposal.blockSize; ++j)
      { BitString poss;
        curr.mutate(&poss, i, j);
        double pfitness = evaluate(poss);
        if (pfitness >= nfitness)
        { nxt = poss;
          nfitness = pfitness;
        }
      }
  } while (nxt != curr);
  return curr;
}

// return yes if the new proposal is better than the old one, not whether
// it's acceptable in absolute terms
bool Individual::isAnImprovement( const BitString &proposal,
                                  const BitString &over )
{ return evaluate(proposal) >= evaluate(over);
}
