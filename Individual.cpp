/*
 * Individual represents a person in the consensus dynamics.
 *
 * Each individual has a fitness function, can evaluate a proposal,
 * and can propose alternatives to any proposal.
 */
#include "Individual.h"

Individual::Individual()
  : fitnesslandscape(0.5)
{}

double Individual::evaluate(const BitString &proposal)
{ return fitnesslandscape.fitness(proposal);
}

BitString Individual::makeProposal(const BitString &proposal)
{ BitString curr = proposal, nxt = proposal;
  do
  { curr = nxt;
    double nfitness = evaluate(nxt);
    for (int i = 0; i < BitString::nBlocks; ++i)
      for (int j = 0; j < BitString::blockSize; ++j)
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
