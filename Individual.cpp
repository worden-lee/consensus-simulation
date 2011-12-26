/*
 * Individual represents a person in the consensus dynamics.
 *
 * Each individual has a fitness function, can evaluate a proposal,
 * and can propose alternatives to any proposal.
 */
#include "Individual.h"
#include "rand.h"
#include <stdlib.h>
#include <sys/time.h>

BlockFitnessLandscape Individual::sharedlandscape("shared seed", 0.5);

Individual::Individual()
  : fitnesslandscape(individual_seed(), 0.5)
{ fitnesslandscape.setExtrinsicLandscape(Individual::sharedlandscape)
    .setWeighting(lparameters.weightingForSharedLandscape());
  cout << "New individual, seed = " << fitnesslandscape.seed << endl;
}

string Individual::individual_seed(void)
{ //static const unsigned int nchars = sizeof(unsigned int) / sizeof(char);
  static char buf[3];
  buf[0] = (unsigned char)(' ' + rand_index(80));
  buf[1] = (unsigned char)(' ' + rand_index(80));
  buf[2] = '\0';
  string s = lparameters.hashKey();
  if (s.length() == 0)
  { struct timeval tv;
    struct timezone tz;
    gettimeofday(&tv, &tz);
    s = int_to_string(tv.tv_sec) + " " + int_to_string(tv.tv_usec);
  }
  return s + " " + &buf[0];
}

double Individual::evaluate(const BitString &proposal)
{ return fitnesslandscape.fitness(proposal);
}

bool Individual::acceptable(const BitString &proposal)
{ return evaluate(proposal) >= 0;
}

BitString Individual::makeProposal(const BitString &proposal)
{ string strat = lparameters.individualProposalStrategy();
  //cout << "individual strategy is '" << strat << "'" << endl;
  if (strat == "best")
  { BitString curr = proposal, nxt = proposal;
    do
    { curr = nxt;
      nxt = fitnesslandscape.bestNeighbor(curr);
    } while (nxt != curr);
    return curr;
  }
  else if (strat == "best neighbor")
  { return fitnesslandscape.bestNeighbor(proposal);
  }
  else if (strat == "any improvement")
  { BitString poss;
    double cfit = evaluate(proposal);
    for (unsigned i = 0; i < proposal.nBits(); ++i)
    { unsigned ii = rand_index(proposal.nBits());
      proposal.mutate(&poss, ii);
      if (evaluate(poss) > cfit)
        return poss;
    }
    return proposal;
  }
  cout << "individual strategy '" << strat << "' not known" << endl;
  return proposal;
}

// return yes if willing to accept the new proposal in place of
// the old one, not whether it's acceptable in absolute terms
bool Individual::acceptableReplacement( const BitString &proposal,
                                        const BitString &over )
{ string bs = lparameters.blockStrategy();
  //cout << "block strategy is " << bs << endl;
  if (bs == "if worse")
    return evaluate(proposal) >= evaluate(over);
  else if (bs == "if worse and not acceptable")
    return acceptable(proposal) || evaluate(proposal) >= evaluate(over);
  else
  { cout << "unknown block strategy " << bs << endl;
    return false;
  }
}
