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

bool Individual::acceptable(const BitString &proposal)
{ return evaluate(proposal) >= 0;
}

BitString Individual::makeProposal(const BitString &proposal)
{ string strat = lparameters.individualProposalStrategy();
  cout << "individual strategy is '" << strat << "'" << endl;
  if (strat == "best")
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
  else if (strat == "best neighbor")
  { BitString curr = proposal;
    double cfit = evaluate(proposal);
    for (unsigned i = 0; i < proposal.nBlocks; ++i)
      for (unsigned j = 0; j < proposal.blockSize; ++j)
      { BitString poss;
        proposal.mutate(&poss, i, j);
        double pfit = evaluate(poss);
        if (pfit > cfit)
        { cfit = pfit;
          curr = poss;
        }
      }
    return curr;
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
  cout << "block strategy is " << bs << endl;
  if (bs == "if worse")
    return evaluate(proposal) >= evaluate(over);
  else if (bs == "if worse and not acceptable")
    return acceptable(proposal) || evaluate(proposal) >= evaluate(over);
  else
  { cout << "unknown block strategy " << bs << endl;
    return false;
  }
}
