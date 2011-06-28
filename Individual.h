/* An agent in the agent-based consensus model
   A Collective (subclass of Community) contains multiple Individuals.
   There is one Collective per site.
*/
#include "FitnessLandscape.h"

class Individual
{
public:
  BlockFitnessLandscape fitnesslandscape;
  Individual();
  string individual_seed(void);
  double evaluate(const BitString &proposal);
  BitString makeProposal(const BitString &proposal);
  bool acceptable(const BitString &proposal);
  bool acceptableReplacement(const BitString &proposal, const BitString &over);
};


