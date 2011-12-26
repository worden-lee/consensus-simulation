#include "FitnessLandscape.h"

// iterator used in constructing MinimumLandscape - see below
template<typename IndividualIteratorT, typename LandscapeT>
class individuals_landscape_iterator : 
  public iterator<input_iterator_tag, LandscapeT*, void>
{
public:
  typedef LandscapeT *_Tp;
  IndividualIteratorT curr, last;
  individuals_landscape_iterator(IndividualIteratorT c, IndividualIteratorT l)
    : curr(c), last(l) {}
  LandscapeT* operator*()
  { if (curr == last) return (LandscapeT*)0;
    else return &(curr->fitnesslandscape);
  }
  individuals_landscape_iterator &operator++()
  { ++curr;
    return *this;
  }
  bool operator==(individuals_landscape_iterator<IndividualIteratorT, LandscapeT>&o)
  { return (curr == o.curr);
  }
  bool operator!=(individuals_landscape_iterator<IndividualIteratorT, LandscapeT>&o)
  { return (curr != o.curr);
  }
};

template<typename LandscapeIteratorT>
class MinimumLandscape : public FitnessLandscape
{
public:
  typedef typename LandscapeIteratorT::_Tp LandscapeT;
  typedef vector<LandscapeT> LandscapesT;
  LandscapesT landscapes;
  MinimumLandscape(LandscapeIteratorT start, LandscapeIteratorT finish)
    : landscapes(start, finish) {}
  double fitness(const BitString &x) const
  { double f = HUGE;
    cout << "MinimumLandscape evaluating " << x.hexString() << ": ";
    for (typename LandscapesT::const_iterator li = landscapes.begin(); 
         li != landscapes.end(); ++li)
    { cout << (*li)->fitness(x) << ' ';
      f = min(f, (*li)->fitness(x));
    }
    cout << ": " << f << endl;
    return f;
  }
};
