#ifndef LDA_HH
#define LDA_HH

#include "lda.hh"
#include "ratings.hh"
#include "matrix.hh"
#include "gpbase.hh"

class LDA {
public:
  LDA(Env &env, Ratings &ratings);
  ~LDA();

  void vb();

private:
  void initialize();
  
  Env &_env;
  Ratings &_ratings;

  uint32_t _n;
  uint32_t _m;
  uint32_t _k;
};


#endif
