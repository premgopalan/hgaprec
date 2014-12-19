#ifndef NORMPREC_HH
#define NORMPREC_HH

#include "env.hh"
#include "normbase.hh"
#include "ratings.hh"
#if 0
#include <time.h>
#endif

class NormPRec {
public:
  NormPRec(Env &env, Ratings &ratings);
  ~NormPRec();

  void vb();

private:
  void initialize();
  void approx_log_likelihood();
  void save_model(); 
  void compute_precision(bool); 

  double elbo();

  void get_phi(NormBase<Matrix> &a, uint32_t ai, 
            NormBase<Matrix> &b, uint32_t bi, 
            Array &phi);
  void get_phi(NormBase<Matrix> &a, uint32_t ai, 
	       NormBase<Matrix> &b, uint32_t bi, 
	       double biasa, double biasb,
	       Array &phi);

  void load_validation_and_test_sets();
  void compute_likelihood(bool validation);
  // double log_factorial(uint32_t n)  const;

  double rating_likelihood(uint32_t p, uint32_t q, yval_t y) const;
  // double rating_likelihood_hier(uint32_t p, uint32_t q, yval_t y) const;
  double prediction_score(uint32_t p, uint32_t q) const;
  double score(uint32_t p, uint32_t q) const;
  void gen_ranking_for_users(bool load_model_state);


  uint32_t duration() const;
  bool is_validation(const Rating &r) const;
  void do_on_stop(); 
    

  Env &_env;
  Ratings &_ratings;

  uint32_t _n;
  uint32_t _m;
  uint32_t _k;
  // uint32_t _t; // for time-series model 
  uint32_t _iter;
  
  NormMatrix _theta;
  NormMatrix _beta;

  NormMatrix _thetabias;
  NormMatrix _betabias;

  //GPArray _thetarate;
  //GPArray _betarate;
  
  CountMap _validation_map;
  CountMap _test_map;
  FreqMap _validation_users_of_movie;
  IDMap _leave_one_out;

  UserMap _sampled_users;
  UserMap _sampled_movies;

  uint32_t _start_time;
  gsl_rng *_r;

  FILE *_hf;
  FILE *_vf;
  FILE *_tf;
  FILE *_af;
  FILE *_pf;
  FILE *_df;

  bool _use_rate_as_score;
  uint32_t _topN_by_user;
  uint32_t _maxval, _minval;
  double _prev_h;
  uint32_t _nh;
};

inline uint32_t
NormPRec::duration() const
{
  time_t t = time(0);
  return t - _start_time;
}

inline bool
NormPRec::is_validation(const Rating &r) const
{
  assert (r.first  < _n && r.second < _m);
  CountMap::const_iterator itr = _validation_map.find(r);
  if (itr != _validation_map.end())
    return true;
  return false;
}


#endif 
