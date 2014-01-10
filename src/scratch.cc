
// void
// GAPRec::set_precision_user_sample()
// {
//   UserMap usermap;
//   do {
//     uint32_t user = gsl_rng_uniform_int(_r, _n);
//     const vector<uint32_t> *movies = _ratings.get_movies(user);

//     if (!movies || movies->size() < 10)
//       continue;
//     usermap[user] = true;
    
//     uint32_t msize = movies->size();
//     uArray liked_movies(msize);
//     for (uint32_t j = 0, k = 0; j < msize; ++j)
//       liked_movies[k++] = movies->at(j);
//     gsl_ran_shuffle(_r, (void *)liked_movies.data(), deg, sizeof(uint32_t));
    
//     for (uint32_t k = 0; k < msize * 0.3; k++) {


//     Rating r;
//     get_random_rating(r);
//     _test_ratings.push_back(r);
//     _test_map[r] = true;
//     n++;
//   }
//   FILE *hef = fopen(Env::file_str("/test-ratings.txt").c_str(), "w");  
//   fprintf(hef, "%s\n", ratingslist_s(_test_map).c_str());
//   fclose(hef);
    

//   } while (

// }


// void
// GAPRec::save_model()
// {
//   FILE *tf = fopen(Env::file_str("/theta.txt").c_str(), "w");  
//   double **gd = _Etheta.data();
//   double s = .0;
//   for (uint32_t i = 0; i < _n; ++i) {
//     const IDMap &m = _network.seq2id();
//     IDMap::const_iterator idt = m.find(i);
//     if (idt != m.end()) {
//       fprintf(tf,"%d\t", i);
//       debug("looking up i %d\n", i);
//       fprintf(tf,"%d\t", (*idt).second);
//       for (uint32_t k = 0; k < _k; ++k) {
// 	if (k == _k - 1)
// 	  fprintf(tf,"%.5f\n", gd[i][k]);
// 	else
// 	  fprintf(tf,"%.5f\t", gd[i][k]);
//       }
//     }
//   }
//   fclose(tf);
// }




class GPArrayGR : public GPBase<Array> {
public:
  GPArrayGR(string name, 
	  double a, double b,
	  uint32_t n, gsl_rng **r): 
    GPBase<Array>(name),
    _n(n),
    _sprior(a), // shape 
    _rprior(b), // rate
    _scurr(n), _snext(n),
    _rnext(n), _rcurr(n),
    _Ev(n), _Elogv(n),
    _r(r) { }
  ~GPArrayGR() {}

  uint32_t n() const { return _n;}
  uint32_t k() const { return 0;}

  const Array &shape_curr() const         { return _scurr; }
  const double rate_curr() const          { return _rcurr; }
  const Array &shape_next() const         { return _snext; }
  const double rate_next() const          { return _rnext; }
  const Array &expected_v() const         { return _Ev;    }
  const Array &expected_logv() const      { return _Elogv; }
  
  Array &shape_curr()       { return _scurr; }
  Array &shape_next()       { return _snext; }

  Array &expected_v()       { return _Ev;    }
  Array &expected_logv()    { return _Elogv; }

  const double sprior() const { return _sprior; }
  const double rprior() const { return _rprior; }
  
  uint32_t n() { return _n; }

  void set_to_prior();
  void set_to_prior_curr();

  void update_shape_next(const Array &phi);
  void update_shape_next(uint32_t n, double v);
  void update_rate_next(double v);

  void swap();
  void compute_expectations();
  void initialize();
  void initialize_exp();

  double compute_elbo_term_helper() const;
  void save_state(const IDMap &m) const;
  void load();

private:
  uint32_t _n;
  double _sprior;
  double _rprior;
  
  Array _scurr;
  Array _snext;
  
  double _rcurr;
  double _rnext;

  Array _Ev;   
  Array _Elogv;
  gsl_rng **_r;
};

inline void
GPArrayGR::set_to_prior()
{
  _snext.set_elements(_sprior);
  _rnext = _rprior;
}

inline void
GPArrayGR::update_shape_next(const Array &sphi)
{
  assert (sphi.size() == _n);
  _snext += sphi;
}

inline void
GPArrayGR::update_shape_next(uint32_t n, double v)
{
  _snext[n] += v;
}

inline void
GPArrayGR::update_rate_next(double v)
{
  _rnext += v;
}

inline void
GPArrayGR::swap()
{
  _scurr.swap(_snext);
  _rcurr = _rnext;
  set_to_prior();
}

inline void
GPArrayGR::set_to_prior_curr()
{
  _scurr.set_elements(_sprior);
  _rcurr = _rprior;
}

inline void
GPArrayGR::compute_expectations()
{
  const double * const ad = _scurr.const_data();
  double *vd1 = _Ev.data();
  double *vd2 = _Elogv.data();
  double a = .0, b = .0;
  for (uint32_t i = 0; i < _n; ++i) {
    make_nonzero(ad[i], _rcurr, a, b);
    vd1[i] = a / b;
    vd2[i] = gsl_sf_psi(a) - log(b);
  }
}

inline void
GPArrayGR::initialize()
{
  double *ad = _scurr.data();
  for (uint32_t i = 0; i < _n; ++i)
    ad[i] = _sprior + 0.01 * gsl_rng_uniform(*_r);
  _rcurr = _rprior;
  set_to_prior();
}

inline double
GPArrayGR::compute_elbo_term_helper() const
{
  const double *etheta = _Ev.const_data();
  const double *elogtheta = _Elogv.const_data();
  const double * const ad = shape_curr().const_data();

  double s = .0;
  double a = .0, b = .0;
  for (uint32_t n = 0; n < _n; ++n)  {
    make_nonzero(ad[n], _rcurr, a, b);
    s += _sprior * log(_rprior) + (_sprior - 1) * elogtheta[n];
    s -= _rprior * etheta[n] + gsl_sf_lngamma(_sprior);
    s -= a * log(b) + (a - 1) * elogtheta[n];
    s += b * etheta[n] + gsl_sf_lngamma(a);
  }
  return s;
}

inline void
GPArrayGR::save_state(const IDMap &m) const
{
  string expv_fname = string("/") + name() + ".tsv";
  string shape_fname = string("/") + name() + "_shape.tsv";
  string rate_fname = string("/") + name() + "_rate.tsv";
  _scurr.save(Env::file_str(shape_fname), m);
  _rcurr.save(Env::file_str(rate_fname), m);
  _Ev.save(Env::file_str(expv_fname), m);
}

inline void
GPArrayGR::load()
{
  string shape_fname = name() + "_shape.tsv";
  string rate_fname = name() + "_rate.tsv";
  _scurr.load(shape_fname);
  _rcurr.load(rate_fname);
  compute_expectations();
}

#endif
