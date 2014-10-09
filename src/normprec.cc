#include "normprec.hh"
#include "utils.hh" 

NormPRec::NormPRec(Env &env, Ratings &ratings)
  : _env(env), _ratings(ratings),
    _n(env.n), _m(env.m), _k(env.k),
    _iter(0),
    _start_time(time(0)),
    _theta("theta", 0.3, 0.3, _n,_k,&_r),
    _beta("beta", 0.3, 0.3, _m,_k,&_r),
    _thetabias("thetabias", 0.3, 0.3, _n, 1, &_r),
    _betabias("betabias", 0.3, 0.3, _m, 1, &_r),
    //_thetarate("thetarate", 0.3, 0.3, _n, &_r),
    //_betarate("betarate", 0.3, 0.3, _m, &_r),
    //_save_ranking_file(false),
    _use_rate_as_score(true),
    _topN_by_user(100),
    _maxval(0), _minval(65536),
    _prev_h(.0), _nh(.0)
    {
  gsl_rng_env_setup();
  const gsl_rng_type *T = gsl_rng_default;
  _r = gsl_rng_alloc(T);
  if (_env.seed)
    gsl_rng_set(_r, _env.seed);
  Env::plog("infer n:", _n);

  _hf = fopen(Env::file_str("/heldout.txt").c_str(), "w");
  if (!_hf)  {
    printf("cannot open heldout file:%s\n",  strerror(errno));
    exit(-1);
  }
  _vf = fopen(Env::file_str("/validation.txt").c_str(), "w");
  if (!_vf)  {
    printf("cannot open heldout file:%s\n",  strerror(errno));
    exit(-1);
  }
  _tf = fopen(Env::file_str("/test.txt").c_str(), "w");
  if (!_tf)  {
    printf("cannot open heldout file:%s\n",  strerror(errno));
    exit(-1);
  }
  _af = fopen(Env::file_str("/logl.txt").c_str(), "w");
  if (!_af)  {
    printf("cannot open logl file:%s\n",  strerror(errno));
    exit(-1);
  }
  _pf = fopen(Env::file_str("/precision.txt").c_str(), "w");
  if (!_pf)  {
    printf("cannot open logl file:%s\n",  strerror(errno));
    exit(-1);
  }  
  _df = fopen(Env::file_str("/ndcg.txt").c_str(), "w");
  if (!_df)  {
    printf("cannot open logl file:%s\n",  strerror(errno));
    exit(-1);
  }  
  if (!_env.write_training)
    load_validation_and_test_sets();

  Env::plog("theta mean:", _theta.mprior());
  Env::plog("theta var:", _theta.vprior());
  Env::plog("beta mean:", _beta.mprior());
  Env::plog("beta var:", _beta.vprior());
}

NormPRec::~NormPRec()
{
  fclose(_hf);
  fclose(_vf);
  fclose(_af);
  fclose(_pf);
  fclose(_tf);
}

void
NormPRec::load_validation_and_test_sets()
{
  char buf[4096];
  sprintf(buf, "%s/validation.tsv", _env.datfname.c_str());
  FILE *validf = fopen(buf, "r");
  assert(validf);
  if (_env.dataset == Env::NYT)
    _ratings.read_nyt_train(validf, &_validation_map);
  else
    _ratings.read_generic(validf, &_validation_map);
  fclose(validf);

  for (CountMap::const_iterator i = _validation_map.begin();
       i != _validation_map.end(); ++i) {
    const Rating &r = i->first;
    _validation_users_of_movie[r.second]++;
  }

  sprintf(buf, "%s/test.tsv", _env.datfname.c_str());
  FILE *testf = fopen(buf, "r");
  assert(testf);
  if (_env.dataset == Env::NYT)
    _ratings.read_nyt_train(testf, &_test_map);
  else
    _ratings.read_generic(testf, &_test_map);
  fclose(testf);
  
  // XXX: keeps one heldout test item for each user
  // assumes leave-one-out
  for (CountMap::const_iterator i = _test_map.begin();
       i != _test_map.end(); ++i) {
    const Rating &r = i->first;
    _leave_one_out[r.first] = r.second;
    debug("adding %d -> %d to leave one out", r.first, r.second);
  }

  printf("+ loaded validation and test sets from %s\n", _env.datfname.c_str());
  fflush(stdout);
  Env::plog("test ratings", _test_map.size());
  Env::plog("validation ratings", _validation_map.size());
}

void
NormPRec::initialize()
{
  _beta.initialize();
  _theta.initialize();

  if (_env.bias) {
    // TODO: look over 
    //_thetabias.initialize2(_m);
    //_thetabias.compute_expectations();
    
    //_betabias.initialize2(_n);
    //_betabias.compute_expectations();
  }
}

// perform inference
void
NormPRec::vb_bias() 
{
  lerr("running vb_bias()");
  initialize();

  Array phi(_k+2);
  Matrix phi_m(_m,_k+2); 
  phi.zero(); 
  phi_m.zero(); 
  while (1) {
    for (uint32_t n = 0; n < _n; ++n) { // for every user 

      Array phi_n(_k+2);
      phi_n.zero(); 

      const vector<uint32_t> *movies = _ratings.get_movies(n);
      for (uint32_t j = 0; j < movies->size(); ++j) {
        uint32_t m = (*movies)[j];
        yval_t y = _ratings.r(n,m);

        // TODO: look over 
        //const double **tbias = _thetabias.expected_logv().const_data();
        //const double **bbias = _betabias.expected_logv().const_data();
        const double **tbias = _thetabias.expected_v().const_data();
        const double **bbias = _betabias.expected_v().const_data();

        get_phi(_theta, n, _beta, m, tbias[n][0], bbias[m][0], phi);

        if (y > 1)
          phi.scale(y);

        // TODO: look over 
        //_thetabias.update_mean_next3(n, 0, phi[_k]);
        //_betabias.update_mean_next3(m, 0, phi[_k+1]);

        // TODO: add current phi to phi_n 
        phi_n.add_to(phi); 
        phi_m.add_slice(m,phi);
      }

      _theta.update_mean_next(n, phi_n);
      _theta.update_var_next(n); 
    }


    _theta.swap(); 
    // _theta.compute_expectations(); 

    // I'M HERE!
    Array p(_k+2);
    for (uint32_t m = 0; m < _m; ++m) { // for every item 

      phi_m.slice(0, m, p); 
      _beta.update_mean_next(m, p); 
      _beta.update_var_next(m); 

    }
    _beta.swap(); 

    // TODO: add bias 


    printf("\r iteration %d", _iter);
    fflush(stdout);    
    if (_iter % _env.reportfreq == 0) {
      compute_likelihood(true);
      compute_likelihood(false);
      save_model();
      compute_precision(false);
    }

    if (_env.save_state_now) {
      lerr("Saving state at iteration %d duration %d secs", _iter, duration());
      do_on_stop();
    }

    _iter++;
  }
}

void
NormPRec::save_model()
{
  _beta.save_state(_ratings.seq2movie());
  _theta.save_state(_ratings.seq2user());

  _betabias.save_state(_ratings.seq2movie());
  _thetabias.save_state(_ratings.seq2user());
}

void
NormPRec::compute_precision(bool save_ranking_file)
{
  printf("unimplemented\n"); 
}

void
NormPRec::get_phi(NormBase<Matrix> &a, uint32_t ai, 
		 NormBase<Matrix> &b, uint32_t bi, 
		 double biasa, double biasb,
		 Array &phi)
{
  assert (phi.size() == a.k() + 2 &&
	  phi.size() == b.k() + 2);
  assert (ai < a.n() && bi < b.n());
  const double  **ea = a.expected_v().const_data();
  const double  **eb = b.expected_v().const_data();
  // phi.zero();
  for (uint32_t k = 0; k < _k; ++k)
    phi[k] = ea[ai][k] + eb[bi][k];
  phi[_k] = biasa;
  phi[_k+1] = biasb;
  phi.lognormalize();
}

void
NormPRec::compute_likelihood(bool validation)
{
  uint32_t k = 0, kzeros = 0, kones = 0;
  double s = .0, szeros = 0, sones = 0;
  
  CountMap *mp = NULL;
  FILE *ff = NULL;
  if (validation) {
    mp = &_validation_map;
    ff = _vf;
  } else {
    mp = &_test_map;
    ff = _tf;
  }

  for (CountMap::const_iterator i = mp->begin();
       i != mp->end(); ++i) {
    const Rating &e = i->first;
    uint32_t n = e.first;
    uint32_t m = e.second;

    yval_t r = i->second;
    double u = rating_likelihood(n,m,r);
    s += u;
    k += 1;
  }

  double a = .0;
  info("s = %.5f\n", s);
  fprintf(ff, "%d\t%d\t%.9f\t%d\n", _iter, duration(), s / k, k);
  printf("s/k %f\n", s/k); 
  fflush(ff);
  a = s / k;  

  if (!validation)
    return;
  
  bool stop = false;
  int why = -1;
  if (_iter > 30) {
    if (a > _prev_h && _prev_h != 0 && fabs((a - _prev_h) / _prev_h) < 0.000001) {
      stop = true;
      why = 0;
    } else if (a < _prev_h)
      _nh++;
    else if (a > _prev_h)
      _nh = 0;

    if (_nh > 2) { // be robust to small fluctuations in predictive likelihood
      why = 1;
      stop = true;
    }
  }
  _prev_h = a;
  FILE *f = fopen(Env::file_str("/max.txt").c_str(), "w");
  fprintf(f, "%d\t%d\t%.5f\t%d\n", 
	  _iter, duration(), a, why);
  fclose(f);
  if (stop) {
    do_on_stop();
    exit(0);
  }
}

double
NormPRec::rating_likelihood(uint32_t p, uint32_t q, yval_t y) const
{
  const double **etheta = _theta.expected_v().const_data();
  const double **ebeta = _beta.expected_v().const_data();
  
  double s = .0;
  for (uint32_t k = 0; k < _k; ++k)
    s += etheta[p][k] * ebeta[q][k];
  
  if (_env.bias) {
    const double **ethetabias = _thetabias.expected_v().const_data();
    const double **ebetabias = _betabias.expected_v().const_data();
    s += ethetabias[p][0] + ebetabias[q][0];
  } 
  
  if (s < 1e-30)
    s = 1e-30;
  
  if (_env.binary_data)
    return y == 0 ? -s : log(1 - exp(-s));    
  return y * log(s) - s - log_factorial(y);
}

void
NormPRec::do_on_stop()
{
  save_model();
  // TODO 
  //gen_ranking_for_users(false);
}

