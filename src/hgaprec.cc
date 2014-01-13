#include "hgaprec.hh"

HGAPRec::HGAPRec(Env &env, Ratings &ratings)
  : _env(env), _ratings(ratings),
    _n(env.n), _m(env.m), _k(env.k),
    _iter(0),
    _start_time(time(0)),
    _theta("theta", 0.3, 0.3, _n,_k,&_r),
    _beta("beta", 0.3, 0.3, _m,_k,&_r),
    _thetabias("thetabias", 0.3, 0.3, _n, 1, &_r),
    _betabias("betabias", 0.3, 0.3, _m, 1, &_r),
    _htheta("htheta", 0.3, 0.3, _n, _k, &_r),
    _hbeta("hbeta", 0.3, 0.3, _m, _k, &_r),
    _thetarate("thetarate", 0.3, 0.3, _n, &_r),
    _betarate("betarate", 0.3, 0.3, _m, &_r),
    _prev_h(.0), _nh(.0),
    _save_ranking_file(false),
    _use_rate_as_score(true)
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
  load_validation_and_test_sets();
}

HGAPRec::~HGAPRec()
{
  fclose(_hf);
  fclose(_vf);
  fclose(_af);
  fclose(_pf);
  fclose(_tf);
}

void
HGAPRec::load_validation_and_test_sets()
{
  char buf[4096];
  sprintf(buf, "%s/validation.tsv", _env.datfname.c_str());
  FILE *validf = fopen(buf, "r");
  assert(validf);
  _ratings.read_generic(validf, &_validation_map);
  fclose(validf);

  sprintf(buf, "%s/test.tsv", _env.datfname.c_str());
  FILE *testf = fopen(buf, "r");
  assert(testf);
  _ratings.read_generic(testf, &_test_map);
  fclose(testf);
  printf("+ loaded validation and test sets from %s\n", _env.datfname.c_str());
  fflush(stdout);
  Env::plog("test ratings", _test_map.size());
  Env::plog("validation ratings", _validation_map.size());
}

void
HGAPRec::initialize()
{
  _beta.initialize();
  _theta.initialize();

  _beta.initialize_exp();
  _theta.initialize_exp();

  if (_env.bias) {
    _thetabias.initialize();
    _thetabias.initialize_exp();
    
    _betabias.initialize();
    _betabias.initialize_exp();
  }

  if (_env.hier) {
    _thetarate.set_to_prior_curr();
    _thetarate.set_to_prior();
    
    _betarate.set_to_prior_curr();
    _betarate.set_to_prior();
    
    _hbeta.initialize();
    _hbeta.initialize_exp();
    //_hbeta.initialize_exp(_betarate.expected_v()[0]);
    
    _htheta.initialize();
    _htheta.initialize_exp();
    //_htheta.initialize_exp(_thetarate.expected_v()[0]);
  }
}


void
HGAPRec::get_phi(GPBase<Matrix> &a, uint32_t ai, 
		 GPBase<Matrix> &b, uint32_t bi, 
		 Array &phi)
{
  assert (phi.size() == a.k() &&
	  phi.size() == b.k());
  assert (ai < a.n() && bi < b.n());
  const double  **eloga = a.expected_logv().const_data();
  const double  **elogb = b.expected_logv().const_data();
  phi.zero();
  for (uint32_t k = 0; k < _k; ++k)
    phi[k] = eloga[ai][k] + elogb[bi][k];
  phi.lognormalize();
}

void
HGAPRec::get_phi(GPBase<Matrix> &a, uint32_t ai, 
		 GPBase<Matrix> &b, uint32_t bi, 
		 double biasa, double biasb,
		 Array &phi)
{
  assert (phi.size() == a.k() + 2 &&
	  phi.size() == b.k() + 2);
  assert (ai < a.n() && bi < b.n());
  const double  **eloga = a.expected_logv().const_data();
  const double  **elogb = b.expected_logv().const_data();
  phi.zero();
  for (uint32_t k = 0; k < _k; ++k)
    phi[k] = eloga[ai][k] + elogb[bi][k];
  phi[_k] = biasa;
  phi[_k+1] = biasb;
  phi.lognormalize();
}

void
HGAPRec::vb()
{
  lerr("running vb()");
  initialize();
  // approx_log_likelihood();

  Array phi(_k);
  while (1) {
    for (uint32_t n = 0; n < _n; ++n) {
      
      const vector<uint32_t> *movies = _ratings.get_movies(n);
      for (uint32_t j = 0; j < movies->size(); ++j) {
	uint32_t m = (*movies)[j];
	yval_t y = _ratings.r(n,m);
	
	get_phi(_theta, n, _beta, m, phi);
	if (y > 1)
	  phi.scale(y);
	
	_theta.update_shape_next(n, phi);
	_beta.update_shape_next(m, phi);
      }
    }

    Array betasum(_k);
    _beta.sum_rows(betasum);
    _theta.update_rate_next(betasum);

    _theta.swap();
    _theta.compute_expectations();

    Array thetasum(_k);
    _theta.sum_rows(thetasum);
    _beta.update_rate_next(thetasum);
    
    _beta.swap();
    _beta.compute_expectations();

    printf("\r iteration %d", _iter);
    fflush(stdout);
    if (_iter % _env.reportfreq == 0) {
      // approx_log_likelihood();
      compute_likelihood(true);
      compute_likelihood(false);
      save_model();
      compute_precision(false);
      if (_env.logl)
	logl();
    }

    if (_env.save_state_now) {
      lerr("Saving state at iteration %d duration %d secs", _iter, duration());
      do_on_stop();
      exit(0);
    }

    _iter++;
  }
}

void
HGAPRec::vb_bias()
{
  lerr("running vb_bias()");
  initialize();

  Array phi(_k+2);
  while (1) {
    for (uint32_t n = 0; n < _n; ++n) {
      
      const vector<uint32_t> *movies = _ratings.get_movies(n);
      for (uint32_t j = 0; j < movies->size(); ++j) {
	uint32_t m = (*movies)[j];
	yval_t y = _ratings.r(n,m);

	const double **tbias = _thetabias.expected_logv().const_data();
	const double **bbias = _betabias.expected_logv().const_data();

	get_phi(_theta, n, _beta, m, tbias[n][0], bbias[m][0], phi);
	
	if (y > 1)
	  phi.scale(y);

	_theta.update_shape_next(n, phi);
	_beta.update_shape_next(m, phi);
	
	_thetabias.update_shape_next(n, 0, phi[_k]);
	_betabias.update_shape_next(m, 0, phi[_k+1]);
      }
    }
    
    Array betasum(_k);
    _beta.sum_rows(betasum);
    _theta.update_rate_next(betasum);

    _theta.swap();
    _theta.compute_expectations();

    Array thetasum(_k);
    _theta.sum_rows(thetasum);
    _beta.update_rate_next(thetasum);
    
    _beta.swap();
    _beta.compute_expectations();

    _thetabias.update_rate_next(0, _m);
    _thetabias.swap();
    _thetabias.compute_expectations();

    _betabias.update_rate_next(0, _n);
    _betabias.swap();
    _betabias.compute_expectations();

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
      exit(0);
    }

    _iter++;
  }
}

void
HGAPRec::vb_hier()
{
  initialize();

  uint32_t x;
  if (_env.bias)
    x = _k+2;
  else
    x = _k;

  Array phi(x);
  while (1) {
    for (uint32_t n = 0; n < _n; ++n) {
      
      const vector<uint32_t> *movies = _ratings.get_movies(n);
      for (uint32_t j = 0; j < movies->size(); ++j) {
	uint32_t m = (*movies)[j];
	yval_t y = _ratings.r(n,m);

	if (_env.bias) {
	  const double **tbias = _thetabias.expected_logv().const_data();
	  const double **bbias = _betabias.expected_logv().const_data();
	  
	  get_phi(_htheta, n, _hbeta, m, tbias[n][0], bbias[m][0], phi);
	} else
	  get_phi(_htheta, n, _hbeta, m, phi);
	
	if (y > 1)
	  phi.scale(y);
	
	_htheta.update_shape_next(n, phi);
	_hbeta.update_shape_next(m, phi);

	if (_env.bias) {
	  _thetabias.update_shape_next(n, 0, phi[_k]);
	  _betabias.update_shape_next(m, 0, phi[_k+1]);
	}
      }
    }
    
    Array betarowsum(_k);
    _hbeta.sum_rows(betarowsum);
    _htheta.set_prior_rate(_thetarate.expected_v(), 
			   _thetarate.expected_logv());
    _htheta.update_rate_next(betarowsum);
    _htheta.swap();
    _htheta.compute_expectations();

    Array thetarowsum(_k);
    _htheta.sum_rows(thetarowsum);
    _hbeta.set_prior_rate(_betarate.expected_v(), 
			  _betarate.expected_logv());
    _hbeta.update_rate_next(thetarowsum);
    _hbeta.swap();
    _hbeta.compute_expectations();

    if (_env.bias) {
      _thetabias.update_rate_next(0, _m);
      _thetabias.swap();
      _thetabias.compute_expectations();
      
      _betabias.update_rate_next(0, _n);
      _betabias.swap();
      _betabias.compute_expectations();
    }

    Array thetacolsum(_n);
    _htheta.sum_cols(thetacolsum);
    _thetarate.update_rate_next(thetacolsum);

    _thetarate.swap();
    _thetarate.compute_expectations();

    Array betacolsum(_m);
    _hbeta.sum_cols(betacolsum);
    _betarate.update_rate_next(betacolsum);

    _betarate.swap();
    _betarate.compute_expectations();
    
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
      exit(0);
    }
    _iter++;
  }
}


void
HGAPRec::compute_likelihood(bool validation)
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
    double u = _env.hier ? rating_likelihood_hier(n,m,r) : rating_likelihood(n,m,r);
    s += u;
    k += 1;
  }

  double a = .0;
  info("s = %.5f\n", s);
  fprintf(ff, "%d\t%d\t%.9f\t%d\n", _iter, duration(), s / k, k);
  fflush(ff);
  a = s / k;  
  
  if (!validation)
    return;
  
  bool stop = false;
  int why = -1;
  if (_iter > 10) {
    if (a > _prev_h && _prev_h != 0 && fabs((a - _prev_h) / _prev_h) < 0.00001) {
      stop = true;
      why = 0;
    } else if (a < _prev_h)
      _nh++;
    else if (a > _prev_h)
      _nh = 0;

    if (_nh > 3) { // be robust to small fluctuations in predictive likelihood
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
HGAPRec::rating_likelihood(uint32_t p, uint32_t q, yval_t y) const
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

double
HGAPRec::rating_likelihood_hier(uint32_t p, uint32_t q, yval_t y) const
{
  const double **etheta = _htheta.expected_v().const_data();
  const double **ebeta = _hbeta.expected_v().const_data();
  
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


double
HGAPRec::log_factorial(uint32_t n)  const
{ 
  double v = log(1);
  for (uint32_t i = 2; i <= n; ++i)
    v += log(i);
  return v;
} 

void
HGAPRec::do_on_stop()
{
  save_model();
  gen_ranking_for_users(false);
}

void
HGAPRec::compute_precision(bool save_ranking_file)
{
  double mhits10 = 0, mhits100 = 0;
  uint32_t total_users = 0;
  FILE *f = 0;
  if (save_ranking_file)
    f = fopen(Env::file_str("/ranking.tsv").c_str(), "w");
  
  if (!save_ranking_file) {
    _sampled_users.clear();
    do {
      uint32_t n = gsl_rng_uniform_int(_r, _n);
      _sampled_users[n] = true;
    } while (_sampled_users.size() < 1000 && _sampled_users.size() < _n / 2);
  } else {
    _sampled_users.clear();
    for (uint32_t n = 0; n < _n; ++n)
      _sampled_users[n] = true;
  }
  
  KVArray mlist(_m);
  for (UserMap::const_iterator itr = _sampled_users.begin();
       itr != _sampled_users.end(); ++itr) {
    uint32_t n = itr->first;
    
    for (uint32_t m = 0; m < _m; ++m) {
      if (_ratings.r(n,m) > 0) { // skip training
	mlist[m].first = m;
	mlist[m].second = .0;
	continue;
      }
      double u = _env.hier ? prediction_score_hier(n, m) : prediction_score(n, m);
      mlist[m].first = m;
      mlist[m].second = u;
    }
    uint32_t hits10 = 0, hits100 = 0;
    mlist.sort_by_value();
    for (uint32_t j = 0; j < mlist.size() && j < _topN_by_user; ++j) {
      KV &kv = mlist[j];
      uint32_t m = kv.first;
      double pred = kv.second;
      Rating r(n, m);

      uint32_t m2 = 0, n2 = 0;
      if (save_ranking_file) {
	IDMap::const_iterator it = _ratings.seq2user().find(n);
	assert (it != _ratings.seq2user().end());
	
	IDMap::const_iterator mt = _ratings.seq2movie().find(m);
	if (mt == _ratings.seq2movie().end())
	  continue;
      
	m2 = mt->second;
	n2 = it->second;
      }

      CountMap::const_iterator itr = _test_map.find(r);
      if (itr != _test_map.end()) {
	int v = itr->second;
	v = _ratings.rating_class(v);
	assert(v > 0);
	if (save_ranking_file) {
	  if (_ratings.r(n, m) == .0) // skip training
	    fprintf(f, "%d\t%d\t%.5f\t%d\n", n2, m2, pred, v);
	}
	
	if (j < 10) {
	  hits10++;
	  hits100++;
	} else if (j < 100) {
	  hits100++;
	}
      } else {
	if (save_ranking_file) {
	  if (_ratings.r(n, m) == .0) // skip training
	    fprintf(f, "%d\t%d\t%.5f\t%d\n", n2, m2, pred, 0);
	}
      }
    }
    mhits10 += (double)hits10 / 10;
    mhits100 += (double)hits100 / 100;
    total_users++;
  }
  if (save_ranking_file)
    fclose(f);
  fprintf(_pf, "%d\t%.5f\t%.5f\n", 
	  total_users,
	  (double)mhits10 / total_users, 
	  (double)mhits100 / total_users);
  fflush(_pf);
}

double
HGAPRec::prediction_score(uint32_t user, uint32_t movie) const
{
  const double **etheta = _theta.expected_v().const_data();
  const double **ebeta = _beta.expected_v().const_data();
  double s = .0;
  for (uint32_t k = 0; k < _k; ++k)
    s += etheta[user][k] * ebeta[movie][k];

  if (_env.bias) {
    const double **ethetabias = _thetabias.expected_v().const_data();
    const double **ebetabias = _betabias.expected_v().const_data();
    s += ethetabias[user][0] + ebetabias[movie][0];
  }
  
  if (_use_rate_as_score)
    return s;
  
  if (s < 1e-30)
    s = 1e-30;
  double prob_zero = exp(-s);
  return 1 - prob_zero;
}

double
HGAPRec::prediction_score_hier(uint32_t user, uint32_t movie) const
{
  const double **etheta = _htheta.expected_v().const_data();
  const double **ebeta = _hbeta.expected_v().const_data();
  double s = .0;
  for (uint32_t k = 0; k < _k; ++k)
    s += etheta[user][k] * ebeta[movie][k];

  if (_env.bias) {
    const double **ethetabias = _thetabias.expected_v().const_data();
    const double **ebetabias = _betabias.expected_v().const_data();
    s += ethetabias[user][0] + ebetabias[movie][0];
  }
  
  if (_use_rate_as_score)
    return s;
  
  if (s < 1e-30)
    s = 1e-30;
  double prob_zero = exp(-s);
  return 1 - prob_zero;
}


void
HGAPRec::gen_ranking_for_users(bool load)
{
  if (load)
    load_beta_and_theta();

  char buf[4096];
  sprintf(buf, "%s/test_users.tsv", _env.datfname.c_str());
  FILE *f = fopen(buf, "r");
  if (!f) { 
    lerr("cannot open %s", buf);
    return;
  }
  //assert(f);
  _sampled_users.clear();
  _ratings.read_test_users(f, &_sampled_users);
  fclose(f);
  compute_precision(true);
  printf("DONE writing ranking.tsv in output directory\n");
  fflush(stdout);
}


void
HGAPRec::load_beta_and_theta()
{
  _beta.load();
  _theta.load();
  if (_env.bias) {
    _betabias.load();
    _thetabias.load();
  }
}

void
HGAPRec::save_model()
{
  if (_env.hier) {
    _hbeta.save_state(_ratings.seq2movie());
    _betarate.save_state(_ratings.seq2movie());
    _htheta.save_state(_ratings.seq2user());
    _thetarate.save_state(_ratings.seq2user());
  } else {
    _beta.save_state(_ratings.seq2movie());
    _theta.save_state(_ratings.seq2user());
  }

  if (_env.bias) {
    _betabias.save_state(_ratings.seq2movie());
    _thetabias.save_state(_ratings.seq2user());
  }
}

void
HGAPRec::logl()
{
  Array phi(_k);
  double s = .0;
  
  const double  **etheta = _theta.expected_v().const_data();
  const double  **ebeta = _beta.expected_v().const_data();
  const double  **elogtheta = _theta.expected_logv().const_data();
  const double  **elogbeta = _beta.expected_logv().const_data();

  for (uint32_t n = 0; n < _n; ++n) {
    const vector<uint32_t> *movies = _ratings.get_movies(n);
    for (uint32_t j = 0; j < movies->size(); ++j) {
      uint32_t m = (*movies)[j];
      yval_t y = _ratings.r(n,m);
      
      get_phi(_theta, n, _beta, m, phi);
      if (y > 1)
	phi.scale(y);

      double v = .0;
      for (uint32_t k = 0; k < _k; ++k)
	v += y * phi[k] * (elogtheta[n][k] + elogbeta[m][k] - log(phi[k]));
      s += v;
      
      for (uint32_t k = 0; k < _k; ++k)
	s -= etheta[n][k] * ebeta[m][k];
    }
  }

  s += _theta.compute_elbo_term();
  s += _beta.compute_elbo_term();
  fprintf(_af, "%.5f\n", s);
  fflush(_af);
}
  

