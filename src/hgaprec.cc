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
    _nmf_theta(_n, _k),
    _nmf_beta(_m, _k),
    _prev_h(.0), _nh(.0),
    _save_ranking_file(false),
    _use_rate_as_score(true),
    _topN_by_user(100)
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

  for (CountMap::const_iterator i = _validation_map.begin();
       i != _validation_map.end(); ++i) {
    const Rating &r = i->first;
    _validation_users_of_movie[r.second]++;
  }

  sprintf(buf, "%s/test.tsv", _env.datfname.c_str());
  FILE *testf = fopen(buf, "r");
  assert(testf);
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
HGAPRec::initialize()
{
  if (!_env.hier) {
    _beta.initialize();
    _theta.initialize();

    _beta.initialize_exp();
    _theta.initialize_exp();

    //_theta.initialize2(_m);
    //_theta.compute_expectations();

    //_beta.initialize2(_n);
    //_beta.compute_expectations();

  } else {

    //_thetarate.set_to_prior_curr();
    //_thetarate.set_to_prior();

    _thetarate.initialize2(_k);
    _thetarate.compute_expectations();

    _betarate.initialize2(_k);
    _betarate.compute_expectations();
    
    //_betarate.set_to_prior_curr();
    //_betarate.set_to_prior();
    
    //_hbeta.initialize2(_n);
    //_hbeta.compute_expectations();
    
    _hbeta.initialize();
    _hbeta.initialize_exp();

    //_hbeta.initialize_exp(_betarate.expected_v()[0]);

    //_htheta.initialize2(_m);
    //_htheta.compute_expectations();
    
    _htheta.initialize();
    _htheta.initialize_exp();

    //_htheta.initialize_exp(_thetarate.expected_v()[0]);
  }
  
  if (_env.bias) {
    _thetabias.initialize2(_m);
    _thetabias.compute_expectations();
    
    _betabias.initialize2(_n);
    _betabias.compute_expectations();
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
HGAPRec::write_training_matrix()
{
  FILE *f = fopen(Env::file_str("/trainm.tsv").c_str(), "w");
  for (uint32_t n = 0; n < _n; ++n) {
    const vector<uint32_t> *movies = _ratings.get_movies(n);
    if (movies->size() == 0)
      continue;
    
    uint32_t c = 0;
    for (uint32_t m = 0; m < _m; ++m)
      if (_ratings.r(n,m) == 0)
	c++;

    if (c == _m - 1) // null row
      continue;

    for (uint32_t m = 0; m < _m; ++m)
      fprintf(f, "%d\t", _ratings.r(n,m));
    fprintf(f, "\n");
  }
  fclose(f);
}

void
HGAPRec::load_nmf_beta_and_theta()
{
  char buf[4096];
  _nmf_theta.load("theta.tsv", 0);
  _nmf_beta.load("beta.tsv", 0);

  IDMap m;
  _nmf_theta.save(Env::file_str("/nmf_theta.tsv").c_str(), m);
  _nmf_beta.save(Env::file_str("/nmf_beta.tsv").c_str(), m);

  sprintf(buf, "%s/test_users.tsv", _env.datfname.c_str());
  FILE *f = fopen(buf, "r");
  if (!f) { 
    lerr("cannot open %s", buf);
    return;
  }
  _sampled_users.clear();
  _ratings.read_test_users(f, &_sampled_users);
  fclose(f);
  compute_precision(true);
  printf("DONE writing ranking.tsv in output directory\n");
  fflush(stdout);
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
      //gen_ranking_for_users(false);
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
    
    if (_env.vb) {
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

      debug("thetabias: %s", _thetabias.expected_v().s().c_str());
      debug("betabias: %s", _betabias.expected_v().s().c_str());

    } else {

      Array betasum(_k);
      _beta.sum_rows(betasum);
      _theta.update_rate_next(betasum);
      Array thetasum(_k);
      _theta.sum_rows(thetasum);
      _beta.update_rate_next(thetasum);

      _thetabias.update_rate_next(0, _m);
      _betabias.update_rate_next(0, _n);
      
      _theta.swap();
      _beta.swap();
      _thetabias.swap();
      _betabias.swap();

      _theta.compute_expectations();
      _beta.compute_expectations();
      _thetabias.compute_expectations();
      _betabias.compute_expectations();
    }

    printf("\r iteration %d", _iter);
    fflush(stdout);    
    if (_iter % _env.reportfreq == 0) {
      compute_likelihood(true);
      compute_likelihood(false);
      save_model();
      compute_precision(false);
      //gen_ranking_for_users(false);
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
      //gen_ranking_for_users(false);
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
  if (_iter > 30) {
    if (a > _prev_h && _prev_h != 0 && fabs((a - _prev_h) / _prev_h) < 0.000001) {
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
  }
  
  KVArray mlist(_m);
  KVIArray ndcglist(_m);
  for (UserMap::const_iterator itr = _sampled_users.begin();
       itr != _sampled_users.end(); ++itr) {
    uint32_t n = itr->first;
    
    //_topN_by_user = _m;
    //for (uint32_t n = 0; n < _n; n++) {
    
    for (uint32_t m = 0; m < _m; ++m) {
      //if (_ratings.r(n,m) > 0) { // skip training
      Rating r(n,m);
      if (_ratings.r(n,m) > 0 || is_validation(r)) { // skip training and validation
	mlist[m].first = m;
	mlist[m].second = .0;
    ndcglist[m].first = m;
    ndcglist[m].second = 0;
	continue;
      }
      double u = .0;
      if (_env.nmf)
	u = prediction_score_nmf(n, m);
      else
	u = _env.hier ? prediction_score_hier(n, m) : prediction_score(n, m);
      mlist[m].first = m;
      mlist[m].second = u;
      Rating r(n,m); 
      ndcglist[m].first = m;
      CountMap::const_iterator itr = _test_map.find(r);
      if (itr != _test_map.end()) {
          ndcglist[m].second = itr->second;
      } else { 
          ndcglist[m].second = 0;
      }
    }
    uint32_t hits10 = 0, hits100 = 0;
    double   dcg10 = .0, dcg100 = .0; 
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
    int v_ = itr->second;
    int v;
	if (v_ >= 1)
	  v = 1;
	else
	  v = 0;
	if (save_ranking_file) {
	  if (_ratings.r(n, m) == .0)  {
	    double hol = _env.hier ? rating_likelihood_hier(n,m,v) : rating_likelihood(n,m,v);
	    fprintf(f, "%d\t%d\t%.5f\t%d\n", n2, m2, pred, v);
	  }
	}
	
	if (v > 0)
	  if (j < 10) {
	    hits10++;
	    hits100++;
        dcg10 += (pow(2.,v_) - 1)/log(j+2);
        dcg100 += (pow(2.,v_) - 1)/log(j+2);
	  } else if (j < 100) {
	    hits100++;
        dcg100 += (pow(2.,v_) - 1)/log(j+2);
	  }
      } else {
	if (save_ranking_file) {
	  if (_ratings.r(n, m) == .0) {
	    double hol = _env.hier ? rating_likelihood_hier(n,m,0) : rating_likelihood(n,m,0);
	    fprintf(f, "%d\t%d\t%.5f\t%d\n", n2, m2, pred, 0);
	  }
	}
      }
    }
    mhits10 += (double)hits10 / 10;
    mhits100 += (double)hits100 / 100;
    total_users++;
    // DCG normalizer
    double dcg10_gt = 0, dcg100_gt = 0;
    bool user_has_test_ratings = true; 
    ndcglist.sort_by_value();
    for (uint32_t j = 0; j < ndcglist.size() && j < _topN_by_user; ++j) {
        int v = ndcglist[j].second; 
        if(v==0) { // all subsequent docs are irrelevant
          if(j==0)
              user_has_test_ratings = false; 
          break;
        }

        if (j < 10) { 
            dcg10_gt += (pow(2.,v) - 1)/log(j+2);
            dcg100_gt += (pow(2.,v) - 1)/log(j+2);
        } else if (j < 100) {
            dcg100_gt += (pow(2.,v) - 1)/log(j+2);
        }
    }
    if(user_has_test_ratings) { 
        cumndcg10 += dcg10/dcg10_gt;
        cumndcg100 += dcg100/dcg100_gt;
    } 
  }
  if (save_ranking_file)
    fclose(f);
  fprintf(_pf, "%d\t%.5f\t%.5f\n", 
	  total_users,
	  (double)mhits10 / total_users, 
	  (double)mhits100 / total_users);
  fflush(_pf);
  fprintf(_df, "%.5f\t%.5f\n", 
  	  cumndcg10 / total_users, 
  	  cumndcg100 / total_users);
  fflush(_df);
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
HGAPRec::prediction_score_nmf(uint32_t user, uint32_t movie) const
{
  const double **etheta = _nmf_theta.const_data();
  const double **ebeta = _nmf_beta.const_data();
  double s = .0;
  for (uint32_t k = 0; k < _k; ++k)
    s += etheta[user][k] * ebeta[movie][k];
  return s;
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
HGAPRec::gen_msr_csv()
{
  load_beta_and_theta();

  FILE *f = fopen(Env::file_str("/pred.csv").c_str(), "w");
  if (!f) {
    lerr("cannot open pred.csv");
    return;
  }
  
  fprintf(f, "User\tHeldOutItem\tHeldOutItemIndex\tUserNegatives\tUserCount\tItemCount\n");
  for (uint32_t n = 0; n < _n; ++n) {
    // User
    IDMap::const_iterator it = _ratings.seq2user().find(n);
    assert (it != _ratings.seq2user().end());
    fprintf(f, "%d\t", it->second);
    debug("user: %d", it->second);

    // id of heldout item
    IDMap::const_iterator ct = _leave_one_out.find(n);
    debug("heldout item for user %d (%d)",   it->second, n);
    assert (ct != _leave_one_out.end());

    uint32_t test_item_seq = ct->second;
    IDMap::const_iterator pt = _ratings.seq2movie().find(test_item_seq);
    assert (pt != _ratings.seq2movie().end());
    fprintf(f, "%d\t", pt->second);
    debug("test item: %d", pt->second);

    // rank of heldout test item
    uint32_t training = 0, negatives = 0;
    KVArray mlist(_m);
    for (uint32_t m = 0; m < _m-1; ++m) {
      Rating r(n, m);
      if (_ratings.r(n,m) > 0 || is_validation(r)) { // skip training non-zero rating
	mlist[m].first = m;
	mlist[m].second = .0;
	training++;
	continue;
      }
      double u = .0;
      if (_env.nmf)
	u = prediction_score_nmf(n, m);
      else
	u = _env.hier ? prediction_score_hier(n, m) : prediction_score(n, m);
      mlist[m].first = m;
      mlist[m].second = u;
      negatives++;
    }
    mlist.sort_by_value();

    uint32_t c = 0, rank = 0;
    for (uint32_t j = 0; j < mlist.size(); ++j) {
      KV &kv = mlist[j];
      uint32_t m = kv.first;
      double pred = kv.second;
      if (m == test_item_seq) {
	rank = c;
	break;
      }
      c++;
    }
    fprintf(f, "%d\t", rank);
    debug("rank: %d", rank);
    
    // UserNegatives
    fprintf(f, "%d\t", negatives);
    debug("user negatives: %d", negatives);

    // UserCount
    fprintf(f, "%d\t", training);
    debug("user count: %d", training);
    
    const vector<uint32_t> *q = _ratings.get_users(test_item_seq);
    uint32_t ntraining_users = 0;
    if (q)
      ntraining_users = q->size();
    
    uint32_t nvalid_users = 0;
    FreqMap::const_iterator itr = _validation_users_of_movie.find(test_item_seq);
    if (itr != _validation_users_of_movie.end())
      nvalid_users = itr->second;
    
    // ItemCount
    fprintf(f, "%d\n", nvalid_users + ntraining_users);
    debug("item count: %d", nvalid_users + ntraining_users);

    printf("\r user %d", n);
    fflush(stdout);
  }
  
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
  if (!_env.hier) {
    _beta.load();
    _theta.load();
  } else {
    _thetarate.load();
    _betarate.load();

    _hbeta.load();
    _htheta.load();
  }
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
  

