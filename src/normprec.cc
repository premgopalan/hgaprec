#include "normprec.hh"
#include "utils.hh" 

NormPRec::NormPRec(Env &env, Ratings &ratings)
  : _env(env), _ratings(ratings),
    _n(env.n), _m(env.m), _k(env.k),
    _iter(0),
    _start_time(time(0)),
    _theta("theta", 1e-1, 0.5, _n,_k,&_r),
    _beta("beta", 1e-1, 0.5, _m,_k,&_r),
    _thetabias("thetabias", 1e-5, 0.0075, _n, 1, &_r),
    _betabias("betabias", 1e-5, 0.0075, _m, 1, &_r),
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

double
NormPRec::elbo()
{

  if(_env.bias) {
    printf("elbo not implemented with bias\n");
    assert(0==1);
  }

  uint32_t x;
  if (_env.bias)
    x = _k+2;
  else
    x = _k;

  Array phi(x);
  double s = .0;

  const double  **etheta = _theta.mean_curr().const_data();
  const double  **ebeta = _beta.mean_curr().const_data();
  const double  **eu = NULL;
  const double  **ei = NULL;

  if (_env.bias) {
    eu = _thetabias.mean_curr().const_data();
    ei = _betabias.mean_curr().const_data();
  }

  Array betaexpsum(_k, true);
  betaexpsum.zero();
  _beta.compute_expectations(); 
  _beta.sum_eexp_rows(betaexpsum);
  Array thetaexpsum(_k, true);
  _theta.compute_expectations(); 
  thetaexpsum.zero();
  _theta.sum_eexp_rows(thetaexpsum);

  double train_ll = 0.;
  uint32_t train_ll_k = 0;

  for (uint32_t n = 0; n < _n; ++n) {
    const vector<uint32_t> *movies = _ratings.get_movies(n);
    for (uint32_t j = 0; j < movies->size(); ++j) {
      uint32_t m = (*movies)[j];
      yval_t y = _ratings.r(n,m);

      if (!_env.bias)
        get_phi(_theta, n, _beta, m, phi);
      else {
        const double **tbias = _thetabias.expected_expv().const_data();
        const double **bbias = _betabias.expected_expv().const_data();
        get_phi(_theta, n, _beta, m, tbias[n][0], bbias[m][0], phi);
      }

      s -= log_factorial(y);

      double v = .0;
      for (uint32_t k = 0; k < _k; ++k)
        s += y * phi[k] * (etheta[n][k] + ebeta[m][k] - log(phi[k]));

      if (_env.bias) {
        //s += y * phi[_k] * (elogu[n][0] - log(phi[_k]));
        //s += y * phi[_k+1] * (elogi[m][0] - log(phi[_k+1]));
      }

      if (_env.bias) {
        s -= eu[n][0];
        s -= ei[m][0];
      }

      // train ll
      train_ll += rating_likelihood(n,m,y);
      train_ll_k++;
      train_ll += rating_likelihood(n, gsl_rng_uniform_int(_r, _m), 0);
      train_ll_k += 1;

    }
  }

  for (uint32_t k = 0; k < _k; ++k)
    s -= thetaexpsum[k] * betaexpsum[k];

  s += _theta.elbo();
  s += _beta.elbo();

  if (_env.bias) {
    //s += _thetabias.elbo();
    //s += _betabias.elbo();
  }

  printf("\t\t\tELBO: %f (train ll: %f)\n", s, train_ll/train_ll_k);

  fprintf(_af, "%.5f\n", s);
  fflush(_af);

  return s;

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
  if(_env.gmf_init) {
    _theta.load_from_gmf(_env.datfname, _k);
    _theta.set_to_prior();
    _theta.compute_expectations();
    
    _beta.load_from_gmf(_env.datfname, _k);
    _beta.set_to_prior();
    _beta.compute_expectations();



  } else {
    _beta.initialize();
    _theta.initialize();
    _beta.compute_expectations();
    _theta.compute_expectations();

    if (_env.bias) {
      _thetabias.initialize2(_m);
      //_thetabias.compute_expectations();
        
      _betabias.initialize2(_n);
      //_betabias.compute_expectations();
    }
  }
}

// perform inference
void
NormPRec::vb() 
{
  lerr("running vb()");
  initialize();
  printf("+ initialized\n"); 

  Array thetaexpsum(_k);
  Array betaexpsum(_k);

  Array *phi; 
  //Matrix *phi_m; 
  Array *phi_m; 
  if (_env.bias) {
      phi = new Array(_k+2);
      //phi_m = new Matrix(_m,_k+2); 
      phi_m = new Array(_k+2); 
  } else {
      phi = new Array(_k);
      //phi_m = new Matrix(_m,_k); 
      phi_m = new Array(_k); 
  }

  Array *phi_n, *p; 
  if (_env.bias) {
    phi_n = new Array(_k+2);
    p = new Array(_k+2);
  } else { 
    phi_n = new Array(_k);
    p = new Array(_k);
  }


  while (1) {
    #if TIME
    clock_t start = clock(), diff;
    #endif
  
    //phi_m->zero(); 

    betaexpsum.zero();
    _beta.sum_eexp_rows(betaexpsum);

    for (uint32_t n = 0; n < _n; ++n) { // for every user 
      if(n % 1000 == 0) { 
        printf("+ iter %d\tfor all users %d/%d\t\t\t\r", _iter, n, _n); 
        fflush(stdout); 
      }
      phi_n->zero(); 
      const vector<uint32_t> *movies = _ratings.get_movies(n);
      int nb_rats = movies->size();
      for (uint32_t j = 0; j < movies->size(); ++j) {
        uint32_t m = (*movies)[j];
        yval_t y = _ratings.r(n,m);

        if (!_env.bias)
            get_phi(_theta, n, _beta, m, *phi);
        else { 
            const double **tbias = _thetabias.expected_v().const_data();
            const double **bbias = _betabias.expected_v().const_data();
            get_phi(_theta, n, _beta, m, tbias[n][0], bbias[m][0], *phi);
        }

        if (y > 1)
          phi->scale(y);

        phi_n->add_to(*phi); 
        //phi_m->add_slice(m,*phi);
      }

      _theta.update_var_next(n, betaexpsum); 
      _theta.update_mean_next(n, *phi_n, betaexpsum);
      #if 0
      printf("\n-elbo theta before (user has %d ratings)\n", nb_rats);
      double e = elbo();
      _theta.swap(); _theta.compute_expectations();
      printf("\n elbo after\n");
      double f = elbo();
      if(e > f) {
        printf("error %f >= %f \n", e, f);
        exit(-1);
      }
      _theta.swap(); _theta.compute_expectations();
      #endif

      if (_env.bias) {
        // TODO 
        //_thetabias.update_mean_next(n, (*phi_n)[_k], betasum.const_data());
        //_thetabias.update_var_next(n); 
      }
    }

    _theta.swap(); 
    _theta.compute_expectations(); 
    _theta.set_next_to_zero();

    thetaexpsum.zero();
    _theta.sum_eexp_rows(thetaexpsum);

    for (uint32_t m = 0; m < _m; ++m) { // for every item 
      if(m % 1000 == 0) { 
          printf("+ iter %d\tfor all items %d/%d\t\t\t\r", _iter, m, _m); 
          fflush(stdout); 
        }

      phi_m->zero(); 

      const vector<uint32_t> *users = _ratings.get_users(m);
      int nb_rats = users->size();
      for (uint32_t j = 0; j < users->size(); ++j) {
        uint32_t n = (*users)[j];
        yval_t y = _ratings.r(n,m);

        if (!_env.bias)
            get_phi(_theta, n, _beta, m, *phi);
        else { 
            const double **tbias = _thetabias.expected_v().const_data();
            const double **bbias = _betabias.expected_v().const_data();
            get_phi(_theta, n, _beta, m, tbias[n][0], bbias[m][0], *phi);
        }

        if (y > 1)
          phi->scale(y);

        phi_m->add_to(*phi); 
      }

      _beta.update_mean_next(m, *phi_m, thetaexpsum); 
      _beta.update_var_next(m, thetaexpsum); 
      if (_env.bias) {
        // TODO 
        //_betabias.update_mean_next(m, (*p)[_k+1]);
        //_betabias.update_var_next(m); 
      }

    }
    //for(uint32_t k=0;k<_k;++k)
    //  val += - thetaexpsum[k] * betaexpsum[k];
    //printf("\nf: %f\ne: %f\n", val, elbo()); 
     
    _beta.swap(); 
    _beta.compute_expectations(); 
    _beta.set_next_to_zero();

    if (_env.bias) { 
        _thetabias.swap();
        _betabias.swap();
    }

    fflush(stdout);    
    if (_iter % _env.reportfreq == 0) {
      compute_likelihood(true);
      compute_likelihood(false);
      save_model();
      elbo();
      compute_precision(false);
    } 

    if (_env.save_state_now) {
      lerr("Saving state at iteration %d duration %d secs", _iter, duration());
      do_on_stop();
    }

    _iter++;
    #if TIME
    diff = clock() - start;
    int msec = diff * 1000 / CLOCKS_PER_SEC;
    printf("\t(%d.%ds)\n", msec/1001, msec%1000);
    #endif

  }
  delete phi_n; 
  delete p; 
  delete phi_m;
  delete phi;
}

void
NormPRec::save_model()
{
  _beta.save_state(_ratings.seq2movie());
  _theta.save_state(_ratings.seq2user());

  if(_env.bias) {
    _betabias.save_state(_ratings.seq2movie());
    _thetabias.save_state(_ratings.seq2user());
  }
}

void
NormPRec::compute_precision(bool save_ranking_file)
{
    double mhits10 = 0, mhits100 = 0;
    double cumndcg10 = 0, cumndcg100 = 0;
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
        for (uint32_t m = 0; m < _m; ++m) {
            Rating r(n,m);
            if (_ratings.r(n,m) > 0 || is_validation(r)) { // skip training and validation
                mlist[m].first = m;
                mlist[m].second = .0;
                ndcglist[m].first = m;
                ndcglist[m].second = 0; 
                continue;
            }
            double u = .0;
            u = prediction_score(n, m);
            mlist[m].first = m;
            mlist[m].second = u;
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
                if (_ratings.test_hit(v_))
                    v = 1;
                else
                    v = 0;

                if (j < 10) {
                    if (v > 0) { //hit
                        hits10++;
                        hits100++;
                    }
                    if (v_ > 0) { //has non-zero relevance
                        dcg10 += (pow(2.,v_) - 1)/log(j+2);
                        dcg100 += (pow(2.,v_) - 1)/log(j+2);
                    }
                } else if (j < 100) {
                    if (v > 0)
                        hits100++;
                    if (v_ > 0)
                        dcg100 += (pow(2.,v_) - 1)/log(j+2);
                }

                if (save_ranking_file) {
                    if (_ratings.r(n, m) == .0)  {
                        double hol = rating_likelihood(n,m,v);
                        //fprintf(f, "%d\t%d\t%.5f\t%d\t%.5f\n", n2, m2, pred, v,
                        //(pow(2.,v_) - 1)/log(j+2));
                        fprintf(f, "%d\t%d\t%.5f\t%d\n", n2, m2, pred, v);
                    }
                }
            } else {
                if (save_ranking_file) {
                    if (_ratings.r(n, m) == .0) {
                        double hol = rating_likelihood(n,m,0);
                        //fprintf(f, "%d\t%d\t%.5f\t%d\t%.5f\n", n2, m2, pred, 0, .0);
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
            if(v==0) { //all subsequent docs are irrelevant
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

void
NormPRec::get_phi(NormBase<Matrix> &a, uint32_t ai, 
		 NormBase<Matrix> &b, uint32_t bi, 
		 Array &phi)
{
  assert (phi.size() == a.k() &&
	  phi.size() == b.k());
  assert (ai < a.n() && bi < b.n());
  const double  **eloga = a.expected_v().const_data();
  const double  **elogb = b.expected_v().const_data();
  phi.zero();
  for (uint32_t k = 0; k < _k; ++k)
    phi[k] = eloga[ai][k] + elogb[bi][k];
  phi.lognormalize();
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
  phi.zero();
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

    // add a fake 0 for that user
    //s += rating_likelihood(n, gsl_rng_uniform_int(_r, _m), 0); 
    //k += 1; 
  }

  double a = .0;
  a = s / k;  
  info("s = %.5f\n", s);
  fprintf(ff, "%d\t%d\t%.9f\t%d\n", _iter, duration(), a, k);
  fflush(ff);

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
NormPRec::score(uint32_t p, uint32_t q) const
{
  const double **etheta = _theta.expected_expv().const_data();
  const double **ebeta = _beta.expected_expv().const_data();
  
  double s = .0;
  for (uint32_t k = 0; k < _k; ++k)
    s += etheta[p][k] * ebeta[q][k];
  
  // TODO: biases!!! 
  if (_env.bias) {
    const double **ethetabias = _thetabias.expected_expv().const_data();
    const double **ebetabias = _betabias.expected_expv().const_data();
    s += ethetabias[p][0] + ebetabias[q][0];
  } 
  
  if (s < 1e-30)
    s = 1e-30;

  return s; 

}


double
NormPRec::rating_likelihood(uint32_t p, uint32_t q, yval_t y) const
{

  double s = score(p,q);

  if (_env.binary_data)
    return y == 0 ? -s : log(1 - exp(-s));    
  return y * log(s) - s - log_factorial(y);
}

void
NormPRec::do_on_stop()
{
  save_model();
  gen_ranking_for_users(false);
}

void
NormPRec::gen_ranking_for_users(bool load)
{
  if (load) { 
    printf("unimplemented\n"); exit(-1); 
  } 

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

double
NormPRec::prediction_score(uint32_t p, uint32_t q) const
{
  double s = score(p, q);

  if (_use_rate_as_score)
    return s;
  
  if (s < 1e-30)
    s = 1e-30;
  double prob_zero = exp(-s);
  return 1 - prob_zero;
}
