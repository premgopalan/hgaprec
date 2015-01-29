#ifndef NORMBASE_HH
#define NORMBASE_HH

#include <gsl/gsl_multimin.h>
#include <gsl/gsl_rng.h>
#include <cmath> 

#define CHECK_GRADIENTS 0


template <class T>
class NormBase { 

public:
  NormBase(string name = ""): _name(name) { }
  virtual ~NormBase() { }
  virtual const T &expected_v() const = 0;
  virtual uint32_t n() const = 0;
  virtual uint32_t k() const = 0;
  virtual void save_state(const IDMap &m) const = 0;
  virtual void load()  = 0;
  string name() const { return _name; }
  double elbo() const;
private:
  string _name;
};

const uint32_t _cg_max_iter = 500; // TODO
const float _cg_convergence  = 1e-5; // TODO

class NormMatrix : public NormBase<Matrix> {

public:
    NormMatrix(string name, double c, double d,
	   uint32_t n, uint32_t k,
	   gsl_rng **r): 
    NormBase<Matrix>(name),
    _n(n), _k(k),
    _mprior(c), // mean 
    _vprior(d), // variance
    _mcurr(n,k),
    _mnext(n,k),
    _vnext(n,k),
    _vcurr(n,k),
    _Eexpv(n,k),
    _r(r) { } 
  virtual ~NormMatrix() {} 

  void vprint(const gsl_vector *x);
  double vnorm(const gsl_vector *x);

  uint32_t n() const { return _n;}
  uint32_t k() const { return _k;}

  const Matrix &mean_curr() const         { return _mcurr; }
  const Matrix &var_curr() const          { return _vcurr; }
  const Matrix &mean_next() const         { return _mnext; }
  const Matrix &var_next() const          { return _vnext; }
  const Matrix &expected_v() const         { return mean_curr(); }
  const Matrix &expected_expv() const  { return _Eexpv; }

  Matrix &mean_curr()       { return _mcurr; }
  Matrix &var_curr()        { return _vcurr; }
  Matrix &mean_next()       { return _mnext; }
  Matrix &var_next()        { return _vnext; }
  Matrix &expected_v()       { return _mcurr;    }

  const double mprior() const { return _mprior; }
  const double vprior() const { return _vprior; }
  void set_to_prior();
  void set_next_to_zero();
  void initialize(); 
  void initialize2(double v);
  float update_mean_next(uint32_t n, const Array &phi, const Array &other);
  float update_var_next(uint32_t n, const Array &other);
  void sum_rows(Array &);
  void sum_eexp_rows(Array &);
  double elbo();
  void compute_expectations();

  void save_state(const IDMap &m) const;
  void load() {}
  void load_from_gmf(string, uint32_t, string);

  double f_user(void *params);
  double f_mean(const gsl_vector * p, void * params); 
  void df_mean(const gsl_vector * x, void * params, gsl_vector * g); 
  void fdf_mean(const gsl_vector * x, void * params, double * f, gsl_vector *g); 
  double f_var(const gsl_vector * p, void * params); 
  void df_var(const gsl_vector * x, void * params, gsl_vector * g); 
  void fdf_var(const gsl_vector * x, void * params, double * f, gsl_vector *g); 

  static double wrap_f_mean(const gsl_vector *x, void *params);
  static void wrap_df_mean(const gsl_vector *x, void *params, gsl_vector *g);
  static void wrap_fdf_mean(const gsl_vector *x, void *params, double *f, gsl_vector *g);
  static double wrap_f_var(const gsl_vector *x, void *params);
  static void wrap_df_var(const gsl_vector *x, void *params, gsl_vector *g);
  static void wrap_fdf_var(const gsl_vector *x, void *params, double *f, gsl_vector *g);

  void swap(); 

private:
  uint32_t _n;
  uint32_t _k;	
  gsl_rng **_r;
  double _mprior;
  double _vprior;

  Matrix _mcurr;      // current variational mean posterior 
  Matrix _mnext;      // to help compute gradient update
  Matrix _vcurr;      // current variational variance posterior
  Matrix _vnext;      // help compute gradient update
  Matrix _Eexpv;      // expected exp weights under variational
        		      // distribution
};

typedef struct bundle {
    NormMatrix *NormMatrixObj; 
    const Array *phi, *other;
    uint32_t id; 
} bundle;

inline void
NormMatrix::load_from_gmf(string dir, uint32_t K, string name="")
{ 
  char buf_mean[1024], buf_var[1024];
  
  if(name=="")
    name = this->name();
  
  sprintf(buf_mean, "%s/gmf-fits/%s_mean-k%d.tsv", dir.c_str(), name.c_str(), K);
  lerr("loading from %s", buf_mean);
  sprintf(buf_var, "%s/gmf-fits/%s_var-k%d.tsv", dir.c_str(), name.c_str(), K);
  lerr("loading from %s", buf_var);
  
  _mcurr.load(buf_mean, 0);
  // gmf keeps a single variance for all entities (users/items) 
  Array var(1); 
  var.load(buf_var);

  double **md = _mcurr.data();
  double **vd = _vcurr.data();
  double **ed = _Eexpv.data();

  // go from normal to log-normal parameters
  //double mean2, var;
  double tmp2;
  for (uint32_t i = 0; i < _n; ++i) {
    for (uint32_t k = 0; k < _k; ++k) {

      //md[i][k] = exp(md[i][k] + var[0]*var[0]/2);
      //vd[i][k] = md[i][k]*md[i][k] * (exp(var[0]*var[0]) -1);
      tmp2 = md[i][k]*md[i][k]; 
      md[i][k] = log(md[i][k] / sqrt(1 + var[0]/tmp2));
      vd[i][k] = sqrt(log(1+ var[0]/tmp2));

    }
  }

  this->compute_expectations();

  IDMap m;
  string expv_fname = string("/") + name + "_exp.tsv";
  _Eexpv.save(Env::file_str(expv_fname), m);
}


inline void
NormMatrix::save_state(const IDMap &m) const
{
  string mean_fname = string("/") + name() + "_mean.tsv";
  string var_fname = string("/") + name() + "_variance.tsv";
  string Eexp_fname = string("/") + name() + ".tsv";
  _mcurr.save(Env::file_str(mean_fname), m);
  _vcurr.save(Env::file_str(var_fname), m);
  _Eexpv.save(Env::file_str(Eexp_fname), m);

}

inline void
NormMatrix::compute_expectations()
{ 
  // compute expectation at the point estimates of the var. distribution
  const double ** const md = _mcurr.const_data();
  const double ** const vd = _vcurr.const_data();
  double **ed1 = _Eexpv.data();
  double a = .0, b = .0;
  for (uint32_t i = 0; i < _mcurr.m(); ++i)
    for (uint32_t j = 0; j < _mcurr.n(); ++j) {
      ed1[i][j] = exp(md[i][j] + vd[i][j]/2);
    }
}


inline void
NormMatrix::sum_rows(Array &v)
{
    const double **ev = _mcurr.const_data();
    for (uint32_t i = 0; i < _n; ++i)
        for (uint32_t k = 0; k < _k; ++k)
            v[k] += ev[i][k];
} 

inline void
NormMatrix::sum_eexp_rows(Array &v)
{
    const double **eexpv = _Eexpv.const_data();
    for (uint32_t i = 0; i < _n; ++i)
        for (uint32_t k = 0; k < _k; ++k)
            v[k] += eexpv[i][k];
}

inline double
NormMatrix::elbo()
{
  double s = 0.;
  double vsquared;
  double prior_mean;

  const double **vd = _vcurr.const_data();
  const double **md = _mcurr.const_data();

  // A, B, C, D
  for(int n=0; n<_n; ++n) {

    for(int k=0; k<_k; ++k)  {
      vsquared = vd[n][k];

      // sigma in derivation -> _vprior here.
      // -log(d\sqrt(2*pi)) -d^2\nu^2
      //s += - log(_vprior * sqrt(2*M_PI)) - _vprior*_vprior * vsquared;
      s += - _vprior * vsquared;

      // - (mu - mu_t-1)^2 / (2*d^2)
      s -= std::pow(md[n][k] - _mprior, 2) / (2 * _vprior);

      // 1/2 * ( log \nu^2 + log 2\pi + 1)
      //s += 0.5 * ( log(vsquared) + log(2*M_PI) + 1);
      s += 0.5 * ( log(vsquared) );

    } // k

  } // n

  return s;
}


// static wrapper-function to be able to callback the member function Display()
inline double 
NormMatrix::wrap_f_mean(const gsl_vector *x, void * params)
{
   NormMatrix& obj = (*((bundle *) params)->NormMatrixObj);
   return obj.f_mean(x, params);
}

inline void 
NormMatrix::wrap_df_mean(const gsl_vector *x, void *params, gsl_vector *g)
{
   NormMatrix& obj = (*((bundle *) params)->NormMatrixObj);
   obj.df_mean(x, params, g);
}

inline void 
NormMatrix::wrap_fdf_mean(const gsl_vector *x, void * params, double *f, gsl_vector *g)
{
   //NormMatrix* me = (NormMatrix*) NormMatrixObj;
   NormMatrix& obj = (*((bundle *) params)->NormMatrixObj);
   obj.fdf_mean(x, params, f, g);
}

// static wrapper-function to be able to callback the member function Display()
inline double 
NormMatrix::wrap_f_var(const gsl_vector *x, void * params)
{
   NormMatrix& obj = (*((bundle *) params)->NormMatrixObj);
   return obj.f_var(x, params);
}

inline void 
NormMatrix::wrap_df_var(const gsl_vector *x, void *params, gsl_vector *g)
{
   NormMatrix& obj = (*((bundle *) params)->NormMatrixObj);
   obj.df_var(x, params, g);
}

inline void 
NormMatrix::wrap_fdf_var(const gsl_vector *x, void * params, double *f, gsl_vector *g)
{
   NormMatrix& obj = (*((bundle *) params)->NormMatrixObj);
   obj.fdf_var(x, params, f, g);
}

inline double
NormMatrix::f_user(void *params)
{
    bundle &b = (bundle &) (*(bundle *) (params));
    NormMatrix& obj = (*b.NormMatrixObj);
    uint32_t id = ((bundle *) params)->id;

    const double ** const md = obj._mcurr.const_data();

    double f = .0;

    gsl_vector *tmp = gsl_vector_alloc(obj._k);

    for(uint32_t k=0; k<_k; ++k) {
      gsl_vector_set(tmp,k,md[id][k]);
      //gsl_vector_set(tmp,k,vd[id][k]);
    }

    // - (x_k - c)^2 / 2*d^2
    //gsl_vector_memcpy(tmp, x);
    gsl_vector_add_constant(tmp, -obj._mprior);
    gsl_vector_mul(tmp, tmp);
    gsl_vector_scale(tmp, 1/(2*obj._vprior));  // OK

    const double ** const vd = obj._vcurr.const_data();
    const double * const pd = b.phi->const_data();
    const double * const od = b.other->const_data();

    double x_k;
    for(uint32_t k=0; k<obj._k; ++k) {
        f -= gsl_vector_get(tmp, k);
        f -= exp(md[id][k] + vd[id][k]/2) * od[k]; // - exp(x + nu^2/2) * exp(other + nu^2_other/2)
        f += pd[k] * md[id][k]; // \phi_mn * x
        f -= vd[id][k]*(obj._vprior);
        f += 0.5*log(vd[id][k]);
    }

    gsl_vector_free(tmp); 
    return f;
}

inline double 
NormMatrix::f_mean(const gsl_vector *x, void * params)
{
    bundle &b = (bundle &) (*(bundle *) (params));
    NormMatrix& obj = (*b.NormMatrixObj);
    uint32_t id = ((bundle *) params)->id;

    double f = .0; 

    //f = -(x_k - c)^2 / 2*d^2 - exp(x + nu^2/2) + \phi_mn * x

    gsl_vector *tmp = gsl_vector_alloc(obj._k);

    // - (x_k - c)^2 / 2*d^2
    gsl_vector_memcpy(tmp, x);
    gsl_vector_add_constant(tmp, -obj._mprior);
    gsl_vector_mul(tmp, tmp);
    gsl_vector_scale(tmp, 1/(2*obj._vprior));

    const double ** const vd = obj._vcurr.const_data();
    const double * const pd = b.phi->const_data();
    const double * const od = b.other->const_data();

    double x_k; 
    for(uint32_t k=0; k<obj._k; ++k) { 
        f -= gsl_vector_get(tmp, k); 
        x_k = gsl_vector_get(x, k); 
        f -= exp(x_k + vd[id][k]/2) * od[k]; // - exp(x + nu^2/2) * exp(other + nu^2_other/2)
        f += pd[k] * x_k; // \phi_mn * x
    }

    gsl_vector_free(tmp); 
    return -f; // maximize
}


inline void 
NormMatrix::df_mean(const gsl_vector * x, void * params, gsl_vector * df)
{
    bundle &b = (bundle &) (*(bundle *) (params));
    NormMatrix& obj = (*b.NormMatrixObj);
    uint32_t id = ((bundle *) params)->id;

    double f = .0;

    //f = x_k / d^2 - exp(x + nu^2/2) + \phi_mn

    // - (x_k-c) / d^2 
    gsl_vector_memcpy(df, x);
    gsl_vector_add_constant(df, -obj._mprior);
    gsl_vector_scale(df, -1/(obj._vprior));

    const double ** const vd = obj._vcurr.const_data();
    const double * const pd = b.phi->const_data(); 
    const double * const od = b.other->const_data(); 

    double x_k, df_k; 
    for(uint32_t k=0; k<obj._k; ++k) { 
        x_k = gsl_vector_get(x, k); 
        gsl_vector_set(df, k, gsl_vector_get(df,k) -exp(x_k + vd[id][k]/2)*od[k] + pd[k]);
    }

    gsl_vector_scale(df, -1); // maximize
}

inline void 
NormMatrix::fdf_mean(const gsl_vector * x, void * params, double * f, gsl_vector * df)
{
    *f = f_mean(x, params);
    df_mean(x, params, df);
}

inline void
NormMatrix::set_next_to_zero()
{
  _mnext.zero();
  _vnext.zero();
}


inline void
NormMatrix::set_to_prior()
{
  _mnext.set_elements(0.);
  _vnext.set_elements(0.);
}

inline void
NormMatrix::swap()
{
  _mcurr.swap(_mnext);
  _vcurr.swap(_vnext);
  set_to_prior();
}


inline void
NormMatrix::initialize()
{
  // TODO: initialize from Gaussian(_mprior, vprior)
  // ideally if we could init from a multi-variate gaussian it may be
  // faster 

  double **ad = _mcurr.data();
  double **bd = _vcurr.data();
  for (uint32_t i = 0; i < _n; ++i) 
    for (uint32_t k = 0; k < _k; ++k) 
      ad[i][k] = _mprior + gsl_rng_uniform(*_r);

      //gsl_ran_gaussian_ziggurat(*_r, _vprior)

  for (uint32_t i = 0; i < _n; ++i)
    for (uint32_t k = 0; k < _k; ++k)
      bd[i][k] = _vprior + 0.1 * gsl_rng_uniform(*_r);
  set_next_to_zero();
}

inline void
NormMatrix::initialize2(double v)
{
  double **ad = _mcurr.data();
  double **bd = _vcurr.data();
  for (uint32_t i = 0; i < _n; ++i) {
    for (uint32_t k = 0; k < _k; ++k) {
      ad[i][k] = _mprior + 0.01 * gsl_rng_uniform(*_r);
      bd[i][k] = _vprior + 0.01 * gsl_rng_uniform(*_r);
    }
  }
  set_to_prior();
}

inline double
NormMatrix::vnorm(const gsl_vector *x) {
  double n=0.;
  for (uint32_t k=0; k<x->size; ++k)
    n += gsl_vector_get(x,k) * gsl_vector_get(x,k);
  return sqrt(n);
}

inline void
NormMatrix::vprint(const gsl_vector *x) {
    for (uint32_t k=0; k<x->size; ++k)
        printf("%f ", gsl_vector_get(x, k));
    printf("\n");
}

inline double 
NormMatrix::f_var(const gsl_vector *x, void * params)
{
    NormMatrix& obj = (*((bundle *) params)->NormMatrixObj);
    uint32_t id = ((bundle *) params)->id;
    const double * const od = (((bundle *) params)->other)->const_data();
    const double ** const md = obj._mcurr.const_data();

    double f = .0; 
    double exp_x_k; 
    for(uint32_t k=0; k<obj._k; ++k) { 

        exp_x_k = exp(gsl_vector_get(x, k)); // variance must be positive
        // -d^2 \nu_k^2 + 0.5 * log(\nu^2) - exp(x + nu^2/2)
        f -= exp_x_k*obj._vprior;
        //f += 0.5*gsl_vector_get(x, k);
        f += 0.5*log(exp_x_k);
        f -= exp(md[id][k] + exp_x_k/2) * od[k];
    }

    return -f; // maximize
}


inline void 
NormMatrix::df_var(const gsl_vector * x, void * params, gsl_vector * df)
{
    NormMatrix& obj = (*((bundle *) params)->NormMatrixObj);
    uint32_t id = ((bundle *) params)->id;

    //df = -d^2 + 1/(2nu_k^2) - 1/2exp(x + nu^2/2)

    const double ** const md = obj._mcurr.const_data();
    const double * const od = (((bundle *) params)->other)->const_data();

    double exp_x_k; 
    for(uint32_t k=0; k<obj._k; ++k) { 
        exp_x_k = exp(gsl_vector_get(x, k)); 
        gsl_vector_set(df, k, -exp_x_k*obj._vprior
                              + 0.5 
                              - 0.5*exp(md[id][k] + gsl_vector_get(x,k) + exp_x_k/2)*od[k]);
                              //- exp_x_k*0.5*exp(md[id][k] + exp_x_k/2)*od[k]);
    }
    gsl_vector_scale(df, -1.); // maximize
}

inline void 
NormMatrix::fdf_var(const gsl_vector * x, void * params, double * f, gsl_vector * df)
{
    *f = f_var(x, params);
    df_var(x, params, df);
}


inline float
NormMatrix::update_mean_next(uint32_t n, const Array &phi, const Array &other)
{
        const double * const pd = phi.const_data(); 

    // update the K-vector mean of a particular user/item
    // involves 
    // 1) phi 
    // 2) priors (c,d)
    // 3) nu
    // keep pointers to all in a bundle 

    gsl_multimin_fdfminimizer * s;
    gsl_multimin_function_fdf mu_obj;

    int iter = 0, i, j;
    int status;
    double f_old, converged;

    bundle b;
    b.NormMatrixObj = this; 
    b.phi = &phi;
    b.other = &other; 
    b.id = n; 

    //printf("optimizing mean %d/%d\n", n, _n); 
    
    mu_obj.f = &NormMatrix::wrap_f_mean;
    mu_obj.df = &NormMatrix::wrap_df_mean;
    mu_obj.fdf = &NormMatrix::wrap_fdf_mean;
    mu_obj.n = _k;
    mu_obj.params = (void *)&b;

    // starting value
    const gsl_multimin_fdfminimizer_type * T;
    //T = gsl_multimin_fdfminimizer_conjugate_fr;
    T = gsl_multimin_fdfminimizer_vector_bfgs2;
    s = gsl_multimin_fdfminimizer_alloc(T, _k);

    gsl_vector* x = gsl_vector_calloc(_k);
    const double ** md = _mcurr.const_data();
    for (uint32_t k = 0; k < _k; ++k) gsl_vector_set(x, k, md[n][k]);
    // tmp 
    #if 0
    gsl_vector* xt = gsl_vector_calloc(_k);
    gsl_vector_memcpy(xt, x); 
    #endif

    #if CHECK_GRADIENTS
    gsl_vector *xg = gsl_vector_alloc(_k);
    gsl_vector *dfh  = gsl_vector_alloc(_k);
    double eps = 1e-5;
    double f1, f2;
    for(uint32_t k=0; k < _k; ++k) {
      gsl_vector_memcpy(xg, x);
      gsl_vector_set(xg, k, gsl_vector_get(xg, k) + eps);
      f1 = wrap_f_mean(xg, (void *)&b);
      gsl_vector_set(xg, k, gsl_vector_get(xg, k) - 2*eps);
      f2 = wrap_f_mean(xg, (void *)&b);
      gsl_vector_set(dfh, k, (f1-f2)/(2*eps));
    }
    gsl_vector* df = gsl_vector_calloc(_k);
    wrap_df_mean(x, (void *)&b, df);
    gsl_vector *t1 = gsl_vector_calloc(_k);
    gsl_vector *t2 = gsl_vector_calloc(_k);
    gsl_vector_memcpy(t1,dfh);
    gsl_vector_sub(t1,df);
    gsl_vector_memcpy(t2,dfh);
    gsl_vector_add(t2,df);
   
    //norm(dh-dy)/norm(dh+dy);
    if(vnorm(t1)/vnorm(t2) > 1e-5) {
        vprint(x);
        printf("mean gradient doesn't match %e", vnorm(t1)/vnorm(t2));
        for(uint32_t k=0; k<_k;++k)
          printf(" (%.10f,%.10f)", gsl_vector_get(df,k), gsl_vector_get(dfh,k));
        printf("\n");
    }
    gsl_vector_free(t1);
    gsl_vector_free(t2);

    #endif

    //gsl_multimin_fdfminimizer_set(s, &nu_obj, x, 0.01, 1e-3);
    gsl_multimin_fdfminimizer_set(s, &mu_obj, x, 0.01, 0.1);

    do
    {
        iter++;
        f_old = s->f;
        status = gsl_multimin_fdfminimizer_iterate(s);
        converged = fabs((f_old - s->f) / f_old);
        //printf("f(mu) = %5.17e ; conv = %5.17e\n", s->f, converged);
        if (status) break;
        status = gsl_multimin_test_gradient(s->gradient, _cg_convergence);
    }
    while ((status == GSL_CONTINUE) && (iter < _cg_max_iter));
    // while ((converged > PARAMS.cg_convergence) &&
    // ((PARAMS.cg_max_iter < 0) || (iter < PARAMS.cg_max_iter)));
    if (iter == _cg_max_iter) {
        printf("warning: cg didn't converge (mu) after %d iterations \n\t(status: %s)\n", _cg_max_iter, gsl_strerror(status));
        printf("x\n");
        vprint(s->x);
        printf("nu^2\n"); 
        const double ** const vd = _vcurr.const_data();
        for (uint32_t k=0; k<_k; ++k)
            printf("%f ", vd[n][k]); 
        printf("\n"); 
        
        const double * const pd = phi.const_data(); 
        printf("phi\n"); 
        for (uint32_t k=0; k<_k; ++k)
            printf("%f ", pd[k]); 
        printf("\n"); 
        //exit(-1); 
    }

    #if 0
    if(iter < 2 && status != GSL_SUCCESS && status != GSL_CONTINUE) {
        printf("\nmean opt. iter %d, status %d, id %d\n", iter, status, n);
        printf("mean begin\n");
        vprint(xt);
        printf("mean\n");
        vprint(s->x);
        const double *const od = other.const_data();
        printf("od\n");
        for(int k=0; k<_k; ++k)
            printf("%f ", od[k]);
        printf("\n");
        const double *const pd = phi.const_data();
        printf("phi\n");
        for(int k=0; k<_k; ++k)
            printf("%f ", pd[k]);
        printf("\n");
        gsl_vector* df = gsl_vector_calloc(_k);
        wrap_df_var(x, (void *)&b, df);
        printf("df\n");
        vprint(df);
        exit(-1);
    }
    #endif

    //printf("x mean after optimization:\n");
    //vprint(s->x, _k);

    Array mean(_k);
    for (uint32_t k = 0; k < _k; ++k)
        mean[k] = gsl_vector_get(s->x, k);
    _mnext.add_slice(n, mean);

    #if 0
    if(wrap_f_mean(xt, (void *)&b) < wrap_f_mean(s->x, (void *)&b)) {
      printf("iter %d, status %d\n", iter, status);
      printf("f_mean before: %f\t", wrap_f_mean(xt, (void *)&b));
      printf("f_mean after: %f\n", wrap_f_mean(s->x, (void *)&b));
    }
    #endif

    gsl_multimin_fdfminimizer_free(s);
    gsl_vector_free(x);

    return 0.; 

}

inline float
NormMatrix::update_var_next(uint32_t n, const Array &other)
{

    int max_iters = -1; 

    // update the K-vector mean of a particular user/item
    // involves 
    // 2) priors (c,d)
    // 3) nu
    // keep pointers to all in a bundle 

    gsl_multimin_fdfminimizer * s;
    gsl_multimin_function_fdf nu_obj;

    int iter = 0, i, j;
    int status;
    double f_old, converged;

    bundle b;
    b.NormMatrixObj = this; 
    b.id = n; 
    b.other = &other; 

    nu_obj.f = &NormMatrix::wrap_f_var;
    nu_obj.df = &NormMatrix::wrap_df_var;
    nu_obj.fdf = &NormMatrix::wrap_fdf_var;
    nu_obj.n = _k;
    nu_obj.params = (void *)&b;

    const gsl_multimin_fdfminimizer_type * T;
    // bfgs2 seems (empirically) as good and faster
    //T = gsl_multimin_fdfminimizer_conjugate_fr;
    T = gsl_multimin_fdfminimizer_vector_bfgs2;
    s = gsl_multimin_fdfminimizer_alloc(T, this->_k);

    gsl_vector* x = gsl_vector_calloc(_k);
    const double **vd = _vcurr.const_data();
    for (uint32_t k = 0; k < _k; ++k) gsl_vector_set(x, k, log(vd[n][k]));

    gsl_multimin_fdfminimizer_set(s, &nu_obj, x, 0.01, 1e-3);
    //gsl_multimin_fdfminimizer_set(s, &nu_obj, x, 0.01, 0.1);

    #if 1
    gsl_vector* xt = gsl_vector_calloc(_k);
    gsl_vector_memcpy(xt, x);
    #endif

    #if CHECK_GRADIENTS
    gsl_vector *xg = gsl_vector_alloc(_k);
    gsl_vector *dfh  = gsl_vector_alloc(_k);
    double eps = 1e-5;
    double f1, f2;
    for(uint32_t k=0; k < _k; ++k) {
      gsl_vector_memcpy(xg, x);
      gsl_vector_set(xg, k, gsl_vector_get(xg, k) + eps/2);
      f1 = wrap_f_var(xg, (void *)&b);
      gsl_vector_set(xg, k, gsl_vector_get(xg, k) - eps);
      f2 = wrap_f_var(xg, (void *)&b);
      gsl_vector_set(dfh, k, (f1-f2)/eps);
    }
    gsl_vector* df = gsl_vector_calloc(_k);
    wrap_df_var(x, (void *)&b, df);
    gsl_vector *t1 = gsl_vector_calloc(_k); 
    gsl_vector *t2 = gsl_vector_calloc(_k); 
    gsl_vector_memcpy(t1,dfh); 
    gsl_vector_sub(t1,df); 
    gsl_vector_memcpy(t2,dfh); 
    gsl_vector_add(t2,df); 
    
    //norm(dh-dy)/norm(dh+dy);
    if(vnorm(t1)/vnorm(t2) > 1e-4) {
        vprint(x);
        printf("var gradient doesn't match %e", vnorm(t1)/vnorm(t2)); 
        for(uint32_t k=0; k<_k;++k) 
          printf(" (%.10f,%.10f)", gsl_vector_get(df,k), gsl_vector_get(dfh,k));
        printf("\n"); 
    }
    gsl_vector_free(t1); 
    gsl_vector_free(t2); 
    #endif
    ///

    do
    {
        iter++;
        f_old = s->f;
        status = gsl_multimin_fdfminimizer_iterate(s);
        converged = fabs((f_old - s->f) / f_old);
        //printf("f(nu) = %5.17e ; conv = %5.17e\n", s->f, converged);
        //vprint(s->x, _k);
        if (status) break;
        status = gsl_multimin_test_gradient(s->gradient, _cg_convergence);
    }
    while ((status == GSL_CONTINUE) && (iter < _cg_max_iter) && (iter < max_iters || max_iters == -1));
    // while ((converged > PARAMS.cg_convergence) &&
    // ((PARAMS.cg_max_iter < 0) || (iter < PARAMS.cg_max_iter)));
    if (iter == _cg_max_iter) {
        printf("warning: cg didn't converge (nu) %d\n", _cg_max_iter);
        printf("x\n");
        vprint(s->x);
        printf("mu\n");
        const double ** const md = _mcurr.const_data();
        for (uint32_t k=0; k<_k; ++k)
            printf("%f ", md[n][k]); 
        printf("\n"); 
        
        printf("\n"); 
        exit(-1); 
    } 
    //if (iter < 2 && status != GSL_SUCCESS && status != GSL_CONTINUE) { 
    //    printf("problem at iteration %d: gsl return code %d\n", iter, status);
    //}

    Array var(_k); 
    for (uint32_t k = 0; k < _k; ++k)
        var[k] = exp(gsl_vector_get(s->x, k));
    _vnext.add_slice(n, var);

    #if 1
    if(wrap_f_var(xt, (void *)&b) < wrap_f_var(s->x, (void *)&b)) {
      printf("f_var before: %.10f\t", wrap_f_var(xt, (void *)&b));
      printf("f_var after: %.10f\n", wrap_f_var(s->x, (void *)&b));
    }
    #endif

    gsl_multimin_fdfminimizer_free(s);
    gsl_vector_free(x);

    return 0.; 

}

#endif 

