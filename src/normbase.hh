#ifndef NORMBASE_HH
#define NORMBASE_HH

// #include "env.hh"

#include <gsl/gsl_multimin.h>
#include <gsl/gsl_rng.h>
#include <cmath> 

template <class T>
class NormBase { 

public:
  NormBase(string name = ""): _name(name) { }
  virtual ~NormBase() { }
  virtual const T &expected_v() const = 0;
  virtual uint32_t n() const = 0;
  virtual uint32_t k() const = 0;
  virtual double compute_elbo_term_helper() const = 0;
  virtual void save_state(const IDMap &m) const = 0;
  virtual void load()  = 0;
  string name() const { return _name; }
  double compute_elbo_term() const;
private:
  string _name;
};


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
    // _Ev(n,k),
    _r(r) { } 
  virtual ~NormMatrix() {} 

  uint32_t _cg_max_iter = 500; // TODO 
  float _cg_convergence  = 1e-5; // TODO
  void pvec(const gsl_vector *x, uint32_t n);

  uint32_t n() const { return _n;}
  uint32_t k() const { return _k;}

  const Matrix &mean_curr() const         { return _mcurr; }
  const Matrix &var_curr() const          { return _vcurr; }
  const Matrix &mean_next() const         { return _mnext; }
  const Matrix &var_next() const          { return _vnext; }
  const Matrix &expected_v() const         { return _mcurr;    }

  Matrix &mean_curr()       { return _mcurr; }
  Matrix &var_curr()        { return _vcurr; }
  Matrix &mean_next()       { return _mnext; }
  Matrix &var_next()        { return _vnext; }
  Matrix &expected_v()       { return _mcurr;    }

  const double mprior() const { return _mprior; }
  const double vprior() const { return _vprior; }
  void set_to_prior();
  void initialize(); 
  void initialize2(double v);
  void update_mean_next(uint32_t n, const Array &phi); 
  void update_var_next(uint32_t n); 
  // compute_elbo_term_helper
  // TODO: look these over 
  double compute_elbo_term_helper() const { return -1.0;}
  void save_state(const IDMap &m) const {}
  void load() {}

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
  //Matrix _Ev;         // expected weights under variational
        		      // distribution
};

typedef struct bundle {
    NormMatrix *NormMatrixObj; 
    const Array *phi; 
    uint32_t id; 
} bundle;

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
NormMatrix::f_mean(const gsl_vector *x, void * params)
{
    NormMatrix& obj = (*((bundle *) params)->NormMatrixObj);
    uint32_t id = ((bundle *) params)->id;

    double f = .0; 

    //f = (x_k - c)^2 / 2*d^2 - exp(x + nu^2/2) + \phi_mn * x

    gsl_vector *tmp = gsl_vector_alloc(obj._k);

    // (x_k - c)^2 / 2*d^2
    gsl_vector_memcpy(tmp, x);
    gsl_vector_add_constant(tmp, -obj._mprior);
    gsl_vector_mul(tmp, tmp);
    gsl_vector_scale(tmp, 1/(2*(obj._vprior*obj._vprior))); 

    const double ** const vd = obj._vcurr.const_data();
    const double * const pd = (*((bundle *) params)->phi).const_data(); 

    double x_k; 
    for(uint32_t k=0; k<obj._k; ++k) { 

        f -= gsl_vector_get(tmp, k); 
        
        x_k = gsl_vector_get(x, k); 
        f -= exp(x_k + (vd[id][k])/2);

        f += pd[k] * x_k; 

    }

    gsl_vector_free(tmp); 
    return -f; // maximize
}


inline void 
NormMatrix::df_mean(const gsl_vector * x, void * params, gsl_vector * df)
{
    NormMatrix& obj = (*((bundle *) params)->NormMatrixObj);
    uint32_t id = ((bundle *) params)->id;

    double f = .0; 

    //f = x_k / d^2 - exp(x + nu^2/2) + \phi_mn

    //gsl_vector *df = gsl_vector_alloc(obj._k);

    // (x_k / d^2) 
    gsl_vector_memcpy(df, x);
    gsl_vector_scale(df, 1/(obj._vprior*obj._vprior)); 

    const double ** const vd = obj._vcurr.const_data();
    const double * const pd = (*((bundle *) params)->phi).const_data(); 

    gsl_vector_add_constant(df, obj._mprior/(obj._vprior*obj._vprior));

    double x_k, nu_k; 
    for(uint32_t k=0; k<obj._k; ++k) { 

        // f -= gsl_vector_get(tmp, k); 
        
        x_k = gsl_vector_get(x, k); 
        gsl_vector_set(df, k, exp(x_k + (vd[id][k])/2) + pd[k]);
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
NormMatrix::set_to_prior()
{
  _mnext.set_elements(_mprior);
  _vnext.set_elements(_vprior);
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
      ad[i][k] = _mprior + 0.01 * gsl_rng_uniform(*_r);
      //gsl_ran_gaussian_ziggurat(*_r, _vprior)

  for (uint32_t k = 0; k < _k; ++k)
    bd[0][k] = _vprior + 0.1 * gsl_rng_uniform(*_r);
  
  for (uint32_t i = 0; i < _n; ++i)
    for (uint32_t k = 0; k < _k; ++k)
      bd[i][k] = bd[0][k];
  set_to_prior();
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

inline void 
NormMatrix::pvec(const gsl_vector *x, uint32_t n) { 
    //printf("vector: "); 
    for (uint32_t k=0; k<n; ++k)
        printf("%f ", gsl_vector_get(x, k));
    printf("\n"); 
}

inline double 
NormMatrix::f_var(const gsl_vector *x, void * params)
{
    NormMatrix& obj = (*((bundle *) params)->NormMatrixObj);
    uint32_t id = ((bundle *) params)->id;

    double f = .0; 

    //f = - d^2 nu_k^2 + 1/2(log \nu^2) - exp(x + nu^2/2)

    // - d^2 nu_k^2
    //gsl_vector_memcpy(tmp, x);
    //gsl_vector_scale(tmp, (obj._vprior*obj._vprior)); 

    const double ** const md = obj._mcurr.const_data();

    double x_k; 
    for(uint32_t k=0; k<obj._k; ++k) { 

        // -d^2 \nu_k^2 + 1/2(log \nu^2) - exp(x + nu^2/2)
        x_k = exp(gsl_vector_get(x, k)); 
        f -= x_k*(obj._vprior*obj._vprior);
        f += 0.5 * log(x_k); 
        if(isnan(f)) 
            printf("x_k: %f; log(x_k): %f \n", x_k, log(x_k)); 

        f -= exp(md[id][k] + x_k/2);
    }

    //printf("f: %f\n", f); 

    return -f; // maximize
}


inline void 
NormMatrix::df_var(const gsl_vector * x, void * params, gsl_vector * df)
{
    NormMatrix& obj = (*((bundle *) params)->NormMatrixObj);
    uint32_t id = ((bundle *) params)->id;

   //df = -d^2 + 1/(2nu_k^2) - 1/2exp(x + nu^2/2) 

    const double ** const md = obj._mcurr.const_data();

    double x_k, nu_k; 
    for(uint32_t k=0; k<obj._k; ++k) { 
        x_k = exp(gsl_vector_get(x, k)); 
        //gsl_vector_set(df, k, -obj._vprior*obj._vprior + 0.5/x_k - 0.5*exp(md[id][k] + x_k/2));
        gsl_vector_set(df, k, -x_k*obj._vprior*obj._vprior + 0.5 - x_k*0.5*exp(md[id][k] + x_k/2));
    }

    gsl_vector_scale(df, -1.0); // maximize
}

inline void 
NormMatrix::fdf_var(const gsl_vector * x, void * params, double * f, gsl_vector * df)
{
    *f = f_var(x, params);
    df_var(x, params, df);
}


inline void
NormMatrix::update_mean_next(uint32_t n, const Array &phi)
{

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
    b.id = n; 

    //printf("optimizing mean %d/%d\n", n, _n); 
    
    mu_obj.f = &NormMatrix::wrap_f_mean;
    mu_obj.df = &NormMatrix::wrap_df_mean;
    mu_obj.fdf = &NormMatrix::wrap_fdf_mean;
    mu_obj.n = _k;
    mu_obj.params = (void *)&b;

    // starting value
    const gsl_multimin_fdfminimizer_type * T;
    // T = gsl_multimin_fdfminimizer_vector_bfgs;
    T = gsl_multimin_fdfminimizer_conjugate_fr;
    // T = gsl_multimin_fdfminimizer_steepest_descent;
    s = gsl_multimin_fdfminimizer_alloc(T, _k);

    gsl_vector* x = gsl_vector_calloc(_k);

    //printf("x before optimization:\n"); 
    //pvec(x, _k); 

    gsl_multimin_fdfminimizer_set(s, &mu_obj, x, 0.01, 1e-3);

    do
    {
        iter++;
        f_old = s->f;
        status = gsl_multimin_fdfminimizer_iterate(s);
        converged = fabs((f_old - s->f) / f_old);
        //printf("f(mu) = %5.17e ; conv = %5.17e\n", s->f, converged);
        if (status) break;
        status = gsl_multimin_test_gradient(s->gradient, this->_cg_convergence);
    }
    while ((status == GSL_CONTINUE) && (iter < this->_cg_max_iter));
    // while ((converged > PARAMS.cg_convergence) &&
    // ((PARAMS.cg_max_iter < 0) || (iter < PARAMS.cg_max_iter)));
    if (iter == this->_cg_max_iter) { 
        printf("warning: cg didn't converge (mu) %d\n", this->_cg_max_iter);
        printf("x\n"); 
        pvec(s->x, _k); 
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
        exit(-1); 
    }

    //printf("x mean after optimization:\n"); 
    //pvec(s->x, _k); 

    Array mean(_k); 
    for (uint32_t k = 0; k < _k; ++k)
        mean[k] = gsl_vector_get(s->x, k);
    _mnext.add_slice(n, mean);

    gsl_multimin_fdfminimizer_free(s);
    gsl_vector_free(x);

}

inline void
NormMatrix::update_var_next(uint32_t n)
{

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

    //printf("optimizing var %d/%d\n", n, _n); 
    
    nu_obj.f = &NormMatrix::wrap_f_var;
    nu_obj.df = &NormMatrix::wrap_df_var;
    nu_obj.fdf = &NormMatrix::wrap_fdf_var;
    nu_obj.n = _k;
    nu_obj.params = (void *)&b;

    // starting value
    const gsl_multimin_fdfminimizer_type * T;
    // T = gsl_multimin_fdfminimizer_vector_bfgs;
    T = gsl_multimin_fdfminimizer_conjugate_fr;
    // T = gsl_multimin_fdfminimizer_steepest_descent;
    s = gsl_multimin_fdfminimizer_alloc(T, this->_k);

    double **bd = _vcurr.data();

    gsl_vector* x = gsl_vector_calloc(_k);
    for (uint32_t k=0; k<_k; ++k)
        gsl_vector_set(x, k, bd[n][k]);

    //printf("x before optimization:\n"); 
    //pvec(x, _k); 
    
    gsl_multimin_fdfminimizer_set(s, &nu_obj, x, 0.01, 1e-3);

    do
    {
        iter++;
        f_old = s->f;
        status = gsl_multimin_fdfminimizer_iterate(s);
        converged = fabs((f_old - s->f) / f_old);
        //printf("f(nu) = %5.17e ; conv = %5.17e\n", s->f, converged);
        //pvec(s->x, _k); 
        if (status) break;
        status = gsl_multimin_test_gradient(s->gradient, this->_cg_convergence);
    }
    while ((status == GSL_CONTINUE) && (iter < this->_cg_max_iter));
    // while ((converged > PARAMS.cg_convergence) &&
    // ((PARAMS.cg_max_iter < 0) || (iter < PARAMS.cg_max_iter)));
    if (iter == this->_cg_max_iter) { 
        printf("warning: cg didn't converge (nu) %d\n", this->_cg_max_iter);
        printf("x\n"); 
        pvec(s->x, _k); 
        printf("mu\n"); 
        const double ** const md = _mcurr.const_data();
        for (uint32_t k=0; k<_k; ++k)
            printf("%f ", md[n][k]); 
        printf("\n"); 
        
        printf("\n"); 
        exit(-1); 
    } 


    //printf("x after optimization:\n"); 
    //pvec(s->x, _k); 

    Array var(_k); 
    for (uint32_t k = 0; k < _k; ++k)
        var[k] = gsl_vector_get(s->x, k);
    _vnext.add_slice(n, var);

    gsl_multimin_fdfminimizer_free(s);
    gsl_vector_free(x);

}


#endif 

