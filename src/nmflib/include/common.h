#ifndef UNIVIE_NMFLIB_COMMON_H
#define UNIVIE_NMFLIB_COMMON_H


//Definition of different options
//-------------------------------------------------------------------

//Configuration options
//---------------------


//Definition of a zero_threshold - smaller values are set to zero
//if not defined, machine epsilon calculated by DLAMCH is used

#define ZERO_THRESHOLD 0.0


//DEBUG_LEVEL defines the level of verbosity
// 2	output everything
// 1	output intermediate results
// 0	output only the final result
#define DEBUG_LEVEL 0

//Turn off error checking by not defining ERROR_CHECKING for possible speedup
#define ERROR_CHECKING

//Definition of different profiles deciding whether or not to print timings
/**
#define PROFILE_LOAD_MATRIX			//timing for loadMatrix()
#define PROFILE_STORE_MATRIX			//timing for storeMatrix()
#define PROFILE_CALCULATE_NORM			//timing for calculateNorm()
#define PROFILE_CALCULATE_MAXCHANGE		//timing for calculateMaxchange()
#define PROFILE_NMF_MU				//timing for nmf_mu()
#define PROFILE_NMF_ALS				//timing for nmf_als()
#define PROFILE_NMF_NEALS			//timing for nmf_neals()
#define PROFILE_NMF_ALSPG			//timing for nmf_alspg()
#define PROFILE_NMF_PG				//timing for nmf_pg()
#define PROFILE_PG_SUBPROB			//timing for pg_subprob()
#define PROFILE_CHECK_ARGUMENTS			//timing for checkArguments()
#define PROFILE_CHECK_MATRICES			//timing for checkMatrices()
#define PROFILE_GENERATE_MATRIX			//timing for generateMatrix()
#define PROFILE_NMF_DRIVER			//timing for nmfDriver()

#define PROFILE_MATLAB_COMPARISON		//timing for a matlab comparison = T(nmfDriver) without loading/storing matrices
*/







//Definition of some datatypes used in nmf_driver and its subroutines
//-------------------------------------------------------------------


// struct idx_double
//------------------
// indexed double structure for storing an index with every double value
// used to resort factorisation results using the qsort() function of the c standard library
typedef struct idx_double {
  double val;
  int idx;
} idx_double;


// enum-types for more intuitive calls
//------------------------------------

// alg_t
//------
// used to define the algorithm to use
// mu		multiplicative update
// als		alternating least squares
// neals 	normal equation alternating least squares
// alspg	alternating least squares using a projected gradient approach
// pg		direct projected gradient approach
typedef enum alg_t {mu, als, neals, alspg, pg} alg_t;	

// init_t
//-------
// used to define the method for matrix-initialisation
// ran		random initialisation
// nndsvd	non negative double singular value decomposition
typedef enum init_t {ran, nndsvd} init_t;	



// options_t
//----------
// used to bundle all possible options for calling nmfDriver in one structure
typedef struct options_t {
  int rep;					//Number of repetitions with new starting matrices
  init_t init;					//Method to use for initialising the starting matrices
  int min_init;					//minimal value for random initialisation
  int max_init;					//maximal value for random initialisaton
  const char* w_out;				//Filename to write final matrix w to
  const char* h_out;				//Filename to write final matrix h to
  double TolX;					//tolerance value for convergence check of maxChange
  double TolFun;				//tolerance value for convergence check of dnorm
  int nndsvd_maxiter;				//maximum iterations for SVD in ddsvd initialisation || -1 -> maxiter = 4 * (m+n)
  int nndsvd_blocksize;				//blocksize for SVD in ddsvd initialisation
  double nndsvd_tol;				//tolerance for SVD in ddsvd initialisation
  int nndsvd_ncv;				//largest number of basis vectors in the Implicitly Restarted Arnoldi Process
} options_t;


#ifdef __cplusplus
extern "C" {
#endif


//defining an inline SWAP-Function
/** swap - swaps two pointers to double
 *
 * Purpose:
 *		Used for matrix "assignments" when no real copy is needed
 *
 * Arguments:
 *
 * a		in & out, 	pointer to pointer to double to be swapped
 * b		in & out, 	pointer to pointer to double to be swapped
 */

static inline void swap(double ** a, double ** b) {
  double* temp = *a;
  *a = *b;
  *b = temp;
}

static inline void swap_singleprec(float ** a, float ** b) {
  float* temp = *a;
  *a = *b;
  *b = temp;
}

 #ifdef __cplusplus
 }
 #endif





#endif
