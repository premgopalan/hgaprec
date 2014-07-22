//blaslapack.h
//version 0.95
//------------

//Purpose:	



#ifndef UNIVIE_NMFLIB_BLASLAPACK_H
#define UNIVIE_NMFLIB_BLASLAPACK_H


#ifdef __cplusplus
extern "C" {
#endif


//declaration of the BLAS/LAPACK routines to be used
//Fortran functions are wrapped in a new c-function
//--------------------------------------------------

// dtrsm - Blas level routine
//----------------------------------
// used to solve a linear system A * x = B with A beeing triangular
static inline double dtrsm(char side, char uplo, char transa, char diag, int m, int n, double alpha, double * a, int lda, double * b, int ldb) {
  extern double dtrsm_(char* side, char* uplo, char* transa, char* diag, int* m, int* n, double* alpha, double * a, int* lda, double * b, int* ldb);
  return dtrsm_(&side, &uplo, &transa, &diag, &m, &n, &alpha, a, &lda, b, &ldb);
}

// dlange - Lapack auxiliary routine
//----------------------------------
// used to calculate the absolute maximum of a matrix and to calculate the frobenius norm
static inline double dlange(char norm, int m, int n, double* a, int lda, double * work) {
  extern double dlange_(char* norm, int* m, int* n, double* a, int* lda, double* work);
  return dlange_(&norm, &m, &n, a, &lda, work);
}

// slange - LAPACK auxiliary routine
//-----------------------------------
// used to calculate the absolut maximum of a matrix and to calculate the frobenius norm
static inline float slange(char norm, int m, int n, float* a, int lda, float * work) {
  extern float slange_(char* norm, int* m, int* n, float* a, int* lda, float* work);
  return slange_(&norm, &m, &n, a, &lda, work);
}


// dlamch - Lapack auxiliary routine
//----------------------------------
// used to calculate machine precision epsilon
static inline double dlamch(char cmach) {
  extern double dlamch_(char* cmach);
  return dlamch_(&cmach);
}

// slamch - Lapack auxiliary routine - single precision !!!
//---------------------------------------------------------
// used to calculate machine precision epsilon
static inline float slamch(char cmach) {
  extern float slamch_(char* cmach);
  return slamch_(&cmach);
}


// dgemm - BLAS routine
//---------------------
// used to calculate general matrix-matrix multiplications of the form
// C = alpha * A * B + beta * C
static inline int dgemm(char transa, char transb, int m, int n, int k, double alpha, double* a, int lda, double* b, int ldb, double beta, double* c, int ldc) {
  extern int dgemm_(char* transa, char* transb, int* m, int* n, int* k, double* alpha, double* a, int* lda, double* b, int* ldb, double* beta, double* c, int* ldc);
  return dgemm_(&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
}

// sgemm - BLAS routine - single precision !!!
//---------------------------------------------
// used to calculate general matrix-matrix multiplications of the form
// C = alpha * A * B + beta * C
static inline int sgemm(char transa, char transb, int m, int n, int k, float alpha, float* a,  int lda, float* b, int ldb, float beta, float* c, int ldc) {
  extern int sgemm_(char* transa, char* transb, int* m, int* n, int* k, float* alpha, float* a, int* lda, float* b, int* ldb, float* beta, float* c, int* ldc);
  return sgemm_(&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
}



// dgeqp3 - Lapack routine
//------------------------
// used to calculate a QR factorization with column pivoting
static inline int dgeqp3(int M, int N, double *A, int LDA, int * JPVT, double * TAU, double * WORK, int LWORK, int *INFO ) {
  extern int dgeqp3_(int *M, int *N, double *A, int *LDA, int *JPVT, double *TAU, double *WORK, int * LWORK, int *INFO);
  return dgeqp3_(&M, &N, A, &LDA, JPVT, TAU, WORK, &LWORK, INFO);
}


// dormqr - Lapack routine
//-----------------------
// used to calculate a matrix-matrix-multiplication of a factor matrix q (result of dgeqp3) and another matrix
static inline int dormqr(char SIDE, char TRANS, int M, int N, int K, double * A, int LDA, double * TAU, double * C, int LDC, double * WORK, int LWORK, int *INFO ) {
  extern int dormqr_(char * SIDE, char * TRANS, int * M, int * N, int * K, double * A, int * LDA, double * TAU, double * C, int * LDC, double * WORK, int * LWORK, int * INFO);
  return dormqr_(&SIDE, &TRANS, &M, &N, &K, A, &LDA, TAU, C, &LDC, WORK, &LWORK, INFO);
}



// dtrtrs - Lapack routine
//------------------------
// used to solve a triangular system of the form a * x = b
static inline int dtrtrs(char UPLO, char TRANS, char DIAG, int N, int NRHS, double * A, int LDA, double * B, int LDB, int *INFO) {
  extern int dtrtrs_(char *UPLO, char *TRANS, char *DIAG, int *N, int *NRHS, double * A, int *LDA, double * B, int *LDB, int *INFO);
  return dtrtrs_(&UPLO, &TRANS, &DIAG, &N, &NRHS, A, &LDA, B, &LDB, INFO);
}


// dgesv - Lapack routine
//-----------------------
// used to calculate a LU-Factorisation of a matrix
static inline int dgesv(int N, int NRHS, double * A, int LDA, int * IPIV, double * B, int LDB, int * INFO) {
  extern int dgesv_(int * N, int * NRHS, double * A, int * LDA, int * IPIV, double * B, int * LDB, int * INFO);
  return dgesv_(&N, &NRHS, A, &LDA, IPIV, B, &LDB, INFO);
}


// dlacpy - Lapack routine
//-------------------------
// used to copy a matrix A or part of it to a matrix B
static inline int dlacpy(char UPLO, int M, int N, double * A, int LDA, double * B, int LDB) {
  extern int dlacpy_(char * UPLO, int * M, int * N, double * A, int * LDA, double * B, int * LDB);
  return dlacpy_(&UPLO, &M, &N, A, &LDA, B, &LDB);
}

//slacpy - LAPACK routine
//-----------------------
// used to copy a matrix A or part of it to a matrix B
static inline int slacpy(char UPLO, int M, int N, float * A, int LDA, float * B, int LDB) {
  extern int slacpy_(char * UPLO, int* M, int* N, float* A, int* LDA, float* B, int* LDB);
  return slacpy_(&UPLO, &M, &N, A, &LDA, B, &LDB);
}



// dscal - Blas routine
//---------------------
// used for multiplying a vector by a constant
static inline int dscal(int N, double DA, double * DX, int INCX) {
  extern int dscal_(int *N, double * DA, double * DX, int * INCX);
  return dscal_(&N, &DA, DX, &INCX);
}


// dcopy - Blas routine
//---------------------
// used for copying a vector
static inline int dcopy(int N, double * DX, int INCX, double * DY, int INCY) {
  extern int dcopy_(int * N, double * DX, int * INCX, double *DY, int * INCY);
  return dcopy_(&N, DX, &INCX, DY, &INCY);
}

// daxpy - Blas level 1 routine
//-----------------------------
// used for calculating dy = da * dx + dy
static inline int daxpy(int N, double DA, double * DX, int INCX, double * DY, int INCY ) {
  extern int daxpy_(int * N, double * DA, double * DX, int * INCX, double * DY, int *INCY);
  return daxpy_(& N, & DA, DX, &INCX, DY, &INCY);
}


// saxpy - BLAS level 1 routine
//-----------------------------
// used for calculating dy = da * dx + dy
static inline int saxpy(int N, float DA, float* DX, int INCX, float* DY, int INCY) {
  extern int saxpy_(int * N, float * DA, float* DX, int* INCX, float* DY, int* INCY);
  return saxpy_(&N, &DA, DX, &INCX, DY, &INCY);  
}


// dorgqr - Lapack routine
//------------------------
// used for calculating the explicit matrix Q out of a matrix factorized and overwritten by dgeqp3
static inline int dorgqr(int M, int N, int K, double * A, int LDA, double * TAU, double * WORK, int LWORK, int * INFO) {
  extern int dorgqr_(int * M, int * N, int * K, double * A, int * LDA, double * TAU, double * WORK, int * LWORK, int * INFO);
  return dorgqr_(&M, &N, &K, A, &LDA, TAU, WORK, &LWORK, INFO);
}



// PROPACK SVD routine
static inline int dlansvd(char jobu, char jobv, int m, int n, int k, int kmax, int (*aprod)(char*, int*, int*, double*, double*, double*, int*), double* U, int ldu, double * Sigma, double* bnd, double* V, int ldv, double tolin, double * work, int lwork, int* iwork, int liwork, double * doption, int * ioption, int * info, double * dparm, int * iparm) {
    extern int dlansvd_(char * jobu, char * jobv, int* m, int* n, int * k, int* kmax, int (*aprod)(char*, int*, int*, double*, double*, double*, int*),double * U, int * ldu, double * Sigma, double * bnd, double * V, int * ldv,double * tolin, double * work, int * lwork, int *iwork, int * liwork,double * doption,int * ioption,int * info,double * dparm, int * iparm);
    return dlansvd_(&jobu, &jobv, &m, &n, &k, &kmax, aprod, U, &ldu, Sigma, bnd, V, &ldv, &tolin, work, &lwork, iwork, &liwork, doption, ioption, info, dparm, iparm); 
}

// dgemv - Lapack routine
//-----------------------
// used to calculate one of the matrix-vector-operations y := alpha * A * x + beta * y _or_ y := alpha * A' * x + beta * y
static inline int dgemv(char TRANS, int M, int N, double ALPHA, double * A, int LDA, double * X, int INCX, double BETA, double * Y, int INCY) {
  extern int dgemv_(char * TRANS, int * M, int * N, double * ALPHA, double * A, int * LDA, double * X, int * INCX, double * BETA, double * Y, int * INCY);
  return dgemv_(&TRANS, &M, &N, &ALPHA, A, &LDA, X, &INCX, &BETA, Y, &INCY);
}


// dlaset - Lapack routine
//------------------------
// used to initialize a m by n matrix to a fixed scalar
static inline int dlaset(char uplo, int m, int n, double alpha, double beta, double * A, int lda) {
  extern int dlaset_(char * uplo, int * m, int * n, double * alpha, double * beta, double * A, int * lda);
  return dlaset_(&uplo, &m, &n, &alpha, &beta, A, &lda);
}


// dsaupd  ARPACK reverse communication interface routine.
//--------------------------------------------------------
static inline int dsaupd(int * IDO, unsigned char BMAT, int N, unsigned char * WHICH, int NEV, double TOL, double * RESID, int NCV, double * V, int LDV, int * IPARAM, int * IPNTR, double * WORKD, double * WORKL, int LWORKL, int * INFO) {
  extern int dsaupd_(int * IDO, unsigned char * BMAT, int * N, unsigned char * WHICH, int * NEV, double * TOL, double * RESID, int * NCV, double * V, int * LDV, int * IPARAM, int * IPNTR, double * WORKD, double * WORKL, int * LWORKL, int * INFO);
  return dsaupd_(IDO, &BMAT, &N, WHICH, &NEV, &TOL, RESID, &NCV, V, &LDV, IPARAM, IPNTR, WORKD, WORKL, &LWORKL, INFO );
}

// dseupd  ARPACK routine that returns Ritz values and (optionally) Ritz vectors
//----------------------------------------------------------------
static inline int dseupd(int RVEC, char HOWMNY, int * SELECT, double * D, double * Z, int LDZ, double SIGMA, unsigned char BMAT, int N, unsigned char * WHICH, int NEV, double TOL, double * RESID, int NCV, double * V, int LDV, int * IPARAM, int * IPNTR, double * WORKD, double *  WORKL, int LWORKL, int * INFO ) {
  extern int dseupd_(int * RVEC, char * HOWMNY, int * SELECT, double * D, double * Z, int * LDZ, double * SIGMA, unsigned char *BMAT, int *N, unsigned char * WHICH, int * NEV, double * TOL, double * RESID, int * NCV, double * V, int *LDV, int * IPARAM, int * IPNTR, double * WORKD, double *  WORKL, int * LWORKL, int * INFO);
  return dseupd_(&RVEC, &HOWMNY, SELECT, D, Z, &LDZ, &SIGMA, &BMAT, &N, WHICH, &NEV, &TOL, RESID, &NCV, V, &LDV, IPARAM, IPNTR, WORKD, WORKL, &LWORKL, INFO);
}  


// dnrm2 - Blas routine
//---------------------
// used for calculating the euclidean norm of a vector x
static inline double dnrm2(int N, double * X, int INCX) {
  extern double dnrm2_(int *N, double * X, int * INCX);
  return dnrm2_(&N, X, &INCX);
}

// dswap BLAS level 1 routine to swap to vectors
//-----------------------------------------------
static inline int dswap(int N, double * DX, int INCX, double * DY, int INCY) {
  extern int dswap_(int * N, double * DX, int * INCX, double * DY, int *INCY);
  return dswap_(&N, DX, &INCX, DY, &INCY);
} 


#ifdef __cplusplus
}
#endif


#endif

