#ifndef UNIVIE_LSI_LIB_CALCULATESVD_H
#define UNIVIE_LSI_LIB_CALCULATESVD_H

/**
 * Calculates the Singular value decomposition for a matrix a
 *
 * m: in, m dimension of a
 * n: in, n dimension of a, n >= m
 * a: in, matrix of size lda * n
 * lda: in, leading dimension of a, lda >= m
 * u: out, matrix of size ldu * m
 * ldu: in, leading dimension of u, ldu >= m
 * vt: out, matrix of size ldvt * n
 * ldvt: in, leading dimension of vt, ldvt >= n
 * s: out, vector of length m
 */
int calculateSVD(double * a, double * u, double * s, double * v, int m, int n, int k, int maxiter, double tol, int ncv, int order);

#endif
