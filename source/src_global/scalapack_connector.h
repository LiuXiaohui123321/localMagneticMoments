#ifndef SCALAPACK_CONNECTOR_H
#define SCALAPACK_CONNECTOR_H

extern "C"
{
	void blacs_gridinit_( int *ictxt, const char *order, const int *nprow, const int *npcol );
	void blacs_gridinfo_( const int *ictxt, int *nprow, int *npcol, int *myprow, int *mypcol );
	int numroc_( const int *n, const int *nb, const int *iproc, const int *srcproc, const int *nprocs );
	void descinit_( 
		int *desc, 
		const int *m, const int *n, const int *mb, const int *nb, const int *irsrc, const int *icsrc, 
		const int *ictxt, const int *lld, int *info);
	// C = a * A.? * B.? + b * C
	void pdgemm_(
		const char *transa, const char *transb,
		const int *M, const int *N, const int *K,
		const double *alpha,
		const double *A, const int *IA, const int *JA, const int *DESCA,
		const double *B, const int *IB, const int *JB, const int *DESCB,
		const double *beta,
		double *C, const int *IC, const int *JC, const int *DESCC);
	void pzgemm_(
		const char *transa, const char *transb,
		const int *M, const int *N, const int *K,
		const double *alpha,
		const complex<double> *A, const int *IA, const int *JA, const int *DESCA,
		const complex<double> *B, const int *IB, const int *JB, const int *DESCB,
		const double *beta,
		complex<double> *C, const int *IC, const int *JC, const int *DESCC);
                void pzgemv_(
		const char *transa,
		const int *M, const int *N,
		const double *alpha,
		const complex<double> *A, const int *IA, const int *JA, const int *DESCA,
		const complex<double> *B, const int *IB, const int *JB, const int *DESCB, const int *K, 
                                           const double *beta,
		complex<double> *C, const int *IC, const int *JC, const int *DESCC,const int *L);
                              void pdgemv_(
		const char *transa,
		const int *M, const int *N,
		const double *alpha,
		const double *A, const int *IA, const int *JA, const int *DESCA,
		const double *B, const int *IB, const int *JB, const int *DESCB, const int *K, 
                                           const double *beta,
		double *C, const int *IC, const int *JC, const int *DESCC,const int *L);

	void pzgetrf_(
		const int *M, const int *N, 
		complex<double> *A, const int *IA, const int *JA, const int *DESCA,
		int *ipiv,  int *info);
}

/*
class ScalapackConnector
{
public:
	static void transpose_desc( int desc_T[9], const int desc[9] )
	{
		desc_T[0] = desc[0];
		desc_T[1] = desc[1];
		desc_T[2] = desc[3];	desc_T[3] = desc[2];
		desc_T[4] = desc[5];	desc_T[5] = desc[4];
		desc_T[6] = desc[6];	desc_T[7] = desc[7];
		desc_T[8] = desc[8];
	}

	static void blacs_gridinit( int &ictxt, const char order, const int nprow, const int npcol )
	{
		blacs_gridinit_(&ictxt, &order, &nprow, &npcol);
	}
	
	static void blacs_gridinfo( const int &ictxt, int &nprow, int &npcol, int &myprow, int &mypcol )
	{
		blacs_gridinfo_( &ictxt, &nprow, &npcol, &myprow, &mypcol );
	}
	
	static int numroc( const int n, const int nb, const int iproc, const int srcproc, const int nprocs )
	{
		return numroc_(&n, &nb, &iproc, &srcproc, &nprocs);
	}
	
	static void descinit( 
		int *desc, 
		const int m, const int n, const int mb, const int nb, const int irsrc, const int icsrc, 
		const int ictxt, const int lld, int &info )
	{
		descinit_(desc, &m, &n, &mb, &nb, &irsrc, &icsrc, &ictxt, &lld, &info);
//		descinit_(desc, &n, &m, &nb, &mb, &irsrc, &icsrc, &ictxt, &lld, &info);
	}
	
	// C = a * A.? * B.? + b * C
	static void pgemm(
		const char transa, const char transb,
		const int M, const int N, const int K,
		const double alpha,
		const double *A, const int IA, const int JA, const int *DESCA,
		const double *B, const int IB, const int JB, const int *DESCB,
		const double beta,
		double *C, const int IC, const int JC, const int *DESCC)
	{
//		int DESCA_T[9], DESCB_T[9], DESCC_T[9];
//		transpose_desc( DESCA_T, DESCA );
//		transpose_desc( DESCB_T, DESCB );
//		transpose_desc( DESCC_T, DESCC );
//		pdgemm_(
//			&transb, &transa,
//			&N, &M, &K,
//			&alpha,
//			B, &JB, &IB, DESCB_T,
//			A, &JA, &IA, DESCA_T,
//			&beta,
//			C, &JC, &IC, DESCC_T);
		pdgemm_(
			&transa, &transb,
			&M, &N, &K,
			&alpha,
			A, &JA, &IA, DESCA,
			B, &JB, &IB, DESCB,
			&beta,
			C, &JC, &IC, DESCC);
	}
};
*/

#endif