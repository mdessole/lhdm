#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h> /* needed for memcpy() */
#include <time.h>

#include "mex.h"
#include "lapack.h"
#include "blas.h"

#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))

#define LAPACK_ROW_MAJOR               101
#define LAPACK_COL_MAJOR               102

typedef const long int diff_t;

clock_t start, end;
clock_t startloc, endloc;
double cpu_time_used;
double cpu_time_usedloc;

char *chn = "N";
char *cht = "T";
char *chl = "L";
char *chf = "F";
char *chc = "C";

typedef struct
{ double val;
  double nrm;
  int idx;
}wrapStruct;

int cmpStruct (const void * a, const void * b)
{

  wrapStruct *wrapStructA = (wrapStruct *)a;
  wrapStruct *wrapStructB = (wrapStruct *)b;

  double valA = wrapStructA->val;
  double valB = wrapStructB->val;

  if (valA > valB)
    return -1;
  else if (valA < valB)
    return 1;
  else
    return 0;

}

int cmpStruct_nrm (const void * a, const void * b)
{

  wrapStruct *wrapStructA = (wrapStruct *)a;
  wrapStruct *wrapStructB = (wrapStruct *)b;

  double nrmA = wrapStructA->nrm;
  double nrmB = wrapStructB->nrm;

  if (nrmA > nrmB)
    return -1;
  else if (nrmA < nrmB)
    return 1;
  else
    return 0;

}


void permute_marked(int m, int ncols,
            double *a, int lda,
            int *jpvt, int n, double *w,
            int *marked,
            int *idxI, int nb_idxI){
  int jb = 0, jblda;
  int jt = ncols-1, jtlda;
  int nscambi = 0;
  int i, j, jc, itmp, jclda;
  double dtmp;

  for (j=0;j<nb_idxI;j++){
    jc = idxI[j];

    while ((jb<jt) && ((marked[jt] == 1) || (marked[jb]==2))){
      itmp = marked[jt];
      marked[jt] = marked[jb];
      marked[jb] = itmp;
      itmp = jpvt[jt];
      jpvt[jt] = jpvt[jb];
      jpvt[jb] = itmp;
      dtmp = w[jt];
      w[jt] = w[jb];
      w[jb] = dtmp;
      jblda = jb*lda;
      jtlda = jt*lda;
      for (i=0;i<m;i++){
    dtmp = a[jtlda];
    a[jtlda] = a[jblda];
    a[jblda] = dtmp;
    jblda++;
    jtlda++;
      }
      nscambi +=1;
      while((jb<ncols) && (marked[jb]!=0) && (marked[jb]!=2))
    jb+=1;
      while((jt>=0) && (marked[jt]!=0) && (marked[jt]!=1))
    jt-=1;
    }

    if (marked[jc]==1){
      while((jb<ncols) && (marked[jb]!=0) && (marked[jb]!=2))
    jb+=1;
      if ((jc<=jb) || (jc< nb_idxI))
    continue;
      if (marked[jb]==0){
    itmp = marked[jc];
    marked[jc] = marked[jb];
    marked[jb] = itmp;
    itmp = jpvt[jc];
    jpvt[jc] = jpvt[jb];
    jpvt[jb] = itmp;
    dtmp = w[jc];
    w[jc] = w[jb];
    w[jb] = dtmp;
    jblda = jb*lda;
    jclda = jc*lda;
    for (i=0;i<m;i++){
      dtmp = a[jclda];
      a[jclda] = a[jblda];
      a[jblda] = dtmp;
      jblda++;
      jclda++;
    }
    nscambi +=1;
    jb+=1;
      }
    }else if (marked[jc]==2){
      while((jt>=0) && (marked[jt]!=0) && (marked[jt]!=1))
    jt-=1;
      if ((jc>=jt) || (jc >= ncols))
    continue;
      if (marked[jt]==0){
    itmp = marked[jc];
    marked[jc] = marked[jt];
    marked[jt] = itmp;
    itmp = jpvt[jc];
    jpvt[jc] = jpvt[jt];
    jpvt[jt] = itmp;
    dtmp = w[jc];
    w[jc] = w[jt];
    w[jt] = dtmp;
    jclda = jc*lda;
    jtlda = jt*lda;
        for (i=0;i<m;i++){
      dtmp = a[jclda];
      a[jclda] = a[jtlda];
      a[jtlda] = dtmp;
      jclda++;
      jtlda++;
    }
    nscambi +=1;
    jt-=1;
      }
    }

  }


  return;
}





void DM_perm(int mfull, int nfull, int ncols, int rank,
         double *w,
         double *a, int lda,
         int *jpvt,
         double thresnrm,
         double thresw,
         double threscos_low,
         int nb_idxMax, int ncMax,
         int *nb_idxI,
         int *workint, double *workdouble, wrapStruct *tmpStr,
         int *maxcw, int *maxcn, int *verbose){
  /*
    This function computes a column ordering based on deviation maximization
    input:
    mfull = numero di righe della matrice originale
    nfull = numero di colonne della matrice originale
    ncols = numero di colonne rimaste da processare
    rank = numeri di colonne già processate
    a = matrix
    lda = leading dimension = nb of rows
    jpvt = vettore che contiene gli indici permutati delle colonne di a
    thresnrm = threshold on the norms
    threscos_low = threshold on the cosines
    nb_idxMax = maximum number of columns selected by DM  <= ncols
    ncMax = maximum number of candidate columns >= nb_idxMax
    nb_idxI = nb of columns selected by DM to be triangularized
    workint = integer workspace
    workdouble = integer workspace
    tmpStr = workspace for sorting

    **auxiliary memory**
    ncols wrapstruct
    ncols+nb_idxMax+nb_idxMax int (mark, idxtmp, idxI )
    (m*nc)+(nc*nc) double (as, cosmat )
  */

    int idxj,i,j,nc, ni, startpos1, startpos2;

    double *as, *cosmat;
    int *idxI, *idxtmp, *mark;
    double one = 1.0, zero = 0.0, maxval, cc, min_nrm, maxnrm;
    int m = mfull - rank; // numero di righe della sottomatrice da fattorizzare
    int sel; // numero di colonne selezionate da DM >=1

    double dtmp;
    int col;
    int tmp;

    //startloc = clock();

    //idxI = auxiliary vector of nb_idxMax
    idxI = workint;
    //mark = auxiliary vector of size ncols
    mark = &workint[nb_idxMax];
    //idxtmp = auxiliary vector of size nb_idxMax
    idxtmp = &workint[nb_idxMax+ncols];

    if (*verbose) printf("Initialize tmpStr, verbose %d \n", *verbose);
    // initialization
    for (i=0;i<ncols;i++){
        // the structure tmpStr is used as wrapper to find the ordering ORD 
        // that sorts the vector norm with quicksort
        tmpStr[i].val = w[rank+i];
        tmpStr[i].idx = i;
        // at the end mark[i]=1 iff the i-th column is selected by DM
        mark[i] = 0;
    }



    //as =  auxiliary vector of size m*ncMax, will contain the matrix a rescaled by the norms
    as = workdouble;
    //cosmat = auxiliary vector of size ncMax*ncMax, will contain the cosine
    cosmat = &workdouble[m*ncMax];

    if (*verbose) printf("sort according to dual variable \n");
    // quicksort call
    // restituisce una struttura ordinata in maniera non crescente a seconda dei valori in val
    // qui ci interessano gli indici idx
    qsort(tmpStr, ncols, sizeof(wrapStruct),cmpStruct);

    if (*verbose) printf("select first columns \n");

    // DM procedure: first add column with largest w component
    *nb_idxI = 1;
    sel = 1;
    idxI[0] = tmpStr[0].idx;
    mark[tmpStr[0].idx] = 1;

    if (ncMax > 1){
        if (*verbose) printf("compute partial norms of candidate columns \n");
        // find the number of candidates based on w
        nc = 0;
        maxval = thresw*tmpStr[0].val;
        maxnrm = 0.0;
        while((nc<ncols) && (nc<ncMax) && (tmpStr[nc].val > maxval)) {
            col = rank + tmpStr[nc].idx;
            dtmp = 0;
            for (int k=rank; k<mfull; k++) {
                dtmp = dtmp + a[col*lda + k]*a[col*lda + k];
            }
            tmpStr[nc].nrm = sqrt(dtmp);
            if (tmpStr[nc].nrm > maxnrm)
                maxnrm = tmpStr[nc].nrm;
            nc++;
        }
        maxnrm = maxnrm*thresnrm;

        if (nc > *maxcw) {*maxcw = nc;}

        if (*verbose){
            printf("candidates %d thres_cos %f \n", nc, threscos_low);
            for ( j=0; j<nc; j++)
                printf("%f\t", tmpStr[j].val );
            printf("\n");
            for (j = 0; j<nc; j++)
                printf("%f\t", tmpStr[j].nrm );
            printf("\n");
        }

        if (0){
          /* sort according to norm value*/
          if (*verbose) {
              printf("\t\t\t\t nc_w = %d \t", nc);
          }
          // order candidates using norms
          qsort(tmpStr, nc, sizeof(wrapStruct), cmpStruct_nrm);
          tmp = 0;
          min_nrm = thresnrm*tmpStr[0].nrm;
          while ((tmp < nc) && (tmpStr[tmp].nrm > min_nrm)) {
              tmp = tmp + 1;
          }
          nc = tmp;
          if (nc > *maxcn) {*maxcn = nc;}
          if (*verbose) {
              printf("nc_nrm = %d \n", nc);
          }
        }

        if (*verbose) printf("normalize columns \n");
        // set as = a(:,ORD)*diag(c), c_j = cc
        for(j=0;j<nc;j++){
            idxj = tmpStr[j].idx;
            cc = 1.0/tmpStr[j].nrm;
            startpos1 = j*m;
            startpos2 = (idxj + rank)*lda + rank;
            for(i=0;i<m;i++) {
                as[startpos1+i] = cc*a[startpos2+i];
            }
        }

        const long int nnc = nc, mm = m;

        if (*verbose) printf("compute cosine matrix nc %d m %d ncMax %d\n", nc, m, ncMax );
        
        // construct the matrix cosmat = as.T*as
        dgemm(cht, chn,
            &nnc, &nnc, &mm,
            &one, as, &mm, as, &mm,
            &zero, cosmat, &nnc);

        if (*verbose){
            for (i = 0; i<nc; i++){
                for (j = 0; j<nc; j++){
                    printf("%f\t", cosmat[i+j*nc] );
                }
                printf("\n");
            }
        }


        // select nb_idxI<=nc columns to be triangularized
        // and save their indices in idxI
        idxtmp[0] = 0;
        i = 1;  // index in ORD of current candidate column
        if (*verbose) printf("DM procedure \n");
        while(i<nc){
            // maximum value of the cosine of the angle between candidate and previously selected columns
            if (tmpStr[i].nrm < maxnrm){
                i++;
                continue;
            }
            maxval = 0.0;
            for(j=0; j<sel; j++) {
                if(maxval < fabs(cosmat[i*nc+idxtmp[j]]))
                    maxval = fabs(cosmat[i*nc+idxtmp[j]]);
            }
            // if the cosine is small enough -> angle large enough
            // then add i to selected columns
            if ((maxval < threscos_low) && (sel < nb_idxMax)) {
                if (*verbose) printf("%f \t", maxval);
                idxtmp[sel] = i;
                idxI[sel] = tmpStr[i].idx;
                mark[tmpStr[i].idx] = 1;
                sel++;
            }
            i++;
        }
        *nb_idxI = sel;
        if (*verbose) printf("\n");
    }

    if (*verbose) printf("permute selected columns \n");
    //permute the selected columns to the leftmost position
    permute_marked(mfull, ncols,
        &a[rank*lda], lda,
        &jpvt[rank], ncols, &w[rank],
        mark,
        idxI, *nb_idxI);

    //endloc = clock();
    //cpu_time_usedloc += ((double) (endloc - startloc)) / CLOCKS_PER_SEC;


    return;
}



double d_sign_DM(double a, double b)
{
  double x;
  x = (a >= 0 ? a : - a);
  return( b >= 0 ? x : -x);
}




void G1_DM(const double *A, const double *B, double *cc, double *ss, double *sig) {
    double XR;
    double YR;

    if (fabs(*A) > fabs(*B)) {
        XR = (*B)/(*A);
        YR = sqrt(1 + XR*XR);
        *cc = d_sign_DM(1/YR, *A);
        *ss = (*cc)*XR;
        *sig = fabs(*A)*YR;
    }
    else if (*B != 0) {
        XR = (*A)/(*B);
        YR = sqrt(1 + XR*XR);
        *ss = d_sign_DM(1/YR, *B);
        *cc = (*ss)*XR;
        *sig = fabs(*B)*YR;
    }
    else {
        *sig = 0;
        *cc = 0;
        *ss = 1;
    }
    return;
}



void solve_triangular(double *a, double *zz, int mda, int nsetp) {
    
    int IP;
    
    for (int L=1; L<=nsetp; L++) {
        IP = nsetp - L;
        if (L != 1) {
            for (int ii=0; ii<=IP; ii++) {
                zz[ii] = zz[ii] - a[(IP+1)*mda + ii]*zz[IP+1];
            }
        }
        zz[IP] = zz[IP]/a[IP*mda + IP];
    }
    //printf("zz[nsetp-1] = %g \n", zz[nsetp-1]);
    return;
}





int nnls_DM(a, mda_, m_, n_, b,
            x, tau_c, tau_w, tau_u, nb_idxMax_, ncMax_, iter_vec)
            double *a;
            int *mda_, *m_, *n_, *nb_idxMax_, *ncMax_;
            double *b, *x, *tau_c, *tau_w, *tau_u;
            int *iter_vec;

    {
  /*  Purpose */
  /*  ======= */

  /*  LHDM solver a NNLS problem with Lawson-Hanson with Deviation Maximization (LHDM) */
  /*  NNLS problem: min ||A*X-B||, s.t. X>=0. */

  /*  Arguments */
  /*  ========= */

  /*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
  /*          On entry, the M-by-N matrix A. */
  /*          On exit, the upper triangle of the array contains the */
  /*          min(M,N)-by-N upper trapezoidal matrix R */

  /*  MDA     (input) INT */
  /*          The leading dimension of the array A. MDA >= max(1,M). */

  /*  M       (input) INT */
  /*          The number of rows of the matrix A. M >= 0. */

  /*  N       (input) INT */
  /*          The number of columns of the matrix A.  N >= 0. */

  /*  B       (input/output) DOUBLE PRECISION array, dimension (M)  */
  /*          On entry, the right-hand side B. */
  /*          On exit, tthe vector Q^T*B */
        
  /*  X       (output) DOUBLE PRECISION array, dimension (N)  */
  /*          Solution vector of NNLS problem. */
  
  /*  TAU_C   (input) DOUBLE */
  /*          Deviation Maximization parameter 0=<tau_c<=1*/

  /*  TAU_W   (input) DOUBLE */
  /*          Deviation Maximization parameter 0=<tau_w<=1*/

  /*  TAU_U   (input) DOUBLE */
  /*          Deviation Maximization parameter 0=<tau_u<=1*/      
        
  /* NB_IDXMAX (input) INT */   
  /*          Maximum number of indices the deviation-maximization can select. */
 
  /* NC_MAX   (input) INT */   
  /*          Maximum number of candidate indices of */
  /*          the deviation-maximization can select. */      
       
  /* ITER_VEC (output) INT array, dimension (3*max(M,N))*/   
  /*          The i-th entry contains the cardinality if the  */
  /*          passive set at the i'th algorithmic step. */    
         
  /*  THRES   (input) DOUBLE PRECISION array, dimension (2) */
  /*          Deviation Maximization parameters:   */
  /*          THRES[0] contains the value of 0=<delta<=1*/
  /*          THRES[1] contains the value of 0=<tau_<=1*/

  
  /*  INFO    (output) INT */
  /*          = 0: successful exit. */
  /*          < 0: if INFO = -i, the i-th argument had an illegal value. */


  /*  Based on */
  /*  "The Lawson-Hanson algorithm with Deviation Maximization: finite convergence and sparse recovery" */
  /*  M. Dessole, M. Dell'Orto, F. Marcuzzi */
  /*  Dipartimento di Matematica "Tullio Levi Civita" */
  /*  University of Padova, Italy */
  
  /*  ===================================================================== */

        
        int verbose = 0;
        int zeroint = 0;   // verbose for DM function
        int counter = 0;
        int tmpint;     
        
        int mda = *mda_;
        int m = *m_;
        int n = *n_;
        int maxmn = max(m,n);
        const long int mm = m, nn = n, mmda = mda, i_2 = 1;
        long int i_1, i_3, nnb_idxMax, nnb_idxI;

        double rnorm[1], *w, *zz;
        int *index, mode[1], iter[1];

        w =  (double*) malloc(n*sizeof(double));
        zz =  (double*) malloc(maxmn*sizeof(double));
        index = (int*) malloc(n*sizeof(int));                
        

        int nb_idxMax = *nb_idxMax_; // maximum number of columns selected by DM
        nnb_idxMax = nb_idxMax;
        int ncMax = *ncMax_;  // maximum number of candidate columns
        int nb_idxI = 0; // nb of columns selected by DM to be triangularized

        int maxcw = 0;
        int maxcn = 0;
        
        int maxcounter = 0;
        
        
        
        // arrays used by DM

        int* workint;
        workint = (int*) malloc((2*nb_idxMax + n)*sizeof(int));
        //printf("workint %d \n", m*ncMax + 2*ncMax);


        double* workdouble;
        workdouble = (double*) malloc((m*ncMax + ncMax*ncMax)*sizeof(double));


        wrapStruct* tmpStr;
        tmpStr = (wrapStruct*) malloc((n)*sizeof(wrapStruct));
        
        
        // arrays used for the triangularization with LAPACKE
        double* tau_h;
        tau_h = (double*) malloc((nb_idxMax)*sizeof(double));

        double* v;
        v = (double*) malloc((m*nb_idxMax)*sizeof(double));

        double* work;
        work = (double*) malloc((max(m,n))*sizeof(double));
   
        double* TU;
        TU = (double*) malloc((nb_idxMax*nb_idxMax)*sizeof(double));

        const long int ldwork = n;//nb_idxMax;//
        double* work_lapacke;
        work_lapacke = (double*) malloc((n*nb_idxMax)*sizeof(double));
        
        
        int wmax_pos;
        int flag;
        
        int i__1, i__2, i__3;
        i__2 = 1;
        
        int nsetpLH;
        
        int i, j, k, jj, ii, tmp_jj, itemp;
        double sm, T, alpha, temp;
        
        
        double cc, ss;
        double unorm;
        double tolerance = 0.0;
        double eps = dlamch_("E");

          
        
        if (verbose) printf("tolerance %g eps %g m %d n %d mda %d\n", tolerance, eps, m, n, mda);
        // set tolerance equal to 10*m*n*eps*||A||_1
        for (jj=0; jj<n;jj++){
            unorm = 0.0;
            for (ii=0; ii<m;ii++){
                unorm += fabs(a[ii+jj*mda]);
            }
            if (unorm > tolerance)
                tolerance = unorm;
        }

        tolerance = tolerance*m*n*eps*10;
        if (verbose) printf("tolerance %g eps %g\n", tolerance, eps);
        *mode = 1;
        if (m<=0 || n<=0) {*mode = 2; return 0;} // da controllare
        
              
        
        int ITER = 0;
        int ITMAX = 2*min(m,n); //min(m,n); //3*min(m,n)


        if (verbose) printf("initializing ...\n");
        // initialize arrays
        for(int i=0; i<n; i++) {
            x[i] = 0.0; index[i] = i; w[i] = 0.0;
        }
        if (verbose) printf("x, index, w initialized. \n");
        for(int i=0; i<maxmn; i++) {
            zz[i] = 0.0;
        }
        if (verbose) printf("zz initialized. \n");
        
        //int startz = 0;
        int nsetp = 0;
        
        
        
        // compute components of the dual vector w
        for (int i=0; i<n; i++) {
            sm = 0;
            for (int j=0; j<m; j++) {
                sm = sm + a[j + i*mda]*b[j];
            }
            w[i] = sm;
        }

        // find out if max(w) > tolerance
        wmax_pos = 0;
        for (int i=0; i<n; i++) {
            if (w[i] > tolerance) {
                wmax_pos = 1;
                break;
            }
        }
        
        
        // OUTER LOOP
        while ((nsetp<n) && wmax_pos && (ITER<ITMAX)) {
            
                     
        
            ITER = ITER + 1;
            
            // CALL DM_perm
            if (ITER<2){
                tmpint = 1;
                DM_perm(m, n, n-nsetp, nsetp,
                w, a, mda,
                index,
                *tau_u,
                *tau_w,
                *tau_c,
                nb_idxMax,
                tmpint,
                &nb_idxI,
                workint, workdouble, tmpStr,
                &maxcw, &maxcn, &zeroint);
            }
            else {
                DM_perm(m, n, n-nsetp, nsetp,
                w, a, mda,
                index,
                *tau_u,
                *tau_w,
                *tau_c,
                nb_idxMax,
                ncMax,
                &nb_idxI,
                workint, workdouble, tmpStr,
                &maxcw, &maxcn, &zeroint);
            }
            
           
            //printf("sel = %d\n", nb_idxI);
            
            //printf("iter = %d \t\t j = %d \t\t w = %g \n", ITER, index[nsetp], w[nsetp]);

                

            if (verbose)
                printf("dm selected = %d, zeroint  %d \n", nb_idxI, zeroint);
            
            
            if (verbose)
                printf("Compute reflectors \n");
            
            for (int k=1; k<=nb_idxI; k++) {
                // compute householder transf for k-th new column
                // possibile problema per k = nb_idxI e nsetp+nb_idxI = n ?
                i__1 = m-nsetp-k+1;
                i__3 = nb_idxI-k;
                i_1 = m-nsetp-k+1;
                i_3 = nb_idxI-k;
                dlarfg_(&i_1,
                    &a[(nsetp+k-1)*mda + nsetp+k-1],
                    &a[(nsetp+k-1)*mda + nsetp+k], &i_2,
                    &tau_h[k-1]);


	            //save the reflector in the k-th column of v
	            if (1){
                    v[(k-1)*mda + k-1] = 1;
                    for (int i=k; i<m-nsetp; i++) {
                        v[(k-1)*mda + i] = a[(nsetp+k-1)*mda +nsetp+i];
	                }
	            }

	            // apply this transf to the other selected columns
	            if (k<nb_idxI)
                    dlarfx_(chl,
                    &i_1, &i_3,
                    &v[(k-1)*mda + k-1],
                    &tau_h[k-1], &a[(nsetp+k)*mda + nsetp+k-1], &mmda,
		            work);
            } // end for k

            i__1 = m-nsetp;
            i_1 = m-nsetp;
            i__3 = n-nsetp-nb_idxI;
            i_3 = n-nsetp-nb_idxI;
            //i__3 = nb_idxI-k;
            nnb_idxI = nb_idxI;

            if (verbose)
                printf("Compute triangular factor \n");

            // compute the triangular factor TU
            dlarft_(chf, chc, &i_1, &nnb_idxI, v, &mmda, tau_h, TU, &nnb_idxMax);

            if (verbose)
                printf("apply WY \n");

            //apply the new transf to the other columns of a
            dlarfb_(chl, cht, chf, chc,
                &i_1, &i_3,&nnb_idxI,
                v, &mmda, TU, &nnb_idxMax,
                &a[(nsetp+nb_idxI)*mda + nsetp], &mmda, work_lapacke, &ldwork);

            //apply the new transf to the rhs
            dlarfb_(chl, cht, chf, chc,
                &i_1, &i_2, &nnb_idxI, v, &mmda,
                TU, &nnb_idxMax, &b[nsetp], &mmda, work_lapacke, &ldwork);


            // zero subdiagonal elements
            for (int col=0; col<nb_idxI; col++) {
                for (int row=nsetp+col+1; row<m; row++) {
                    a[(nsetp+col)*mda + row] = 0;
	            }
            }
            
            
            

            // UPDATE INDICES
            
            nsetpLH = nsetp+1;
            
            nsetp = nsetp + nb_idxI;

            iter_vec[ITER-1] = nsetp;

            counter = 0;
            alpha = 0;

            
            while(alpha==0) {
                counter = counter + 1;
                //printf("counter = %d\n", counter);
            
            // copy b[1,...,nsetp] in zz
            for (int k=0; k<nsetp; k++) {
                zz[k] = b[k];
            }
            
            // solve the triangular system
            // store the solution temporarily in zz
            
            solve_triangular(a, zz, mda, nsetp);
            
            
            // see if all new constrained coeffs are feasible
            // if not compute alpha

            
            
            alpha = 2;
            for (int i=0; i<nsetp; i++) {
                if (zz[i] <= 0.0) {
                    T = - x[i]/(zz[i] - x[i]);
                    //printf("i = %d\t\t T = %g\n", i, T);
                    if (alpha > T) {
                        alpha = T;
                        jj = i;
                    }
                }
            }
            if (verbose)
                printf("alpha = %g\n", alpha);
            if (alpha==0) {
                 
                nsetp = nsetp - (int)ceil(nb_idxI/pow(2, counter));
                if (nsetp<nsetpLH) {
                    nsetp = nsetpLH;
                }
                //printf("nsetp = %d\n", nsetp);
                if (nsetp<=0 || counter>100) {
                    printf("Exiting on alpha inner loop \n");
                    goto L10;}
            }    
                
            }// end while alpha 
            //printf("iter = %d \t\t counter = %d\n", ITER, counter);
            
            if (counter>maxcounter) {
                maxcounter = counter;
            }
            
            // if all new constrained coeffs are feasible the alpha will still = 2
            // if so do not enter the secondary loop
            
            // INNER LOOP
            while ((alpha != 2) && (ITER<ITMAX)) {
                
                
                
                ITER = ITER + 1;
                //printf("ENTER INNER LOOP \n");
                if (verbose) printf("ENTER INNER LOOP \n");
                counter = 0;
                
                for (int i=0; i<nsetp; i++) {
                    x[i] = x[i] + alpha*(zz[i] - x[i]);
                    if (x[i]<=0.0) counter++;
                }
                
         //       printf("\n x  (%d)\n", nsetp);
         //       for (int i=45; i<nsetp; i++) {
         //           printf("%g \n", x[i]);
         //       }
                
                
                if (verbose) printf("counter non positive entries %d, alpha %g, jj %d \n", counter, alpha, jj);

                
                // modify a and b and the index arrays to move coefficient jj from set p to set z
                
                tmp_jj = jj;
                
                flag = 1;
                while (flag) {
                    flag = 0;
                        
                    jj = tmp_jj;
                    //printf("jj = %d\n", jj);
                        
                    x[jj] = 0;
    
                    if (jj < nsetp - 1) {
                        for (int j=jj+1; j<nsetp; j++) {
                            G1_DM(&a[j*mda + j-1], &a[j*mda + j], &cc, &ss, &a[j*mda + j-1]);
                            a[j*mda + j] = 0;
                            for (int L=j+1; L<n; L++) {
                                // apply procedure G2
                                temp = a[L*mda + j-1];
                                a[L*mda + j-1] = cc*temp + ss*a[L*mda + j];
                                a[L*mda + j] = -ss*temp + cc*a[L*mda + j];
                            }
    
                            // apply procedure G2 to b
                            temp = b[j-1];
                            b[j-1] = cc*temp + ss*b[j];
                            b[j] = -ss*temp + cc*b[j];
                        }
                    }
    
                    // permute columns
                    for (int row=0; row<m; row++) {  // sposta a sinistra una riga alla volta
                        temp = a[jj*mda + row];
                        for (int col=jj+1; col<nsetp; col++) {
                            a[(col-1)*mda + row] = a[col*mda + row];
                        }
                        a[(nsetp-1)*mda + row] = temp;
                    }
    
                    // permute components of x, w and index
                    // x
                    for (int k=jj+1; k<nsetp; k++) {
                        x[k-1] = x[k];
                    }
                    x[nsetp-1] = 0;
    
                    // w
                    temp = w[jj];
                    for (int k=jj+1; k<nsetp; k++) {
                        w[k-1] = w[k];
                    }
                    w[nsetp-1] = temp;
    
    
                    // index
                   itemp = index[jj];
                   for (int k=jj+1; k<nsetp; k++) {
                       index[k-1] = index[k];
                   }
                   index[nsetp-1] = itemp;
                        
                        
                   nsetp = nsetp-1;  
                        
                        
                   // see if the remaining coeffs in set p are feasible.
                   // they should be because of the way alpha was determined.
                   // if any are infeasible it is due to round-off error.
                   // any that are nonpositive will be set to zero and moved from set p to set z
    
                   for (int jj=0; jj<nsetp; jj++) {
                       if (x[jj] <= 0.0) {
                           //printf("jj = %d \t\t x[jj] = %g\n", jj, x[jj]);
                           tmp_jj = jj;
                           flag = 1; 
                           break;
                       }
                   }
          
               } // end while flag
            
                
            // copy b into zz; then solve again and loop back
            for (int i=0; i<m; i++) {zz[i] = b[i];}
                
            solve_triangular(a, zz, mda, nsetp);  // save the result in zz
                
            alpha = 2;
            for (int i=0; i<nsetp; i++) {
                if (zz[i] <= 0.0) {
                    T = - x[i]/(zz[i] - x[i]);
                    if (alpha > T) {
                        alpha = T;
                        jj = i;
                    }
                }
            }
                
            iter_vec[ITER-1] = nsetp;    
                
            } // end while interno
            
                
            // update x
            for (int i=0; i<nsetp; i++) {
                x[i] = zz[i];
            }
            
                
            // compute components of the dual vector w
            for (int i=0; i<n; i++) {
                sm = 0;
                for (int j=nsetp; j<m; j++) {
                    sm = sm + a[j + i*mda]*b[j];
                }
                w[i] = sm;
            }

            // find out if max(w) > tolerance
            wmax_pos = 0;
            for (int i=nsetp; i<n; i++) {
                if (w[i] > tolerance) {
                    wmax_pos = 1;
                    break;
                }
            }    
                
            // all new coeffs are positive, loop back to beginning    
        
            //printf("\t\t\t\t\t\t nsetp = %d\n", nsetp);
        } // end while esterno
        // end of main loop    
        
        // compute the norm of the final residual vector
        sm = 0;
        if (nsetp < m) {
            for (int i=nsetp; i<m; i++) {
                sm = sm + b[i]*b[i];
            }
        }
        else {
            for (int j=0;j<n;j++) {
                w[j] = 0;
            }
        }
        *rnorm = sqrt(sm);     // *rnorm è posta uguale a zero se nsetp=m
        
        // copy x in zz
        for (int k=0; k<n; k++) {
            zz[k] = x[k];
        }

        // copy zz in x permuting its components
        for (int k=0; k<n; k++) {
            j = index[k];
            x[j] = zz[k];
        }

        *iter = ITER;
            
            
        if (ITER==ITMAX) {             // problema se termina con iter=itmax
            printf("NNLS quitting on iteration count. \n");
            *mode = 3;
        }
        
        if (verbose) {
            printf("maxcw = %d \n", maxcw);
            printf("maxcn = %d \n", maxcn);
        }
        L10:    
        if (verbose) printf("freeing...\n");
        free(w);
        if (verbose) printf("w done\n");
        free(zz);
        if (verbose) printf("zz done\n");
        free(index);
        if (verbose) printf("index done\n");
        free(tmpStr);
        if (verbose) printf("tmpStr done\n");
        free(workint);
        if (verbose) printf("workint done\n");
        free(workdouble);
        if (verbose) printf("workdouble done\n");
        free(tau_h);
        if (verbose) printf("tau_h done\n");
        free(v);
        if (verbose) printf("v done\n");
        free(work);
        if (verbose) printf("work done\n");
        free(TU);
        if (verbose) printf("TU done\n");
        free(work_lapacke);
        if (verbose) printf("work_lapacke done\n");
        if (verbose) printf("returning...\n");

        end = clock();
        cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
        
        iter_vec[2*min(m,n)] = maxcounter;
        

        if (verbose) {
            printf("time = %f \n", cpu_time_used);
            printf("loc = %f \n", cpu_time_usedloc);
            printf("passive set = %d\n", nsetp);
        }    
        
    return 0;    
    }
        
        
        
        
        




void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int verbose = 0;
  double *A, *b, *x, *delta, *tau, *thresnrm, *thres; // *rnorm, *w, *zz,
  double *A2, *b2;
  mxArray *Awork, *bwork;
  int *iter_vec;
  int nb_idxMax, ncMax; //*index, *mode, *iter,
  double *blocksize;
  int info;
  ptrdiff_t m,n,mda;      /* matrix dimensions */

 if (verbose) printf("start mexFunction \n");

  if(0){

    if (nrhs < 4) {
      printf("nrhs < 4\n");
      //return;
    }

    /* Check for proper number of input and output arguments */
    if (nrhs < 4) {
      //mexErrMsgIdAndTxt( "MATLAB:nnls_dm:minrhs",
      //		       "At least four input arguments required.");
      printf("MATLAB:nnls_dm:minrhs\nAt least four input arguments required.\n");
      return;
    }
    else if(nlhs > 2){
      //    mexErrMsgIdAndTxt( "Matlab:nnls_dm:maxlhs",
      //		       "Too many output arguments.");

      printf("Matlab:nnls_dm:maxlhs\nToo many output arguments.\n");
      return;
    }
    else
      printf("Number of inputs/outputs is OK.\n");

    return;
  }
    
  if (verbose) printf("declare variables \n");
  A = mxGetPr(prhs[0]); /* first input: matrix A */
  b = mxGetPr(prhs[1]); /* second input: vector b*/
  thres = mxGetPr(prhs[2]);
  delta = &thres[0];
  tau = &thres[1];
  thresnrm = &thres[2];
  blocksize = mxGetPr(prhs[3]);
  nb_idxMax =  (int) blocksize[0];
  ncMax = (int) blocksize[1];

  /* dimensions of input matrices */
  m = mxGetM(prhs[0]);
  n = mxGetN(prhs[0]);
  mda = m;

  if (verbose) printf("copy inputs m %d n %d  \n", m, n);
  /*NNLS works in-place, so we copy the inputs first*/
  Awork = mxCreateDoubleMatrix(m, n, mxREAL);
  A2 = mxGetPr(Awork);
  memcpy(A2, A, m*n*mxGetElementSize(prhs[0]));
  if (verbose) printf("A ok \n");
  bwork = mxCreateDoubleMatrix(m, 1, mxREAL);
  b2 = mxGetPr(bwork);
  memcpy(b2, b, m*1*mxGetElementSize(prhs[1]));
  if (verbose) printf("b ok \n");

  if (verbose) printf("create outputs \n");
  /* create outputs x, iter_vec */
  plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL); //x
  x = mxGetPr(plhs[0]);

  plhs[1] = (int*) mxCreateNumericMatrix(2*min(m,n)+1, 1, mxINT32_CLASS, mxREAL); //iter_vec
  iter_vec = mxGetPr(plhs[1]);

    

  info = nnls_DM(A2, &mda, &m, &n, b2, x,
                 delta, tau, thresnrm, &nb_idxMax, &ncMax, iter_vec);
 if (verbose) printf("free workspace \n");
  mxDestroyArray(Awork);
  mxDestroyArray(bwork);

  return;
}
