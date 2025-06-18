/****************************************************/
/* disloctest_mex.c - MEX interface to disloctest.c */
/* written by Z. Liu, on May 2005                   */
/****************************************************/

/* 01-2011. AMB. Output U, D, S, change calling conventions, add some error
   checking. */

#include <mex.h>
#ifdef TV
#include <sys/time.h>

static double difftime(struct timeval* t1, struct timeval* t2)
{
  static const double us = 1.0e6;
  return (t2->tv_sec*us + t2->tv_usec - t1->tv_sec*us - t1->tv_usec)/us;
}
#endif

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  char buffer[256];
  double *poutput;
  double *obs_xyz, *mvert_xyz, *bbb_bc, *displist;
  int numobs, numvert, numdispl;
  double mu, nu;
  int m, n, nrows, i, j, k;
  int calc_displ, calc_strain, calc_stress;
  mxArray *mxpoutput, *mxbbb_bc;

#ifdef TV
  struct timeval t1, t2, t3, t4;
  double d1, d2;
  gettimeofday(&t1, 0);
#endif

  /* check argument syntax */
  if (nrhs != 6 || nlhs > 3) {
    mexPrintf("tridisloc3d 0.0\n"
	      "Usage: [U D S] = disloctest_mex(xyz,nd,el,comp,mu,nu)\n");
    return;
  }

  /* check station coordinate matrix */
  obs_xyz = mxGetPr(prhs[0]);
  numobs = mxGetN(prhs[0]);
  n = mxGetM(prhs[0]);
  if (n != 3)
    mexErrMsgTxt("First argument must be a 3xn matrix containing n station "
		 "coordinates.");

  /* check mesh vertex */
  mvert_xyz = mxGetPr(prhs[1]);
  numvert = mxGetN(prhs[1]);
  n = mxGetM(prhs[1]);
  if (n !=3)
    mexErrMsgTxt("2nd argument must be a 3xn matrix containing n vertex "
		 "coordinates.");

  /* check connectivity list */
  displist = mxGetPr(prhs[2]);
  numdispl = mxGetN(prhs[2]);
  n = mxGetM(prhs[2]);
  if (n !=3)
    mexErrMsgTxt("3rd argument must be a 3xn matrix containing n triangle "
		 "connectivity lists");

  /* boundary condition */
  bbb_bc = mxGetPr(prhs[3]);
  m = mxGetN(prhs[3]);
  n = mxGetM(prhs[3]);
  if (m != numdispl)
    mexErrMsgTxt("4th argument and 3rd arguments should have the same "
		 "dimension.");
  if (n !=3)
    mexErrMsgTxt("4th argument must be a 3xn matrix containing [SS DS OP] "
		 "for each triangle.");

  /* shear modulus */
  m = mxGetM(prhs[4]);
  n = mxGetN(prhs[4]);
  if (m != 1 || n != 1)
    mexErrMsgTxt("5th argument must be a scalar (shear modulus).");
  mu = mxGetScalar(prhs[4]);

  /* poison ratio */
  m = mxGetM(prhs[5]);
  n = mxGetN(prhs[5]);
  if (m != 1 || n != 1)
    mexErrMsgTxt("6th argument must be a scalar (poisson ratio); if stress is "
		 "not requested, set mu to an arbitrary scalar.");
  nu = mxGetScalar(prhs[5]);

  /* AMB. Run through some tests of the input. */
  for (j = 0; j < numobs; j++)
    if (obs_xyz[3*j + 2] > 0.0) {
      sprintf(buffer, "Observation point %d is above free surface.", j+1);
      mexWarnMsgTxt(buffer);
    }
  for (j = 0; j < numvert; j++)
    if (mvert_xyz[3*j + 2] > 0.0) {
      sprintf(buffer, "Triangle vertex %d is above free surface.", j+1);
      mexWarnMsgTxt(buffer);
    }

  /* AMB. Swap rows 1 and 2 in bbb_bc so that the input to this function is
     in the order [SS DS OP]. disloctest uses the order [DS SS OP]. */
  mxbbb_bc = mxDuplicateArray(prhs[3]);
  bbb_bc = mxGetPr(mxbbb_bc);
  for (j = 0; j < numdispl; j++) {
    double tmp = bbb_bc[3*j];
    bbb_bc[3*j] = bbb_bc[3*j + 1];
    bbb_bc[3*j + 1] = tmp;
  }

  /* Determine outputs. */
  nrows = 3;
  calc_displ = 1;
  calc_strain = 0;
  calc_stress = 0;
  if (nlhs >= 2) {
    nrows += 9;
    calc_strain = 1;
  }
  if (nlhs >= 3) {
    nrows += 6;
    calc_stress = 1;
  }
  mxpoutput = mxCreateDoubleMatrix(nrows, numobs, mxREAL);
  poutput = mxGetPr(mxpoutput);

#ifdef TV
  gettimeofday(&t2, 0);
#endif
  disloctest(poutput, obs_xyz, mvert_xyz, displist, bbb_bc, numobs, numvert,
	     numdispl, mu, nu, calc_displ, calc_strain, calc_stress);
#ifdef TV
  gettimeofday(&t3, 0);
#endif

  /* Copy to output arrays. A little inefficient, but no way around it. Time
     testing shows that the call to disloctest consumes >= 99.5% of the time
     spent in this mex file for two triangular dislocations and one observation
     point; for more complex inputs, the time spent here is even more
     negligible. */
  plhs[0] = mxCreateDoubleMatrix(3, numobs, mxREAL);
  if (nlhs >= 2) plhs[1] = mxCreateDoubleMatrix(9, numobs, mxREAL);
  if (nlhs >= 3) plhs[2] = mxCreateDoubleMatrix(6, numobs, mxREAL);
  for (j = 0, k = 0; j < numobs; j++) {
    for (i = 0; i < 3; i++, k++)
      mxGetPr(plhs[0])[3*j + i] = poutput[k];
    if (nlhs >= 2)
      for (i = 0; i < 9; i++, k++)
	mxGetPr(plhs[1])[9*j + i] = poutput[k];
    if (nlhs >= 3)
      for (i = 0; i < 6; i++, k++)
	mxGetPr(plhs[2])[6*j + i] = poutput[k];
  }

  mxDestroyArray(mxbbb_bc);
  mxDestroyArray(mxpoutput);

#ifdef TV
  gettimeofday(&t4, 0);
  d1 = difftime(&t1, &t4);
  d2 = difftime(&t2, &t3);
  mexPrintf("tt = %f  dt = %f  \% = %e\n", d1, d2, (d1 - d2)/d1);
#endif
}

