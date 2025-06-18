Installation
------------

mex for glnx64:

mex tridisloc3d.c poly3d_tridisloc3d.c matrixpoly3d1.c safetan.c infcoeff.c -o tridisloc3d -lm


Use
---

Type 'help tridisloc3d' at Matlab command line.

IF YOU USE THIS CODE YOU MUST CITE
Thomas, A. L. (1993), Poly3D: A three-dimensional, polygonal element, 
 displacement discontinuity boundary element computer program with 
 applications to fractures, faults, and cavities in the Earth's crust, 
 M.S. thesis, Dep. of Geol. and Environ. Sci., Stanford Univ., Stanford, Calif.



Revision history
----------------

Original. June 1993. Andrew L. Thomas 

Unknown evolution.

x.05-2005. Z. Liu.
/* disloctest_mec - MEX interface to disloctest.c */
/* written by Z. Liu, on May 2005                   */
modified on 10-27-2006. remove features that cause memory leaking.

01-2011. AMB (ambrad@cs.stanford.edu)
 - Renamed disloctest_mex to tridisloc3d and changed a few things to make inputs
   and outputs consistent with disloc3d.
 - Renamed disloctest.c to poly3d_tridisloc3d.c to reflect lineage.
 - Modified disloctest.c and infcoeff.c so that we can output U, D, S, where D
   and S are optional. Computation time is faster if only U is requested
   (in this case, computation time is slower by far less than 1% than in
   disloctest_mex; the small increase is due to error checking, flipping some
   values to follow disloc3d conventions, and some bits of logic to distinguish
   between U and U,D,S calls).
   - U is the 3xn matrix of displacements, as previously.
   - D is the 9xn matrix of displacement derivatives.
   - S is the 6xn matrix of stresses.
 - Take mu as an additional input.
 - Some cosmetic changes.
 - Ran valgrind and fixed a memory bug and freed up some unfreed memory.
 - Comparing tridisloc3d with disloc3d on two triangles and one rectangle shows
   that tridisloc3d is 100-300 times slower than disloc3d. Here are some time
   breakdowns:
   - The error checking, memory allocation, and swapping of values in
     tridisloc3d.c require negligible time, far less than 1% of total time.
   - In poly3d_tridisloc3d.c, >= 99% of time is spent in displ_strain or lower
     in the call stack. Hence the linked-list code, memory allocation and
     cleanup, and so on are negligible.
