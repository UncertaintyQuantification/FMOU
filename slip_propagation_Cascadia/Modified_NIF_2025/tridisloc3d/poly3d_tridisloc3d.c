/*****  modified by ZL on 05-2005, major changes include:
        a. argument passing
        b. remove remote stress b.c., and stress, strain calculation
        c. add freeup memory operation, otherwise, mex file crashes due to
           memory leaking!
        d. remove lots of file reading operations to speed things up.
        e. still keep dynamic linklist, but this should be changed eventually!
        f. two new functions added: disloc_test, freeup
*****************************************************************************/

/**** modified on 10-27-2006
       * remove features that cause memory leaking ****/

/* 01-2011. AMB. See README.txt for changes. In short, add strain and stress
   outputs and fix some errors. Renamed from disloctest.c. */

/******************************** Includes **********************************/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "matrixpoly3d.h"
#include "safetan.h"
#include "infcoeff.h"

/* For CPU timing, define TV. */
#ifdef TV
#include <sys/time.h>
static double difftime(struct timeval* t1, struct timeval* t2)
{
  static const double us = 1.0e6;
  return (t2->tv_sec*us + t2->tv_usec - t1->tv_sec*us - t1->tv_usec)/us;
}
#endif

/******************************** Constants *********************************/
/*-------------
 MISCELLANEOUS
--------------*/
#define MAXFILE	 			256			/* Max length of file names			*/
#define MAX_ERROR_MSG		256			/* Max length of error messages		*/
#define MAXWORDS			50			/* Max # words on getwords() line	*/
#define GLOBAL_NAME			"global"	/* Global coord system name			*/
#define ELT_CSYS_NAME		"elocal"	/* Element coord system name		*/
#define ERROR				-1			/* Return value	for function errors	*/
#define FALSE				0			/* False flag						*/
#define TRUE				1			/* True flag						*/
#define BVECTOR_BC			0			/* Burgers vector BC component flag	*/
#define TRACTION_BC			1			/* Traction BC component flag		*/


/*--------------------------------------
 NUMERICAL LIMITS USED WHEN CHECKING...
---------------------------------------*/
#define SWAP_TINY			1.0e-10		/* ...if vertices must be swapped	*/
#define TINY_ANGLE			0.5*PI/180.	/* ...if elt coord sys can be calc	*/
#define BVERT_TINY			1.0e-10		/* ...if point lies below a vertex	*/
#define COPLANAR_LIMIT		30.			/* ...if elt vertices are co-planar	*/

/*-------------------------------------
 PRINT OPTION ARRAY SIZE AND POSITIONS
--------------------------------------*/
#define NUM_PR_OPTS			5
#define DISPL				0
#define	STRAIN				1
#define PSTRAIN				2
#define STRESS				3
#define PSTRESS				4

/*--------------------------------------------------------------
 CHARS USED IN INPUT FILE PRINT STRINGS TO ENABLE PRINT OPTIONS
---------------------------------------------------------------*/
#define DISPL_CHAR			'd'
#define STRESS_CHAR			's'
#define STRAIN_CHAR			'e'
#define TRACTION_CHAR		't'
#define BVECTOR_CHAR		'b'
#define PRINCIPAL_CHAR		'p'

/*----------------------------------------
 INPUT FILE FORMAT FOR DEFINING CONSTANTS
-----------------------------------------*/
#define CONST_NAME_POS		0
#define CONST_VALUE_POS		2
#define CONST_NUM_PARAMS	3

/*------------------------------------------------------
 INPUT FILE FORMAT FOR DEFINING USER COORDINATE SYSTEMS
-------------------------------------------------------*/
#define CS_NAME_POS			0
#define CS_PARENT_POS		1
#define CS_ORIGIN_POS		2
#define CS_ROT_POS			5
#define CS_ROT_ORDER_POS	8
#define CS_NUM_PARAMS		9

/*------------------------------------------------
 INPUT FILE FORMAT FOR DEFINING OBSERVATION GRIDS
-------------------------------------------------*/
#define OG_NAME_POS			0
#define OG_DIMEN_POS		1
#define OG_PRINT_OPS_POS	2
#define OG_INPUT_CSYS_POS	3
#define OG_OBSPT_CSYS_POS	4
#define OG_DATA_CSYS_POS	5
#define OG_BEGIN_POS		6
#define OG_END_POS			9
#define OG_NUMPTS_POS		12
#define OG_MIN_NUM_PARAMS	9

/*---------------------------------------
 INPUT FILE FORMAT FOR DEFINING VERTICES
----------------------------------------*/
#define V_CHAR			    'v'
#define V_CHAR_POS			0
#define V_NAME_POS			1
#define V_CSYS_POS			2
#define V_X_POS				3
#define V_NUM_PARAMS		6

/*--------------------------------------
 INPUT FILE FORMAT FOR DEFINING OBJECTS
---------------------------------------*/
#define OBJ_CHAR			'o'
#define OBJ_CHAR_POS		0
#define OBJ_NAME_POS		1
#define	OBJ_PRINT_OPS_POS	2
#define OBJ_POS_CSYS_POS	3
#define OBJ_MIN_NUM_PARAMS	2

/*---------------------------------------
 INPUT FILE FORMAT FOR DEFINING ELEMENTS
----------------------------------------*/
#define E_CHAR			   'e'
#define E_CHAR_POS			0
#define E_NUM_VERT_POS		1
#define E_BC_CSYS_POS		2
#define E_BC_TYPE_POS		3
#define E_BC_POS			4
#define E_VERTEX_POS		7
#define E_MIN_NUM_PARAMS	10

/********************************** Macros **********************************/
#define RADIANS(A) ((A)*PI/180.)			/* Convert degrees to radians		*/
#define MAX(A,B) (((A) > (B)) ? (A):(B))		/* MAX macro				*/

/******************************** Structures ********************************/
struct csys_s {						/* -- COORDINATE SYSTEM STRUCTURE -	*/
	char			*name;			/* Coordinate system name		*/
	double			origin[3];		/* Coord sys origin	(global)	*/
	double			local_rot[3][3];	/* (To) global rotation matrix		*/
	struct csys_s	*next;				/* Ptr to next c.s. in linked list	*/
};
typedef struct csys_s csys_t;

struct obs_grid_s {					/* -- OBSERVATION GRID STRUCTURE --	*/
	char			*name;			/* Observation grid name		*/
	int			dimension;		/* Dimension of observation grid	*/
	double			begin[3];		/* Obs grid beginning coords		*/
	double			end[3];			/* Obs grid ending coords		*/
	int			numpts[3];		/* No of obs pts along x1,x2,x3		*/
	int			print[NUM_PR_OPTS];	/* Print options array			*/
	csys_t			*endpt_csys;		/* Input coordinate system		*/
	csys_t			*obspt_csys;		/* Observation grid coord system	*/
	csys_t			*outp_csys;		/* Output coord sys for obs grid	*/
	struct obs_grid_s 	*next;			/* Ptr to next o.l. in linked list	*/
};
typedef struct obs_grid_s obs_grid_t;

struct vert_s {						/* ------- VERTEX STRUCTURE -------	*/
	char			*name;			/* Vertex name				*/
	double			x[3];			/* Vertex coordinates (global)		*/
	csys_t			*csys;			/* Coordinate system for vertex		*/
	struct vert_s		*next;			/* Ptr to next vertex in linked list*/
};
typedef struct vert_s vert_t;

struct disloc_seg_s {					/* --- DISLOC SEGMENT STRUCTURE ---	*/
	double 			elt_b[3][3];		/* Proj of element b to segment b	*/
	double			trend;			/* Strike of plunging leg of d.s.	*/
	double			plunge;			/* Plunge of plunging leg of d.s.	*/
	double			local_rot[3][3];	/* Local-to-global rotation matrix	*/
	vert_t			*vert[2];		/* Dislocation segment vertices		*/
};
typedef struct disloc_seg_s disloc_seg_t;

struct elt_s {						/* ------ ELEMENT STRUCTURE -------	*/
	int			num_vertices;		/* Number of vertices			*/
	int			bc_type[3];		/* Boundary condition type array	*/
	double			bc[3];			/* Boundary condition magnitudes	*/
	csys_t			elt_csys;		/* Element-local coordinate system	*/
	double			*b[3];			/* Burgers vector array			*/
	disloc_seg_t	 	*disloc_seg;		/* Dislocation segment array		*/
	csys_t			*bc_csys;		/* Ptr to coord sys for element BCs	*/
	struct elt_s	*next;				/* Ptr to next elt in linked list	*/
};
typedef struct elt_s elt_t;

struct obj_s {						/* ------ OBJECT STRUCTURE -------- */
	char			*name;			/* Object type				*/
	int			print[NUM_PR_OPTS];	/* Print options			*/
	csys_t			*pos_csys;		/* Position coordinate system		*/
	elt_t			*first_elt;		/* Pointer to first element		*/
	elt_t			*last_elt;		/* Pointer to last element		*/
	struct obj_s		*next;			/* Ptr to next obj in linked list	*/
};
typedef struct obj_s obj_t;

/*************************** External Variables *****************************/
int			half_space_E = TRUE;		/* Half/whole space flag		*/
int			check_cond_num_E = TRUE;	/* Check matrix condition num flag	*/
double		cond_num_E = -1.0;			/* Matrix condition number		*/
char		infile_E[MAXFILE];			/* Input file name			*/
char		outfile_E[MAXFILE];			/* Output file name			*/
int			linenum_E = 0;			/* Current line # in input file		*/
int			num_elts_E = 0;			/* Number of elements			*/
int			below_vertex_E = FALSE;		/* Below vertex flag				*/
double		null_value_E = -999.0;		/* Null output value				*/

/*********************************************************************************************************/
/************************* NEW: added 98-12-09 to reflect the problem of observation point near vertices */
double			coef_exclu_E   = 0.0;		/* coef exclusion value */
int			    near_vertex_E = FALSE;		/* nearnest vertex flag for observation point */
/*       EXPLANATIONS:
         Let's obs_pt be an observation point.
         Put flag near_vertex_E to FALSE.
         For each vertex in current project:
         Let's d be the mean length of segments containing v
         Let's l be the distance from v to a choosen observation point
         Then if l<d*coef_exclu => near_vertex_E = TRUE

         Then, before printing computed values for this observation point:
         if near_vertex_E=TRUE => print null_value_E
         otherwise print computed values

    So coef_exclu_E = 1.0, means 100% of the mean length of segments containing v.
*/
/*********************************************************************************************************/

double		shear_mod_E		= -1.0;		/* Shear modulus					*/
double		psn_ratio_E		= -1.0;		/* Poisson's ratio					*/
double		youngs_mod_E	= -1.0;		/* Young's modulus					*/
double		bulk_mod_E		= -1.0;		/* Bulk modulus						*/
double		lame_lambda_E	= -1.0;		/* Lame's lambda					*/

int			print_elt_geom_E = FALSE;	/* Print element geometry flag		*/
char		*elt_geom_csys_name_E =NULL;/* Element geometry coord sys name	*/

int			rem_stress_bc_E = TRUE;		/* Remote stress vs strain bc flag	*/
double		rem_stress_E[3][3];			/* Remote stress tensor				*/
double		rem_strain_E[3][3];			/* Remote strain tensor				*/

double		*b_vector_E;				/* Burger's vector array			*/
double		**ic_matrix_E;				/* Influence coeff matrix			*/

csys_t		*first_csys_E     = NULL;	/* 1st memb of csys linked list		*/
obs_grid_t	*first_obs_grid_E = NULL;	/* 1st memb of obs grid linked list	*/
obj_t		*first_obj_E      = NULL;	/* 1st memb of obj linked list		*/
elt_t		*first_elt_E      = NULL;	/* 1st memb of elt linked list		*/
vert_t		*first_vert_E     = NULL;	/* 1st memb of vert linked list	*/


/************************* ANSI Function Declarations ***********************/
void	determine_burgers_vectors(void);
void	displ_strain(int calc_displ, int calc_strain, double x[3],
		double displ[3], double strain[3][3], elt_t *omit_elt, int elt_num);
void	displ_strain_poly_elt(int calc_displ, int calc_strain,
		elt_t *current_elt, double x[3], double displ[3],
		double strain[3][3], elt_t *omit_elt, int under, int elt_num, int ielem);
csys_t	*find_csys(const char *name);
vert_t	*find_vert(const char *name);
void	get_elt_info(elt_t **current_elt, obj_t *current_obj, double *displist, double *bbb_bc,
		int numdispl);
void	displ_strain_ics_poly_elt(int calc_displ, int calc_strain,
		elt_t *current_elt, double x[3],
		double displ_ic[3][3], double strain_ic[3][3][3],
		elt_t *omit_elt, int elt_num, int ielem);
void	print_obs_grid_data(double *poutput, int numobs, double mu);
void	p_error(const char *error_msg, const char *line);
void	display_msg(char *_msg);
void	print_obs_pt_data(obs_grid_t *current_obs_grid, double x[3],
		double displ[3], double strain[3][3]);
int		read_objs_elts_verts(double *mvert_xyz, double *displist, double *bbb_bc, int numvert, int numdispl);
int		read_observation_grids(double *obs_xyz, int numobs, char *calc_code);
void	setup_global_coords(void);
void    disloctest(double *poutput, double *obs_xyz, double *mvert_xyz, double *displist,
		   double *bbb_bc, int numobs, int numvert, int numdispl, double mu, double nu,
		   int calc_displ, int calc_strain, int calc_stress);
void    freeup(void);

/******************************* function disloctest ************************************/
/* AMB. Include strain and stress. */
 void disloctest(double *poutput, double *obs_xyz, double *mvert_xyz, double *displist,
		 double *bbb_bc, int numobs, int numvert, int numdispl, double mu, double nu,
		 int calc_displ, int calc_strain, int calc_stress)
{
  /* AMB */
  char calc_code[4];
  int k;
  /* Determine what needs to be computed. */
  k = 0;
  if (calc_displ) {
    calc_code[k] = DISPL_CHAR;
    k++;
  }
  if (calc_strain) {
    calc_code[k] = STRAIN_CHAR;
    k++;
  }
  if (calc_stress) {
    calc_code[k] = STRESS_CHAR;
    k++;
  }
  calc_code[k] = '\0';

	 /* constants */
        psn_ratio_E = nu;
	half_space_E = TRUE;

	/* Defines the global coordinate system, making it the first member in the linked list
           of coordinate systems (first_csys_E).
        ----------------------------------------------------*/
	setup_global_coords();

	/* Read observation grids
	-------------------------*/
	read_observation_grids(obs_xyz, numobs, calc_code);

	/* Read elements and vertices
	-----------------------------*/
	read_objs_elts_verts(mvert_xyz, displist, bbb_bc, numvert, numdispl);

	/* Solve for burger's vector for each element
	---------------------------------------------*/
	determine_burgers_vectors();

	/* Calculate displacements and stresses along observation grids
	---------------------------------------------------------------*/
	print_obs_grid_data(poutput, numobs, mu);

       /* free up allocated memory */
        freeup();
 }


/************************ Function freeup ***********************************/
/*** free dynamically allocated memory, necessary, otherwise,
     segmentation fault error occurs !

     AMB: free ->name and a few other fields.
*****************************************************************************/
void freeup()
{

	obs_grid_t  *current_obs_grid;
	elt_t       *current_elt;
	vert_t      *current_vert;
	obj_t       *current_obj;
	csys_t      *current_csys;
        /* observation grid */
	current_obs_grid = first_obs_grid_E;
	while (current_obs_grid != NULL)
	{
		first_obs_grid_E = current_obs_grid->next;
		free(current_obs_grid->name);
		free(current_obs_grid);
		current_obs_grid = first_obs_grid_E;
	}
        /* element list */
	current_elt = first_elt_E;
	while (current_elt != NULL)
	{
		first_elt_E = current_elt->next;
		free(current_elt->elt_csys.name);
		free(current_elt->disloc_seg);
		free(current_elt);
		current_elt = first_elt_E;
	}
        /* vertex list */
	current_vert = first_vert_E;
	while (current_vert != NULL)
	{
		first_vert_E = current_vert->next;
		free(current_vert->name);
		free(current_vert);
		current_vert = first_vert_E;
	}
        /* object list */
	current_obj = first_obj_E;
	while (current_obj != NULL)
	{
		first_obj_E = current_obj->next;
		free(current_obj->name);
		free(current_obj);
		current_obj = first_obj_E;
	}
        /* coordinate system */
	current_csys = first_csys_E;
	free(first_csys_E->name);
	while (current_csys != NULL)
	{
		first_csys_E = current_csys->next;
		free(current_csys);
		current_csys = first_csys_E;
	}
}


/************************ Function: displ_strain ****************************
* Calculates the total displacement and/or strain at a point due to ALL
* elements.
*
* In:	   calc_displ	- calculate displacements flag
*		   calc_strain  - calculate strains flag
*		   x			- coords (global) of pt at which to calc displ & strain
*		   omit_elt	    - element to omit when calculating displs (NULL = none)
*
* Out:	   displ  		- displacement vector (global coords)
*		   strain 		- strain tensor (global coords)
*
*****************************************************************************/
#if defined(__STDC__) || defined(ANSI) /* ANSI */
void     displ_strain(int calc_displ, int calc_strain, double x[3],double displ[3], double strain[3][3], elt_t *omit_elt, int elt_num)
#else
void displ_strain(calc_displ, calc_strain, x, displ, strain, omit_elt, elt_num)
int				calc_displ;
int				calc_strain;
double			x[3];
double			displ[3];
double			strain[3][3];
elt_t		    *omit_elt;
int elt_num;
#endif

{
	elt_t	*current_elt;
	double	elt_displ[3];
	double	elt_strain[3][3];
/*Declarations for correction of the "shadow effect" */
         int i;
         double orient;
      	double under_plane;
      	int under;
      	double normal [3];
      	double data1 [3];
      	double data [3];
      	vert_t *verta;
      	vert_t *vertb;
      	vert_t *vertc;
      	double seg1 [3];
      	double seg2 [3];
      	double inside;
      	double inside_vector [3];
      	int inside_test;
      	int inside_test_fin;
      	double x3_global [3];
      	disloc_seg_t *disloc_seg;
/* end of the declaration */

/*Declarations for coef_exclu */
         double dist;
/* end of the declaration */

/* Element number contributing to current element (elt_num) - for debugging*/
        int ielem;

	initialize_vector(displ,0.0);
	initialize_matrix(strain,0.0);

	/* Loop over each element
	-------------------------*/
	current_elt = first_elt_E;
   near_vertex_E = FALSE;
        ielem=0;
	while (current_elt != NULL)
   {
         ielem++;
		/* Test UNDER to determine if the data point is inside the "shadow zone"
		---------------------------------------------------------------------*/
		under = 0;
		inside_test_fin = 1;

		/* Determine if the element has a positive side up or down (orient)
		   and determine if the data is under the plane
			defined by the element (under_plane)
		-----------------------------------------------------------------*/
		disloc_seg = current_elt->disloc_seg;
		x3_global [0] = 0;
		x3_global [1] = 0;
		x3_global [2] = 1;
		verta = disloc_seg [0].vert[0];
		vertb = disloc_seg [0].vert[1];
		vertc = disloc_seg [1].vert[1];
		subtract_vectors (vertb->x, verta->x, seg1);
		subtract_vectors (vertc->x, verta->x, seg2);
		cross_product (seg1, seg2, normal);
		orient = dot_product (normal, x3_global);
		subtract_vectors (x, verta->x, data1);
		under_plane = dot_product (normal, data1);

      /* Determine if obs_pt is nearnest vertex verta, vertb or vertc
         1) compute d = mean length of the 3 disloc_seg
         2) compute distance from obs_pt to verta and see if this diance is < d*coef_exclu
         3) compute distance from obs_pt to vertb and see if this diance is < d*coef_exclu
         4) compute distance from obs_pt to vertc and see if this diance is < d*coef_exclu
         if condition is TRUE, break the loop over elements and set the flag near_vertex_E=TRUE
            and return.

         We use a new function "distance" define in matrix.h & matrix.c
      */
      dist = 0.0;
      dist  = distance(verta->x,vertb->x);
      dist += distance(vertb->x,vertc->x);
      dist += distance(vertc->x,verta->x);
      dist *= coef_exclu_E/3.0;

      if (distance(x,verta->x)<=dist || distance(x,vertb->x)<=dist || distance(x,vertc->x)<=dist)
      {
         near_vertex_E = TRUE;
         break;
      }

		/* Determine if the data point is in the "rigid body"
			(cf. explanation of the bug)
		---------------------------------------------------*/

		for (i = 0; i < current_elt->num_vertices; i++)
		{
			inside = 0;
			verta = disloc_seg [i].vert[0];
			vertb = disloc_seg [i].vert[1];
			subtract_vectors (vertb->x, verta->x, seg1);
		   subtract_vectors (x, verta->x, data);
			cross_product (seg1, data, inside_vector);
			inside = dot_product (x3_global, inside_vector);

			if (orient > 0)
			{
               if (inside > 0 && under_plane < 0)
                  inside_test = 1;
				   else
                  inside_test = 0;
			}
			if (orient < 0)
			{
            if (inside < 0 && under_plane >  0)
               inside_test = 1;
				else
               inside_test = 0;
			}
			if (orient == 0)
				under = 0;

			inside_test_fin *= inside_test;
		}

		/*Gives a value 1 to under for data points under the element
		and positive side up
		Gives a value 2  to under for data points under the element
		and positive side down
		Gives a value 0  to under for data points not under the element
		---------------------------------------------------------------*/
/* PRL CHANGE: Take into account when point x is in the plane, and plane is dipping,
               by checking for very small values of under_plane (that can produce
               false positive inside_test_fin).

               under_plane=sqrt(under_plane*under_plane);
               if (under_plane < BVERT_TINY) inside_test_fin = 0;
   END PRL CHANGE */

		if (inside_test_fin > 0)
      {
			if (orient > 0)
            under = 1;
		   if (orient < 0)
            under = 2;
                        /* PRL CHANGE: Take into account when point x is in the plane, and plane is dipping,
                                       by checking for very small values of under_plane (that can produce
                                       false positive inside_test_fin). */

                        under_plane=sqrt(under_plane*under_plane);
                        if (under_plane < BVERT_TINY)
            under = 0;
                        /* END PRL CHANGE */
		}
      else
         under = 0;

		/* Calculate displacement and strain due to this element
		--------------------------------------------------------*/
		displ_strain_poly_elt(calc_displ,calc_strain,current_elt,
			x,elt_displ,elt_strain,omit_elt,under,elt_num,ielem);

		/* Add the contribution of this element to the total displ & strain
		-----------------------------------------------------------------*/
		add_vectors(displ,elt_displ,displ);
		add_matrices(strain,elt_strain,strain);
		current_elt = current_elt->next;
	}

	/* Adjust strain for remote strain
	----------------------------------*/
   if (near_vertex_E!=TRUE)
	   add_matrices(strain,rem_strain_E,strain);
}


/****************** Function: displ_strain_poly_elt ********************
* Calculates the displacement and/or strain at a point due to a single
* polygonal element.
*
* In:	calc_displ	- calculate displacements flag
*		calc_strain - calculate strains flag
*		current_elt	- element being considered
*		x			- coords (global) of pt at which to calc displ & strain
*		omit_elt	- element to omit when calculating displs (NULL = none)
* Out:	elt_displ  	- displacement vector (global coords) due to this elt
*		elt_strain	- strain tensor (global coords) due to this elt
****************************************************************************/
void     displ_strain_poly_elt(int calc_displ, int calc_strain,elt_t *t_strain[3][3], elt_t *omit_elt, int under, int elt_num, int ielem)
{current_elt, double x[3], double elt_displ[3],double el
	double	displ_ic[3][3];
	double	strain_ic[3][3][3];
	int		i, j, k;
	/*Declarations for shadow effect correction*/
	double bglobal[3];
	double e[3][3];
	double eg[3][3];
	double bg[3][3];

	double value=0.0;
	/* Initialize displacement vector and strain tensor
	---------------------------------------------------*/
	initialize_vector(elt_displ,0.0);
	initialize_matrix(elt_strain,0.0);

	/* Calculate the displacement and strain influence coefficients
	---------------------------------------------------------------*/
	displ_strain_ics_poly_elt(calc_displ,calc_strain,current_elt,x,
		displ_ic,strain_ic,omit_elt,elt_num,ielem);

/*        if (elt_num == 77) fprintf(stderr," %lg %lg %lg  %lg\n %lg %lg %lg  %lg\n %lg %lg %lg  %lg\n",displ_ic[0][0],displ_ic[0][1],displ_ic[0][2],(*current_elt->b[0]),displ_ic[1][0],displ_ic[1][1],displ_ic[1][2],(*current_elt->b[1]),displ_ic[2][0],displ_ic[2][1],displ_ic[2][2],(*current_elt->b[2]));
*/

	/* Superpose the contribution from each burger's vector component
	-----------------------------------------------------------------*/
	for (i=0; i < 3; i++)
   {
		for (j=0; j < 3; j++)
      {
		 	elt_displ[i] += displ_ic[j][i] * (*current_elt->b[j]);
         	for (k=0; k < 3; k++)
				elt_strain[i][j] += strain_ic[k][i][j] * (*current_elt->b[k]);
		}
	}

 /* Shadow effect correction in case of data point under the element
---------------------------------------------------------------------------------*/

	if (under > 0)
   {
	   	bglobal[0]=0;
      	bglobal[1]=0;
      	bglobal[2]=0;

		e[0][0]=1;e[0][1]=0;e[0][2]=0;
		e[1][0]=0;e[1][1]=1;e[1][2]=0;
		e[2][0]=0;e[2][1]=0;e[2][2]=1;

		eg[0][0]=1;eg[0][1]=0;eg[0][2]=0;
		eg[1][0]=0;eg[1][1]=1;eg[1][2]=0;
		eg[2][0]=0;eg[2][1]=0;eg[2][2]=1;

		for (i=0; i < 3; i++)
      {
			/* Transform the burger vector's component in the global coordinate system
			--------------------------------------------------------------------------- */
			rotate_vector(INVERSE_ROT,current_elt->elt_csys.local_rot,e[i]);
			scalar_vector_mult (*current_elt->b[i], e[i],bg[i]);
		}
		for (i=0; i < 3; i++)
      {
		   for (j=0; j<3; j++)
            		bglobal[i]+=dot_product(bg[j], eg[i]);

			/*Corrects the displacement by the corresponding burger's vector (global coordsys)
			--------------------------------------------------------------------------------*/
			if (under == 1)
			   elt_displ[i] -= bglobal[i];
         else
			if (under == 2)
				elt_displ[i] += bglobal[i];
		}

	}/* end of if (under>0)*/
}



/***************** Function: displ_strain_ics_poly_elt **********************
* Calculates the displacement and/or strain influence coefficients at a point
* due to a polygonal element.
*
* In:	calc_displ	- calculate displacements flag
*		calc_strain - calculate strains flag
*		current_elt	- element being considered
*		x			- coords (global) of pt at which to calc displ & strain
*		omit_elt	- element to omit when calculating displs (NULL = none)
*
* Out:	displ_ic[i][j]		- the jth component of displ (global coords)
*                        	  due to a unit ith Burgers vector component
*                        	  (bc_coord_sys coords)
*		strain_ic[i][j][k]	- the jk component of strain (global coords)
*                        	  due to a unit ith Burgers vector component
*                        	  (bc_coord_sys coords)
****************************************************************************/
#if defined(__STDC__) || defined(ANSI) /* ANSI */
void     displ_strain_ics_poly_elt(int calc_displ, int calc_strain,elt_t *current_elt, double x[3], double displ_ic[3][3],double strain_ic[3][3][3], elt_t *omit_elt, int elt_num,int ielem)
#else
void displ_strain_ics_poly_elt(calc_displ, calc_strain, current_elt, x,
	displ_ic, strain_ic, omit_elt, elt_num,ielem)
int		calc_displ;
int		calc_strain;
int elt_num;
int ielem;
elt_t	*current_elt;
double	x[3];
double	displ_ic[3][3];
double	strain_ic[3][3][3];
elt_t	*omit_elt;
#endif

{
	int		i, j, k, l;
	int		seg;
	int		swap;
	vert_t	 *vert1;
	vert_t	 *vert2;
	double	depth1;
	double	depth2;
	double	beta;
	double	temp_double;
	double	temp_vector[3];
	double	r;
	double	z3;
	double	y1[3];
	double	y2[3];
	double	displ_ic1[3][3];
	double	displ_ic2[3][3];
	double	displ_ic3[3][3];
	double	displ_ic4[3][3];
	double	strain_ic1[3][3][3];
	double	strain_ic2[3][3][3];
	double	strain_ic3[3][3][3];
	double	strain_ic4[3][3][3];
	disloc_seg_t	 *disloc_seg;

        /* PRL change, dbug, if == 1, prints some stuff in comninou_displ_ics */
        int dbug;

        dbug=0;
        /*if (elt_num == 4 && (ielem == 1 || ielem == 15)) dbug=1; */

	/* Initialize the influence coefficients
	----------------------------------------*/
	for (i=0; i < 3; i++) {
		initialize_vector(displ_ic[i],0.0);
		initialize_matrix(strain_ic[i],0.0);
	}

	/* Loop over the element's dislocation segments
	-----------------------------------------------*/
	disloc_seg = current_elt->disloc_seg;
	for (seg = 0; seg < current_elt->num_vertices; seg++) {

                if (dbug==1) fprintf(stderr," Seg %i\n",seg);

		/* Determine the segment vertices
		---------------------------------*/
		vert1 = disloc_seg[seg].vert[0];
		vert2 = disloc_seg[seg].vert[1];

		/* Compute the vectors from the segments vertices to the
		  	data point
		---------------------------------------------------------*/
		subtract_vectors(x,vert1->x,y1);
		subtract_vectors(x,vert2->x,y2);

		/* Rotate vert-to-obs_point vectors to segment-local coords
		-----------------------------------------------------------*/
		rotate_vector(FORWARD_ROT,disloc_seg[seg].local_rot,y1);
		rotate_vector(FORWARD_ROT,disloc_seg[seg].local_rot,y2);

		depth1 = -vert1->x[2];
		depth2 = -vert2->x[2];
		beta   = PI/2.0 - disloc_seg[seg].plunge;

		if (((sqrt(y1[0]*y1[0]+y1[1]*y1[1]) < BVERT_TINY) && y1[2] >= 0.0) ||
			((sqrt(y2[0]*y2[0]+y2[1]*y2[1]) < BVERT_TINY) && y2[2] >= 0.0)) {
			below_vertex_E = TRUE;
			return;
		}

		/* If x lies along dipping leg of angular dislocations, swap
		   the vertex order, so singularity will be avoided
		------------------------------------------------------------*/
		swap = FALSE;
		z3 = y1[0]*sin(beta) + y1[2]*cos(beta);
		r  = vector_magnitude(y1);
		if ((r - z3) < SWAP_TINY) {
			swap = TRUE;
			copy_vector(y2,temp_vector);
			copy_vector(y1,y2);
			copy_vector(temp_vector,y1);
			for (i=0; i < 2; i++) {
				y1[i] *= -1.0;
				y2[i] *= -1.0;
			}
			temp_double = depth2;
			depth2 = depth1;
			depth1 = temp_double;
			beta = PI - beta;
                        if (elt_num==77) fprintf(stderr," swap is TRUE, beta = %lg\n",beta);
		}

		/* Be careful for the calculation of vertical elements
		------------------------------------------------------*/
		if (beta==0.0) {
			beta = 1.0e-14;
                        /*if (elt_num == 77) fprintf(stderr," beta set to: %lg\n",beta);*/
                }


		/* Calculate displacement influence coeffs
		------------------------------------------*/
		if (calc_displ) {

			/* Avoid displacement discontinuity when calculating displ
			   inf coeff of an element on itself
			----------------------------------------------------------*/
			if  (current_elt != omit_elt) {

				/*comninou_displ_ics(y1,depth1,beta,psn_ratio_E,
					half_space_E,displ_ic1,dbug);*/
                 comninou_displ_ics(y1,depth1,beta,psn_ratio_E,
					half_space_E,displ_ic1);
				/*comninou_displ_ics(y2,depth2,beta,psn_ratio_E,
					half_space_E,displ_ic2,dbug);*/
                 comninou_displ_ics(y2,depth2,beta,psn_ratio_E,
					half_space_E,displ_ic2);
                                /*if (elt_num == 80) {
                                     fprintf(stderr," y1: %lg %lg %lg,  depth1: %lg\n", y1[0],y1[1],y1[2],depth1);
                                     fprintf(stderr," displ_ic1: %lg %lg %lg\n      %lg %lg %lg\n %lg %lg %lg\n", displ_ic1[0][0], displ_ic1[0][1], displ_ic1[0][2], displ_ic1[1][0], displ_ic1[1][1], displ_ic1[1][2], displ_ic1[2][0], displ_ic1[2][1], displ_ic1[2][2]);
                                     fprintf(stderr," y2: %lg %lg %lg,  depth2: %lg\n", y2[0],y2[1],y2[2],depth2);
                                     fprintf(stderr," displ_ic2: %lg %lg %lg\n      %lg %lg %lg\n %lg %lg %lg\n", displ_ic2[0][0], displ_ic2[0][1], displ_ic2[0][2], displ_ic2[1][0], displ_ic2[1][1], displ_ic2[1][2], displ_ic2[2][0], displ_ic2[2][1], displ_ic2[2][2]);

                                }*/
				/* Superpose the angular dislocation influence coeffs into
		   	   	a dislocation segment influence coeff
				----------------------------------------------------------*/
				subtract_matrices(displ_ic1,displ_ic2,
					displ_ic3);
                                /*if (dbug == 1) {
                                     fprintf(stderr," displ_ic3: \t%lg \t%lg \t%lg\n \t\t%lg \t%lg \t%lg\n \t\t%lg \t%lg \t%lg\n", displ_ic3[0][0], displ_ic3[0][1], displ_ic3[0][2], displ_ic3[1][0], displ_ic3[1][1], displ_ic3[1][2], displ_ic3[2][0], displ_ic3[2][1], displ_ic3[2][2]);
                                }*/

				/* Swap the vertices back to proper order (if necessary)
				--------------------------------------------------------*/
				if (swap) {
					scalar_vector_mult(-1.0,displ_ic3[2],
						displ_ic3[2]);
					for (i=0; i < 3; i++) {
						displ_ic3[i][0] *= -1.0;
						displ_ic3[i][1] *= -1.0;
					}
				}
                                if (dbug == 1) {
                                     fprintf(stderr," displ_ic3: \t%lg \t%lg \t%lg\n \t\t%lg \t%lg \t%lg\n \t\t%lg \t%lg \t%lg\n", displ_ic3[0][0], displ_ic3[0][1], displ_ic3[0][2], displ_ic3[1][0], displ_ic3[1][1], displ_ic3[1][2], displ_ic3[2][0], displ_ic3[2][1], displ_ic3[2][2]);
                                     fprintf(stderr," disloc_seg[%i].elt_b: \t%lg \t%lg \t%lg\n \t\t\t%lg \t%lg \t%lg\n \t\t\t%lg \t%lg \t%lg\n", seg,disloc_seg[seg].elt_b[0][0], disloc_seg[seg].elt_b[0][1], disloc_seg[seg].elt_b[0][2], disloc_seg[seg].elt_b[1][0], disloc_seg[seg].elt_b[1][1], disloc_seg[seg].elt_b[1][2], disloc_seg[seg].elt_b[2][0], disloc_seg[seg].elt_b[2][1], disloc_seg[seg].elt_b[2][2]);
                                }

				/* Transform from disloc segment to element influence coeffs
				------------------------------------------------------------*/
				for (i=0; i < 3; i++) {
					initialize_vector(displ_ic4[i],0.0);
					for (j=0; j < 3; j++) {
						for (k=0; k < 3; k++) {
							displ_ic4[i][j] +=
								disloc_seg[seg].elt_b[i][k] *
								displ_ic3[k][j];
						}
					}
				}
                                if (dbug == 1) {
                                     fprintf(stderr," displ_ic4: \t%lg \t%lg \t%lg\n \t\t%lg \t%lg \t%lg\n \t\t%lg \t%lg \t%lg\n", displ_ic4[0][0], displ_ic4[0][1], displ_ic4[0][2], displ_ic4[1][0], displ_ic4[1][1], displ_ic4[1][2], displ_ic4[2][0], displ_ic4[2][1], displ_ic4[2][2]);
                                }


				/* Rotate from C&D to global coordinates
				----------------------------------------*/
				for (i=0; i < 3; i++) {
					rotate_vector(INVERSE_ROT,disloc_seg[seg].local_rot,
						displ_ic4[i]);
				}
                                if (dbug == 1) {
                                     fprintf(stderr," G displ_ic4: \t%lg \t%lg \t%lg\n \t\t%lg \t%lg \t%lg\n \t\t%lg \t%lg \t%lg\n", displ_ic4[0][0], displ_ic4[0][1], displ_ic4[0][2], displ_ic4[1][0], displ_ic4[1][1], displ_ic4[1][2], displ_ic4[2][0], displ_ic4[2][1], displ_ic4[2][2]);
                                }


				/* Superpose the contribution of this dislocation segment
				---------------------------------------------------------*/
				for (i=0; i < 3; i++) {
					add_vectors(displ_ic[i],displ_ic4[i],
						displ_ic[i]);
				}
                                if (dbug == 1) {
                                     fprintf(stderr," displ_ic: \t%lg \t%lg \t%lg\n \t\t%lg \t%lg \t%lg\n \t\t%lg \t%lg \t%lg\n", displ_ic[0][0], displ_ic[0][1], displ_ic[0][2], displ_ic[1][0], displ_ic[1][1], displ_ic[1][2], displ_ic[2][0], displ_ic[2][1], displ_ic[2][2]);
                                }


			}
		}

		/* Calculate strain influence coeffs
		------------------------------------*/
		if (calc_strain) {

			comninou_strain_ics(y1,depth1,beta,psn_ratio_E,
				half_space_E,strain_ic1);
			comninou_strain_ics(y2,depth2,beta,psn_ratio_E,
				half_space_E,strain_ic2);

			/* Superpose the angular dislocation influence coeffs into
		   	   a dislocation segment influence coeff
			----------------------------------------------------------*/
			for (i=0; i < 3; i++) {
				subtract_matrices(strain_ic1[i],strain_ic2[i],
					strain_ic3[i]);
			}

			/* Swap the vertices back to proper order (if necessary)
			--------------------------------------------------------*/
			if (swap) {
				scalar_matrix_mult(-1.0,strain_ic3[2],
					strain_ic3[2]);
				for (i=0; i < 3; i++) {
					strain_ic3[i][0][2] *= -1;
					strain_ic3[i][2][0] *= -1;
					strain_ic3[i][1][2] *= -1;
					strain_ic3[i][2][1] *= -1;
				}
			}


			for (i=0; i < 3; i++) {
				initialize_matrix(strain_ic4[i],0.0);
				for (j=0; j < 3; j++) {
					for (k=0; k < 3; k++) {
						for (l=0; l < 3; l++) {
							strain_ic4[i][j][k] +=
								disloc_seg[seg].elt_b[i][l] *
								strain_ic3[l][j][k];
						}
					}
				}
			}

			/* Rotate strain inf coeffs to global coords
			--------------------------------------------*/
			for (i=0; i < 3; i++) {
				rotate_tensor(INVERSE_ROT,disloc_seg[seg].local_rot, strain_ic4[i]);
			}

			/* Superpose the contribution of this dislocation segment
			---------------------------------------------------------*/
			for (i=0; i < 3; i++) {
				add_matrices(strain_ic[i],strain_ic4[i],
					strain_ic[i]);
			}

		}

	} /* loop over dislocation segments */

}


/********************** Function: print_obs_grid_data ******************
* Loops through linked list of observation grids, calculating and printing
* the requested displacement, strain, and stress data for each to the
* output file.
************************************************************************/
/* AMB. Include strain and stress. However, PSTRAIN and PSTRESS are still not
   supported. */
void print_obs_grid_data(double *poutput, int numobs, double mu)
{
  obs_grid_t	*current_obs_grid;
  double	x[3];
  double	dx[3];
  double	displ[3];
  double	strain[3][3];
  double	*begin;
  double	*end;
  int		*numpts;
  char	        error_msg[MAX_ERROR_MSG];
  int		*print;
  int		calc_displ;
  int		calc_strain;
  int		calc_stress;
  int           izero;
  int           i, j, k, nrows, os;
  double        lambda, theta, nu;
  
  /* local counter */
  int count=0;

#ifdef TV
  struct timeval t1, t2;
  double ac = 0.0;
#endif

  /* for stress */
  nu = psn_ratio_E;
  lambda = 2.0*mu*nu/(1.0 - 2.0*nu);

  izero=0;
  /* Loop over each observation grid
     ----------------------------------*/
  current_obs_grid = first_obs_grid_E;
  while (current_obs_grid != NULL) {
    
    /* Determine data required by obs grid print options
       ----------------------------------------------------*/
    print = current_obs_grid->print;
    calc_displ = print[DISPL];
    calc_stress = print[STRESS] || print[PSTRESS];
    calc_strain = calc_stress || print[STRAIN] || print[PSTRAIN];
    
    /* Remenber that begin_pt & end_pt are in GLOBAL coordinate system.*/
    begin  = current_obs_grid->begin;
    end    = current_obs_grid->end;
    numpts = current_obs_grid->numpts;
    subtract_vectors(end,begin,dx);
    
    /* Process the observation grid
       -------------------------------*/
    switch (current_obs_grid->dimension) {
      
    case 0:
      copy_vector(begin,x);
      /* Transform current observation point into global CSys
	 ----------------------------------------------------*/
      transform_position_vector(INVERSE_TRANS,
				current_obs_grid->endpt_csys->origin,
                                      current_obs_grid->endpt_csys->local_rot,
				x);
#ifdef TV
      gettimeofday(&t1, 0);
#endif
      displ_strain(calc_displ,calc_strain,x,displ,strain,
		   NULL, izero);
#ifdef TV
      gettimeofday(&t2, 0);
      ac += difftime(&t1, &t2);
#endif
      /*** in displ_strain, near_vertex_E is set ***/

      nrows = 0;
      if (calc_displ)  nrows += 3;
      if (calc_strain) nrows += 9;
      if (calc_stress) nrows += 6;
      os = 0;
      if (calc_displ) {
	if (near_vertex_E == 0) {
	  for (i = 0; i < 3; i++) poutput[count*nrows + i] = displ[i];
	} else {
	  for (i = 0; i < 3; i++) poutput[count*nrows + i] = null_value_E;
	}
	os += 3;
      }
      if (calc_strain) {
	if (near_vertex_E == 0) {
	  for (i = 0, k = 0; i < 3; i++)
	    for (j = 0; j < 3; j++, k++)
	      poutput[count*nrows + os + k] = strain[i][j];
	} else {
	  for (i = 0; i < 9; i++) poutput[count*nrows + os + i] = null_value_E;
	}
	os += 9;
      }
      if (calc_stress) {
	if (near_vertex_E == 0) {
	  theta = 0.0;
	  for (i = 0; i < 3; i++) theta += strain[i][i];
	  poutput[count*nrows + os + 0] = lambda*theta + 2.0*mu*strain[0][0];
	  poutput[count*nrows + os + 1] = mu*(strain[1][0] + strain[0][1]);
	  poutput[count*nrows + os + 2] = mu*(strain[2][0] + strain[0][2]);
	  poutput[count*nrows + os + 3] = lambda*theta + 2.0*mu*strain[1][1];
	  poutput[count*nrows + os + 4] = mu*(strain[2][1] + strain[1][2]);
	  poutput[count*nrows + os + 5] = lambda*theta + 2.0*mu*strain[2][2];
	} else {
	  for (i = 0; i < 6; i++) poutput[count*nrows + os + i] = null_value_E;
	}
      }
      break;
      
    default:
      sprintf(error_msg,
	      "Invalid dimension (%d) for observation grid",
	      current_obs_grid->dimension);
      p_error(error_msg,NULL);
      
    }
    
    
    current_obs_grid = current_obs_grid->next;
    count++;
    near_vertex_E= FALSE;
    
  } /*while*/

#ifdef TV
  printf("ac = %f\n", ac);
#endif
}



/************************** Function: find_csys *************************
* Returns a pointer to the coordinate system named by name, or NULL if no
* such coordinate system exists.
*
* In:	name	- name of the coordinate system to find
*****************************************************************************/
#if defined(__STDC__) || defined(ANSI) /* ANSI */
csys_t*  find_csys(const char *name)
#else
csys_t *find_csys(name)
const char *name;
#endif

{
	csys_t *current_csys;

	current_csys = first_csys_E;

	while (current_csys != NULL) {
		if (!strcmp(name,current_csys->name))
			break;
		current_csys = current_csys->next;
	}

	return(current_csys);
}


/***************************** Function: find_vert ***************************
* Return a pointer to the vertex named by name, or NULL if no such vertex
* exists.
*
* In:	name	- name of the vertex to find
*****************************************************************************/
#if defined(__STDC__) || defined(ANSI) /* ANSI */
vert_t*  find_vert(const char *name)
#else
vert_t *find_vert(name)
const char *name;
#endif

{
	vert_t	*current_vert;

	current_vert = first_vert_E;

	while (current_vert != NULL) {
		if (!strcmp(name,current_vert->name))
			break;
		current_vert = current_vert->next;
	}

	return(current_vert);
}


/***************************** Function: p_error *****************************
* Prints an error message to stderr and calls exit().  If line != NULL, the
* line and line number (from the input file) on which the error occured are
* printed as well.
*
* In:	error_msg	- error message to print
* 		line		- input file line number at which the error occurred
*					  (NULL = N/A)
*****************************************************************************/
#if defined(__STDC__) || defined(ANSI) /* ANSI */
void     p_error(const char *error_msg, const char *line)
#else
void p_error(error_msg, line)
const char	*error_msg;
const char	*line;
#endif

{
	fprintf(stderr,"\nerror: %s",error_msg);
	if (line == NULL) {
		fprintf(stderr,"\n");
	} else {
		fprintf(stderr," (%s, line %d)\n",infile_E,linenum_E);
		fprintf(stderr,"       %s\n",line);
	}
	exit(1);

}

/***************************** Function: display_msg *****************************
* Prints an message to stderr.
*
* In:	error_msg	- error message to print
*****************************************************************************/
#if defined(__STDC__) || defined(ANSI) /* ANSI */
void     display_msg(char *_msg)
#else
void display_msg(_msg)
char	*_msg;
#endif

{
	fprintf(stderr,"%s\n",_msg);
}


/************************ Function: setup_global_coords **********************
* Defines the global coordinate system, making it the first member in the
* linked list of coordinate systems (first_csys_E).
*****************************************************************************/
#if defined(__STDC__) || defined(ANSI) /* ANSI */
void     setup_global_coords(void)
#else
void setup_global_coords()
#endif

{

	int		i;

	/* Set first coord system to global coordinates
	-----------------------------------------------*/
	first_csys_E = (csys_t *) calloc((size_t) 1,
		sizeof(csys_t));
	if (!first_csys_E)
		p_error("Cannot allocate memory (calloc) for global coord sys",
			NULL);

	first_csys_E->name = (char *) malloc((size_t)
					     strlen(GLOBAL_NAME)+1);
	if (!first_csys_E->name)
		p_error( "Cannot allocate memory for global coord system name",
			NULL);
	strcpy(first_csys_E->name,GLOBAL_NAME);

	for (i=0; i < 3; i++) {
		first_csys_E->origin[i] = 0;
	}
	initialize_matrix(first_csys_E->local_rot,0.0);
	for (i=0; i < 3; i++) {
		first_csys_E->local_rot[i][i] = 1.0;
	}

}


/********************** Function: read_observation_grids ********************/
/* AMB. Include strain and stress. */
int read_observation_grids(double *obs_xyz, int numobs, char *calc_code)
{
	int		numwords;
	const char	*word[MAXWORDS];
	int		i;
	int     j;
	int		dimension;
	char	error_msg[MAX_ERROR_MSG];
	obs_grid_t	 *current_obs_grid;
	int		num_ones;
	int		numpts;
	char	temp_char;

	/* Read in observation grids
	----------------------------*/
	for (j=0; j < numobs ; j++) {

               /* Get grid dimension
		---------------------*/
		dimension = 0;
		numwords = 9;
		word[0]="ObsGrid";
		word[1]="0";
		word[2] = calc_code;
		word[3]="global";
		word[4]="global";
		word[5]="global";

		/* Allocate memory for observation grid
		---------------------------------------*/
		if (first_obs_grid_E == NULL) {
			first_obs_grid_E = (obs_grid_t *)
				calloc((size_t) 1,sizeof(obs_grid_t));
			if (!first_obs_grid_E)
				p_error("Cannot allocate memory (calloc) for obs grid",
					NULL);
			current_obs_grid = first_obs_grid_E;
		} else {
			current_obs_grid->next = (obs_grid_t *)
				calloc((size_t) 1,sizeof(obs_grid_t));
			if (!current_obs_grid->next)
				p_error("Cannot allocate memory (calloc) for obs grid",
					NULL);
			current_obs_grid = current_obs_grid->next;
		}

		/* Set the observation grid dimension
		-------------------------------------*/
		current_obs_grid->dimension = dimension;

		/* Get observation grid name
		----------------------------*/
		current_obs_grid->name = (char *) malloc((size_t)
			strlen(word[OG_NAME_POS])+1);
		if (!current_obs_grid->name)
			p_error("Cannot allocate memory for observation grid name",
			 NULL);
		strcpy(current_obs_grid->name,word[OG_NAME_POS]);

		/* Get the print options
		------------------------*/
		i = 0;
		current_obs_grid->print[DISPL]   = FALSE;
		current_obs_grid->print[STRAIN]  = FALSE;
		current_obs_grid->print[STRESS]  = FALSE;
		current_obs_grid->print[PSTRAIN] = FALSE;
		current_obs_grid->print[PSTRESS] = FALSE;
		while ((temp_char = word[OG_PRINT_OPS_POS][i]) != '\0') {
			switch (temp_char) {
				case DISPL_CHAR:
					current_obs_grid->print[DISPL] = TRUE;
					break;
				case STRAIN_CHAR:
					current_obs_grid->print[STRAIN] = TRUE;
					break;
				case STRESS_CHAR:
					current_obs_grid->print[STRESS] = TRUE;
					break;
				case PRINCIPAL_CHAR:
					i++;
					switch (word[OG_PRINT_OPS_POS][i]) {
						case STRAIN_CHAR:
							current_obs_grid->print[PSTRAIN] = TRUE;
							break;
						case STRESS_CHAR:
							current_obs_grid->print[PSTRESS] = TRUE;
							break;
						default:
							p_error("Invalid observation grid print option", NULL);
					}
					break;
				default:
					p_error("Invalid observation grid print option",NULL);
			}
			i++;
		}


		/* Get the input coordinate system
		----------------------------------*/
		if ((current_obs_grid->endpt_csys =
			find_csys(word[OG_INPUT_CSYS_POS])) == NULL) {
			p_error("Undefined coordinate system",NULL);
		}

		/* Get the observation point coordinate system
		----------------------------------------------*/
		if ((current_obs_grid->obspt_csys =
			find_csys(word[OG_OBSPT_CSYS_POS])) == NULL) {
			p_error("Undefined coordinate system",NULL);
		}

		/* Get the output coodinate system
		-----------------------------------*/
		if ((current_obs_grid->outp_csys =
			find_csys(word[OG_DATA_CSYS_POS])) == NULL) {
			p_error("Undefined coordinate system",NULL);
		}

		/* Get the beginning & ending coordinates
		-----------------------------------------*/
		for (i=0; i < 3; i++) {
                        current_obs_grid->begin[i] = obs_xyz[j*3 + i];
		}

		/* Get number of points along each coord axis in grid
		----------------------------------------------------*/
		num_ones = 0;
		for (i=0; i < ((dimension == 2) ? 3 : dimension); i++) {
			current_obs_grid->numpts[i] = numpts =
				atoi(word[OG_NUMPTS_POS+i]);
			if (numpts < 2) {
				if (dimension == 2 && numpts == 1) {
						num_ones++;
				} else {
					sprintf(error_msg,
						"%dD observation grid axes require %d or more points",
						dimension,((dimension == 2) ? 1 : 2));
					p_error(error_msg,NULL);
				}
			}
		}
		if (dimension == 2 && num_ones != 1) {
			p_error("1 (& only 1) 2D observation grid axis requires 1 point",NULL);
		}

      /*             98-12-09
		   Do not convert the begin & end pts to global coords !!!!!
         So, no transformation in global coordinate system is required
		*/

      /*
		transform_position_vector(INVERSE_TRANS,
			current_obs_grid->endpt_csys->origin,
			current_obs_grid->endpt_csys->local_rot,
			current_obs_grid->begin);
		transform_position_vector(INVERSE_TRANS,
			current_obs_grid->endpt_csys->origin,
			current_obs_grid->endpt_csys->local_rot,
			current_obs_grid->end);
      */

	}

	return(0);
}


/* function read_objs_elts_verts */
int  read_objs_elts_verts(double *mvert_xyz, double *displist, double *bbb_bc, int numvert, int numdispl)
{
	int     j,i;
	int     numwords;
	char    buffer[10];
	vert_t *current_vert;
	obj_t  *current_obj;
	elt_t  *current_elt;

	/* Read in vertices
	--------------------------------*/
	for (j=0; j<numvert; j++) {

	       sprintf(buffer,"%d",j+1);
	     /* Allocate memory for vertex*/
		if (first_vert_E == NULL) {
			first_vert_E = (vert_t *) calloc((size_t) 1,sizeof(vert_t));
			if (!first_vert_E)
				p_error("Cannot allocate memory (calloc) for vertex",NULL);
			current_vert = first_vert_E;
		} else {
			current_vert->next = (vert_t *)
			calloc((size_t) 1,sizeof(vert_t));
			if (!current_vert->next)
			p_error("Cannot allocate memory (calloc) for vertex",buffer);
			current_vert = current_vert->next;
		}
		/* Get vertex name
		------------------*/
		current_vert->name = (char *) malloc((size_t) strlen(buffer)+1);
		if (!current_vert->name)
			p_error("Cannot allocate memory for vertex name",
			buffer);
		strcpy(current_vert->name,buffer);

	  /* Get the coodinate system
		----------------------------*/
		sprintf(buffer,"%s","global");
		if ((current_vert->csys = find_csys(buffer))
			== NULL) {
			p_error("Undefined coordinate system", buffer);
		}
		/* Read the vertex position
		---------------------------*/
		for (i=0; i < 3; i++)
		current_vert->x[i] = mvert_xyz[j*3 + i];
		/* Transform vertex position vector to global coords
		----------------------------------------------------*/
		transform_position_vector(INVERSE_TRANS,
		current_vert->csys->origin,
		current_vert->csys->local_rot,
		current_vert->x);
	}
	/* Read object info */
	/* Allocate memory for object */
	if (first_obj_E == NULL) {
		first_obj_E = (obj_t *) calloc((size_t) 1,sizeof(obj_t));
		if (!first_obj_E)
			p_error("Cannot allocate memory (calloc) for object",
			NULL);
		current_obj = first_obj_E;
	} else {
		current_obj->next = (obj_t *)calloc((size_t) 1,sizeof(obj_t));
		if (!current_obj->next)
			p_error("Cannot allocate memory (calloc) for object",
			NULL);
		current_obj = current_obj->next;
	}

	/* Get object name
	------------------*/
	sprintf(buffer,"%s","fault");
	current_obj->name = (char *) malloc((size_t)strlen(buffer)+1);
	if (!current_obj->name)
		p_error("Cannot allocate memory for object name",NULL);
	strcpy(current_obj->name,buffer);


		/* Get the coordinate system for element positions
		--------------------------------------------------*/
	sprintf(buffer,"%s","global");
		if ((current_obj->pos_csys =
			find_csys(buffer)) == NULL) {
			p_error("Undefined coordinate system",NULL);
		}
	current_obj->first_elt = NULL;

	/* Read element info */
   	get_elt_info(&current_elt,current_obj, displist, bbb_bc, numdispl);


	/* Make sure last object contains element(s)
	--------------------------------------------*/
	if (current_obj->first_elt == NULL)
		p_error("Last object contains no elements",NULL);

	/*return(0);*/
}

/************************* Function: get_elt_info ************************
* Reads element info and sets up new member in linked list of elements.
* In/Out:	current_elt		- the current (last) element in the linked list
*			current_obj		- the object to which current_elt belongs
*****************************************************************************/
void     get_elt_info(elt_t **current_elt, obj_t *current_obj, double *displist,
					  double *bbb_bc, int numdispl)
{
	int		num_vertices;
	int		i, j;
	int     k;
	const char   *word[MAXWORDS];
	char   buffer1[10],buffer2[10],buffer3[10];
	vert_t	 *vert1;
	vert_t	 *vert2;
	double  dx[3];
	double	trend;
	double	plunge;
	double	rot2[3][3];
	double	rot1[3][3];
	double	trac_bc_adjust[3];
	int		seg;
	vert_t	*vert[3];
	double	vert_x[3];
	double		vector1[3];
	double		vector2[3];
	double	normal_vector[3];
	double		x[3][3];
	double		global_x[3][3];
	disloc_seg_t	*disloc_seg;

        int temp;


	for(k=0; k<numdispl; k++)
	{
	    /* Increment num_elts_E */
		num_elts_E++;
		num_vertices = 3;
		/* Allocate memory for this element
		-----------------------------------*/
		if (first_elt_E == NULL) {
			first_elt_E = (elt_t *) malloc(sizeof(elt_t));
			if (!first_elt_E)
				p_error("Cannot allocate memory for element",NULL);
			*current_elt = first_elt_E;
			(*current_elt)->next = NULL;
		} else {
			(*current_elt)->next = (elt_t *) malloc(sizeof(elt_t));
			if (!(*current_elt)->next)
				p_error("Cannot allocate memory for element",NULL);
			*current_elt = (*current_elt)->next;
			(*current_elt)->next = NULL;
		}

		/* Set object pointers to first and last elements */
		if (first_obj_E == NULL)
			p_error("No objects defined. Element must be part of an object",NULL);
		if (current_obj->first_elt == NULL)
			current_obj->first_elt = *current_elt;
		current_obj->last_elt = *current_elt;

		/* Set the number of vertices */
		(*current_elt)->num_vertices = num_vertices;

		/* Allocate memory for dislocation segment array */
		(*current_elt)->disloc_seg = (disloc_seg_t *)
		calloc((size_t) num_vertices, sizeof(disloc_seg_t));
		if (!(*current_elt)->disloc_seg)
			p_error("Cannot allocate memory for dislocation segment array",	NULL);

		/* Get dislocation segment vertices*/
		for (i = 0; i < num_vertices; i++) {
			/********  setup word **********/
                        temp = displist[k*3];
			sprintf(buffer1,"%d",temp);
			word[E_VERTEX_POS] = buffer1;
                        /*printf("%f %s\n",displist[k*3],buffer1); */
                        temp = displist[k*3+1];
			sprintf(buffer2,"%d",temp);
			word[E_VERTEX_POS+1] = buffer2;
                        /*printf("%f %s\n",displist[k*3+1],buffer2);*/
                        temp = displist[k*3+2];
                        sprintf(buffer3,"%d",temp);
			word[E_VERTEX_POS+2] = buffer3;
                       /*printf("%f %s\n",displist[k*3+2],buffer3);*/

            /**********************************/
			j = (i == (num_vertices-1)) ? 0 : i+1;
			if (i != 0) {
				(*current_elt)->disloc_seg[i].vert[0] =
					(*current_elt)->disloc_seg[i-1].vert[1];
			} else {
				if (((*current_elt)->disloc_seg[i].vert[0] =
					find_vert(word[E_VERTEX_POS+i])) == NULL) {
					p_error("Undefined vertex",NULL);
				}
			}
			if (((*current_elt)->disloc_seg[i].vert[1] =
				find_vert(word[E_VERTEX_POS+j])) == NULL) {
				p_error("Undefined vertex",NULL);
			}
		}

		/* Inititalize rotation matrices
		--------------------------------*/
		initialize_matrix(rot2,0.0);
		initialize_matrix(rot1,0.0);
		initialize_vector((*current_elt)->elt_csys.origin,0.0);

		/* Loop over the dislocation segments */
		disloc_seg = (*current_elt)->disloc_seg;
		for (seg = 0; seg < (*current_elt)->num_vertices; seg++) {
		/* Calculate this dislocation segment's first vertex's contribution
		   to the element center */
			for (i=0; i < 3; i++) {
				(*current_elt)->elt_csys.origin[i] +=
					disloc_seg[seg].vert[0]->x[i] / (*current_elt)->num_vertices;
			}

		/* Determine the two vertices for this dislocation segment
		----------------------------------------------------------*/
		vert1 = disloc_seg[seg].vert[0];
		vert2 = disloc_seg[seg].vert[1];

		/* Calculate trend and plunge of this dislocation segment
		--------------------------------------------------------*/
		subtract_vectors(vert2->x,vert1->x,dx);
		trend = PI/2.0 - safe_atan2(dx[1],dx[0]);
		disloc_seg[seg].trend = trend;
		plunge = -safe_atan(dx[2],sqrt(dx[0]*dx[0] + dx[1]*dx[1]));
		disloc_seg[seg].plunge = plunge;


		/* Calculate the segment-local (Comninou & Dunders) to global
		   coordinates rotation matrix
		------------------------------------------------------------*/
		rot2[0][0] = 1.0;
		rot2[1][1] = rot2[2][2] = -1.0;
		rot1[0][0] = rot1[1][1] = sin(trend);
		rot1[1][0] = -(rot1[0][1] = cos(trend));
		rot1[2][2] = 1.0;
		matrix_mult(rot2,rot1,disloc_seg[seg].local_rot);

	} /*for dislocation segments */

	/* Compute element's local coordinate system
	--------------------------------------------*/

	/* AMB: Handle this name like all the others. */
	(*current_elt)->elt_csys.name =
	  (char *) malloc((size_t) strlen(ELT_CSYS_NAME)+1);
	if (!(*current_elt)->elt_csys.name)
	  p_error("Cannot allocate memory for elt-local coord sys name",
		  NULL);
	strcpy((*current_elt)->elt_csys.name, ELT_CSYS_NAME);

	for (i=0; i < 3; i++) {
		for (j=0; j < 3; j++) {
			global_x[i][j] = (i == j) ? 1.0 : 0.0;
		}
	}
	for (i=0; i < 3; i++) {
		vert[i] = disloc_seg[(i*(*current_elt)->num_vertices)/3].vert[0];
	}
	subtract_vectors(vert[1]->x,vert[0]->x,vector1);
	subtract_vectors(vert[2]->x,vert[0]->x,vector2);
	normalize_vector(vector1);
	normalize_vector(vector2);

	cross_product(vector1,vector2,x[2]);
	if (vector_magnitude(x[2]) < TINY_ANGLE)
		p_error("Cannot calc element normal. Elt must have a very odd shape.",NULL);
	normalize_vector(x[2]);
	cross_product(global_x[2],x[2],x[1]);
	if (vector_magnitude(x[1]) < TINY_ANGLE)
		copy_vector(global_x[1],x[1]);
	normalize_vector(x[1]);
	cross_product(x[1],x[2],x[0]);
	normalize_vector(x[0]);

	for (i=0; i < 3; i++) {
		for (j=0; j < 3; j++) {
			(*current_elt)->elt_csys.local_rot[i][j] =
			dot_product(x[i],global_x[j]);
		}
	}

	/* Check that all vertices are co-planar
	----------------------------------------*/
	for (seg=0; seg < num_vertices; seg++) {
		copy_vector((*current_elt)->disloc_seg[seg].vert[0]->x,vert_x);
		transform_position_vector(FORWARD_TRANS,
			(*current_elt)->elt_csys.origin,
			(*current_elt)->elt_csys.local_rot,vert_x);
		if (fabs(vert_x[2]) >
			fabs(sqrt(vert_x[0]*vert_x[0]+vert_x[1]*vert_x[1])/COPLANAR_LIMIT))
			p_error("Vertices are not co-planar",NULL);
	}

	/* Get the bc coodinate system
	------------------------------*/
	word[E_BC_CSYS_POS] = "elocal";
	if (!strcmp(word[E_BC_CSYS_POS],ELT_CSYS_NAME)) {
		(*current_elt)->bc_csys = &((*current_elt)->elt_csys);
	} else if (((*current_elt)->bc_csys =
		find_csys(word[E_BC_CSYS_POS])) == NULL) {
		p_error("Undefined coordinate system", NULL);
	}

	/* Read the boundary condition types
	------------------------------------*/
	word[E_BC_TYPE_POS]="bbb";

	for (i=0; i < 3; i++) {
		if (word[E_BC_TYPE_POS][i] == BVECTOR_CHAR)
			(*current_elt)->bc_type[i] = BVECTOR_BC;
		else if (word[E_BC_TYPE_POS][i] == TRACTION_CHAR)
			(*current_elt)->bc_type[i] = TRACTION_BC;
		else {
			p_error("Invalid boundary condition type",
				NULL);
		}
	}

	/* Calculate adjustment to traction BCs due to remote stresses
	--------------------------------------------------------------*/
	scalar_vector_mult(-1.0,x[2],normal_vector);
	cauchy(rem_stress_E,normal_vector,trac_bc_adjust);
	rotate_vector(FORWARD_ROT,(*current_elt)->bc_csys->local_rot,trac_bc_adjust);

	/* Read the boundary condition values
	-------------------------------------*/
	for (i=0; i < 3; i++) {
		(*current_elt)->bc[i] = bbb_bc[k*3+i];

		/* Adjust traction BC components for remote stress
		--------------------------------------------------*/
		if ((*current_elt)->bc_type[i] == TRACTION_BC)
			(*current_elt)->bc[i] -= trac_bc_adjust[i];
	}

	/* Calculate the projection of each unit component of the element
	   burgers vector into a segment-local coordinates burgers vector
	   (for calculating influence coefficients)
	-----------------------------------------------------------------*/
	for (seg = 0; seg < (*current_elt)->num_vertices; seg++) {
		for (i=0; i < 3; i++) {
			initialize_vector(disloc_seg[seg].elt_b[i],0.0);
			disloc_seg[seg].elt_b[i][i] = -1.0;
			rotate_vector(INVERSE_ROT,
				(*current_elt)->bc_csys->local_rot,
				disloc_seg[seg].elt_b[i]);
			rotate_vector(FORWARD_ROT,disloc_seg[seg].local_rot,
				disloc_seg[seg].elt_b[i]);
		}
	}
} /*end of for on numdispl */
}


/******************** Function: determine_burgers_vectors ********************
* Sets up a system of linear equations to solve for unknown Burgers vector
* components.  Each traction boundary condition component given in the input
* file leads to one equation and one unknown.  Solves the system of equations
* using the functions d_ludcmp() and d_lubksb() adapted from the book
* "Numerical Recipes in C, 2nd Ed." (Press et al., 1992).
*
* Modifies:
*	b_vector[]          - static array holding (previously) unknown Burgers
*                         vector components
*   "current_elt"->b[i] - pointer to ith component of Burgers vector.
*****************************************************************************/
#if defined(__STDC__) || defined(ANSI) /* ANSI */
void     determine_burgers_vectors()
#else
void determine_burgers_vectors()
#endif

{
	elt_t	*elt1;
	int		num_eqns = 0;
	int		eqn_num = 0;
	int		i;


	/* Determine number of simultaneous equations
	---------------------------------------------*/
	elt1 = first_elt_E;
	while (elt1 != NULL) {
		for (i=0; i < 3; i++) {
			if (elt1->bc_type[i] == TRACTION_BC)
				num_eqns++;
		}
		elt1 = elt1->next;
	}

	/* printf( "\n\n  -->nbr eq= %i", num_eqns ); */

	/* Set up vectors
	-----------------*/
	eqn_num = 1;
	elt1 = first_elt_E;
	while (elt1 != NULL) {
		for (i=0; i < 3; i++) {
			switch (elt1->bc_type[i]) {
				case BVECTOR_BC:
					elt1->b[i] = &elt1->bc[i];
					break;
			}
		}
		elt1 = elt1->next;
	}
	/* printf( "\n done with determining burgers_vectors ...\n"); */
}
