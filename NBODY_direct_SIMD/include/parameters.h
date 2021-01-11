#ifndef __PARAMETERS_H__
#define __PARAMETERS_H__

#include "agreements.h"



#ifdef _SINGLE_PREC_

/* Warning: this may be not safe to use different types for COORDINATES_T, VALUES_T and REAL_T... */
#define COORDINATES_T float
#define COORDINATES_T_MIN FLT_MIN
#define VALUES_T float
#define REAL_T float

/* See IO.h: */
#define FMB_IO_FREAD_COORDINATE()  FMB_IO_Fread_float()
#define FMB_IO_FWRITE_COORDINATE(x) FMB_IO_Fwrite_float(x)
#define FMB_IO_FREAD_REAL()  FMB_IO_Fread_float()
#define FMB_IO_FWRITE_REAL(x) FMB_IO_Fwrite_float(x)
#define FMB_IO_FREAD_VALUE()  FMB_IO_Fread_float()
#define FMB_IO_FWRITE_VALUE(x) FMB_IO_Fwrite_float(x)

#else /* #ifdef _SINGLE_PREC_ */

/* Warning: this may be not safe to use different types for COORDINATES_T, VALUES_T and REAL_T... */
#define COORDINATES_T double
#define COORDINATES_T_MIN DBL_MIN
#define VALUES_T double
#define REAL_T double

#define FMB_IO_FREAD_COORDINATE()  FMB_IO_Fread_double()
#define FMB_IO_FWRITE_COORDINATE(x) FMB_IO_Fwrite_double(x)
#define FMB_IO_FREAD_REAL()  FMB_IO_Fread_double()
#define FMB_IO_FWRITE_REAL(x) FMB_IO_Fwrite_double(x)
#define FMB_IO_FREAD_VALUE()  FMB_IO_Fread_double()
#define FMB_IO_FWRITE_VALUE(x) FMB_IO_Fwrite_double(x)

#endif /* #ifdef _SINGLE_PREC_ */



/************ Indexes: ************/
typedef short height_ind_t;

/* Type of the indexes of the expansions: */
typedef int exp_ind_t;


/* Thread indexes: */
typedef int thread_ind_t; /* so that it can be stored in a 'void *' */

/* MPI Processes indexes: */
typedef int proc_ind_t; /* so that it can be stored in a 'void *' */

/* bodies_ind_t: */
typedef long bodies_ind_t; /* Can NOT be unsigned (see bodies_it_t for example). */
#define FORMAT_BODIES_IND_T "%li"
#define FMB_IO_FREAD_BODIES_IND()   FMB_IO_Fread_long()
#define FMB_IO_FWRITE_BODIES_IND(x) FMB_IO_Fwrite_long(x)






/************ FMB_Info: ***********/
typedef struct{

  REAL_T dt;  /* 0.001 */
  REAL_T tend ; /*0.001 */ 

  long nb_bodies; /* for this process only */
  long total_nb_bodies; /* sum of 'FMB_Info.nb_bodies' over all MPI processes */


  REAL_T eps_soft; /* Softening parameter ('epsilon') in gravitational computations. */
  REAL_T eps_soft_square; /* eps_soft_square = eps_soft . eps_soft */ 
  
  bool save; /* If set to TRUE, we save position, mass, force and/or potential of all particles.
	      * Default value: FALSE. */

  bool sum;  /* If set to TRUE, we compute and display the sum of the forces and/or potential over all particles.
	      * Default value: FALSE. */

} FMB_Info_t;

extern FMB_Info_t FMB_Info;




/* Some default values for 'FMB_Info': */
#define DFLT_NB_LESS_S_MIN 0
#define DFLT_H_MIN_BLAS_WITHOUT_RECOPIES 5








/***** Equality between floating point numbers: (for position_Are_too_close()) *****/
/*** 1st method: ***
 * #define GAP_TOLERATED_FOR_POSITIONS 0.0000000001 /\* 10^(-10) *\/ 
 * And then test like: (FMB_FABS(position_Get_x(p1) - position_Get_x(p2)) < GAP_TOLERATED_FOR_POSITIONS)   */
/*** 2nd method: ***
 * WARNING: this is not optimal (see "Comparing floating point numbers", Bruce Dawson
 * (http://www.cygnus-software.com/papers/comparingfloats/comparingfloats.htm)). */
#define COORDINATES_T_MAX_RELATIVE_ERROR 0.00001 /* default value for 'maxRelativeError' */
FMB_INLINE bool COORDINATES_T__ARE_ALMOST_EQUAL(COORDINATES_T A, COORDINATES_T B, COORDINATES_T maxRelativeError){
  COORDINATES_T relativeError;
  if (FMB_FABS(A - B) < COORDINATES_T_MIN){
    return TRUE;
  }
  if (FMB_FABS(B) > FMB_FABS(A)){
    relativeError = FMB_FABS((A - B) / B);
  }
  else {
    relativeError = FMB_FABS((A - B) / A);
  }
  if (relativeError <= maxRelativeError){
    return TRUE;
  }
  return FALSE;
}







/* Numbering of the exceptions for the nearest neighbors indexes: */
#define NOT_SET -1
#define BORDER -2 /* Only for "free-space boundaries conditions". */
#define NO_NEIGHBOR -3 /* Only for the non uniform case. */



/* BOUNDARY CONDITIONS: */
#define FREE_SPACE_BOUNDARY_CONDITIONS 1
#define PERIODIC_BOUNDARY_CONDITIONS 2
#define BOUNDARY_CONDITIONS FREE_SPACE_BOUNDARY_CONDITIONS

/* Direct computation between two bodies: */
#define NO_MUTUAL 0
#define MUTUAL 1


/***** To sort the "results" files: *****/
#ifdef __AIX__
/* The maximum number of fields of the 'sort' command on AIX 5.1 is 10, 
 * which is not enough: the 'sort command' fails for large files ... 
 * Use: ~/bin/gnu_sort.sh */
#define SORT_COMMANDE_FORMAT "sort -o %s -n -k 7,7 %s "
#endif /* __AIX__ */

#ifdef __LINUX__
/* The sorting fails with '-n' on pollux ... */
#define SORT_COMMANDE_FORMAT "sort -o %s -k 7,7 %s "
#endif /* __LINUX__ */ 






#endif 




