#ifndef __POSITION_H__
#define __POSITION_H__

#include "parameters.h"



/*********************************************************************************************
**********************************************************************************************

   POSITION_T

**********************************************************************************************
*********************************************************************************************/

/* The position fields are stored as an array 
 * so that they are contiguous in memory: 
 * this MIGHT be necessary for the direct computation with BLAS, 
 * but this MIGHT be slower for the direct computation without BLAS!
 * We therefore let the two possibilities ...  
 *
 * After testing: the position fields stored as arrays are not slower 
 * for the direct computation without BLAS! */

#ifdef _POSITION_T_STORED_WITH_ARRAY_ 

#define POS_X 0
#define POS_Y 1
#define POS_Z 2
typedef struct {
  COORDINATES_T dat[3];
} position_t;

#else /* #ifdef _POSITION_T_STORED_WITH_ARRAY_ */

typedef struct {
  COORDINATES_T x;
  COORDINATES_T y;
  COORDINATES_T z;
} position_t;

#endif /* #ifdef _POSITION_T_STORED_WITH_ARRAY_ */








/*! For the spherical coordinates we use the "physical" notation \n
  for theta and phi. That is to say that phi is the longitudinal
  coordinate whereas theta is the colatitudinal coordinate
  (the colatitudinal coordinate "starts" from the z-axis 
  whereas the latitudinal coordinate "starts" from the x-y plan 
  (i.e. the equator for the earth). */
typedef struct { /* (r, theta, phi) */
  COORDINATES_T r;
  COORDINATES_T th; /* for "theta": theta is in the range [0, Pi] */
  COORDINATES_T ph; /* for "phi": phi is in the range ]-Pi, Pi] */
} spherical_position_t;







#ifdef _POSITION_T_STORED_WITH_ARRAY_
FMB_INLINE void position_Set_x(position_t *FMB_RESTRICT p, COORDINATES_T v){  
  p->dat[POS_X] = v;  
} 
FMB_INLINE COORDINATES_T position_Get_x(FMB_CONST position_t *FMB_RESTRICT p){ 
  return p->dat[POS_X]; 
} 
FMB_INLINE COORDINATES_T *position_Get_p_x(position_t *FMB_RESTRICT p){ 
  return p->dat + POS_X; 
} 

FMB_INLINE void position_Set_y(position_t *FMB_RESTRICT p, COORDINATES_T v){  
  p->dat[POS_Y] = v;  
}  
FMB_INLINE COORDINATES_T position_Get_y(FMB_CONST position_t *FMB_RESTRICT p){ 
  return p->dat[POS_Y]; 
} 
FMB_INLINE COORDINATES_T *position_Get_p_y(position_t *FMB_RESTRICT p){ 
  return p->dat + POS_Y; 
} 
#else /* #ifdef _POSITION_T_STORED_WITH_ARRAY_ */
FMB_INLINE void position_Set_x(position_t *FMB_RESTRICT p, COORDINATES_T v){  
  p->x = v;  
} 
FMB_INLINE COORDINATES_T position_Get_x(FMB_CONST position_t *FMB_RESTRICT p){ 
  return p->x; 
} 
FMB_INLINE COORDINATES_T *position_Get_p_x(position_t *FMB_RESTRICT p){ 
  return &(p->x); 
} 

FMB_INLINE void position_Set_y(position_t *FMB_RESTRICT p, COORDINATES_T v){  
  p->y = v;  
}  
FMB_INLINE COORDINATES_T position_Get_y(FMB_CONST position_t *FMB_RESTRICT p){ 
  return p->y; 
} 
FMB_INLINE COORDINATES_T *position_Get_p_y(position_t *FMB_RESTRICT p){ 
  return &(p->y); 
} 
#endif /* #ifdef _POSITION_T_STORED_WITH_ARRAY_ */


/*********************************************************************************************
   3D
*********************************************************************************************/

#ifdef _POSITION_T_STORED_WITH_ARRAY_
/* For force_Direct_computation_[with_potential_]mutual() in 3D: */
#define POSITION_GET_X(p) ((p)->dat[POS_X]) 
#define POSITION_GET_Y(p) ((p)->dat[POS_Y]) 
#define POSITION_GET_Z(p) ((p)->dat[POS_Z]) 
FMB_INLINE void position_Set_z(position_t *FMB_RESTRICT p, COORDINATES_T v){  
  p->dat[POS_Z] = v;  
} 
FMB_INLINE COORDINATES_T position_Get_z(FMB_CONST position_t *FMB_RESTRICT p){ 
  return p->dat[POS_Z]; 
} 
FMB_INLINE COORDINATES_T *position_Get_p_z(position_t *FMB_RESTRICT p){ 
  return p->dat + POS_Z; 
} 
#else /* #ifdef _POSITION_T_STORED_WITH_ARRAY_ */
/* For force_Direct_computation_[with_potential_]mutual() in 3D: */
#define POSITION_GET_X(p) ((p)->x) 
#define POSITION_GET_Y(p) ((p)->y) 
#define POSITION_GET_Z(p) ((p)->z) 
FMB_INLINE void position_Set_z(position_t *FMB_RESTRICT p, COORDINATES_T v){  
  p->z = v;  
} 
FMB_INLINE COORDINATES_T position_Get_z(FMB_CONST position_t *FMB_RESTRICT p){ 
  return p->z; 
} 
FMB_INLINE COORDINATES_T *position_Get_p_z(position_t *FMB_RESTRICT p){ 
  return &(p->z); 
} 
#endif /* #ifdef _POSITION_T_STORED_WITH_ARRAY_ */


FMB_INLINE void position_Initialize(position_t *FMB_RESTRICT p){ 
  position_Set_x(p, (COORDINATES_T)0); 
  position_Set_y(p, (COORDINATES_T)0); 
  position_Set_z(p, (COORDINATES_T)0); 
} 

FMB_INLINE void position_Affect(position_t *FMB_RESTRICT p_to_affect, 
				FMB_CONST position_t *FMB_RESTRICT p_new_values){ 
  position_Set_x(p_to_affect, position_Get_x(p_new_values)); 
  position_Set_y(p_to_affect, position_Get_y(p_new_values));  
  position_Set_z(p_to_affect, position_Get_z(p_new_values)); 
} 

FMB_INLINE char position_Are_too_close(FMB_CONST position_t *FMB_RESTRICT p1, 
				       FMB_CONST position_t *FMB_RESTRICT p2){ 
  return ((COORDINATES_T__ARE_ALMOST_EQUAL(position_Get_x(p1), position_Get_x(p2), COORDINATES_T_MAX_RELATIVE_ERROR)  
 	   && COORDINATES_T__ARE_ALMOST_EQUAL(position_Get_y(p1), position_Get_y(p2), COORDINATES_T_MAX_RELATIVE_ERROR) 
 	   && COORDINATES_T__ARE_ALMOST_EQUAL(position_Get_z(p1), position_Get_z(p2), COORDINATES_T_MAX_RELATIVE_ERROR)) ?  
 	  TRUE :  
 	  FALSE); 
} 

FMB_INLINE void position_Negate(position_t *FMB_RESTRICT p){ 
  position_Set_x(p, - position_Get_x(p)); 
  position_Set_y(p, - position_Get_y(p)); 
  position_Set_z(p, - position_Get_z(p)); 
} 

/* Because of the FMB_RESTRICT, we can not have: p_target==p_src. */ 
FMB_INLINE void position_Add(position_t *FMB_RESTRICT p_target, 
			     FMB_CONST position_t *FMB_RESTRICT p_src){ 
#ifdef _DEBUG_
  if (p_target == p_src){
    FMB_error("ERROR: in position_Add(p_target, p_src), p_target == p_src.\n");
  }
#endif /* #ifdef _DEBUG_ */
  position_Set_x(p_target, position_Get_x(p_target) + position_Get_x(p_src)); 
  position_Set_y(p_target, position_Get_y(p_target) + position_Get_y(p_src));  
  position_Set_z(p_target, position_Get_z(p_target) + position_Get_z(p_src)); 
} 

/* Because of the FMB_RESTRICT, we can not have: p_target==p_src. */ 
FMB_INLINE void position_Substract(position_t *FMB_RESTRICT p_target, 
				   FMB_CONST position_t *FMB_RESTRICT p_src){ 
#ifdef _DEBUG_
  if (p_target == p_src){
    FMB_error("ERROR: in position_Substract(p_target, p_src), p_target == p_src.\n");
  }
#endif /* #ifdef _DEBUG_ */
  position_Set_x(p_target, position_Get_x(p_target) - position_Get_x(p_src)); 
  position_Set_y(p_target, position_Get_y(p_target) - position_Get_y(p_src));  
  position_Set_z(p_target, position_Get_z(p_target) - position_Get_z(p_src)); 
} 

FMB_INLINE void position_Mult_by_real(position_t *FMB_RESTRICT p_pos, 
				      COORDINATES_T real){ 
  position_Set_x(p_pos, position_Get_x(p_pos) * real);  
  position_Set_y(p_pos, position_Get_y(p_pos) * real);  
  position_Set_z(p_pos, position_Get_z(p_pos) * real); 
} 



FMB_INLINE COORDINATES_T position_Compute_square_distance(FMB_CONST position_t *FMB_RESTRICT p1, 
							  FMB_CONST position_t *FMB_RESTRICT p2){ 
  return ((position_Get_x(p1) - position_Get_x(p2)) * (position_Get_x(p1) - position_Get_x(p2)))  
    + ((position_Get_y(p1) - position_Get_y(p2)) * (position_Get_y(p1) - position_Get_y(p2)))  
    + ((position_Get_z(p1) - position_Get_z(p2)) * (position_Get_z(p1) - position_Get_z(p2))) ;  
} 


/* WARNING: the way this function is written may cause imbalance problems in terms
   of computation between leafs : there should be a better way of comparing positions 
   such that all leafs do roughly the same amount of computation in the direct computation
   phase ... */
/* return: *p1 "<" *p2 */
FMB_INLINE bool position_Is_lower(FMB_CONST position_t *FMB_RESTRICT p1, 
				  FMB_CONST position_t *FMB_RESTRICT p2){
  if (FMB_FABS(position_Get_x(p1) - position_Get_x(p2)) < EPSILON)
    if (FMB_FABS(position_Get_y(p1) - position_Get_y(p2)) < EPSILON)
      if (FMB_FABS(position_Get_z(p1) - position_Get_z(p2)) < EPSILON)
	return FALSE;
      else
	return position_Get_z(p1) < position_Get_z(p2);
    else 
      return position_Get_y(p1) < position_Get_y(p2);
  else 
    return position_Get_x(p1) < position_Get_x(p2);
}








/*! The "prec" parameter should be equal to "low" when we want to display positions 
  of bodies, and equal to "high" when we want to display computed force vectors of bodies. */
FMB_INLINE void position_Display(FMB_CONST position_t *FMB_RESTRICT p, 
				 FILE *f, 
				 precision_double_t prec){ 
  fprintf(f, "(");

  if (prec == low){
    fprintf(f, LOW_PRECISION_DOUBLE_FPRINTF, position_Get_x(p));
    fprintf(f, ", ");
    fprintf(f, LOW_PRECISION_DOUBLE_FPRINTF, position_Get_y(p));
    fprintf(f, ", ");
    fprintf(f, LOW_PRECISION_DOUBLE_FPRINTF, position_Get_z(p));
  }
  else {
    fprintf(f, HIGH_PRECISION_DOUBLE_FPRINTF, position_Get_x(p));
    fprintf(f, ", ");
    fprintf(f, HIGH_PRECISION_DOUBLE_FPRINTF, position_Get_y(p));
    fprintf(f, ", ");
    fprintf(f, HIGH_PRECISION_DOUBLE_FPRINTF, position_Get_z(p));
  }
  
  fprintf(f, ")");
} 





/*! The "prec" parameter should be equal to "low" when we want to display positions 
  of bodies, and equal to "high" when we want to display computed force vectors of bodies. */
FMB_INLINE void pos_xyz_Display(FMB_CONST COORDINATES_T pos_x, 
				FMB_CONST COORDINATES_T pos_y, 
				FMB_CONST COORDINATES_T pos_z, 
				FILE *f, 
				precision_double_t prec){ 
  fprintf(f, "(");

  if (prec == low){
    fprintf(f, LOW_PRECISION_DOUBLE_FPRINTF, pos_x);
    fprintf(f, ", ");
    fprintf(f, LOW_PRECISION_DOUBLE_FPRINTF, pos_y);
    fprintf(f, ", ");
    fprintf(f, LOW_PRECISION_DOUBLE_FPRINTF, pos_z);
  }
  else {
    fprintf(f, HIGH_PRECISION_DOUBLE_FPRINTF, pos_x);
    fprintf(f, ", ");
    fprintf(f, HIGH_PRECISION_DOUBLE_FPRINTF, pos_y);
    fprintf(f, ", ");
    fprintf(f, HIGH_PRECISION_DOUBLE_FPRINTF, pos_z);
  }
  
  fprintf(f, ")");
} 




FMB_INLINE bool position_Are_relative_corner_positions_correct(position_t *p_corner0_position,
							       position_t *p_corner1_position){
  if ((position_Get_x(p_corner0_position) > position_Get_x(p_corner1_position))
      || (position_Get_y(p_corner0_position) > position_Get_y(p_corner1_position))
      || (position_Get_z(p_corner0_position) > position_Get_z(p_corner1_position)))
    return FALSE;
  else
    return TRUE;
}










#endif







