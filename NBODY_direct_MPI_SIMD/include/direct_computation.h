#ifndef __DIRECT_COMPUTATION_H__
#define __DIRECT_COMPUTATION_H__

#include "bodies.h"






/*********************************************************************************************
**********************************************************************************************

   Compute the interaction (forces and/or potentials)

**********************************************************************************************
*********************************************************************************************/



/*! The mutual interaction principle is always used in 'bodies_Compute_own_interaction()'. */
void bodies_Compute_own_interaction(bodies_t *FMB_RESTRICT p_b, COORDINATES_T eps_carre);


/*void bodies_Compute_own_interaction_MPI_V1(bodies_t *FMB_RESTRICT p_b_fixe, COORDINATES_T* posx_tmp, COORDINATES_T* posy_tmp, COORDINATES_T* posz_tmp, VALUES_T* values_tmp);*/

/*void bodies_Compute_own_interaction_MPI(COORDINATES_T* p_posx, COORDINATES_T* p_posy, COORDINATES_T* p_posz, COORDINATES_T* p_fx, COORDINATES_T* p_fy, COORDINATES_T* p_fz, VALUES_T* p_values, bodies_ind_t nb_bodies_loc, COORDINATES_T* posx_tmp, COORDINATES_T* posy_tmp, COORDINATES_T* posz_tmp, VALUES_T* values_tmp);*/

void bodies_Compute_own_interaction_MPI(bodies_t *FMB_RESTRICT p_b_fixe, COORDINATES_T* posx_tmp, COORDINATES_T* posy_tmp, COORDINATES_T* posz_tmp, VALUES_T* values_tmp, COORDINATES_T eps_carre);

#endif /* #ifdef __DIRECT_COMPUTATION_H__ */ 
