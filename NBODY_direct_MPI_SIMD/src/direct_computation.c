#include "direct_computation.h"
#include <omp.h>
#include <mpi.h>



/*********************************************************************************************
**********************************************************************************************

   bodies_Compute_own_interaction

**********************************************************************************************
*********************************************************************************************/

void bodies_Compute_own_interaction_MPI(bodies_t *FMB_RESTRICT p_b, COORDINATES_T* posx_tmp, COORDINATES_T* posy_tmp, COORDINATES_T* posz_tmp, VALUES_T* values_tmp, COORDINATES_T eps_carre){
	int i;	
  for(i=0;i<(p_b->nb_bodies);i++){
		int j;	

		COORDINATES_T fx = 0;
		COORDINATES_T fy = 0;
		COORDINATES_T fz = 0;
		#pragma omp simd reduction(+:fx,fy,fz)
		for(j=0;j<(p_b->nb_bodies);j++){			
	  	COORDINATES_T Pxjmi = posx_tmp[j] - p_b->p_pos_x[i];
	  	COORDINATES_T Pyjmi = posy_tmp[j] - p_b->p_pos_y[i];
	  	COORDINATES_T Pzjmi = posz_tmp[j] - p_b->p_pos_z[i];

			COORDINATES_T r_carre = Pxjmi*Pxjmi + Pyjmi*Pyjmi + Pzjmi*Pzjmi;			
			COORDINATES_T scal = r_carre + eps_carre;				
			scal *= sqrt(scal);
			
			COORDINATES_T mij = p_b->p_values[i] * values_tmp[j];
			COORDINATES_T temp = mij/scal;
			
			COORDINATES_T Fjx = temp*Pxjmi;
			COORDINATES_T Fjy = temp*Pyjmi;
			COORDINATES_T Fjz = temp*Pzjmi;
			
			fx += Fjx;
			fy += Fjy;
			fz += Fjz;
		}
		p_b->p_fx[i] = fx;
		p_b->p_fy[i] = fy;
		p_b->p_fz[i] = fz;
  }
}

void bodies_Compute_own_interaction(bodies_t *FMB_RESTRICT p_b, COORDINATES_T eps_carre){

  int i;
	for(i=0;i<(p_b->nb_bodies);i++){
		int j;	

		COORDINATES_T fx = 0;
		COORDINATES_T fy = 0;
		COORDINATES_T fz = 0;
		#pragma omp simd reduction(+:fx,fy,fz)
		for(j=0;j<(p_b->nb_bodies);j++){			
	  	COORDINATES_T Pxjmi = p_b->p_pos_x[j] - p_b->p_pos_x[i];
	  	COORDINATES_T Pyjmi = p_b->p_pos_y[j] - p_b->p_pos_y[i];
	  	COORDINATES_T Pzjmi = p_b->p_pos_z[j] - p_b->p_pos_z[i];

			COORDINATES_T r_carre = Pxjmi*Pxjmi + Pyjmi*Pyjmi + Pzjmi*Pzjmi;			
			COORDINATES_T scal = r_carre + eps_carre;				
			scal *= sqrt(scal);
			
			COORDINATES_T mij = p_b->p_values[i] * p_b->p_values[j];
			COORDINATES_T temp = mij/scal;
			
			COORDINATES_T Fjx = temp*Pxjmi;
			COORDINATES_T Fjy = temp*Pyjmi;
			COORDINATES_T Fjz = temp*Pzjmi;
			
			fx += Fjx;
			fy += Fjy;
			fz += Fjz;
		}
		p_b->p_fx[i] = fx;
		p_b->p_fy[i] = fy;
		p_b->p_fz[i] = fz;
	} 
}
