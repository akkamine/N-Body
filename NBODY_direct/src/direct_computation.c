#include "direct_computation.h"
#include <omp.h>
#include <mpi.h>



/*********************************************************************************************
**********************************************************************************************

   bodies_Compute_own_interaction

**********************************************************************************************
*********************************************************************************************/

void bodies_Compute_own_interaction(bodies_t *FMB_RESTRICT p_b){
  int i,j;
  COORDINATES_T Pxjmi,Pyjmi,Pzjmi,temp,scal, r_carre, mij,Fjx,Fjy,Fjz;

/* SANS MUTUAL SEQUENTIEL */

  for(i=0;i<(p_b->nb_bodies);i++){
  	bodies_Set_fx(p_b,i,0);
  	bodies_Set_fy(p_b,i,0);
  	bodies_Set_fz(p_b,i,0);
    for(j=0;j<p_b->nb_bodies;j++){
      if(i!=j){
      	Pxjmi=bodies_Get_pos_x(p_b, j)-bodies_Get_pos_x(p_b, i);
				Pyjmi=bodies_Get_pos_y(p_b, j)-bodies_Get_pos_y(p_b, i);
				Pzjmi=bodies_Get_pos_z(p_b, j)-bodies_Get_pos_z(p_b, i);

				r_carre = Pxjmi*Pxjmi + Pyjmi*Pyjmi + Pzjmi*Pzjmi;			
				scal = r_carre + pow(FMB_Info.eps_soft,2);				
				scal *= sqrt(scal);
				mij = bodies_Get_value(p_b, i)*bodies_Get_value(p_b, j);
				temp = mij/scal;
				Fjx = temp*Pxjmi;
				Fjy = temp*Pyjmi;
				Fjz = temp*Pzjmi;
				bodies_Set_fx(p_b, i, bodies_Get_fx(p_b,i)+Fjx);
				bodies_Set_fy(p_b, i, bodies_Get_fy(p_b,i)+Fjy);
				bodies_Set_fz(p_b, i, bodies_Get_fz(p_b,i)+Fjz);
      }
    }
  }

/* AVEC MUTUAL SEQUENTIEL */

/*	COORDINATES_T Fjx,Fjy,Fjz;
	for(i=0;i<(p_b->nb_bodies);i++){
		bodies_Set_fx(p_b,i,0);
		bodies_Set_fy(p_b,i,0);
		bodies_Set_fz(p_b,i,0);
	}

	for(i=0;i<(p_b->nb_bodies);i++){
		for(j=0;j<i;j++){
			Pxjmi=bodies_Get_pos_x(p_b, j)-bodies_Get_pos_x(p_b, i);
			Pyjmi=bodies_Get_pos_y(p_b, j)-bodies_Get_pos_y(p_b, i);
			Pzjmi=bodies_Get_pos_z(p_b, j)-bodies_Get_pos_z(p_b, i);
			      		
			r_carre = Pxjmi*Pxjmi + Pyjmi*Pyjmi + Pzjmi*Pzjmi;			
			scal = r_carre + FMB_Info.eps_soft*FMB_Info.eps_soft;				
			scal *= sqrt(scal);
			mij = bodies_Get_value(p_b, i)*bodies_Get_value(p_b, j);
			temp = mij/scal;

			Fjx=temp*Pxjmi;
			Fjy=temp*Pyjmi;
			Fjz=temp*Pzjmi;

			bodies_Set_fx(p_b, i, bodies_Get_fx(p_b,i)+Fjx);
			bodies_Set_fy(p_b, i, bodies_Get_fy(p_b,i)+Fjy);
			bodies_Set_fz(p_b, i, bodies_Get_fz(p_b,i)+Fjz);

			bodies_Set_fx(p_b, j, bodies_Get_fx(p_b,j)-Fjx);
			bodies_Set_fy(p_b, j, bodies_Get_fy(p_b,j)-Fjy);
			bodies_Set_fz(p_b, j, bodies_Get_fz(p_b,j)-Fjz);
		
		}		
	}*/


} 


/* Mutual opti */
/*
  COORDINATES_T Fjx,Fjy,Fjz;
  int biais;
  position_t* lforce;
#pragma omp parallel
  {
#pragma omp single
    lforce = malloc(sizeof(position_t)*p_b->nb_bodies*omp_get_num_threads());
#pragma omp for schedule(static)
    for(i=0;i<p_b->nb_bodies*omp_get_num_threads();i++){
      lforce[i].x=0;
      lforce[i].y=0;
      lforce[i].z=0;
    }
#pragma omp for private(j,Bi,Bj,Pximj,Pyimj,Pzimj,Fjx,Fjy,Fjz,scal,temp,biais) schedule(dynamic)
    for(i=0;i<p_b->nb_bodies;i++){
      Bi=p_b->p_bodies[i];
      for(j=0;j<i;j++){
				Bj=p_b->p_bodies[j];
				Pximj=Bi.position.x-Bj.position.x;
				Pyimj=Bi.position.y-Bj.position.y;
				Pzimj=Bi.position.z-Bj.position.z;
				scal=Pximj*Pximj+Pyimj*Pyimj+Pzimj*Pzimj+FMB_Info.eps_soft*FMB_Info.eps_soft;
				temp=(Bi.v*Bj.v)/(scal*sqrt(scal));
				Fjx=temp*Pximj;
				Fjy=temp*Pyimj;
				Fjz=temp*Pzimj;
				biais=omp_get_thread_num()*p_b->nb_bodies;
				lforce[i+biais].x-=Fjx;
				lforce[i+biais].y-=Fjy;
				lforce[i+biais].z-=Fjz;
				lforce[j+biais].x+=Fjx;
				lforce[j+biais].y+=Fjy;
				lforce[j+biais].z+=Fjz; 
      }
    }
#pragma omp for private(j)
    for(i=0;i<p_b->nb_bodies;i++){
      p_b->p_bodies[i].force_vector.x=0;
      p_b->p_bodies[i].force_vector.y=0;
      p_b->p_bodies[i].force_vector.z=0;
      for(j=0;j<omp_get_num_threads();j++){
	p_b->p_bodies[i].force_vector.x+=lforce[j*p_b->nb_bodies+i].x;
	p_b->p_bodies[i].force_vector.y+=lforce[j*p_b->nb_bodies+i].y;
	p_b->p_bodies[i].force_vector.z+=lforce[j*p_b->nb_bodies+i].z;
      }
    }
  }
*/
//}



//#endif /* #ifdef _BODIES_SPLIT_DATA_ */
