#include "direct_computation.h"
#include <omp.h>
#include <mpi.h>



/*********************************************************************************************
**********************************************************************************************

   bodies_Compute_own_interaction

**********************************************************************************************
*********************************************************************************************/


void bodies_Compute_own_interaction_MPI(bodies_t *FMB_RESTRICT p_b_fixe, COORDINATES_T* posx_tmp, COORDINATES_T* posy_tmp, COORDINATES_T* posz_tmp, VALUES_T* values_tmp){
	
	int i,j;
	COORDINATES_T Pxjmi,Pyjmi,Pzjmi,temp,scal;
	COORDINATES_T r_carre, mij,Fjx,Fjy,Fjz;
	
  for(i=0;i<(p_b_fixe->nb_bodies);i++){
    for(j=0;j<(p_b_fixe->nb_bodies);j++){
      Pxjmi=posx_tmp[j]-bodies_Get_pos_x(p_b_fixe, i);
			Pyjmi=posy_tmp[j]-bodies_Get_pos_y(p_b_fixe, i);
			Pzjmi=posz_tmp[j]-bodies_Get_pos_z(p_b_fixe, i);

			r_carre = Pxjmi*Pxjmi + Pyjmi*Pyjmi + Pzjmi*Pzjmi;			
			scal = r_carre + FMB_Info.eps_soft*FMB_Info.eps_soft;				
			scal *= sqrt(scal);
			
			mij = bodies_Get_value(p_b_fixe, i)*values_tmp[j];
			temp = mij/scal;
			Fjx = temp*Pxjmi;
			Fjy = temp*Pyjmi;
			Fjz = temp*Pzjmi;
			bodies_Set_fx(p_b_fixe, i, bodies_Get_fx(p_b_fixe,i)+Fjx);
			bodies_Set_fy(p_b_fixe, i, bodies_Get_fy(p_b_fixe,i)+Fjy);
			bodies_Set_fz(p_b_fixe, i, bodies_Get_fz(p_b_fixe,i)+Fjz);
    }
  }
}

/*void bodies_Compute_own_interaction_MPI_odd(COORDINATES_T* p_posx, COORDINATES_T* p_posy, COORDINATES_T* p_posz, COORDINATES_T* p_fx, COORDINATES_T* p_fy, COORDINATES_T* p_fz, VALUES_T* p_values, bodies_ind_t nb_bodies_loc, COORDINATES_T* posx_tmp, COORDINATES_T* posy_tmp, COORDINATES_T* posz_tmp, VALUES_T* values_tmp){*/
/*	*/
/*	int i,j;*/
/*	COORDINATES_T Pxjmi,Pyjmi,Pzjmi,temp,scal;*/
/*	COORDINATES_T r_carre, mij,Fjx,Fjy,Fjz;*/
/*	*/
/*  for(i=0;i<nb_bodies_loc;i++){*/
/*    for(j=0;j<nb_bodies_loc;j++){*/
/*      Pxjmi=posx_tmp[j]-p_posx[i];*/
/*			Pyjmi=posy_tmp[j]-p_posy[i];*/
/*			Pzjmi=posz_tmp[j]-p_posz[i];*/

/*			r_carre = Pxjmi*Pxjmi + Pyjmi*Pyjmi + Pzjmi*Pzjmi;			*/
/*			scal = r_carre + FMB_Info.eps_soft*FMB_Info.eps_soft;				*/
/*			scal *= sqrt(scal);*/
/*			*/
/*			mij = p_values[i]*values_tmp[j];*/
/*			temp = mij/scal;*/
/*				*/
/*			p_fx[i] += temp*Pxjmi;*/
/*			p_fy[i] += temp*Pyjmi;*/
/*			p_fz[i] += temp*Pzjmi;*/
/*    }*/
/*  }*/
/*}*/



void bodies_Compute_own_interaction(bodies_t *FMB_RESTRICT p_b){

  int i,j;
  COORDINATES_T Pxjmi,Pyjmi,Pzjmi,temp,scal, r_carre, mij;

#ifdef _BODIES_SPLIT_DATA_  


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
				bodies_Set_fx(p_b, i, bodies_Get_fx(p_b,i)+temp*Pxjmi);
				bodies_Set_fy(p_b, i, bodies_Get_fy(p_b,i)+temp*Pyjmi);
				bodies_Set_fz(p_b, i, bodies_Get_fz(p_b,i)+temp*Pzjmi);
      }
    }
  }


/* OK - Mutual Sequentiel */
/*
	COORDINATES_T Fjx,Fjy,Fjz;
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
	}  
*/

/* MPI Sans Mutual 
	int my_rank, nb_proc;
	MPI_Comm_size(MPI_COMM_WORLD,&nb_proc);
	MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
	
*/

/* AVEC MUTUAL AVEC MPI*/
/*
	int rang,p;
	MPI_Status status;
	MPI_Request request;
	
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&p);
	MPI_Comm_rank(MPI_COMM_WORLD,&rang);

	if(rang!=0)
		lforce = calloc(bodies_Nb_bodies,sizeof(position_t));
	
  for(i=0;i<p_b->nb_bodies;i++){
  	if(i%p==rang){
  
			for(j=0;j<i;j++){
				Pxjmi=bodies_Get_pos_x(p_b, j)-bodies_Get_pos_x(p_b, i);
				Pyjmi=bodies_Get_pos_y(p_b, j)-bodies_Get_pos_y(p_b, i);
				Pzjmi=bodies_Get_pos_z(p_b, j)-bodies_Get_pos_z(p_b, i);
				
				scal=(Pxjmi*Pxjmi+Pyjmi*Pyjmi+Pzjmi*Pzjmi)+FMB_Info.eps_soft*FMB_Info.eps_soft;
			
		   	temp=(bodies_Get_value(p_b, i)*bodies_Get_value(p_b, j))/(scal*sqrt(scal));
			
				Fx=temp*Pxjmi;
				Fy=temp*Pyjmi;
				Fz=temp*Pzjmi;
				
				if(rang!=0){
					lforce[i].x+=Fx;
					lforce[i].y+=Fy;
					lforce[i].z+=Fz;
					lforce[j].x-=Fx;
					lforce[j].y-=Fy;
					lforce[j].z-=Fz;
				}
				else{
					bodies_Set_fx(p_b, i, bodies_Get_fx(p_b,i)+temp*Pxjmi);
					bodies_Set_fy(p_b, i, bodies_Get_fy(p_b,i)+temp*Pyjmi);
					bodies_Set_fz(p_b, i, bodies_Get_fz(p_b,i)+temp*Pzjmi); 
					bodies_Set_fx(p_b, j, bodies_Get_fx(p_b,j)-temp*Pxjmi);
					bodies_Set_fy(p_b, j, bodies_Get_fy(p_b,j)-temp*Pyjmi);
					bodies_Set_fz(p_b, j, bodies_Get_fz(p_b,j)-temp*Pzjmi); 
		  }
		}
  }
  
  //Faire la somme des buffer par pair jusqu'à n'avoir plus qu'un buffer de taille nb_bodies

*/

}  /* Fermeture de la fonction
#else /* #ifdef _BODIES_SPLIT_DATA_ */

/* OK - Mutual avec critical */
/*
  COORDINATES_T Fjx,Fjy,Fjz;
#pragma omp parallel for
  for(i=0;i<(p_b->nb_bodies);i++){
    p_b->p_bodies[i].force_vector.x=0;
    p_b->p_bodies[i].force_vector.y=0;
    p_b->p_bodies[i].force_vector.z=0;
}
#pragma omp parallel for private(j,Bi,Bj,Pximj,Pyimj,Pzimj,Fjx,Fjy,Fjz,scal,temp) schedule(dynamic)
  for(i=0;i<(p_b->nb_bodies);i++){
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
#pragma omp critical
{
      p_b->p_bodies[i].force_vector.x-=Fjx;
      p_b->p_bodies[i].force_vector.y-=Fjy;
      p_b->p_bodies[i].force_vector.z-=Fjz;
      p_b->p_bodies[j].force_vector.x+=Fjx;
      p_b->p_bodies[j].force_vector.y+=Fjy;
      p_b->p_bodies[j].force_vector.z+=Fjz;
}
    }     
  }
*/

/* OK - Sans Mutual */
/*
#pragma omp parallel for private (j,Bi,Bj,Pximj,Pyimj,Pzimj,scal,temp) schedule (static)
  for(i=0;i<(p_b->nb_bodies);i++){
    Bi=p_b->p_bodies[i];
    Bi.force_vector.x=0;
    Bi.force_vector.y=0;
    Bi.force_vector.z=0;
    for(j=0;j<p_b->nb_bodies;j++){
      if(i!=j){
      Bj=p_b->p_bodies[j];
      Pximj=Bi.position.x-Bj.position.x;
      Pyimj=Bi.position.y-Bj.position.y;
      Pzimj=Bi.position.z-Bj.position.z;
      scal=Pximj*Pximj+Pyimj*Pyimj+Pzimj*Pzimj+FMB_Info.eps_soft*FMB_Info.eps_soft;
      temp=(Bi.v*Bj.v)/(scal*sqrt(scal));
      Bi.force_vector.x-=temp*Pximj;
      Bi.force_vector.y-=temp*Pyimj;
      Bi.force_vector.z-=temp*Pzimj;
      }
    }
    p_b->p_bodies[i]=Bi;
  }

*/
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



#endif /* #ifdef _BODIES_SPLIT_DATA_ */
