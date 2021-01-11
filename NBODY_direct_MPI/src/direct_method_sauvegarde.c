#include "direct_method.h"
#include "IO.h" 
#include "bodies.h"
#include <mpi.h>
#include <omp.h>




/* Here are the initialization of the global variables: */
bodies_t bodies;
char *Direct_data_file;
bool Direct_are_data_bzipped2 = FALSE; 
position_t center;
COORDINATES_T half_side;

FMB_Info_t FMB_Info;


/* See definition in 'FMB.c'. */
void bunzip2_file(const char *filename);
void bzip2_file(const char *filename);












/*********************************************************************************************
**********************************************************************************************

   Direct_method_Init

**********************************************************************************************
*********************************************************************************************/
void Direct_method_Init(){

  /* Checking: */
  if (f_output == NULL){
    FMB_error("'f_output' must be set.\n");
  }


  /************************************ eps_soft_square: **********************************************/
  fprintf(f_output, "Softening parameter: %.1e\n", FMB_Info.eps_soft); 
  FMB_Info.eps_soft_square = FMB_Info.eps_soft * FMB_Info.eps_soft;

  /* Clear 'center' and 'half_side': */
  position_Initialize(&center);
  half_side = (COORDINATES_T) 0.0;

}






/*********************************************************************************************
**********************************************************************************************

   Direct_method_Data

**********************************************************************************************
*********************************************************************************************/
void Direct_method_Data(char *data_file){
  bodies_ind_t k;
  bodies_ind_t nb_bodies; 

  if (INFO_DISPLAY(2)){
    fprintf(f_output, "Opening data file \'%s\' for direct computation... \n", data_file); 
  }

  /* Initialize Input operations: */    
  FMB_IO_InitI(data_file);
  
  FMB_IO_Scan_header(&nb_bodies, &center, &half_side);

  if (INFO_DISPLAY(1)){
    fprintf(f_output, "Bodies number: ");
    fprintf(f_output, FORMAT_BODIES_IND_T, nb_bodies);
    fprintf(f_output, "\n"); 
    fflush(f_output);
  }


  bodies_Initialize(&bodies, nb_bodies);

  for (k=0; k<nb_bodies ;++k){
    /* We have to use 'bodies_Add()'! */
    body_t body_tmp;
    body_Initialize(&body_tmp);

    if (FMB_IO_Scan_body(&body_tmp) != 1){
      FMB_error("In Direct_method_Data(): FMB_IO_Scan_body() failed for body #%i\n", k);
    }

/*     if (k<100){ body_Display(&body_tmp, f_output); }  */
    
    bodies_Add(&bodies, &body_tmp);
  }

  bodies_ClearFP(&bodies);

  /* Terminate Input operations: */
  FMB_IO_TerminateI();

}




/*********************************************************************************************
 ********************************************************************************************
**********************************************************************************************

 Direct_method_Data_bodies

**********************************************************************************************
*********************************************************************************************/
 /* Same as Direct_method_Data() but we use the position and values
  * of all bodies stored in 'p_b' (instead of the bodies stored
  * in the file "data_file" in Direct_method_Data()). */
void Direct_method_Data_bodies(bodies_t *p_b){
  
  bodies_it_t it;

  bodies_Initialize(&bodies, bodies_Nb_bodies(p_b));

  for (bodies_it_Initialize(&it, p_b);
       bodies_it_Is_valid(&it);
       bodies_it_Go2Next(&it)){
    /* We have to use 'bodies_Add()'! */
    body_t body_tmp;
    bodies_it_Get_body(&it, &body_tmp);
    bodies_Add(&bodies, &body_tmp);
  }

  bodies_ClearFP(&bodies);

}





/*********************************************************************************************
**********************************************************************************************

   Direct_method_Compute

**********************************************************************************************
*********************************************************************************************/
void Direct_method_Compute(){

    /********************* With reciprocity: **********************************************/
    /* Compute the force and the potential: */
    bodies_Compute_own_interaction(&bodies);        


    /**************** Possible scaling with CONSTANT_INTERACTION_FACTOR: ********************/
    /* We can also use CONSTANT_INTERACTION_FACTOR only for the total potential energy ... */
#ifdef _USE_CONSTANT_INTERACTION_FACTOR_
    bodies_Scale_with_CONSTANT_INTERACTION_FACTOR(&bodies);
#endif /* #ifdef _USE_CONSTANT_INTERACTION_FACTOR_ */


}














/*********************************************************************************************
**********************************************************************************************

   Direct_method_Terminate

**********************************************************************************************
*********************************************************************************************/
void Direct_method_Terminate(){

  bodies_Free(&bodies);

  if (Direct_are_data_bzipped2){
    /* We recompress the data file: */
    bzip2_file(Direct_data_file);
  }
  FMB_free(Direct_data_file);

}



























/*********************************************************************************************
**********************************************************************************************

   sum

**********************************************************************************************
*********************************************************************************************/
void Direct_method_Sum(char *results_file,
		       unsigned long step_number_value,
		       bodies_t *p_bodies, 
		       VALUES_T total_potential_energy,double* bufx, double* bufy, double* bufz){

/* Sans MPI */
/*  int i;
  COORDINATES_T sumx=0,sumy=0,sumz=0;

  for(i=0;i<bodies_Nb_bodies(p_bodies);i++){
    sumx+=bodies_Get_fx(p_bodies, i);
    sumy+=bodies_Get_fy(p_bodies, i); 
    sumz+=bodies_Get_fz(p_bodies, i); 
  }
  FILE *f;
  if((f=fopen("somme_forces.txt","w"))==NULL){
    printf("Erreur à l'ouverture du fichier");
    return;
  }
  fprintf(f,"Somme des forces :\nSur x : %1.15e \nSur y : %1.15e \nSur z : %1.15e\n",sumx,sumy,sumz);
  fclose(f);
*/

/* AVEC MPI */
	int my_rank,nb_proc;
	MPI_Status status;
	MPI_Request request;
  int i,k;
  COORDINATES_T tmpx, tmpy, tmpz;
  COORDINATES_T sumx=0,sumy=0,sumz=0;
  /*COORDINATES_T bufx[bodies_Nb_bodies(p_bodies)/nb_proc];
  COORDINATES_T bufy[bodies_Nb_bodies(p_bodies)/nb_proc];
  COORDINATES_T bufz[bodies_Nb_bodies(p_bodies)/nb_proc];*/
  int i_2 = 1;

	/*Initialisation MPI*/
	//MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&nb_proc);
	MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
	printf("en debut \n");
/*	
	MPI_Scatter(p_bodies->p_fx,bodies_Nb_bodies(p_bodies)/nb_proc,MPI_DOUBLE,bufx,bodies_Nb_bodies(p_bodies)/nb_proc,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Scatter(p_bodies->p_fy,bodies_Nb_bodies(p_bodies)/nb_proc,MPI_DOUBLE,bufy,bodies_Nb_bodies(p_bodies)/nb_proc,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Scatter(p_bodies->p_fz,bodies_Nb_bodies(p_bodies)/nb_proc,MPI_DOUBLE,bufz,bodies_Nb_bodies(p_bodies)/nb_proc,MPI_DOUBLE,0,MPI_COMM_WORLD);
*/
  for(i=0;i<bodies_Nb_bodies(p_bodies)/nb_proc;i++){
    sumx+=bufx[i];
    sumy+=bufy[i]; 
    sumz+=bufz[i];
  }
  for(i=0;i<(int)(log(nb_proc));i++){ //Il faut que le nombre de processeurs soit une puissance de 2
    if(my_rank%(i_2*2)==1){
      MPI_Send(&sumx,1,MPI_DOUBLE,my_rank-i_2,0,MPI_COMM_WORLD);
      MPI_Send(&sumy,1,MPI_DOUBLE,my_rank-i_2,0,MPI_COMM_WORLD);
      MPI_Send(&sumz,1,MPI_DOUBLE,my_rank-i_2,0,MPI_COMM_WORLD);
    }
    else{
    	if((my_rank+i_2)%(i_2*2)==1){
    		MPI_Recv(&tmpx,1,MPI_DOUBLE,my_rank+i_2,0,MPI_COMM_WORLD,&status);
				sumx += tmpx;
    		MPI_Recv(&tmpy,1,MPI_DOUBLE,my_rank+i_2,0,MPI_COMM_WORLD,&status);
				sumy += tmpy;
    		MPI_Recv(&tmpz,1,MPI_DOUBLE,my_rank+i_2,0,MPI_COMM_WORLD,&status);
				sumz += tmpz;
			}
		}
		i_2*2;
	}
	
	if(my_rank==0){
		printf("en fin\n");
		FILE *f;
		if((f=fopen("somme_forces.txt","w"))==NULL){
		  printf("Erreur à l'ouverture du fichier");
		  return;
		}
		fprintf(f,"Somme des forces :\nSur x : %1.15e \nSur y : %1.15e \nSur z : %1.15e\n",sumx,sumy,sumz);
		fclose(f);
  }
}








/*********************************************************************************************
**********************************************************************************************

   save 

**********************************************************************************************
*********************************************************************************************/
void Direct_method_Dump_bodies(char *results_filename,
			       unsigned long step_number_value,
			       bodies_t *p_bodies){
  bodies_it_t it;

  /* Initialize Ouput operations: */    
  FMB_IO_InitO(results_filename);
  
  if (FMB_IO_Info.output_format != NEMO_format){
    
    /********** FMB file format: **********/
    if (FMB_IO_Info.output_format == FMB_binary_format){
      FMB_error("Unable to write the 'header' for FMB_binary_format in Direct_method_Dump_bodies(). \n");
    }
    FMB_IO_Print_header(step_number_value, FALSE /* only_position_and_value */,
			bodies_Nb_bodies(p_bodies), &center, half_side);
    
    for (bodies_it_Initialize(&it, p_bodies);
	 bodies_it_Is_valid(&it);
	 bodies_it_Go2Next(&it)){ 
      
      FMB_IO_Print_body_from_bodies_it(&it, FALSE /* only_position_and_value */);
    } /* for */
    
  } /* if (FMB_IO_Info.output_format != NEMO_format) */
  else {
    /********** NEMO file format: **********/
    FMB_IO_Print_all_bodies_from_bodies_t(p_bodies);
  } /* else (FMB_IO_Info.output_format != NEMO_format) */
  
  /* Terminate Output operations: */    
  FMB_IO_TerminateO();

}







/*********************************************************************************************
**********************************************************************************************

   integration en temps 

**********************************************************************************************
*********************************************************************************************/

void Intermediaire_Integration(COORDINATES_T prod, bodies_t *p_b, int i){  
  prod=FMB_Info.dt/(2*bodies_Get_value(p_b,i));
	bodies_Set_sx(p_b, i, bodies_Get_sx(p_b,i) + bodies_Get_fx(p_b, i)*prod);
  bodies_Set_sy(p_b, i, bodies_Get_sy(p_b,i) + bodies_Get_fy(p_b, i)*prod);
  bodies_Set_sz(p_b, i, bodies_Get_sz(p_b,i) + bodies_Get_fz(p_b, i)*prod);
}

void Direct_Method_Time_Integration(bodies_t *p_b){
  int i;
  COORDINATES_T prod;
  
  #pragma omp parallel for
  for(i=0;i<(p_b->nb_bodies);i++){
    Intermediaire_Integration(prod,p_b,i);
    bodies_Set_pos_x(p_b,i,bodies_Get_pos_x(p_b,i)+FMB_Info.dt*(bodies_Get_sx(p_b,i)));
    bodies_Set_pos_y(p_b,i,bodies_Get_pos_y(p_b,i)+FMB_Info.dt*(bodies_Get_sy(p_b,i)));
		bodies_Set_pos_z(p_b,i,bodies_Get_pos_z(p_b,i)+FMB_Info.dt*(bodies_Get_sz(p_b,i)));
  }
  Direct_method_Compute();
  #pragma omp parallel for
  for(i=0;i<(p_b->nb_bodies);i++)
    Intermediaire_Integration(prod,p_b,i);
}































