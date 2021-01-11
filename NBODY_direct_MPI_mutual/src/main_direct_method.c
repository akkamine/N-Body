/* #include <stdio.h> */
/* #include <stdlib.h> */
/* #include <string.h> */
#include <sys/stat.h>

#include "direct_method.h"
#include "IO.h" 
#include "direct_computation.h"
//#include "noeud.h"

/* For FMB_Info.save: */
#define RESULTS_DIR "../NBODY_direct_MPI_mutual_results_3520039/"
#define RESULTS_FILE "results_"

/*** For MPI: ***/
#include <mpi.h>

/*** For timers: ***/
#include <sys/time.h>

double my_gettimeofday(){
  struct timeval tmp_time;
  gettimeofday(&tmp_time, NULL);
  return tmp_time.tv_sec + (tmp_time.tv_usec * 1.0e-6L);
}



/* See definition at the end of this file. */
int parse_command(int argc, 
		  char **argv,
		  char **p_data_file,
		  char **p_results_file);

 

/*********************************************************************************************
**********************************************************************************************

   MAIN 

**********************************************************************************************
*********************************************************************************************/

int main(int argc, char **argv){

	int my_rank,nb_proc, i;
	MPI_Status status;
	int nb_bodies_loc;
	int source,dest;
	bodies_t bodies_loc;
	/* Données circulant sur l'anneau */
	COORDINATES_T* posx_tmp;
	COORDINATES_T* posy_tmp;
	COORDINATES_T* posz_tmp;
	VALUES_T* values_tmp;
	/* Permet de faire le send */
	COORDINATES_T* posx_send;
	COORDINATES_T* posy_send;
	COORDINATES_T* posz_send;
	VALUES_T* values_send;
	
	/* Forces circulant sur l'anneau */
	COORDINATES_T* fx_tmp;
	COORDINATES_T* fy_tmp;
	COORDINATES_T* fz_tmp;
	/* Permet de faire le send */
	COORDINATES_T* fx_send;
	COORDINATES_T* fy_send;
	COORDINATES_T* fz_send;	
	
	COORDINATES_T prod;
  //COORDINATES_T sumx_loc, sumy_loc, sumz_loc;
	//COORDINATES_T sumx_glob, sumy_glob, sumz_glob;
	COORDINATES_T sum_loc[3];
	COORDINATES_T sum_glob[3];
    
	long nb_steps = 0;
	REAL_T tstart, tend;
	int tnow ; 
	tstart = 0 ; 
	tnow = tstart ; 
	tend = 10.0 ; 

	char *data_file = NULL;
	char *results_file = NULL;

	/* Timers: */
	double t_start = 0.0, t_end = 0.0;
    
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&nb_proc);
	MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

	/********************************* Options on command line: ***************************************/
	f_output = stdout; /* by default */
	parse_command(argc, argv, &data_file, &results_file);

	if(my_rank==0){
		/******************************** Files and FILE* : ***************************************/
		if (INFO_DISPLAY(1)){
			fprintf(f_output, "*** Compute own interactions of the box defined in \"%s\" ***.\n", 
					data_file);
		}
    
    /********************************* Reading files by processus 0 *******************************/

		Direct_method_Init();
		/***************************** Bodies'positions and masses initialization: ****************/
		Direct_method_Data(data_file);
		/* init leapfrog */
		
		for(i=0;i<bodies.nb_bodies;i++){
			prod=FMB_Info.dt/(2.*bodies_Get_value(&bodies,i));
			bodies_Set_sx(&bodies, i, bodies_Get_sx(&bodies,i) + bodies_Get_fx(&bodies, i)*prod);
			bodies_Set_sy(&bodies, i, bodies_Get_sy(&bodies,i) + bodies_Get_fy(&bodies, i)*prod);
			bodies_Set_sz(&bodies, i, bodies_Get_sz(&bodies,i) + bodies_Get_fz(&bodies, i)*prod);
		}

		if (INFO_DISPLAY(1)){ 
		  fprintf(f_output, "Start Time : %lf \t End Time : %lf \t dt : %lf \n",tstart, tend, FMB_Info.dt);
		  fprintf(f_output, "Number of steps: %lu\n", (unsigned long) ((tend-tstart)/FMB_Info.dt));
		}
		nb_bodies_loc = bodies_Nb_bodies(&bodies)/nb_proc;
	}
	
	/****************************** Répartition des données *************************************/
	
  MPI_Bcast(&nb_bodies_loc,1,MPI_INT,0,MPI_COMM_WORLD);
	
	bodies_loc.nb_bodies = nb_bodies_loc;
	bodies_loc.size_allocated = nb_bodies_loc;
	bodies_loc.p_pos_x = malloc(sizeof(COORDINATES_T)*nb_bodies_loc);	
	bodies_loc.p_pos_y = malloc(sizeof(COORDINATES_T)*nb_bodies_loc);	
	bodies_loc.p_pos_z = malloc(sizeof(COORDINATES_T)*nb_bodies_loc);	
	bodies_loc.p_sx = malloc(sizeof(COORDINATES_T)*nb_bodies_loc);	
	bodies_loc.p_sy = malloc(sizeof(COORDINATES_T)*nb_bodies_loc);	
	bodies_loc.p_sz = malloc(sizeof(COORDINATES_T)*nb_bodies_loc);	
	bodies_loc.p_fx = malloc(sizeof(COORDINATES_T)*nb_bodies_loc);	
	bodies_loc.p_fy = malloc(sizeof(COORDINATES_T)*nb_bodies_loc);	
	bodies_loc.p_fz = malloc(sizeof(COORDINATES_T)*nb_bodies_loc);	
	bodies_loc.p_values = malloc(sizeof(VALUES_T)*nb_bodies_loc);	
	
	/*
	bodies_Initialize(&bodies_loc, nb_bodies_loc);
	
	for(i=0; i<nb_bodies_loc; i++){
		body_t body_tmp;
		body_Initialize(&body_tmp);
		bodies_Add(&bodies_loc, &body_tmp);
	}
	*/
	
	MPI_Scatter(bodies.p_pos_x,nb_bodies_loc,MPI_DOUBLE,bodies_loc.p_pos_x,nb_bodies_loc,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Scatter(bodies.p_pos_y,nb_bodies_loc,MPI_DOUBLE,bodies_loc.p_pos_y,nb_bodies_loc,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Scatter(bodies.p_pos_z,nb_bodies_loc,MPI_DOUBLE,bodies_loc.p_pos_z,nb_bodies_loc,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Scatter(bodies.p_sx,nb_bodies_loc,MPI_DOUBLE,bodies_loc.p_sx,nb_bodies_loc,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Scatter(bodies.p_sy,nb_bodies_loc,MPI_DOUBLE,bodies_loc.p_sy,nb_bodies_loc,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Scatter(bodies.p_sz,nb_bodies_loc,MPI_DOUBLE,bodies_loc.p_sz,nb_bodies_loc,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Scatter(bodies.p_fx,nb_bodies_loc,MPI_DOUBLE,bodies_loc.p_fx,nb_bodies_loc,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Scatter(bodies.p_fy,nb_bodies_loc,MPI_DOUBLE,bodies_loc.p_fy,nb_bodies_loc,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Scatter(bodies.p_fz,nb_bodies_loc,MPI_DOUBLE,bodies_loc.p_fz,nb_bodies_loc,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Scatter(bodies.p_values,nb_bodies_loc,MPI_DOUBLE,bodies_loc.p_values,nb_bodies_loc,MPI_DOUBLE,0,MPI_COMM_WORLD);

	/* Données circulant sur l'anneau */
	posx_tmp = malloc(sizeof(COORDINATES_T)*nb_bodies_loc);
	posy_tmp = malloc(sizeof(COORDINATES_T)*nb_bodies_loc);
	posz_tmp = malloc(sizeof(COORDINATES_T)*nb_bodies_loc);
	values_tmp = malloc(sizeof(VALUES_T)*nb_bodies_loc);
	memcpy(values_tmp,bodies_loc.p_values,sizeof(VALUES_T)*nb_bodies_loc);
	
	/* Permet de faire le send */
	posx_send = malloc(sizeof(COORDINATES_T)*nb_bodies_loc);
	posy_send = malloc(sizeof(COORDINATES_T)*nb_bodies_loc);
	posz_send = malloc(sizeof(COORDINATES_T)*nb_bodies_loc);
	values_send = malloc(sizeof(VALUES_T)*nb_bodies_loc);
	
	/* Permet de faire le MUTUAL */
	fx_tmp = calloc(sizeof(COORDINATES_T),nb_bodies_loc);
	fy_tmp = calloc(sizeof(COORDINATES_T),nb_bodies_loc);
	fz_tmp = calloc(sizeof(COORDINATES_T),nb_bodies_loc);
	fx_send = calloc(sizeof(COORDINATES_T),nb_bodies_loc);
	fy_send = calloc(sizeof(COORDINATES_T),nb_bodies_loc);
	fz_send = calloc(sizeof(COORDINATES_T),nb_bodies_loc);
	
  /******************************************************************************************/
  /********************************** Start of the simulation: ******************************/
  /******************************************************************************************/

		/* Start timer: */
		t_start = my_gettimeofday();
		
		while ( FMB_Info.dt*tnow < FMB_Info.tend ) { 
	
				/********************* Direct method computation: ************************************/
				memcpy(posx_tmp,bodies_loc.p_pos_x,sizeof(COORDINATES_T)*nb_bodies_loc);
				memcpy(posy_tmp,bodies_loc.p_pos_y,sizeof(COORDINATES_T)*nb_bodies_loc);
				memcpy(posz_tmp,bodies_loc.p_pos_z,sizeof(COORDINATES_T)*nb_bodies_loc);
				
				for(i=0;i<nb_bodies_loc;i++){
					fx_tmp[i] = 0;
					fy_tmp[i] = 0;
					fz_tmp[i] = 0;
				}
				
				/* calcul des forces entre les corps locaux */
				bodies_Compute_own_interaction(&bodies_loc);
				
				/* définition des source et destination du processus */
				dest = (my_rank+1)%nb_proc;
				source = (nb_proc+my_rank-1)%nb_proc;
				
				for(i=0;i<nb_proc/2;i++){ 
					if(my_rank%2==0){
						/* envoie des positions */
						MPI_Send(posx_tmp,nb_bodies_loc,MPI_DOUBLE,dest,0,MPI_COMM_WORLD);
						MPI_Send(posy_tmp,nb_bodies_loc,MPI_DOUBLE,dest,0,MPI_COMM_WORLD);
						MPI_Send(posz_tmp,nb_bodies_loc,MPI_DOUBLE,dest,0,MPI_COMM_WORLD);
						/* envoie des valeurs */
						MPI_Send(values_tmp,nb_bodies_loc,MPI_DOUBLE,dest,0,MPI_COMM_WORLD);
						/* envoie des forces */
						MPI_Send(fx_tmp,nb_bodies_loc,MPI_DOUBLE,dest,0,MPI_COMM_WORLD);
						MPI_Send(fy_tmp,nb_bodies_loc,MPI_DOUBLE,dest,0,MPI_COMM_WORLD);
						MPI_Send(fz_tmp,nb_bodies_loc,MPI_DOUBLE,dest,0,MPI_COMM_WORLD);
											
						/* réception des positions */			
						MPI_Recv(posx_tmp,nb_bodies_loc,MPI_DOUBLE,source,0,MPI_COMM_WORLD,&status);	
						MPI_Recv(posy_tmp,nb_bodies_loc,MPI_DOUBLE,source,0,MPI_COMM_WORLD,&status);	
						MPI_Recv(posz_tmp,nb_bodies_loc,MPI_DOUBLE,source,0,MPI_COMM_WORLD,&status);
						/* réception des valeurs */
						MPI_Recv(values_tmp,nb_bodies_loc,MPI_DOUBLE,source,0,MPI_COMM_WORLD,&status);
						/* réception des forces */
						MPI_Recv(fx_tmp,nb_bodies_loc,MPI_DOUBLE,source,0,MPI_COMM_WORLD,&status);	
						MPI_Recv(fy_tmp,nb_bodies_loc,MPI_DOUBLE,source,0,MPI_COMM_WORLD,&status);	
						MPI_Recv(fz_tmp,nb_bodies_loc,MPI_DOUBLE,source,0,MPI_COMM_WORLD,&status);
					}
					else{
						memcpy(posx_send,posx_tmp,sizeof(COORDINATES_T)*nb_bodies_loc);
						memcpy(posy_send,posy_tmp,sizeof(COORDINATES_T)*nb_bodies_loc);
						memcpy(posz_send,posz_tmp,sizeof(COORDINATES_T)*nb_bodies_loc);
						memcpy(values_send,values_tmp,sizeof(VALUES_T)*nb_bodies_loc);
						memcpy(fx_send,fx_tmp,sizeof(COORDINATES_T)*nb_bodies_loc);
						memcpy(fy_send,fy_tmp,sizeof(COORDINATES_T)*nb_bodies_loc);
						memcpy(fz_send,fz_tmp,sizeof(COORDINATES_T)*nb_bodies_loc);
						
						/* réception des positions */			
						MPI_Recv(posx_tmp,nb_bodies_loc,MPI_DOUBLE,source,0,MPI_COMM_WORLD,&status);	
						MPI_Recv(posy_tmp,nb_bodies_loc,MPI_DOUBLE,source,0,MPI_COMM_WORLD,&status);	
						MPI_Recv(posz_tmp,nb_bodies_loc,MPI_DOUBLE,source,0,MPI_COMM_WORLD,&status);
						/* réception des valeurs */
						MPI_Recv(values_tmp,nb_bodies_loc,MPI_DOUBLE,source,0,MPI_COMM_WORLD,&status);
						/* réception des forces */
						MPI_Recv(fx_tmp,nb_bodies_loc,MPI_DOUBLE,source,0,MPI_COMM_WORLD,&status);	
						MPI_Recv(fy_tmp,nb_bodies_loc,MPI_DOUBLE,source,0,MPI_COMM_WORLD,&status);	
						MPI_Recv(fz_tmp,nb_bodies_loc,MPI_DOUBLE,source,0,MPI_COMM_WORLD,&status);

						/* envoie des positions */
						MPI_Send(posx_send,nb_bodies_loc,MPI_DOUBLE,dest,0,MPI_COMM_WORLD);
						MPI_Send(posy_send,nb_bodies_loc,MPI_DOUBLE,dest,0,MPI_COMM_WORLD);
						MPI_Send(posz_send,nb_bodies_loc,MPI_DOUBLE,dest,0,MPI_COMM_WORLD);
						/* envoie des valeurs */
						MPI_Send(values_send,nb_bodies_loc,MPI_DOUBLE,dest,0,MPI_COMM_WORLD);
						/* envoie des forces */
						MPI_Send(fx_send,nb_bodies_loc,MPI_DOUBLE,dest,0,MPI_COMM_WORLD);
						MPI_Send(fy_send,nb_bodies_loc,MPI_DOUBLE,dest,0,MPI_COMM_WORLD);
						MPI_Send(fz_send,nb_bodies_loc,MPI_DOUBLE,dest,0,MPI_COMM_WORLD);
					}
					/* calcul des forces */
					if( (nb_proc%2!=0) || (nb_proc%2==0 && (i<nb_proc/2-1 || my_rank<nb_proc/2)))
						bodies_Compute_MUTUAL_interaction_MPI(&bodies_loc,posx_tmp,posy_tmp,posz_tmp,values_tmp,fx_tmp,fy_tmp,fz_tmp);
				}
				/* envoie des forces aux processus propriétaires */
				/* nouvelles destinations et sources */
				dest = (my_rank+nb_proc/2+(nb_proc%2!=0))%nb_proc;
				source = (my_rank-(nb_proc/2+(nb_proc%2!=0))+nb_proc)%nb_proc;
				if(my_rank < nb_proc/2){
					/* envoie des forces */
					MPI_Send(fx_tmp,nb_bodies_loc,MPI_DOUBLE,dest,0,MPI_COMM_WORLD);
					MPI_Send(fy_tmp,nb_bodies_loc,MPI_DOUBLE,dest,0,MPI_COMM_WORLD);
					MPI_Send(fz_tmp,nb_bodies_loc,MPI_DOUBLE,dest,0,MPI_COMM_WORLD);
					/* reception des forces */
					MPI_Recv(fx_tmp,nb_bodies_loc,MPI_DOUBLE,source,0,MPI_COMM_WORLD,&status);	
					MPI_Recv(fy_tmp,nb_bodies_loc,MPI_DOUBLE,source,0,MPI_COMM_WORLD,&status);
					MPI_Recv(fz_tmp,nb_bodies_loc,MPI_DOUBLE,source,0,MPI_COMM_WORLD,&status);
				}
				else{
					memcpy(fx_send,fx_tmp,sizeof(COORDINATES_T)*nb_bodies_loc);
					memcpy(fy_send,fy_tmp,sizeof(COORDINATES_T)*nb_bodies_loc);
					memcpy(fz_send,fz_tmp,sizeof(COORDINATES_T)*nb_bodies_loc);
					/* reception des forces */
					MPI_Recv(fx_tmp,nb_bodies_loc,MPI_DOUBLE,source,0,MPI_COMM_WORLD,&status);	
					MPI_Recv(fy_tmp,nb_bodies_loc,MPI_DOUBLE,source,0,MPI_COMM_WORLD,&status);
					MPI_Recv(fz_tmp,nb_bodies_loc,MPI_DOUBLE,source,0,MPI_COMM_WORLD,&status);
					
					/* envoie des forces */
					MPI_Send(fx_send,nb_bodies_loc,MPI_DOUBLE,dest,0,MPI_COMM_WORLD);
					MPI_Send(fy_send,nb_bodies_loc,MPI_DOUBLE,dest,0,MPI_COMM_WORLD);
					MPI_Send(fz_send,nb_bodies_loc,MPI_DOUBLE,dest,0,MPI_COMM_WORLD);
				}			
				for(i=0;i<nb_bodies_loc;i++){
					bodies_Set_fx(&bodies_loc, i, bodies_Get_fx(&bodies_loc,i)+fx_tmp[i]);
					bodies_Set_fy(&bodies_loc, i, bodies_Get_fy(&bodies_loc,i)+fy_tmp[i]);
					bodies_Set_fz(&bodies_loc, i, bodies_Get_fz(&bodies_loc,i)+fz_tmp[i]);			
				}						

    		for(i=0;i<nb_bodies_loc;i++){
    			prod = (FMB_Info.dt)/bodies_Get_value(&bodies_loc,i);

					bodies_Set_pos_x(&bodies_loc,i,bodies_Get_pos_x(&bodies_loc,i)+(FMB_Info.dt*(bodies_Get_sx(&bodies_loc,i))));
					bodies_Set_pos_y(&bodies_loc,i,bodies_Get_pos_y(&bodies_loc,i)+(FMB_Info.dt*(bodies_Get_sy(&bodies_loc,i))));
					bodies_Set_pos_z(&bodies_loc,i,bodies_Get_pos_z(&bodies_loc,i)+(FMB_Info.dt*(bodies_Get_sz(&bodies_loc,i))));
					
    			bodies_Set_sx(&bodies_loc, i, bodies_Get_sx(&bodies_loc,i) + bodies_Get_fx(&bodies_loc, i)*prod);
  				bodies_Set_sy(&bodies_loc, i, bodies_Get_sy(&bodies_loc,i) + bodies_Get_fy(&bodies_loc, i)*prod);
  				bodies_Set_sz(&bodies_loc, i, bodies_Get_sz(&bodies_loc,i) + bodies_Get_fz(&bodies_loc, i)*prod);
				}
		
		if(FMB_Info.save && nb_steps%2==0){
			MPI_Gather(bodies_loc.p_pos_x,nb_bodies_loc,MPI_DOUBLE,bodies.p_pos_x,nb_bodies_loc,MPI_DOUBLE,0,MPI_COMM_WORLD);
			MPI_Gather(bodies_loc.p_pos_y,nb_bodies_loc,MPI_DOUBLE,bodies.p_pos_y,nb_bodies_loc,MPI_DOUBLE,0,MPI_COMM_WORLD);
			MPI_Gather(bodies_loc.p_pos_z,nb_bodies_loc,MPI_DOUBLE,bodies.p_pos_z,nb_bodies_loc,MPI_DOUBLE,0,MPI_COMM_WORLD);
			MPI_Gather(bodies_loc.p_values,nb_bodies_loc,MPI_DOUBLE,bodies.p_values,nb_bodies_loc,MPI_DOUBLE,0,MPI_COMM_WORLD);
				
			if(my_rank == 0){	
				/************************* Save the positions and the forces: ***************************/
					if (FMB_Info.save){
						if (results_file == NULL){
							/* The 'results' filename has not been set yet: */
						#define TMP_STRING_LENGTH 10
							char step_number_string[TMP_STRING_LENGTH];
							int  results_file_length = 0;
	
							/* Find the relative filename in 'data_file': */
							char *rel_data_file = strrchr(data_file, '/') + 1 ; /* find last '/' and go to the next character */
	
							results_file_length = strlen(RESULTS_DIR) + strlen(RESULTS_FILE) + strlen(rel_data_file) + TMP_STRING_LENGTH + 1 /* for '\0' */ ; 
							results_file = (char *) FMB_malloc_with_check(results_file_length * sizeof(char));
							strncpy(step_number_string, "", TMP_STRING_LENGTH);
							sprintf(step_number_string, "_%lu", nb_steps ); 
							strncpy(results_file, "", results_file_length);
							strcpy(results_file, RESULTS_DIR); 
							strcat(results_file, RESULTS_FILE); 
							strcat(results_file, rel_data_file); 
							strcat(results_file, step_number_string); 
	//						printf("%s\n",step_number_string);
						#undef TMP_STRING_LENGTH
						}
						
						/* Create directory RESULTS_DIR: */
						{	struct stat filestat;
							if (stat (RESULTS_DIR, &filestat) != 0) {
								/* The directory RESULTS_DIR does not exist, we create it: */
								mkdir(RESULTS_DIR, 0700); 
							}
						}
						Direct_method_Dump_bodies(results_file, nb_steps, &bodies);
						FMB_free(results_file);
						results_file= NULL ; 
						
					} /* if (FMB_Info.save) */
				} /* if(my_rank == 0) */			
		} /* 	if(FMB_Info.save && nb_steps%2==0) */	
			
		      /*********************** Sum of the forces *********************************************/
			if(FMB_Info.sum && nb_steps%2==0){
/*	          Direct_method_Sum_V1(&bodies_loc,&sumx_loc,&sumy_loc,&sumz_loc);*/
/*            MPI_Reduce(&sumx_loc,&sumx_glob,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);*/
/*            MPI_Reduce(&sumy_loc,&sumy_glob,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);*/
/*            MPI_Reduce(&sumz_loc,&sumz_glob,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);*/
				
      	Direct_method_Sum_V2(&bodies_loc,sum_loc);
        MPI_Reduce(sum_loc,sum_glob,3,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
        if(my_rank==0){
        	FILE *force_sum;
          if((force_sum=fopen("somme_forces.txt","a"))==NULL)
        	  printf("Erreur à l'ouverture du fichier");

		      else{
						/* fprintf(force_sum,"Somme des forces à l'étape %lu :\nSur x : %1.15e \nSur y : %1.15e \nSur z : %1.15e\n",nb_steps,sumx_glob,sumy_glob,sumz_glob);*/
		      	fprintf(force_sum,"Somme des forces à l'étape %lu :\nSur x : %1.15e \nSur y : %1.15e \nSur z : %1.15e\n",nb_steps,sum_glob[0],sum_glob[1],sum_glob[2]);
		        fclose(force_sum);
		      }
		    }
			} /* if(FMB_Info.sum && nb_steps%2==0) */
		
		tnow++;
		nb_steps ++ ; 
  	}  /* while ( FMB_Info.dt*tnow <= FMB_Info.tend )  */

	/* End timer: */
	t_end = my_gettimeofday();
				
	printf("Temps d'exécution pour le processus %d : %lf secondes\n",my_rank,t_end-t_start);

  /******************************************************************************************/
  /********************************** End of the simulation: ********************************/
  /******************************************************************************************/

	if(my_rank==0){
		Direct_method_Terminate();

		/********************** Close FILE* and free memory before exiting: ***********************/
		if (argc == 3)
		  if (fclose(f_output) == EOF)
		    perror("fclose(f_output)");
		
		FMB_free(data_file);
		
	}
  /****************************************** EXIT ******************************************/
  //exit(EXIT_SUCCESS);
	MPI_Finalize();
}















/*********************************************************************************************
**********************************************************************************************

   usage

**********************************************************************************************
*********************************************************************************************/

void usage(){
  char mes[300] = "";
  
  sprintf(mes, "Usage : a.out [-h] %s [-o output_filename] --in[r]=data_filename %s \n"
	  , "[--soft value]"
	  , ""
	  );
  
  fprintf(stderr, "%s", mes);


  fprintf(stderr, "\nDescription of the short options:\n"); 
/*   fprintf(stderr, "\t -v \t\t\t Display the version.\n"); */
  fprintf(stderr, "\t -h \t\t\t Display this message.\n"); 
  fprintf(stderr, "\t -i 'level' \t\t Info display level (0, 1 or 2).\n");
  fprintf(stderr, "\t -o 'output_filename' \t Otherwise stdout.\n");


  fprintf(stderr, "\nDescription of the long options:\n");
  fprintf(stderr, "\t --in='filename' \t Input data filename.\n");

  fprintf(stderr, "\t --save \t\t Save position, mass, force and/or potential of all particles.\n");

  /* Unused in this code: */
  /*   fprintf(stderr, "\t --out='filename' \t Output data filename for '--save' option.\n"); */

  fprintf(stderr, "\t --sum  \t\t Compute and display the sum of the forces and/or potential over all particles.\n");
  fprintf(stderr, "\t --soft='value' \t Softening parameter.\n");
  fprintf(stderr, "\t --dt='value' \t\t Leapfrog integration timestep \n");
  fprintf(stderr, "\t --tend='value' \t Time to stop integration \n");

  /* We use only NEMO file format in this code: */
  /*   fprintf(stderr, "\t --it='value' \t\t input  data format ('fma' for FMB ASCII, 'fmb' for FMB binary, 'nemo').\n"); */
  /*   fprintf(stderr, "\t --ot='value' \t\t output data format ('fma' for FMB ASCII, 'fmah' for FMB ASCII human readable, 'fmb' for FMB binary, 'nemo').\n"); */

  fprintf(stderr, "\n");
  exit(EXIT_FAILURE);
}










/*********************************************************************************************
**********************************************************************************************

   parse_command

**********************************************************************************************
*********************************************************************************************/


/* Long option codes for 'val' field of struct option. 
 * Ascii codes 65 -> 90 ('A'->'Z') and 97 -> 122 ('a'->'z') 
 * are reserved for short options */
/* Same code as in main.c: */
#define LONGOPT_CODE_SOFT    14
#define LONGOPT_CODE_SAVE    24
#define LONGOPT_CODE_SUM     25
#define LONGOPT_CODE_IT 34
#define LONGOPT_CODE_OT 35
#define LONGOPT_CODE_IN 41
#define LONGOPT_CODE_OUT 43
#define LONGOPT_CODE_DT 49
#define LONGOPT_CODE_TEND 48 

int parse_command(int argc, 
		  char **argv,
		  char **p_data_file,
		  char **p_results_file){
  char options[]="hi:o:";
  int curr_opt;
  /*   opterr = 0; */

  /* Default values: */
  FMB_Info.dt = 0.1;
  FMB_Info.tend = 10.0 ; 
  FMB_Info.eps_soft = 0.01;

  
  struct option longopts[] = {
			      {"soft",
			       required_argument,
			       NULL, 
			       LONGOPT_CODE_SOFT},
			      {"dt",
			       required_argument,
			       NULL, 
			       LONGOPT_CODE_DT},
			      {"tend",
			       required_argument,
			       NULL, 
			       LONGOPT_CODE_TEND},
			      {"save",
			       no_argument,
             NULL,
			       LONGOPT_CODE_SAVE},
			      {"sum",
			       no_argument,
			       NULL, 
			       LONGOPT_CODE_SUM},
			      {"it",
			       required_argument,
			       NULL, 
			       LONGOPT_CODE_IT},
			      {"ot",
			       required_argument,
			       NULL, 
			       LONGOPT_CODE_OT},
			      {"in",
			       required_argument,
			       NULL, 
			       LONGOPT_CODE_IN},
			      {"out",
			       required_argument,
			       NULL, 
			       LONGOPT_CODE_OUT},
 			      {0}}; /* last element of the array  */
  

  
  if (argc == 1){
    usage();
  }


  /* Default values: see direct_*/
  FMB_Info.eps_soft = 0.0;

  
  curr_opt=getopt_long(argc, argv, options, longopts, NULL);
  while(curr_opt != (int) EOF){

    switch(curr_opt){
    case 'h' : 
      usage();
      break;
    case 'i':
      FMB_IO_Info.info_display_level = atoi(optarg);
      if (FMB_IO_Info.info_display_level != 0 && 
	  FMB_IO_Info.info_display_level != 1 &&
	  FMB_IO_Info.info_display_level != 2 &&
	  FMB_IO_Info.info_display_level != 3){
	FMB_error("Wrong FMB_IO_Info.info_display_level value.\n");
      }
      break;
    case 'o':
      if ((f_output = fopen(optarg, "w")) == NULL){
	perror("fopen(\'output_filename\', \"w\")");
      }	  
      break;
    case '?' : 
      usage();
      break;

      
    case LONGOPT_CODE_SOFT:
      FMB_Info.eps_soft = (REAL_T) atof(optarg);
      break;
    case LONGOPT_CODE_DT:
      FMB_Info.dt =(REAL_T) atof(optarg);
      break;
    case LONGOPT_CODE_TEND:
      FMB_Info.tend =(REAL_T) atof(optarg) ; 
      break;
    case LONGOPT_CODE_SAVE:
      FMB_Info.save = TRUE;
      break;
    case LONGOPT_CODE_SUM:
      FMB_Info.sum  = TRUE;
      break;      
    case LONGOPT_CODE_IT:
      if (strcmp(optarg, "fma") == 0){
	FMB_IO_Info.input_format = FMB_ASCII_format;
      }
      else {
	if (strcmp(optarg, "fmah") == 0){
	  FMB_error("FMB_ASCII_human_format is only for \"output format\", not for \"input format\".\n");
	} 
	else {
	  if (strcmp(optarg, "fmb") == 0){
	    FMB_IO_Info.input_format = FMB_binary_format;
	  }
	  else {
	    if ((strcmp(optarg, "nemo") == 0) || (strcmp(optarg, "NEMO"))){
	      FMB_IO_Info.input_format = NEMO_format;
	    }
	    else {
	      FMB_error("Unknow format for --it option!\n");	    
	    }
	  }
	}
      }
      FMB_IO_Info.input_format_from_cmd_line = TRUE;
      break;
    case LONGOPT_CODE_OT:
      if (strcmp(optarg, "fma") == 0){
	FMB_IO_Info.output_format = FMB_ASCII_format;
      }
      else {
	if (strcmp(optarg, "fmah") == 0){
	  FMB_IO_Info.output_format = FMB_ASCII_human_format;
	}
	else {
	  if (strcmp(optarg, "fmb") == 0){
	    FMB_IO_Info.output_format = FMB_binary_format;
	  }
	  else {
	    if ((strcmp(optarg, "nemo") == 0) || (strcmp(optarg, "NEMO"))){
	      FMB_IO_Info.output_format = NEMO_format;
	    }
	    else {
	      FMB_error("Unknow format for --it option!\n");	    
	    }
	  }
	}
      }
      FMB_IO_Info.output_format_from_cmd_line = TRUE;
      break;
    case LONGOPT_CODE_IN:
      if (*p_data_file != NULL){ FMB_ERROR_BRIEF(); }
      *p_data_file = (char *) FMB_malloc_with_check((strlen(optarg) + 1 /* for '\0' */) * sizeof(char));
      strcpy(*p_data_file, optarg); 
      break;
    case LONGOPT_CODE_OUT:
      if (*p_results_file != NULL){ FMB_ERROR_BRIEF(); }
      *p_results_file = (char *) FMB_malloc_with_check((strlen(optarg) + 1 /* for '\0' */) * sizeof(char));
      strcpy(*p_results_file, optarg); 
      break;
    } /* switch */
    curr_opt=getopt_long(argc, argv, options, longopts, NULL);
  }

  /* Check that an input data filename has been provided: */
  if (*p_data_file == NULL){
    FMB_error("No 'input data filename' provided.\n");
  }

  return 0;
}




















