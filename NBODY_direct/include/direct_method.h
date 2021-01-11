#ifndef __DIRECT_METHOD_H__
#define __DIRECT_METHOD_H__


#include "bodies.h"
#include "memory_allocation.h"



/***************************** Global variables: *******************************/
/* We recall here the declarations of some useful global variables. */

extern bodies_t bodies;
extern FMB_Info_t FMB_Info;







/**************************** Direct Method: ************************************/
void Direct_method_Init();

void Direct_method_Data_bodies(bodies_t *p_b);
void Direct_method_Data(char *data_file);

void Direct_method_Compute();

void Direct_method_Move();
void KnD_Direct_method_Move(REAL_T dt ) ; 
void K_Direct_method_Move(REAL_T dt ) ; 

void Direct_method_Terminate();

void Direct_Method_Time_Integration(bodies_t *p_b);


/**************************** Auxiliary functions: ******************************/

void Direct_method_Sum_V1(bodies_t *bodies_loc, COORDINATES_T* sumx, COORDINATES_T* sumy, COORDINATES_T* sumz);
void Direct_method_Sum_V2(bodies_t *bodies_loc, COORDINATES_T* sum);

void Direct_method_Dump_bodies(char *results_file,
			       unsigned long step_number_value,
			       bodies_t *p_bodies);


void Direct_method_Compute();
void Intermediaire_Integration(bodies_t *p_b, int i);


#endif /* #ifndef __DIRECT_METHOD_H__ */






