#include "noeud.h"

COORDINATES_T distance(position_t position1,position_t position2);
void intermediaire_tree_Compute(noeud* noeud,body_t* corps,COORDINATES_T accuracy);
void tree_Compute(bodies_t *FMB_RESTRICT p_b,COORDINATES_T accuracy);
