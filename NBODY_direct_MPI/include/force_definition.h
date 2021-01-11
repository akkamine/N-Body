#ifndef __FORCE_DEFINITION_H__
#define __FORCE_DEFINITION_H__

/*! \file 
  We describe here all the features of the force field between the bodies
  (electrostatic, gravitational ...) */

/*! Universal gravitational constant: */
#define G ((VALUES_T) 6.67259E-11)
/* In 2D: */
/* #define POTENTIAL_SIGN  */
/* #define UNIT_VECTOR_COMPONENT(tgt_comp, src_comp) ((tgt_comp) - (src_comp))   */
/* In 3D: */
#define POTENTIAL_SIGN -


#define CONSTANT_INTERACTION_FACTOR G 










#endif /* #ifndef FORCE_DEFINITION_H */


