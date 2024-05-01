#ifndef _vectdefs_h
#define _vectdefs_h

#define THREEDIM // dimension of the system
#define RBC // rigid boundary conditions or periodic boundary conditions
#define CR 1 // coefficient of restitution (if CR=1 then elastic collision, if CR=0 then perfectly inelastic collision) 

#if !defined(NDIM) && !defined(TWODIM) && !defined(THREEDIM)
#define THREEDIM                              /* supply default dimensions  */
#endif

#if defined(THREEDIM) || (NDIM==3)
#undef  TWODIM
#define THREEDIM
#define NDIM 3
#endif

#if defined(TWODIM) || (NDIM==2)
#undef  THREEDIM
#define TWODIM
#define NDIM 2
#endif

#if !defined(BC) && !defined(RBC) && !defined(PBC)
#define RBC                              /* supply default boundary conditions  */
#endif

#if defined(RBC)
#undef PBC
#define RBC
#define BC 1
#endif

#if defined(PBC)
#undef RBC
#define PBC
#define BC 2
#endif

//typedef real vector[NDIM];
////typedef real matrix[NDIM][NDIM];

#endif  /* ! _vectdefs_h */

