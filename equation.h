#ifndef __CPU_RICCI_H
#define __CPU_RICCI_H

#include "first_macro_1d.h"     /* Where global parameters are defined */
#include "structs_1d.h"         /* Where structures are defined */


struct mhd_par_1d { 
  struct GRID_PAR *grid_ptr;
  FLOAT mm;      /* mass of field */
  FLOAT sigma;   /* dissipation parameter */
  FLOAT R10;     /* Boundary condition at x=0, m = R0*p */
  FLOAT R11;     /* Boundary condition at x=1, p = R1*m */ 
  FLOAT c;
  FLOAT s; 
  FLOAT a; 




  /* normals to fases, first digit = coord. second = value */

  FLOAT nx_10;
  FLOAT nx_11;

};



#ifdef IMEX

void FF(struct GRID_PAR *grid_ptr, struct field_array *fields_ptr, struct field_array *derivs_ptr, struct FUNCTION_PAR *function_par);
void FS(struct GRID_PAR *grid_ptr, struct field_array *fields_ptr, struct field_array *derivs_ptr, struct FUNCTION_PAR *function_par);
void FI(struct GRID_PAR *grid_ptr, struct field_array *fields_ptr, struct field_array *derivs_ptr, struct FUNCTION_PAR *function_par);

#else

void FF(struct GRID_PAR *grid_ptr, struct field_array *fields_ptr, struct field_array *derivs_ptr, struct FUNCTION_PAR *function_par);

#endif // IMEX


static inline FLOAT Fx(FLOAT *u, struct GRID_PAR *grid_1d_ptr, struct FUNCTION_PAR *function_par, int);
static inline FLOAT mysign_zero(double d);
static inline FLOAT MM3(double a, double b, double c, double weight);

static inline FLOAT Source(FLOAT *u, struct GRID_PAR *grid_1d_ptr, struct FUNCTION_PAR *function_par, int);




#endif /* __CPU_RICCI_H */
