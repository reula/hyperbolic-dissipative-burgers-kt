/*********************************************************************
* Here are the global definitions used by all the subroutines        *
*                                                                    *       
*********************************************************************/
 
#ifdef FIRST_MACRO_1D_H
#else 
#define FIRST_MACRO_1D_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>
/* #include <curses.h> */
/* #include <ncurses.h> */
#include <time.h>       
#include <string.h>    /* functions to get data via web (stadin) */
/* #include <rfftw.h>  */    /*  needed for taking fft in derivs and main */

//#define DEBUG
//#define DEBUG_ADISCO
//#define DEBUG_PYGRAPH


#define DISSIPATION
#undef DISSIPATION  

#define IMAX          /* Mixed implicit explicit stiff solvers need more than one RHS and an explicit inversion */
#undef IMAX

#define FLUX
#ifndef FLUX
#define F_DIFF // Finite DIFFERENCES SBP with Penalties or Periodic
#endif

#define PERIODIC
//#undef PERIODIC

#ifdef FLUX
#ifdef PERIODIC
#define N_Ghost 0
#else
#define N_Ghost 0
//It was previously 2
#endif // PERIODIC
#endif // FLUX

 
 
#ifdef F_DIFF
#define PENALTY // if not it is Olsson, must undef PERIODIC....
//#undef PENALTY
#define N_Ghost 0
#endif // F_DIFF


#define CATANEO
//#undef CATANEO

#define SOURCE
//#undef SOURCE


/* That is the number of variables that enter the system            */


#ifndef CATANEO
#define N_FIELDS 1
#else
#define N_FIELDS 2
#endif

/* The name of the fields */

#ifndef CATANEO
#define U    0  // scalar variable
#else
#define U    0
#define V    1
#endif

#ifdef FLUX
#define N_DERIVS 0
#else 
#ifndef CATANEO
#define N_DERIVS 1
#define UX    0
#else
#define N_DERIVS 2
#define UX    0
#define VX    1
#endif




#endif

#ifndef CATANEO
#define N_FLUXES  1   /* usually the same as N_FIELDS */
#define N_AUX 1
#define FU    0
#else
#define N_FLUXES  2   /* usually the same as N_FIELDS */
#define N_AUX 2
#define FU    0
#define FV    1
#endif



//#define B  1         /* this is for the dissipative operators */
/* That is the number of gridpoints times the number of fields      */

#define N_TOTAL ((N_GRIDPOINTS_1)*(N_FIELDS))


/* The total number of plots (used in struct plot) */

#define N_PLOTS 2
#define UU 0    // el lugar para plotear el vinculo.
#define VV 1



/* The number of points at which we want local values */
#define N_POINTS 0



/* The number PI */

#ifdef PI

#undef PI
#endif
#define PI 3.141592653589793238462643

/* --------------------------------------------------------------*/

/* The size of variables */

#define FLOAT double



/* -------------------- FUNCTIONS -------------------------------*/

/* different types of plotting routines, could be:
   graph(er-3d), map(le)3, map(le)5, mtv, idl, sv, fsv or gen 
   this are changes on the type of file where 
   the routine is defined */

#define SDF  /* use to change the plot structure */
#undef SDF
#undef SV
#undef NETCDF
#define ASCHI
#define PYGRAPH

#ifdef SDF
#include "sdf.h"
#define ADISCO adisco_sdf_1d
#define ADISCO_POINT adisco_txt_point_1d
	//#else
	//#define ADISCO adisco_dummy_1d
#endif

#ifdef ASCHI
#define ADISCO adisco_aschi_1d
#endif

#ifdef NETCDF
#define ADISCO adisco_netcdf_1d
#endif

#ifdef SV
#define ADISCO adisco_fsv_1d   /* #define ADISCO adisco_sv_3d */
#endif



/* -----------------------------------------------------------------*/


/* input output definitions */

#undef WEB_INPUT              /* define only one of these two */
#define FILE_INPUT

#ifdef WEB_INPUT
#define OUT_PUT stdout
#endif
#ifdef FILE_INPUT
#define OUT_PUT file_data_ptr
#endif

/* used in routines to get data from the web */
#define MAX_ENTRIES 400
#define MAX_ENTRY 200

/* To print the values of some macros */
/* Warning: only use where OUT_PUT is properly defined */

#define Strg(x) #x
#define PRINT_MACRO_VALUE(x) fprintf(OUT_PUT,"<li> Using %s = %s  </br>\n",#x,Strg(x));
#define GET_MACRO_VALUE(x)   sprintf(macro_value_strg,"%s",Strg(x))

#define MAIL  /* To send e-mails with output data or saying the run finished */
#undef MAIL


/* ----------------------- structures -------------------------*/

#define GRID_PAR grid_1d         /* grid parameters */
#define PLOT_PAR plot_1d         /* ploting parameters */
#define INI_PAR mhd_ini_par_1d  /* where initial data parameters are stored */
#define FUNCTION_PAR mhd_par_1d /* equation parameters */


/* other structures */

#ifdef WEB_INPUT
#define INPUT_FUNCTION wave_1d_input_web
#endif
#ifdef FILE_INPUT
#define INPUT_FUNCTION wave_1d_input_file
#endif


/* -----------------------------------------------------------------*/

/* different functions for the integrator, here goes
   most of the physical input */

#ifdef IMAX 
#define FF ff_eq_F              /* Non stif part of RHS */
#define FS ff_eq_S              /* Stiff part of RHS */
#define FI ff_eq_I              /* Inversion for implicit part */
#else
#define FF ff_eq
#endif

/* different arguments for the function FF */

#define FF_ARG struct grid_1d *grid_1d, struct field_array  *fields, struct field_array  *derivs, struct FUNCTION_PAR *function_par
#define FF_ARG_DUMMY struct GRID_PAR *, struct field_array  *, struct field_array  *, struct FUNCTION_PAR *

/* different functions to take derivatives derivQ_1d, derivQ_3_1d, derivD_1d, derivQQ_1d, deriv_strand_third_order_boundaries_sixth_interior, deriv_strand_fourth_order_boundaries_eight_interior, etc */

#ifdef PERIODIC

//#define DERIV derivD_Per_1d
//#define DERIV derivQ_Per_1d
#define DERIV deriv_6_Per_1d
//#define DERIV deriv_8_Per_1d

#else 

//#define DERIV derivS_1d
#define DERIV deriv_strand_fourth_order_boundaries_eight_interior_1d

#endif // PERIODIC

#define NON_OPTIMIZED_EIGHT_ORDER_OPERATOR
#undef NON_OPTIMIZED_EIGHT_ORDER_OPERATOR

#define OPTIMIZED_EIGHT_ORDER_OPERATOR
//#undef OPTIMIZED_EIGHT_ORDER_OPERATOR

/* -----------------------------------------------------------------*/

/* different dissipative operators diss_KO_4_D_1d, diss_KO_4_00_D_1d */

#ifdef DISSIPATION
#define DISS diss_KO6_Per_1d
#endif



/* different runge-kutta routines */

#ifdef IMAX
//#define RKX imexrkL343 // for this one define IMAX
#define RKX imexrkssp3 // for this one define IMAX paralell with tvd3

#else
//#define RKX rk3
#define RKX tvd3  
//#define RKX rk4
#endif




/* -----------------------------------------------------------------*/

#endif
