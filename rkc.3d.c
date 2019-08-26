/*********************************************************************
*                                                                    *
*  rkc -- integrates    dy/dt = f(y,t)     using runge-kutta         *
*                                                                    *
* Parameters:                                                        *
*      v_init_ptr     -- pointer to where the y initial is           *
*                        and where the y final goes                  *
*      h              -- time step                                   *
*      f              -- pointer to function which computes f(y,t)   *
*                                                                    *
* Returns: nothing                                                   *
*                                                                    *
* WARNING uses global variables for intermediate values              *
*********************************************************************/

#include "first_macro_3d.h"  /* Where global parameters are defined */
#include "structs_3d.h"      /* Where structures are defined */
#include "grid.h"
#include "derivs_3d.h"       /* Where derivatives functions are defined */
#include "gen_3d.h"
#include "globals.h"
#include "rkc.h"
#include "assert.h" 

void rk4(struct field_array *v_init_ptr,
       struct GRID_PAR *grid_ptr,
       struct FUNCTION_PAR *equation_par,
       void (* FF)(FF_ARG_DUMMY)) {

FLOAT h = (*grid_ptr).time_step;
FLOAT x;      /* initial time */    
FLOAT xh;     /* initial time + hh */
FLOAT hh;     /* hh half the h */
FLOAT h6;     /* one sixth of h */
int n_fields = (*grid_ptr).n_fields;
long int n_var = ((*grid_ptr).final_grid_1-(*grid_ptr).start_grid_1)
               *((*grid_ptr).final_grid_2-(*grid_ptr).start_grid_2)
               *((*grid_ptr).final_grid_3-(*grid_ptr).start_grid_3);

/****************** global variables defined above main *****************/


/* struct field_array dv_temp;  */      /* derivative temporal value         */ 
/* struct field_array dv_sum;   */      /* sum values of deriv               */ 
/* struct field_array v_temp;   */      /* intermediate value of v           */ 

/* struct field_array *dv_temp_ptr = &dv_temp;  */  /* derivative temporal value        */ 
/* struct field_array *dv_sum_ptr = &dv_sum;    */  /* sum          values of deriv     */ 
/* struct field_array *v_temp_ptr = &v_temp;    */  /* intermediate value of v          */ 


/******************* end global variables *******************************/


x = (*v_init_ptr).time;
hh = h*(0.5);
h6 = h/(6.0);
xh = x + hh;


//printf("initial time = %f\n",(*v_init_ptr).time);

 FF(grid_ptr,v_init_ptr,&globals.dv_temp,equation_par);       /* take derivatives at initial time */

/* put derivative in dev_sum = (k1 + 2k2 + 2k3 + k4) */
/* dev_sum = k1 */
#ifdef DEBUG
 printf("after first function evaluation in RKC \n");
#endif

 {int field_ind;
 long int ind;
   for (field_ind = 0; field_ind < n_fields; field_ind++){ 
#pragma omp parallel for collapse(1)
       for (ind = 0; ind < n_var; ++ind){
	   globals.dv_sum.u[field_ind][ind] = globals.dv_temp.u[field_ind][ind];
       }
   }
 }


/* asign values to v_temp (v_temp = v_init + hh*dv_init) */

 {int field_ind;
 long int ind;
   for (field_ind = 0; field_ind < n_fields; field_ind++){ 
#pragma omp parallel for collapse(1)
       for (ind = 0; ind < n_var; ++ind){
   globals.v_temp.u[field_ind][ind] = (*v_init_ptr).u[field_ind][ind] 
                          + globals.dv_temp.u[field_ind][ind]*hh;
       }
   }
 }

/* asign time to v_temp */
globals.v_temp.time = xh;

//printf("intermediate time = %f\n",(*v_temp_ptr).time);

/* take again derivative at time xh = t + h/2 and put them in dv_temp */
 FF(grid_ptr,&globals.v_temp,&globals.dv_temp,equation_par);


/* put derivative in dev_sum = (k1 + 2k2 + 2k3 + k4) */
/* dev_sum = k1 + 2k2 */

 {int field_ind;
 long int ind;
   for (field_ind = 0; field_ind < n_fields; field_ind++){ 
#pragma omp parallel for collapse(1)
       for (ind = 0; ind < n_var; ++ind){
	   globals.dv_sum.u[field_ind][ind] = globals.dv_sum.u[field_ind][ind] 
	       + 2. * globals.dv_temp.u[field_ind][ind];
       }
   }
 }



/* asign values to v_temp = v_init + (h/2)*dv_temp  */
 {int field_ind;
 long int ind;
   for (field_ind = 0; field_ind < n_fields; field_ind++){ 
#pragma omp parallel for collapse(1)
       for (ind = 0; ind < n_var; ++ind){
	   globals.v_temp.u[field_ind][ind] = (*v_init_ptr).u[field_ind][ind] 
	       + globals.dv_temp.u[field_ind][ind]*hh;
       }
   }
 }

 /* time does not need to be updated */

 /* takes derivatives at x+h/2, v_temp (new) and put them in dv_med */
 FF(grid_ptr,&globals.v_temp,&globals.dv_temp,equation_par);


/* put derivative in dev_sum = (k1 + 2k2 + 2k3 + k4) */
/* dev_sum = k1 + 2k2 + 2k3 */

 {int field_ind;
 long int ind;
   for (field_ind = 0; field_ind < n_fields; field_ind++){ 
#pragma omp parallel for collapse(1)
       for (ind = 0; ind < n_var; ++ind){
	   globals.dv_sum.u[field_ind][ind] = globals.dv_sum.u[field_ind][ind] 
	       + 2. * globals.dv_temp.u[field_ind][ind];
       }
   }
 }


/* asign values to v_temp = v_init + h*dv_temp */
 {int field_ind;
 long int ind;   
   for (field_ind = 0; field_ind < n_fields; field_ind++){ 
#pragma omp parallel for collapse(1)
     for (ind = 0; ind < n_var; ++ind){
	 globals.v_temp.u[field_ind][ind] = (*v_init_ptr).u[field_ind][ind] 
	     + globals.dv_temp.u[field_ind][ind]*h;
     }
   }
 }

 /* asigns new time to v_temp = x + h */
 globals.v_temp.time = x + h;

//printf("final time = %f\n",(*v_temp_ptr).time);

/* takes derivatives at x+h with new v_temp and puts them in dv_temp */
 FF(grid_ptr,&globals.v_temp,&globals.dv_temp,equation_par);

/* asign value to v_fina = v_init + h6*(dv_temp + dv_temp)*/
 {int field_ind;
 long int ind;   
 for (field_ind = 0; field_ind < n_fields; field_ind++){ 
#pragma omp parallel for collapse(1)
     for (ind = 0; ind < n_var; ++ind){
	 (*v_init_ptr).u[field_ind][ind] = (*v_init_ptr).u[field_ind][ind] 
	     + h6*(globals.dv_sum.u[field_ind][ind] 
		   +  globals.dv_temp.u[field_ind][ind]);
     }
 }
 }

 /* asign returning time */
(*v_init_ptr).time = x + h;
/*  printf("time at RKC = %f",x+h); */
}


/***********************************************************************************
 *                                                                                 *
 *   Third order Runge Kutta                                                       *
 *                                                                                 *
 **********************************************************************************/


void rk3(struct field_array *v_init_ptr,
       struct GRID_PAR *grid_ptr,
       struct FUNCTION_PAR *equation_par,
       void (* FF)(FF_ARG_DUMMY)) {
FLOAT h = (*grid_ptr).time_step;
FLOAT x;      /* initial time */    
FLOAT xh;     /* initial time + hh */
FLOAT hh;     /* hh half the h */
FLOAT h9;     /* one ninth of h */
FLOAT h34;

int n_fields = (*grid_ptr).n_fields;
long int n_var = ((*grid_ptr).final_grid_1-(*grid_ptr).start_grid_1)
               *((*grid_ptr).final_grid_2-(*grid_ptr).start_grid_2)
               *((*grid_ptr).final_grid_3-(*grid_ptr).start_grid_3);

/****************** global variables defined above main *****************/


/* struct field_array dv_temp;  */      /* derivative temporal value         */ 
/* struct field_array dv_sum;   */      /* sum values of deriv               */ 
/* struct field_array v_temp;   */      /* intermediate value of v           */ 

/* struct field_array *dv_temp_ptr = &dv_temp;  */  /* derivative temporal value        */ 
/* struct field_array *dv_sum_ptr = &dv_sum;    */  /* sum          values of deriv     */ 
/* struct field_array *v_temp_ptr = &v_temp;    */  /* intermediate value of v          */ 


/******************* end global variables *******************************/


x = (*v_init_ptr).time;

hh = h*(0.5);
h9 = h/(9.0);
h34 = h*3.0/4.0;
xh = x + hh;




 FF(grid_ptr,v_init_ptr,&globals.dv_temp,equation_par);       /* take derivatives at initial time */

/* put derivative in dev_sum = (2*k1 + 3*k2 + 4*k3) */
/* dev_sum = 2*k1 */
#ifdef DEBUG
 printf("after first function evaluation in RKC \n");
#endif

 {int field_ind;
 long int ind;
   for (field_ind = 0; field_ind < n_fields; field_ind++){ 
#pragma omp parallel for collapse(1)
       for (ind = 0; ind < n_var; ++ind){
		globals.dv_sum.u[field_ind][ind] = 2.0 * globals.dv_temp.u[field_ind][ind];
       }
   }
 }


/* asign values to v_temp (v_temp = v_init + hh*dv_init) */

 {int field_ind;
 long int ind;
   for (field_ind = 0; field_ind < n_fields; field_ind++){ 
#pragma omp parallel for collapse(1)
       for (ind = 0; ind < n_var; ++ind){
   globals.v_temp.u[field_ind][ind] = (*v_init_ptr).u[field_ind][ind] 
                          + hh * globals.dv_temp.u[field_ind][ind];
       }
   }
 }

/* asign time to v_temp */
globals.v_temp.time = xh;

/* take again derivative at time xh = t + h/2 and put them in dv_temp */
 FF(grid_ptr,&globals.v_temp,&globals.dv_temp,equation_par);


/* put derivative in dev_sum = (2*k1 + 3*k2 + 4*k3) */
/* dev_sum = 2*k1 + 3*k2 */

 {int field_ind;
 long int ind;
   for (field_ind = 0; field_ind < n_fields; field_ind++){ 
#pragma omp parallel for collapse(1)
       for (ind = 0; ind < n_var; ++ind){
	   globals.dv_sum.u[field_ind][ind] = globals.dv_sum.u[field_ind][ind] 
	       + 3.0 * globals.dv_temp.u[field_ind][ind];
       }
   }
 }



/* asign values to v_temp = v_init + (3*h/4)*dv_temp  */
 {int field_ind;
 long int ind;
   for (field_ind = 0; field_ind < n_fields; field_ind++){ 
#pragma omp parallel for collapse(1)
       for (ind = 0; ind < n_var; ++ind){
	   globals.v_temp.u[field_ind][ind] = (*v_init_ptr).u[field_ind][ind] 
	       + h34 * globals.dv_temp.u[field_ind][ind];
       }
   }
 }


/* asign time to v_temp */
globals.v_temp.time = x + h34;

 /* takes derivatives at x+3*h/4, v_temp (new) and put them in dv_med */

 FF(grid_ptr,&globals.v_temp,&globals.dv_temp,equation_par);




/* asign value to v_fina = v_init + h9*(dv_sum + 4*dv_temp)*/
 {int field_ind;
 long int ind;   
 for (field_ind = 0; field_ind < n_fields; field_ind++){ 
#pragma omp parallel for collapse(1)
     for (ind = 0; ind < n_var; ++ind){
	 (*v_init_ptr).u[field_ind][ind] = (*v_init_ptr).u[field_ind][ind] 
	     + h9*(globals.dv_sum.u[field_ind][ind] 
		   +  4.0 * globals.dv_temp.u[field_ind][ind]);
     }
 }
 }

 /* asign returning time */
(*v_init_ptr).time = x + h;
/*  printf("time at RKC = %f",x+h); */
}



void tvd3(struct field_array *v_init_ptr, 
	 struct GRID_PAR *grid_ptr,
	 struct FUNCTION_PAR *equation_par,
	 void (* FF)(FF_ARG_DUMMY)) {

/****************** global variables defined above main *****************/

/* union fields dv_temp;  */       /* derivative temporal value         */
/* union fields dv_sum;   */       /* sum values of deriv               */
/* union fields v_temp;   */       /* intermediate value of v           */


/* union fields *dv_temp_ptr; */        /* derivative temporal value        */
/* union fields *dv_sum_ptr;   */       /* sum          values of deriv     */
/* union fields *v_temp_ptr; */         /* intermediate value of v          */

/* dv_temp_ptr = &dv_temp; */           /* derivative temporal value        */
/* dv_sum_ptr = &dv_sum;  */            /* sum          values of deriv     */
/* v_temp_ptr = &v_temp; */             /* intermediate value of v          */

/******************* end global variables *******************************/


/************************************************************************/
/* The scheme is:                                                       */
/* u^1 = u^n + h*F(u^n,t)                                               */
/* u^2 = 3u^n/4 + u^1/4 + h/4*F(u^1,t+h)                                */
/* u^(n+1) = u^n/3 + 2*u^2/3 * 2h/3*F(u^2,t+h/2)                        */
/*                                                                      */
/* ver paper de Luis del 2008                                           */
/************************************************************************/


FLOAT x = (*v_init_ptr).time;
FLOAT h = (*grid_ptr).time_step;
FLOAT h2 = h/2.;
FLOAT h4 = h*(0.250);
FLOAT h3 = h/3.;

int n_fields = (*grid_ptr).n_fields;
long int n_var = ((*grid_ptr).final_grid_1-(*grid_ptr).start_grid_1)
               *((*grid_ptr).final_grid_2-(*grid_ptr).start_grid_2)
               *((*grid_ptr).final_grid_3-(*grid_ptr).start_grid_3);

 FF(grid_ptr,v_init_ptr,&globals.dv_temp,equation_par);       /* take derivatives at initial time */


/* asign values to v_temp (v_temp = v_init + h*dv_init) */

{int field_ind;
 long int ind;
   for (field_ind = 0; field_ind < n_fields; field_ind++){ 
#pragma omp parallel for collapse(1)
       for (ind = 0; ind < n_var; ++ind){
	   globals.v_temp.u[field_ind][ind] = (*v_init_ptr).u[field_ind][ind]
									    + globals.dv_temp.u[field_ind][ind]*h;
       }
   }
 }

/* asign time to v_temp */
globals.v_temp.time = x+h;



/* take again derivative at time xh = t + h and put them in dv_temp */
 FF(grid_ptr,&globals.v_temp,&globals.dv_temp,equation_par);



/* asign values to v_temp = 3*v_init/4 + v_temp/4 + (h/4)*dv_temp  */
 
 {int field_ind;
 long int ind;
   for (field_ind = 0; field_ind < n_fields; field_ind++){ 
#pragma omp parallel for collapse(1)
       for (ind = 0; ind < n_var; ++ind){
	   globals.v_temp.u[field_ind][ind] = 3. * (*v_init_ptr).u[field_ind][ind] / 4.
		                               + globals.v_temp.u[field_ind][ind] / 4.
									   + globals.dv_temp.u[field_ind][ind]*h4;
       }
   }
 }
 
 
/* asign time to v_temp */
globals.v_temp.time = x + h2;

 /* takes derivatives at x+h/2, v_temp (new) and put them in dv_med */
 FF(grid_ptr,&globals.v_temp,&globals.dv_temp,equation_par);


{int field_ind;
 long int ind;
   for (field_ind = 0; field_ind < n_fields; field_ind++){ 
#pragma omp parallel for collapse(1)
       for (ind = 0; ind < n_var; ++ind){
	   (*v_init_ptr).u[field_ind][ind] = (*v_init_ptr).u[field_ind][ind] / 3.
		                               + 2. * globals.v_temp.u[field_ind][ind] / 3.
									   + 2. * globals.dv_temp.u[field_ind][ind] * h3;
       }
   }
 }

 /* asign returning time */
(*v_init_ptr).time = x + h;
/*  printf("time at RKC = %f",x+h); */




}


