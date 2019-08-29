/*********************************************************************
*                                                                    *
* inidat -- provides initial data                                    *
*                                                                    *
* Parameters:                                                        *
* y_a_ptr        -- pointer where to write the initial data          *
* initial_time   -- initial time                                     *
*                                                                    *
* Returns: pointer to field_array where data was writen              *
*                                                                    *
*********************************************************************/

#include "first_macro_1d.h"  /* Where global parameters are defined */
#include "structs_1d.h"      /* Where structures are defined */
//#include "derivs_1d.h"       /* Where derivatives functions are defined */
#include "inidat.h"


/* struct field_array *inidat(struct field_array *y_a_ptr) { */

void inidat(struct field_array *y_a_ptr, 
	    struct GRID_PAR *grid_1d_ptr,
	    struct INI_PAR *ini_par_ptr){

/* -------> grid parameters <--------------------------*/


    int ni_1 = (*grid_1d_ptr).start_grid + N_Ghost; 
    int nf_1 = (*grid_1d_ptr).final_grid - N_Ghost; 
    //int ni_1 = (*grid_1d_ptr).start_grid; 
    //int nf_1 = (*grid_1d_ptr).final_grid;

    FLOAT xi = (*grid_1d_ptr).initial_x;
    FLOAT xf = (*grid_1d_ptr).final_x;
    FLOAT temp, xl, xr;

/* -------> initial data parameters <------------------*/



#ifdef PERIODIC
FLOAT twoPIdN1 = 2.*PI*(xf-xi)/(FLOAT)(nf_1-ni_1);
FLOAT one_dN1 = (xf-xi)/(FLOAT)(nf_1-ni_1);
#else
 FLOAT twoPIdN1 = 2.*PI*(xf-xi)/(FLOAT)(nf_1-ni_1-1);
 FLOAT one_dN1 = (xf-xi)/(FLOAT)(nf_1-ni_1-1);
#endif // PERIODIC



  /* Parameters from main */

  FLOAT a0 = (*ini_par_ptr).a0;                
  FLOAT k_a_10 = (*ini_par_ptr).k_a_10; 
  FLOAT shift_a0 = (*ini_par_ptr).shift_a0;          

/* Amplitude sin(k_a1* x + shift_a1) in U1 */

  FLOAT a1 = (*ini_par_ptr).a1;                
  FLOAT k_a_11 = (*ini_par_ptr).k_a_11;  
  FLOAT shift_a1 = (*ini_par_ptr).shift_a1;          

/* Amplitude of cos(k_a0* x + shift_a0) in U0 */

  FLOAT c0 = (*ini_par_ptr).c0;                
  FLOAT k_c_10 = (*ini_par_ptr).k_c_10;  
  FLOAT shift_c0 = (*ini_par_ptr).shift_c0;          

/* Amplitude of cos(k_c1* x + shift_c0) in U1 */
  FLOAT c1 = (*ini_par_ptr).c1;                
  FLOAT k_c_11 = (*ini_par_ptr).k_c_11;
  FLOAT shift_c1 = (*ini_par_ptr).shift_c1;          

/* Amplitude of b0*exp(-((x-c0_1)^2 + (y-c0_2)^)/sigma_b0) in U0 */
  FLOAT b0 = (*ini_par_ptr).b0;                
  FLOAT sigma_b0 = (*ini_par_ptr).sigma_b0;
  FLOAT c0_1 = (*ini_par_ptr).c0_1;


/* Amplitude of exp(cos(k_b1*x)^2/sigma_b1) in U1 */
  FLOAT b1 = (*ini_par_ptr).b1;                
  FLOAT sigma_b1 = (*ini_par_ptr).sigma_b1;
  FLOAT c1_1 = (*ini_par_ptr).c1_1;

/* Global wave motion */

  FLOAT v1 = (*ini_par_ptr).v1;
  
  
  int initial_data_type = (*ini_par_ptr).initial_data_type;

#ifdef DEBUG
  printf("m = %f\n",m);
#endif

/*---------> values for different fields <-------------*/

/* union fields y; */

/* first the time */

/* y.time = initial_time; */

(*y_a_ptr).time = (*grid_1d_ptr).initial_time;


{register int g_ind1; 
 FLOAT x;
 xl = -0.5;
 xr = 0.5;
 
switch (initial_data_type){

  case 1: // 
   

/* u[0] */

	//for (g_ind1 = ni_1; g_ind1 < nf_1; ++g_ind1) {
		for (g_ind1 = 0; g_ind1 < nf_1 + N_Ghost; ++g_ind1) {

			x= xi + (FLOAT)g_ind1*one_dN1;
		
			//(*y_a_ptr).u[U][g_ind1] = (fabs(x)-0.5) < 0?  1.0:0.0;
			(*y_a_ptr).u[U][g_ind1] = cos(PI*x) - 0.5;
    	}

  break;


  default: //unknown initial data type
    printf("unknown data type =  %d\n", initial_data_type);
    exit(1);
    break;
    
	}    
}


DERIV(grid_1d_ptr, y_a_ptr->u[U], y_a_ptr->u[V]);

//  printf("<LI>Inidat finished </br></LI>\n"); 
}
