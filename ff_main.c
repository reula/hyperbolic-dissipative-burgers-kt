#define _POSIX_C_SOURCE 200112L
#include <stdlib.h>
#include <sys/time.h>
#include <fenv.h>
	

#include "first_macro_1d.h"  /* Where global parameters are defined */
#include "structs_1d.h"      /* Where structures are defined */
//#include "derivs_1d.h"       /* Where derivatives functions are defined */
//#include "gen_1d.h"
#include "globals.h"
#include "inidat.h"
#include "input.h"
#include "equation.h" 	
#include "integ.h"
#include "rkc.h"
#include "adisco_1d.h" 

// 4M hugepage boundary
#define HUGEPAGE_SIZE (1 << 22)

/***********************   Global variables   ****************************/
struct globals globals;

/***********************   Helper functions   ******************************/

static inline void get_max_V(struct field_array *y, struct GRID_PAR *grid, FLOAT *max);
static inline void get_max_DU(struct field_array *y, struct GRID_PAR *grid, FLOAT *max);
static inline FLOAT norm_Energy(struct GRID_PAR *grid_1d_ptr,
        struct FUNCTION_PAR *function_par_ptr,
		struct field_array  *fields_ptr);

/* Allocate a set of fields at once in huge pages if possible */
void alloc_field_data(size_t nfields, size_t field_elems, FLOAT ** data) {
    FLOAT * memblock;

    size_t alloc_size = nfields * field_elems * sizeof(FLOAT);
    // Request a hugepage-aligned chunk of memory
    if (posix_memalign((void **) &memblock, HUGEPAGE_SIZE, alloc_size) != 0) {
        fprintf(stderr, "out of memory in posix_memaling\n");
        exit(1);
    }

#ifdef MADV_HUGEPAGE
    // Give the OS a hint
    madvise(memblock, alloc_size, MADV_HUGEPAGE);
#endif

    // Copy pointers to each field
    for (size_t i = 0; i < nfields; ++i) {
        data[i] = &memblock[i * field_elems];
    }
}

/* Free a set of fields allocated with alloc_field_data */
void free_field_data(FLOAT ** data) {
    free(*data);
}

void * safe_malloc(size_t size) {
    void * rv = malloc(size);
    if (rv == NULL) {
        fprintf(stderr, "out of memory in safe_malloc\n");
        exit(1);
    }
    return rv;
}

/***********************   Global variables   ****************************/


/*
void norm_L2(struct GRID_PAR *grid_1d_ptr, 
		struct FUNCTION_PAR *function_par_ptr,
		struct field_array  *fields_ptr);
		
void norm_Energy(struct GRID_PAR *grid_1d_ptr,
        struct FUNCTION_PAR *function_par_ptr,
		struct field_array  *fields_ptr);

*/

int main() {

  /* variable declarations */

  struct GRID_PAR grd;
  struct GRID_PAR *grd_ptr = &grd;
  FLOAT h;
  struct field_array y;


  /*   Ploting names */

  struct PLOT_PAR plot; 
  struct PLOT_PAR *plot_ptr = &plot;

  /*   Initial data parameters */

  struct INI_PAR init_parameters;
  struct INI_PAR *init_par_ptr = &init_parameters;
  
  
  /*   Function parameters  */


  struct FUNCTION_PAR equation_parameters;
  struct FUNCTION_PAR *equation_par_ptr = &equation_parameters;

  /* Parameters coming from first_macro */


#ifdef FILE_INPUT
  FILE *file_data_ptr;
#endif

// For measuring time 

struct timeval start, stop, end;
double tiempo;





  /* Get data from web page or data file */
 



INPUT_FUNCTION(grd_ptr, equation_par_ptr, 
		 init_par_ptr, plot_ptr);

 printf("out of input function\n");

file_data_ptr = plot_ptr->input_data_file_ptr;
PRINT_MACRO_VALUE(RKX)
PRINT_MACRO_VALUE(DERIV)
PRINT_MACRO_VALUE(FF)
#ifdef IMAX
PRINT_MACRO_VALUE(FF)
PRINT_MACRO_VALUE(FS)
PRINT_MACRO_VALUE(FI)
#else
PRINT_MACRO_VALUE(FF)
#endif //IMAX
PRINT_MACRO_VALUE(ADISCO)
#ifdef DISSIPATION
PRINT_MACRO_VALUE(DISS)
#endif

    /* ------------------------------------ Allocate memory ------------------------------------ */

    /* Allocation #1:------->  auxiliary fields for rkc */

    alloc_field_data(N_FIELDS, grd.n_grid_pts, globals.dv_temp.u);
    alloc_field_data(N_FIELDS, grd.n_grid_pts, globals.dv_sum.u);
    alloc_field_data(N_FIELDS, grd.n_grid_pts, globals.v_temp.u);
    
    
    /* Allocation #2: ------> The fields */
    alloc_field_data(N_FIELDS, grd.n_grid_pts, y.u);

    /* Allocation #3:-------> their derivatives (used in evaluating the function) */
    alloc_field_data(N_DERIVS, grd.n_grid_pts, globals.dfields.du);

    /* Allocation #4:--------->  other auxiliary fields */
    alloc_field_data(N_AUX, grd.n_grid_pts, globals.auxfields.u_aux);

    /* Allocation #5 :---------> plot vector memory */
//    plot_ptr->plot_field = (double *)safe_malloc((plot_ptr->grid_plot_pts_1) * sizeof(double));

//exit(0);

#ifdef SDF
    plot_ptr->plot_field = safe_malloc((plot_ptr->grid_plot_pts_1) * sizeof(double));
#endif
#ifdef PYGRAPH
	plot_ptr->plot_field_pygraph = safe_malloc((plot_ptr->grid_plot_pts_1) * sizeof(float)*2);
#endif	
    

/* Array of names for plots */

 

/* Relation between fields and plot names */

 sprintf(plot_ptr->name[U],"U");
 sprintf(plot_ptr->name[V],"V");
  
    sprintf(plot_ptr->window_name[U], "%s_U_%d"
	  ,plot_ptr->output_file_name,grd.n_grid_pts);
   sprintf(plot_ptr->window_name[V], "%s_V_%d"
	  ,plot_ptr->output_file_name,grd.n_grid_pts);
	  
  plot_ptr->initial_x= grd.initial_x;
  plot_ptr->final_x = grd.final_x;


  plot_ptr->n_plots = N_PLOTS;

  
/* Relation between fields and plot names */

    plot_ptr->pointers[U] = U;
    plot_ptr->pointers[V] = V;

  plot_ptr->initial_time = grd.initial_time;       
  plot_ptr->final_time = grd.final_time;   
  plot_ptr->time_slice = 0;
  //plot_ptr->grid_plot_factor_1 = 1;

 

  /* Open output file (some times used only for compatibility) */
	printf("opening file\n");


			plot_ptr = ADISCO('O',  &plot, &grd, &y); 
//			plot_ptr = adisco_aschi_1d('O',  &plot, &grd, &y);
			
//			exit(0);
	
//  ADISCO_POINT('O', &plot, &grd, &y);



  /* makes initial time-interval */

 
  h = grd.time_step; // to be used in the equations.

  /* sends input data to file / browser */

fprintf(OUT_PUT,"<li> Total number of Time Steps = %f </br>\n",(double)(grd.data_steps*grd.int_steps));

fprintf(OUT_PUT,"<li> Number of Time Steps per unit time = %f </br>\n",1.0/h);

fprintf(OUT_PUT,"<li> Time_Step / Space_Step_x = (h/(xf-xi)*(n_gridpts-1)) = %f </br>\n",h*(double)grd.n_grid_pts/(grd.final_x-grd.initial_x));
fprintf(OUT_PUT,"</ul>%c</br>",10);



fflush(stdout);

/* send input data to the screen */

#ifdef FILE_INPUT
printf("Total number of Time Steps = %f \n",(double)(grd.data_steps*grd.int_steps));

printf("Number of Time Steps per unit time = %f \n",1.0/h);

printf("Time_Step / Space_Step_x= (h/(xf-xi)*n_gridpts) = %f \n",h*(double)grd.n_grid_pts/(grd.final_x-grd.initial_x));
printf("\n");

fflush(stdout);
#endif

// start counting time -----------------------------------------------------

gettimeofday(&start, NULL);



// *******************************************************
// BIG LOOP START
// *******************************************************

#ifdef BIG_LOOP
int loop, big_loop = 40;
#else
int loop, big_loop = 1;
#endif

for (loop = 1; loop <= big_loop; loop++){

FLOAT V_max = 0.0;
FLOAT dU_max = 0.0;



// *******************************************************

/*     creates initial data                            */

inidat(&y,grd_ptr,&init_parameters);

/* write initial data to file                          */ 

/* plot data */
      
	
 

//    	plot_ptr = ADISCO('P',  &plot, &grd, &y); 

	
//    ADISCO_POINT('P',  &plot, &grd, &y); 
	
		plot_ptr = adisco_aschi_1d('P',  &plot, &grd, &y);

		plot_ptr = adisco_pygraph_1d('W',  &plot, &grd, &y);


  /*     creates potential function */

		// exit(0); 

/* inipot(pot_ptr,&pot_parameters, &grd, plot_ptr); */


#ifdef BUMP_TEST
exit(0);
#endif



  /* Take data_steps */ 
  


  {long int k_outer;
  for (k_outer=1; k_outer<= grd.data_steps; k_outer++) {
/*       printf("h = %f\n",h);  */
/*       printf("time = %f\n",y.a.time); */

#ifdef IMAX
integ(&y,grd_ptr,equation_par_ptr,FF,FS,FI,RKX); 
#else
integ(&y,grd_ptr,equation_par_ptr,FF,RKX); 
#endif //IMAX


//printf("time after integ in main = %f",y.time); 
//printf('time = %f, V_max = % f, dUdx_max = %f \n', time, V_max, dUdx_max);
//fflush(stdout);
/* printf("Out of integ \n");  */

get_max_V(&y, grd_ptr, &V_max);
get_max_DU(&y, grd_ptr, &dU_max);

#ifndef BIG_LOOP
printf("time = %f, Energy = %f, V_max = %f, dU_max = %f \n", y.time, norm_Energy(grd_ptr, equation_par_ptr, &y), V_max, dU_max);
#endif

      //printf("...");
      fflush(stdout);
      /* Do pointwise output */
//      ADISCO_POINT('P',  &plot, &grd, &y); 
      /* Do 1d output */
      if ((k_outer%grd.factor_1d_steps)==0){
			  #ifdef SDF
			  		  plot_ptr = ADISCO('P', &plot, &grd, &y);  
			  #endif
		  plot_ptr->time_slice = k_outer;
#ifndef BIG_LOOP
		  plot_ptr = adisco_pygraph_1d('A',  &plot, &grd, &y);
		  plot_ptr = ADISCO('P', &plot, &grd, &y);
#endif
		
      }

	 


/* printf("�</br>\n");  */ 
//printf("r"); 
/* printf("%c",7); */ 
fflush(stdout);
  }



  }
	
//  plot_ptr = adisco_aschi_1d('P', &plot, &grd, &y); // printing the last value in aschi at .dat files.

//  plot_ptr = ADISCO('P', &plot, &grd, &y);
//		  norm_Energy(&grd, equation_par_ptr,&y);


// ************************************************
// BIG LOOP END
// ************************************************

#ifdef BIG_LOOP
printf(" %f,  %f, %f \n", equation_parameters.s, V_max, dU_max);

equation_parameters.s = equation_parameters.s + 0.015; 

#endif //BIG_LOOP

}
// ************************************************


gettimeofday(&end, NULL);

tiempo = ((double) end.tv_sec + end.tv_usec / 1000000.0) - ((double) start.tv_sec + start.tv_usec / 1000000.0);


fprintf(OUT_PUT,"<ul>%c</br>",10);
fprintf(OUT_PUT,"<li> Execution time = %u secs. Time = %u", (unsigned)(clock()/CLOCKS_PER_SEC), tiempo);
fprintf(OUT_PUT,"</ul>%c</br>",10);

#ifdef FILE_INPUT
printf("\n");
printf("Execution (CPU) time = %u secs. Time = %u", (unsigned)(clock()/CLOCKS_PER_SEC), tiempo);
printf("\n");
#endif

/* close output file */
//    plot_ptr = ADISCO('C',  &plot, &grd, &y); 
//    ADISCO_POINT('C', &plot, &grd, &y);
//	plot_ptr = adisco_aschi_1d('C',  &plot, &grd, &y); 

#ifdef FILE_INPUT
  fclose(plot_ptr->input_data_file_ptr);
#endif
printf("%c",7);
printf("finishing \n");
return(0);
}

/*
static inline void get_max(char inst, struct field_array *y, struct GRID_PAR *grid_1d_ptr, FLOAT *max){
		int ni_1 = (*grid_1d_ptr).start_grid; 
		int nf_1 = (*grid_1d_ptr).final_grid; 
//		FLOAT xi = (*grid_1d_ptr).initial_x;
//		FLOAT xf = (*grid_1d_ptr).final_x;
//		FLOAT dt = (*grid_1d_ptr).time_step;
//		FLOAT h_1 = (FLOAT)(nf_1-ni_1)/(xf-xi);
		int grid_ind1;

	switch (inst) {// switch
		case 'V':{// case 'V' 
			for (grid_ind1 = ni_1; grid_ind1 < nf_1; ++grid_ind1){
					if (fabs(y->u[V][grid_ind1]) > (*max)){
						(*max) = (FLOAT)fabs(y->u[V][grid_ind1]);
					}
			}
		}
		break;
		
		case 'D':{
			for (grid_ind1 = ni_1; grid_ind1< nf_1; ++grid_ind1){
					if (fabs(y->u[U][(grid_ind1+1) % nf_1] - y->u[U][grid_ind1] )  > (*max)){(*max) = fabs(y->u[U][(grid_ind1+1) % nf_1] - y->u[U][grid_ind1] );}
			}
		}
		break;
	}
}
* */

static inline void get_max_V(struct field_array *y, struct GRID_PAR *grid_1d_ptr, FLOAT *max){
		int ni_1 = (*grid_1d_ptr).start_grid; 
		int nf_1 = (*grid_1d_ptr).final_grid; 
//		FLOAT xi = (*grid_1d_ptr).initial_x;
//		FLOAT xf = (*grid_1d_ptr).final_x;
//		FLOAT dt = (*grid_1d_ptr).time_step;
//		FLOAT h_1 = (FLOAT)(nf_1-ni_1)/(xf-xi);
		int grid_ind1;
 
			for (grid_ind1 = ni_1; grid_ind1 < nf_1; ++grid_ind1){
					if (fabs(y->u[V][grid_ind1]) > (*max)){
						(*max) = (FLOAT)fabs(y->u[V][grid_ind1]);
					}
			}
		}

static inline void get_max_DU(struct field_array *y, struct GRID_PAR *grid_1d_ptr, FLOAT *max){
		int ni_1 = (*grid_1d_ptr).start_grid; 
		int nf_1 = (*grid_1d_ptr).final_grid; 
//		FLOAT xi = (*grid_1d_ptr).initial_x;
//		FLOAT xf = (*grid_1d_ptr).final_x;
//		FLOAT dt = (*grid_1d_ptr).time_step;
//		FLOAT h_1 = (FLOAT)(nf_1-ni_1)/(xf-xi);
		int grid_ind1;


			for (grid_ind1 = ni_1; grid_ind1< nf_1; ++grid_ind1){
					if (fabs(y->u[U][(grid_ind1+1) % nf_1] - y->u[U][grid_ind1] )  > (*max)){(*max) = fabs(y->u[U][(grid_ind1+1) % nf_1] - y->u[U][grid_ind1] );}
			}
		}

static inline FLOAT norm_Energy(struct GRID_PAR *grid_1d_ptr,
        struct FUNCTION_PAR *function_par_ptr,
		struct field_array  *y){
			
		int ni_1 = (*grid_1d_ptr).start_grid; 
		int nf_1 = (*grid_1d_ptr).final_grid; 
		FLOAT xi = (*grid_1d_ptr).initial_x;
		FLOAT xf = (*grid_1d_ptr).final_x;
		FLOAT dt = (*grid_1d_ptr).time_step;
		FLOAT h_1 = (FLOAT)(nf_1-ni_1)/(xf-xi);
		int grid_ind1;

		FLOAT E = 0.0;;
		
			for (grid_ind1 = ni_1; grid_ind1< nf_1; ++grid_ind1){
				E = E + function_par_ptr->s * y->u[U][grid_ind1]*y->u[U][grid_ind1] + function_par_ptr->c * y->u[U][grid_ind1]*y->u[U][grid_ind1];
			}
			
		return(E / h_1);
		}
		
