
/*********************************************************************
*                                                                    *
* F -- evaluates the function f(y,t), lots of physical imput on it   *
*                                                                    *
* Parameters:                                                        *
*       fields -- pointer to field_vector from where to extract (y,t)*
*       derivs -- pointer to field_vector where to put derivatives   *
*     wave_par -- pointer to parameter struct                        *
*                                                                    *
* Returns: nothing                                                   *
*                                                                    *
*********************************************************************/
#include "first_macro_1d.h"  /* Where global parameters are defined */
#include "structs_1d.h"      /* Where structures are defined */
#include "derivs_1d.h"       /* Where derivatives functions are defined */
#include "equation.h"


/***********************************************************************/

void ff_eq(struct GRID_PAR *grid_1d_ptr,
		struct field_array *fields_ptr, 
		struct field_array *derivs_ptr,
		struct FUNCTION_PAR *function_par) {




    int ni_1 = (*grid_1d_ptr).start_grid + N_Ghost; 
    int nf_1 = (*grid_1d_ptr).final_grid - N_Ghost; 
    int n_gridpts_1 = nf_1 - ni_1;
#ifdef PERIODIC
    int n_mod = nf_1 + N_Ghost; // for the periodic case
#else
   int n_mod = nf_1 + 2*N_Ghost;
#endif

    FLOAT xi = (*grid_1d_ptr).initial_x;
    FLOAT xf = (*grid_1d_ptr).final_x;
    FLOAT dt = (*grid_1d_ptr).time_step;
    
#ifdef PERIODIC
    FLOAT one_dN1 = (xf-xi)/(FLOAT)(nf_1-ni_1);
    FLOAT h_1 = (FLOAT)(nf_1-ni_1)/(xf-xi);
#else
    FLOAT one_dN1 = (xf-xi)/(FLOAT)(nf_1-ni_1-1);
    FLOAT h_1 = (FLOAT)(nf_1-ni_1-1)/(xf-xi);
#endif
    FLOAT time;
    FLOAT x;
#ifdef FLUX	
    FLOAT u_pp[N_FIELDS], u_pm[N_FIELDS], u_mp[N_FIELDS], u_mm[N_FIELDS];
    FLOAT Dp[N_FIELDS], Dpp[N_FIELDS], Dm[N_FIELDS], Dmm[N_FIELDS];
    FLOAT v_x0[N_FIELDS], v_xp[N_FIELDS], v_xm[N_FIELDS];
    FLOAT a_p, a_m, a_pp, a_pm, a_mp, a_mm;
    FLOAT H_p[N_FIELDS], H_m[N_FIELDS];
	//    FLOAT theta = 1.;   // factor for limiter
	//    FLOAT theta = 2.;   // factor for limiter
    FLOAT theta = 1.5;   // factor for limiter
#endif 

#ifdef SOURCE
	FLOAT u[N_FIELDS];
#endif

#ifdef F_DIFF
	FLOAT u[N_FIELDS], u_ext[N_FIELDS];
	FLOAT RP=1., RM = 0.;

#ifdef PENALTY

	FLOAT factor, L;
  	char macro_value_strg[100];

    GET_MACRO_VALUE(DERIV);
    if (strcmp(macro_value_strg,"derivS_1d")==0) {
      factor = 2.0;
    }
    else if (strcmp(macro_value_strg,"derivQ_1d")==0) {
      factor = 48.0/17.0;
    }
    else if (strcmp(macro_value_strg,"derivQ_3_1d")==0) {
      factor = 11.0/3.0;
    }
    else if (strcmp(macro_value_strg,"deriv_strand_third_order_boundaries_sixth_interior_1d")==0) {
      factor = 43200.0/13649.0;
    }
    else if (strcmp(macro_value_strg,"deriv_strand_fourth_order_boundaries_eight_interior_1d")==0) {
      factor = 5080320.0/1498139.0;
    }
    else {
      factor=2.0; printf("check factor por penalty!!!%s!!!",macro_value_strg);
    }
    
    L=(FLOAT)(nf_1-ni_1-1)/(xf-xi);
  	L=factor*L;
	
	
	FLOAT U0[N_FIELDS], UP[N_FIELDS], UM[N_FIELDS], T0[N_FIELDS], TP[N_FIELDS], TM[N_FIELDS], V0, VP, VM;
   
   	{int i;
   	for (i = 0; i < N_FIELDS; ++i){
		   U0[i]=0.0;
		   UP[i]=0.0;
		   UM[i]=0.0;
		   T0[i]=0.0;
		   TP[i]=0.0;
		   TM[i]=0.0;
	   }
	}

#endif // PERIODIC

#endif


    FLOAT sx,sy,e,bx,by,B2,S2,div;

    
    FLOAT s = (*function_par).s; // Value for By
    FLOAT a = (*function_par).a; // Value for Bx
    FLOAT c = (*function_par).c; // Value for Sx

    FLOAT sigma = (*function_par).sigma;
    

    /* normals */

    FLOAT nx_10 = (*function_par).nx_10;
    FLOAT nx_11 = (*function_par).nx_11;

    FLOAT nx;


    
    
  /* first the time */

  time = (*fields_ptr).time; 
  (*derivs_ptr).time = (*fields_ptr).time;




  /* ghost points */

/*
  {register int i, grid_ind1;
	for (i = 0; i < N_FIELDS; ++i){
		for (grid_ind1 = 0; grid_ind1 < N_Ghost; ++grid_ind1){
			(*fields_ptr).u[i][grid_ind1] = (*fields_ptr).u[i][ni_1];
			(*fields_ptr).u[i][nf_1 + grid_ind1] = (*fields_ptr).u[i][nf_1-1];
		}
	}
  }
  */

 // This is fixed in inidat.c
 /*
	{register grid_ind1;
  		for (grid_ind1 = 0; grid_ind1 < N_Ghost; ++grid_ind1){
			(*fields_ptr).u[BY][grid_ind1] = (*fields_ptr).u[BY][ni_1];
			(*fields_ptr).u[BY][nf_1 + grid_ind1] = (*fields_ptr).u[BY][nf_1-1];
			(*fields_ptr).u[SX][grid_ind1] = 0.20;
			(*fields_ptr).u[SX][nf_1 + grid_ind1] = -0.20;
			(*fields_ptr).u[SY][grid_ind1] = (*fields_ptr).u[SY][ni_1];
			(*fields_ptr).u[SY][nf_1 + grid_ind1] = (*fields_ptr).u[SY][nf_1-1];
		}
	}
	*/

  /* inner points */
#ifdef FLUX    
	{register int grid_ind1, i;
//#pragma omp parallel for 

        // here we go one point inside the ghost zone (the -1 in the second line reaches the first point)
#ifdef PERIODIC
		for (grid_ind1 = ni_1; grid_ind1< nf_1; ++grid_ind1){
#else             
		for (grid_ind1 = ni_1+1; grid_ind1< nf_1-1; ++grid_ind1){ 
#endif
			for (i = 0; i < N_FIELDS; ++i){ 
                
                Dp[i] = (*fields_ptr).u[i][(grid_ind1+1) % n_mod] - (*fields_ptr).u[i][grid_ind1];
//				Dpp[i] = (*fields_ptr).u[i][(grid_ind1+2) % n_mod] - (*fields_ptr).u[i][(grid_ind1+1) % n_mod];
				Dm[i] = (*fields_ptr).u[i][grid_ind1] - (*fields_ptr).u[i][(grid_ind1-1 + n_mod) % n_mod];
//				Dmm[i] = (*fields_ptr).u[i][(grid_ind1-1 + n_mod) % n_mod] - (*fields_ptr).u[i][(grid_ind1-2 + n_mod) % n_mod];
				globals.auxfields.u_aux[i][grid_ind1] = 0.5*h_1*(mysign_zero(Dp[i])+mysign_zero(Dm[i]))*MM3(Dp[i]+Dm[i],Dp[i],Dm[i],theta);
//				v_xp[i] = 0.5*h_1*(mysign_zero(Dpp[i])+mysign_zero(Dp[i]))*MM3(Dpp[i]+Dp[i],Dpp[i],Dp[i],theta);
//				v_xm[i] = 0.5*h_1*(mysign_zero(Dm[i])+mysign_zero(Dmm[i]))*MM3(Dm[i]+Dmm[i],Dm[i],Dmm[i],theta);
            }
        }
    }
                

	{register int grid_ind1, i;
//	#pragma omp parallel for 

#ifdef PERIODIC
		for (grid_ind1 = ni_1; grid_ind1< nf_1; ++grid_ind1){
#else        
		for (grid_ind1 = ni_1+1; grid_ind1< nf_1-1; ++grid_ind1){ 
#endif
			for (i = 0; i < N_FIELDS; ++i){ 

				v_x0[i] = globals.auxfields.u_aux[i][grid_ind1];
                v_xp[i] = globals.auxfields.u_aux[i][(grid_ind1+1) % n_mod];
                v_xm[i] = globals.auxfields.u_aux[i][(grid_ind1 + n_mod -1) % n_mod];

				u_pp[i] = (*fields_ptr).u[i][(grid_ind1+1) % n_mod] - 0.5*one_dN1*v_xp[i];
				u_pm[i] = (*fields_ptr).u[i][grid_ind1]   + 0.5*one_dN1*v_x0[i];
				u_mp[i] = (*fields_ptr).u[i][grid_ind1]   - 0.5*one_dN1*v_x0[i];
				u_mm[i] = (*fields_ptr).u[i][(grid_ind1 + n_mod -1) % n_mod] + 0.5*one_dN1*v_xm[i];
				
#ifdef SOURCE
				u[i] = (*fields_ptr).u[i][grid_ind1];
#endif
			}
			
                

/* 
		// maximum propagation speed is 1
*/
			a_p = 1.0;
			a_m = 1.0;
		
	 
			for (i = 0; i < N_FIELDS; ++i){ 
		
				H_p[i] = 0.5*(Fx(u_pp, grid_1d_ptr, function_par,i) + Fx(u_pm, grid_1d_ptr, function_par,i));
				H_p[i] = H_p[i] - 0.5*a_p*(u_pp[i]-u_pm[i]);
		
				H_m[i] = 0.5*(Fx(u_mp, grid_1d_ptr, function_par,i) + Fx(u_mm, grid_1d_ptr, function_par,i));
				H_m[i] = H_m[i] - 0.5*a_m*(u_mp[i]-u_mm[i]);
		
				(*derivs_ptr).u[i][grid_ind1] = -h_1*(H_p[i]-H_m[i]);
				
#ifdef SOURCE
				(*derivs_ptr).u[i][grid_ind1] = (*derivs_ptr).u[i][grid_ind1] + Source(u,grid_1d_ptr, function_par,i);
#endif				
				
			}

   }
   
   //now we set the first order KT operator at the borders
   #ifndef PERIODIC
   
   //first we average both sides
   
   (*fields_ptr).u[U][ni_1] = ((*fields_ptr).u[U][ni_1] + (*fields_ptr).u[U][nf_1-1])*0.5;
   (*fields_ptr).u[U][nf_1-1] = (*fields_ptr).u[U][ni_1];
   
   //first point
   grid_ind1 = ni_1;
   			a_p = 1.0;
		
   for (i = 0; i < N_FIELDS; ++i){ 
		u_pp[i] = (*fields_ptr).u[i][(grid_ind1+1)];
		u_pm[i] = (*fields_ptr).u[i][grid_ind1];
	}
	for (i = 0; i < N_FIELDS; ++i){ 
		H_p[i] = Fx(u_pp, grid_1d_ptr, function_par,i);
		H_p[i] = H_p[i] - a_p*u_pp[i];
	
		H_m[i] = Fx(u_pm, grid_1d_ptr, function_par,i);
		H_m[i] = H_m[i] - a_p*u_pm[i];
		
		(*derivs_ptr).u[i][grid_ind1] = -h_1*(H_p[i]-H_m[i]);
	}
	
	grid_ind1 = nf_1-1;
    
			a_m = 1.0;
   for (i = 0; i < N_FIELDS; ++i){ 
		u_pp[i] = (*fields_ptr).u[i][(grid_ind1)];
		u_pm[i] = (*fields_ptr).u[i][grid_ind1-1];
	}
	for (i = 0; i < N_FIELDS; ++i){ 
		H_p[i] = Fx(u_pp, grid_1d_ptr, function_par,i);
		H_p[i] = H_p[i] + a_m*u_pp[i];
	
		H_m[i] = Fx(u_pm, grid_1d_ptr, function_par,i);
		H_m[i] = H_m[i] + a_m*u_pm[i];
		
		(*derivs_ptr).u[i][grid_ind1] = -h_1*(H_p[i]-H_m[i]);

		
	}
	
   #endif
   
 }

#endif // FLUX
 
#ifdef F_DIFF

// First upload the fluxes in arrays.
	{register int grid_ind1, i;

		for (grid_ind1 = ni_1; grid_ind1< nf_1; ++grid_ind1){

			for (i = 0; i < N_FIELDS; ++i){ 
				u[i] = (*fields_ptr).u[i][(grid_ind1)];
			}

			for (i = 0; i < N_FLUXES; ++i){                
				globals.auxfields.u_aux[i][grid_ind1] = -Fx(u, grid_1d_ptr, function_par,i);
            }
        }

		// Take derivatives of fluxes

		for (i = 0; i < N_DERIVS; ++i){	
			DERIV(grid_1d_ptr, globals.auxfields.u_aux[i], globals.dfields.du[i]);
		}	
	}	

	{register int grid_ind1, i;
//	#pragma omp parallel for 

		for (grid_ind1 = ni_1; grid_ind1< nf_1; ++grid_ind1){

			for (i = 0; i < N_FIELDS; ++i){ 
				(*derivs_ptr).u[i][grid_ind1] = globals.dfields.du[i][grid_ind1];
				
			}
		}

   }
#endif

}

/*************************************************************************************/

 static inline FLOAT MM3(double a, double b, double c, double weight){	// (2*D0,Dp,Dm)
  weight = weight*2.;
  
  if (fabs(a) <= (weight*fabs(b))) {
  return (fabs(a) <= (weight*fabs(c))) ? fabs(a)*.5 : fabs(c);} 
  else {
  return (fabs(b) <= fabs(c)) ? fabs(b) : fabs(c);
	}
	}
	
 //double MM2(double a, double b, double weight){	
  //if (fabs(a) <= weight*fabs(b)) {
  //return fabs(a) <= weight*fabs(c) ? fabs(a) : fabs(c);} 
  //else {
  //return fabs(b) <= fabs(c) ? fabs(b) : fabs(c);
	//}
	//}

static inline FLOAT mysign_zero_m_one(double d){
    if (d > 0.) {return 1.0;}
	else return -1.0;
}

static inline FLOAT mysign_zero(double d){
    if (d > 0.) {return 1.0;}
	else if (d < 0.) return -1.0;
	else return 0.0;
}



/* -------------- Fluxes -------------------------*/

// WARNING HARDCODED Bx = a.

static inline FLOAT  Fx(double *u, struct GRID_PAR *grid_1d_ptr, struct FUNCTION_PAR *function_par,int i){
	FLOAT c = (*function_par).c;
	FLOAT s = (*function_par).s;
		
//	FLOAT a = (*function_par).a;

		 switch(i)
    {
    case U:    return -(0.5*u[U]*u[U] + s*u[V]); break;
	case V:    return -c*u[U] ; break;
	default:   return printf("out of range in funtion Fx"); exit(0); break;
	}
	return u[U];
	
}


static inline FLOAT  Source(double *u, struct GRID_PAR *grid_1d_ptr, struct FUNCTION_PAR *function_par,int i){
	FLOAT c = (*function_par).c;
//	FLOAT s = (*function_par).s;
		
//	FLOAT a = (*function_par).a;

		 switch(i)
    {
    case U:    return 0.0; break;
	case V:    return -c*u[V] ; break;
	default:   return printf("out of range in funtion Fx"); exit(0); break;
	}
	return u[U];
	
}



