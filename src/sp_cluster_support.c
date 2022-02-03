/**
 * (c) Author:    Woongkyu Jee, woong.jee.16@ucl.ac.uk, wldndrb1@gmail.com
 * Created:   02.06.2019 ~
 * 	
 * University College London, Department of Chemistry
 **/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include"sp_cluster_type.h"

#define SP_SUPPORT_TRUE 1
#define SP_SUPPORT_FALSE -1

#define DEBUG_SUPPORT

/* Get the Lowest energy state */
int sp_cluster_support_get_lowest_state( gsl_vector* v )
{
    int Return = 0;
    double min;

    min = gsl_vector_get(v,0);
    for(int i=0;i<4;i++)
    {
        if( gsl_vector_get(v,i) < min  )
        {
            min = gsl_vector_get(v,i);
            Return = i;
        }
    }
    return Return;
    // Returns the index of element in gsl_vector* eval
}

void sp_cluster_support_sign_gs_eigenvector( void* sp_sys_void )
{
    sp_cluster_system* sp_sys = (sp_cluster_system*)sp_sys_void;
    int low_idx;

    for(int i=0;i<sp_sys->number_of_sp_ion;i++)
    {
        low_idx = sp_cluster_support_get_lowest_state(sp_sys->sp_ion[i].eigen_value);
        if( gsl_matrix_get(sp_sys->sp_ion[i].eigen_vector,0,low_idx) < 0. )          // eigen_vector : matrix(4,4) // eigen_vector_gs : vector(4)
        {
            gsl_matrix_set(sp_sys->sp_ion[i].eigen_vector,0,low_idx,-1.*gsl_matrix_get(sp_sys->sp_ion[i].eigen_vector,0,low_idx));
            gsl_matrix_set(sp_sys->sp_ion[i].eigen_vector,1,low_idx,-1.*gsl_matrix_get(sp_sys->sp_ion[i].eigen_vector,1,low_idx));
            gsl_matrix_set(sp_sys->sp_ion[i].eigen_vector,2,low_idx,-1.*gsl_matrix_get(sp_sys->sp_ion[i].eigen_vector,2,low_idx));
            gsl_matrix_set(sp_sys->sp_ion[i].eigen_vector,3,low_idx,-1.*gsl_matrix_get(sp_sys->sp_ion[i].eigen_vector,3,low_idx));
        }
    }
    return;
}

void sp_cluster_support_load_gs_eigenvector( void* sp_sys_void )
{
    sp_cluster_system* sp_sys = (sp_cluster_system*)sp_sys_void;
    int low_idx;
    for(int i=0;i<sp_sys->number_of_sp_ion;i++)
    {   low_idx = sp_cluster_support_get_lowest_state(sp_sys->sp_ion[i].eigen_value);

        gsl_vector_set(sp_sys->sp_ion[i].eigen_vector_gs,0,gsl_matrix_get(sp_sys->sp_ion[i].eigen_vector,0,low_idx));
        gsl_vector_set(sp_sys->sp_ion[i].eigen_vector_gs,1,gsl_matrix_get(sp_sys->sp_ion[i].eigen_vector,1,low_idx));
        gsl_vector_set(sp_sys->sp_ion[i].eigen_vector_gs,2,gsl_matrix_get(sp_sys->sp_ion[i].eigen_vector,2,low_idx));
        gsl_vector_set(sp_sys->sp_ion[i].eigen_vector_gs,3,gsl_matrix_get(sp_sys->sp_ion[i].eigen_vector,3,low_idx));
    }
    return;
}


/* Matrix Viewer: Print out a N x N matrix on console */
void sp_cluster_support_matrix_view( const gsl_matrix* m )
{
    if( m != NULL )
    {   for(int i=0;i<m->size1;i++)
        {   for(int j=0;j<m->size2;j++)
            {   //printf("%s%12.4lf",gsl_matrix_get(m,i,j)>0.?"+":"",gsl_matrix_get(m,i,j));
                printf("%12.4e",gsl_matrix_get(m,i,j));
            }
        puts("");
        }
    }
    else
        puts("sp_cluster_support_matrix_view input 'm' (gsl_matrix*) is a null pointer ... in SP_Support.c or SP_Support.h");

    return;
}

/* Vector Viewer: Print out a vector on console */
void sp_cluster_support_vector_view( const gsl_vector* v )
{
    if( v != NULL )
    {   for(int i=0;i<v->size;i++)
        {       //printf("%s%12.4f",gsl_vector_get(v,i)>0.?"+":"",gsl_vector_get(v,i));
                printf("%12.4e",gsl_vector_get(v,i));
        }
        puts("");
    }
    else
        puts("sp_cluster_support_vector_view input 'v' (gsl_vector*) is a null pointer ... in SP_Support.c or SP_Support.h");

    return;
}


/* Matrix Viewer: Print out a N x N matrix on file */
void sp_cluster_support_matrix_view_f( FILE* fp, const gsl_matrix* m )
{
    if( m != NULL )
    {   for(int i=0;i<m->size1;i++)
        {   for(int j=0;j<m->size2;j++)
                fprintf(fp,"%s%.18lf\t",gsl_matrix_get(m,i,j)>0.?"+":"",gsl_matrix_get(m,i,j));
        fprintf(fp,"\n");
        }
    }
    else
        fputs("sp_cluster_support_matrix_view input 'm' (gsl_matrix*) is a null pointer ... in SP_Support.c or SP_Support.h",fp);

    fflush(fp);

    return;
}

/* Vector Viewer: Print out a vector on file */
void sp_cluster_support_vector_view_f( FILE* fp, const gsl_vector* v )
{
    if( v != NULL )
    {   for(int i=0;i<v->size;i++)
            fprintf(fp,"%s%.18f\t",gsl_vector_get(v,i)>0.?"+":"",gsl_vector_get(v,i));
        fprintf(fp,"\n");
    }
    else
        fputs("sp_cluster_support_vector_view input 'v' (gsl_vector*) is a null pointer ... in SP_Support.c or SP_Support.h",fp);

    fflush(fp);

    return;
}



/* Transformation Matrix calculator */
/* 
 * this function takes a vector, and the vector is transformed along transformed z-axis
 * the final return is the transformation matrix (rank 2 tensor), and the this data type is
 * gsl_matrix*
 */

gsl_matrix* sp_cluster_support_get_transformation_matrix( const gsl_vector* v )
{
    // gsl_vector_get(v,0) == 'x'
    // gsl_vector_get(v,1) == 'y'
    // gsl_vector_get(v,2) == 'z'
    gsl_matrix* pReturn = NULL;
    
    const double rxy = sqrt(pow(gsl_vector_get(v,0),2.)+pow(gsl_vector_get(v,1),2.));
    const double R   = sqrt(pow(gsl_vector_get(v,0),2.)+pow(gsl_vector_get(v,1),2.)+pow(gsl_vector_get(v,2),2.));
    double n1, n2, tmp; // dummy variables for workspace
    pReturn = gsl_matrix_calloc(4,4);   // Normally a rotation matrix is 3x3 
                                        // For the 1st column or row has an element of '1' at T_11.
 
    if( pReturn != NULL )
    {
        gsl_matrix_set(pReturn,0,0,1.);
        
        if( gsl_vector_get(v,0) == 0. && gsl_vector_get(v,1) == 0. && gsl_vector_get(v,2) > 0. ) // if vector 'v' is on z-axis
        {   gsl_matrix_set(pReturn,1,1,1.); gsl_matrix_set(pReturn,2,2,1.); gsl_matrix_set(pReturn,3,3,1.);
            // set the matrix as I
        }
        else if( gsl_vector_get(v,0) == 0. && gsl_vector_get(v,1) == 0. && gsl_vector_get(v,2) < 0. )   // if vector 'v' is on negative z-axis
        {   gsl_matrix_set(pReturn,1,1,1.); gsl_matrix_set(pReturn,2,2,1.); gsl_matrix_set(pReturn,3,3,-1.);
            // set the matrix has xy-plane reflection
        }
        else
        {   gsl_matrix_set(pReturn,3,1,gsl_vector_get(v,0)/R);
            gsl_matrix_set(pReturn,3,2,gsl_vector_get(v,1)/R);
            gsl_matrix_set(pReturn,3,3,gsl_vector_get(v,2)/R);  // set k' in the transformed (local) symmetry

            gsl_matrix_set(pReturn,2,1,gsl_vector_get(v,2)*gsl_vector_get(v,0)/rxy);
            gsl_matrix_set(pReturn,2,2,gsl_vector_get(v,2)*gsl_vector_get(v,1)/rxy);
            gsl_matrix_set(pReturn,2,3,-R*sqrt(1.-gsl_vector_get(v,2)*gsl_vector_get(v,2)/R/R));

            n1 = 1./sqrt(pow(gsl_matrix_get(pReturn,2,1),2.) 
               + pow(gsl_matrix_get(pReturn,2,2),2.)
               + pow(gsl_matrix_get(pReturn,2,3),2.));

            for(int i=0;i<3;i++)
            {   tmp = gsl_matrix_get(pReturn,2,i+1);
                tmp = tmp*n1;
                gsl_matrix_set(pReturn,2,i+1,tmp);  // set j' in the transformed (local) symmetry
            }

            gsl_matrix_set(pReturn,1,1,n1/R*(pow(gsl_vector_get(v,2),2.)*gsl_vector_get(v,1)/rxy
                        + R*gsl_vector_get(v,1)*sqrt(1.-pow(gsl_vector_get(v,2),2.)/R/R)));
            gsl_matrix_set(pReturn,1,2,n1/R*(-pow(gsl_vector_get(v,2),2.)*gsl_vector_get(v,0)/rxy
                        - R*gsl_vector_get(v,0)*sqrt(1.-pow(gsl_vector_get(v,2),2.)/R/R)));
            gsl_matrix_set(pReturn,1,3,0.);
    
            n2 = 1./sqrt(pow(gsl_matrix_get(pReturn,1,1),2.)
                    + pow(gsl_matrix_get(pReturn,1,2),2.)
                    + pow(gsl_matrix_get(pReturn,1,3),2.));

            for(int i=0;i<3;i++)
            {   tmp = gsl_matrix_get(pReturn,1,i+1);
                tmp = tmp/n2;
                gsl_matrix_set(pReturn,1,i+1,tmp);  // set i' in the transformed (local) symmetry 
            }
        }
    }   
    else
        puts("sp_cluster_get_trans_mat 'pReturn' alloc error ... in SP_support.c or SP_support.h");

    return pReturn;
}   /* SEE THE DETAILS IN THE 'N1' RING BINDER */



double sp_cluster_support_kronecker_delta( int a, int b )
{
    double Return = 0.;

    if( a == b )
        Return = 1.;

    return Return;
}

// norm of vector
double sp_cluster_support_get_norm( double v1, double v2, double v3 )
{   return pow(v1*v1+v2*v2+v3*v3,0.5);
}

// SPLINER

double** sp_cluster_support_get_spline( const double** data, const int knot_stride )
{
    double** pReturn = (double**)malloc((knot_stride-1)*sizeof(double*));
    for(int i=0;i<knot_stride-1;i++)
        pReturn[i] = (double*)calloc(4,sizeof(double));

    // Work Space
    const int n = knot_stride;
    double f[n];  double t[n];
    double a[n];  double b[n];  double c[n];  double d[n];
    double z[n];  double L[n];  double h[n];
    double alpha[n];  double mu[n];

    // Measuring slopes
    const double rslope = 0.;
    const double lslope = (data[3][1]-data[0][1])/(data[3][0]-data[0][0]);

    // Memset
    memset(f,0,n*sizeof(double));   memset(t,0,n*sizeof(double));
    memset(a,0,n*sizeof(double));   memset(b,0,n*sizeof(double));   memset(c,0,n*sizeof(double));   memset(d,0,n*sizeof(double));
    memset(z,0,n*sizeof(double));   memset(L,0,n*sizeof(double));   memset(h,0,n*sizeof(double));
    memset(alpha,0,n*sizeof(double));   memset(mu,0,n*sizeof(double));

    // Make Spline

    for(int i=0;i<n;i++)
    {   f[i] = data[i][1];
        t[i] = data[i][0];      }

    h[0] = t[1] - t[0];
    alpha[0] = 3.*(f[1]-f[0])/h[0]-3.*lslope;
    L[0] = 2.*h[0];
    mu[0] = 0.5;
    z[0] = alpha[0]/L[0];
    b[0] = lslope;
    // End Left Init

    for(int k=1;k<n-1;k++)
    {
        h[k] = t[k+1] - t[k];
        
        alpha[k] = (3./(h[k]*h[k-1]))*(f[k+1]*h[k-1]-f[k]*(h[k]+h[k-1])+f[k-1]*h[k]);
        
        L[k] = 2.*(h[k]+h[k-1])-h[k-1]*mu[k-1];

        mu[k] = h[k]/L[k];

        z[k] = (alpha[k]-h[k-1]*z[k-1])/L[k];
    }

    // Right End init
    alpha[n-1] = 3.*rslope - 3.*(f[n-1]-f[n-2])/h[n-2];
    L[n-1] = h[n-2]*(2.-mu[n-2]);
    z[n-1] = (alpha[n-1]-h[n-2]*z[n-2])/L[n-1];
    c[n-1] = z[n-1];

    for(int j=n-2;j>=0;j--)
    {
        c[j] = z[j]-mu[j]*c[j+1];

        b[j] = (f[j+1]-f[j])/h[j] - h[j]*(c[j+1]+2.*c[j])/3.;

        d[j] = (c[j+1]-c[j])/(3.*h[j]);

        a[j] = f[j];
    }
    // R-E END

    // SAVE DATA
    for(int i=0;i<knot_stride-1;i++)
    {   pReturn[i][0] = d[i];
        pReturn[i][1] = c[i] - 3.*t[i]*d[i];
        pReturn[i][2] = b[i] - 2.*c[i]*t[i] + 3.*t[i]*t[i]*d[i];
        pReturn[i][3] = a[i] - b[i]*t[i] + c[i]*t[i]*t[i] - d[i]*t[i]*t[i]*t[i];
    }

    return pReturn;
}

//int sp_cluster_support_get_atom_number(

void sp_cluster_support_print_xyz( void* sp_sys_void, const double cur_energy, const int rank, const int numtasks )
{
	int offset;
	// cast the type into "sp_cluster_system*"
	sp_cluster_system* sp_sys = (sp_cluster_system*)sp_sys_void;

	int atom_number=0;
	for(int i=0;i<sp_sys->number_of_classic_ion;i++)
	{	if( sp_sys->classic_ion[i].if_shell == SP_SUPPORT_FALSE )
			atom_number++;
		// CNT ONLY WHEN ITS CORE
	}

	if( cur_energy != 0. )
	{	
		//printf("\t%d\n",sp_sys->number_of_classic_ion+sp_sys->number_of_sp_ion);
		printf("\t%d\n",atom_number+sp_sys->number_of_sp_ion);
		printf(" SCF DONE %.6lf\n",cur_energy);
	}

	// CORE POSITION ONLY ... "xyz" Format Compatible 
	for(int n=0;n<sp_sys->number_of_sp_ion+sp_sys->number_of_classic_ion;n++)
	{
    	    offset = n-sp_sys->number_of_classic_ion;
	    if( n < sp_sys->number_of_classic_ion )
	    {
		if( sp_sys->classic_ion[n].if_shell == SP_SUPPORT_FALSE ) // if it is core
		{	fprintf(stdout,"%3s%12.6lf%12.6lf%12.6lf\n", sp_sys->classic_ion[n].atom_name,
			gsl_vector_get(sp_sys->classic_ion[n].core_position,0), gsl_vector_get(sp_sys->classic_ion[n].core_position,1), gsl_vector_get(sp_sys->classic_ion[n].core_position,2));
		}
	    }
	    else // this is for printing sp-ions
	    {
		    offset = n - sp_sys->number_of_classic_ion;
		    fprintf(stdout,"%3s%12.6lf%12.6lf%12.6lf\n", sp_sys->sp_ion[offset].atom_name,
			    gsl_vector_get(sp_sys->sp_ion[offset].core_position,0), gsl_vector_get(sp_sys->sp_ion[offset].core_position,1), gsl_vector_get(sp_sys->sp_ion[offset].core_position,2));
	    }
	}
	
	printf("\n");
	printf(" CONFIGURATION_XYZ_SC_INFO\n");
	printf(" %d\t%d\n",sp_sys->number_of_classic_ion,sp_sys->number_of_sp_ion);

	// SHOW SHELL CORE POSITION BOTH
	for(int n=0;n<sp_sys->number_of_sp_ion+sp_sys->number_of_classic_ion;n++)
	{
    	    offset = n-sp_sys->number_of_classic_ion;
	    if( n < sp_sys->number_of_classic_ion )
	    {
		if( sp_sys->classic_ion[n].if_shell == SP_SUPPORT_FALSE )
		{
			fprintf(stdout,"%3s%3s%12.6lf%12.6lf%12.6lf\n", sp_sys->classic_ion[n].atom_name,"c",
			gsl_vector_get(sp_sys->classic_ion[n].core_position,0), gsl_vector_get(sp_sys->classic_ion[n].core_position,1), gsl_vector_get(sp_sys->classic_ion[n].core_position,2));
		}
		else if( sp_sys->classic_ion[n].if_shell == SP_SUPPORT_TRUE )
		{
			fprintf(stdout,"%3s%3s%12.6lf%12.6lf%12.6lf\n", sp_sys->classic_ion[n].atom_name,"s",
			gsl_vector_get(sp_sys->classic_ion[n].core_position,0), gsl_vector_get(sp_sys->classic_ion[n].core_position,1), gsl_vector_get(sp_sys->classic_ion[n].core_position,2));
		}
	    }
	    else // this is for printing sp-ions
	    {
		    offset = n - sp_sys->number_of_classic_ion;
		    fprintf(stdout,"%3s%15.6lf%12.6lf%12.6lf\n", sp_sys->sp_ion[offset].atom_name,
			    gsl_vector_get(sp_sys->sp_ion[offset].core_position,0), gsl_vector_get(sp_sys->sp_ion[offset].core_position,1), gsl_vector_get(sp_sys->sp_ion[offset].core_position,2));
	    }
	}


	return;
}


/// DIIS_SUPPORT TOOLS

/*
int if_diis;       
int diis_max_depth;
int diis_cur_depth;

gsl_vector** diis_error_vector;                         // SAVE AS CIRCULAR QUEUE  
gsl_vector*  diis_coefficient_vector;
gsl_permutation* diis_ws_p;          
gsl_matrix*  diis_error_matrix;      
gsl_matrix*  diis_error_matrix_inv;  
*/

void sp_cluster_support_diis_error_vector_queue_init( void* sp_sys_void )
{
    sp_cluster_system* sp_sys = (sp_cluster_system*)sp_sys_void;
    sp_sys->diis_error_vector_queue_front = 0;
    sp_sys->diis_error_vector_queue_rear  = 0;
    return;
}

int sp_cluster_support_diis_error_vector_queue_isfull( void* sp_sys_void )
{   sp_cluster_system* sp_sys = (sp_cluster_system*)sp_sys_void;
    if( (sp_sys->diis_error_vector_queue_rear+1)%sp_sys->diis_max_depth == sp_sys->diis_error_vector_queue_front )
        return SP_SUPPORT_TRUE;
    else
        return SP_SUPPORT_FALSE;
}

int sp_cluster_support_diis_error_vector_queue_isempty( void* sp_sys_void )
{   sp_cluster_system* sp_sys = (sp_cluster_system*)sp_sys_void;
    if( sp_sys->diis_error_vector_queue_front == sp_sys->diis_error_vector_queue_rear )
        return SP_SUPPORT_TRUE;
    else
        return SP_SUPPORT_FALSE;
}

void sp_cluster_support_diis_error_vector_enqueue( void* sp_sys_void )
{
    sp_cluster_system* sp_sys = (sp_cluster_system*)sp_sys_void;
    double delta_evec[4];
    int low_idx;

    double delta_evec_sum = 0.;

    //if( !((sp_sys->diis_error_vector_queue_rear+1)%sp_sys->diis_max_depth == sp_sys->diis_error_vector_queue_front) )   // check if queue is full
    if( sp_cluster_support_diis_error_vector_queue_isfull( sp_sys ) == SP_SUPPORT_FALSE )
    {
        sp_sys->diis_error_vector_queue_rear = (sp_sys->diis_error_vector_queue_rear+1)%sp_sys->diis_max_depth;         // rear index ++
        sp_sys->diis_cur_depth++;

        for(int i=0;i<sp_sys->number_of_sp_ion;i++)
        {
            low_idx = sp_cluster_support_get_lowest_state(sp_sys->sp_ion[i].eigen_value);

            // current eigenvector - (loaded) previous eigenvector_gs
            delta_evec[0] = gsl_matrix_get(sp_sys->sp_ion[i].eigen_vector,0,low_idx) - gsl_vector_get(sp_sys->sp_ion[i].eigen_vector_gs,0);
            delta_evec[1] = gsl_matrix_get(sp_sys->sp_ion[i].eigen_vector,1,low_idx) - gsl_vector_get(sp_sys->sp_ion[i].eigen_vector_gs,1);
            delta_evec[2] = gsl_matrix_get(sp_sys->sp_ion[i].eigen_vector,2,low_idx) - gsl_vector_get(sp_sys->sp_ion[i].eigen_vector_gs,2);
            delta_evec[3] = gsl_matrix_get(sp_sys->sp_ion[i].eigen_vector,3,low_idx) - gsl_vector_get(sp_sys->sp_ion[i].eigen_vector_gs,3);

            delta_evec_sum += delta_evec[0]*delta_evec[0] + delta_evec[1]*delta_evec[1] + delta_evec[2]*delta_evec[2] + delta_evec[3]*delta_evec[3];

            gsl_vector_set(sp_sys->diis_error_vector[ sp_sys->diis_error_vector_queue_rear ],i*4+0,delta_evec[0]);
            gsl_vector_set(sp_sys->diis_error_vector[ sp_sys->diis_error_vector_queue_rear ],i*4+1,delta_evec[0]);
            gsl_vector_set(sp_sys->diis_error_vector[ sp_sys->diis_error_vector_queue_rear ],i*4+2,delta_evec[0]);
            gsl_vector_set(sp_sys->diis_error_vector[ sp_sys->diis_error_vector_queue_rear ],i*4+3,delta_evec[0]);
        

            // loading .. previous eigen_vector_gs (note that the function 'sp_cluster_support_load_gs_eigenvector()' must be called before enqueue)
            gsl_vector_set(sp_sys->diis_prev_eigen_vector[ sp_sys->diis_error_vector_queue_rear ],i*4+0, gsl_vector_get(sp_sys->sp_ion[i].eigen_vector_gs,0));
            gsl_vector_set(sp_sys->diis_prev_eigen_vector[ sp_sys->diis_error_vector_queue_rear ],i*4+1, gsl_vector_get(sp_sys->sp_ion[i].eigen_vector_gs,1));
            gsl_vector_set(sp_sys->diis_prev_eigen_vector[ sp_sys->diis_error_vector_queue_rear ],i*4+2, gsl_vector_get(sp_sys->sp_ion[i].eigen_vector_gs,2));
            gsl_vector_set(sp_sys->diis_prev_eigen_vector[ sp_sys->diis_error_vector_queue_rear ],i*4+3, gsl_vector_get(sp_sys->sp_ion[i].eigen_vector_gs,3));

            #ifdef DEBUG_SUPPORT
            printf("\ndelta evec\n");
            printf("%12.6lf%12.6lf%12.6lf%12.6lf\n",delta_evec[0],delta_evec[1],delta_evec[2],delta_evec[3]);
            printf("%12.6lf%12.6lf%12.6lf%12.6lf\n",
            gsl_matrix_get(sp_sys->sp_ion[i].eigen_vector,0,low_idx),
            gsl_matrix_get(sp_sys->sp_ion[i].eigen_vector,1,low_idx),
            gsl_matrix_get(sp_sys->sp_ion[i].eigen_vector,2,low_idx),
            gsl_matrix_get(sp_sys->sp_ion[i].eigen_vector,3,low_idx));
            #endif
/*
            printf("%12.6lf%12.6lf%12.6lf%12.6lf\n",gsl_vector_get(sp_sys->diis_prev_eigen_vector[ sp_sys->diis_error_vector_queue_rear ], i*4 + 0 ),
            gsl_vector_get(sp_sys->diis_prev_eigen_vector[ sp_sys->diis_error_vector_queue_rear ], i*4 + 1 ),
            gsl_vector_get(sp_sys->diis_prev_eigen_vector[ sp_sys->diis_error_vector_queue_rear ], i*4 + 2 ),
            gsl_vector_get(sp_sys->diis_prev_eigen_vector[ sp_sys->diis_error_vector_queue_rear ], i*4 + 3 ));
*/
        /* Note that the both,
        
            Error vector ( sp_sys->diis_error_vector )
            Prev  evecs  ( sp_sys->diis_prev_eigen_vector )

            Follows the form of (gsl_vector)v[i][j]

            v[i]    ... 'i' -> index in the queue // queue is circular type, max dept is 'sp_sys->diis_max_depth'
            v[i][j] ... 'j' -> DATA, s_sp1, px_sp1, py_sp1, pz_sp1 // s_sp2, px_sp2, py_sp2, pz_sp2 // ....
                                i.e., length of number_of_sp_ions * 4 (each lone pair ground-state MO saved with stride of 4)
        */
        }
        printf("delta_evec_sum: %12.8lf\n",sqrt(delta_evec_sum)/sp_sys->number_of_sp_ion);

    }
    return;     // if queue is full do nothing
}     

void sp_cluster_support_diis_error_vector_dequeue( void* sp_sys_void )
{   sp_cluster_system* sp_sys = (sp_cluster_system*)sp_sys_void;
    if( sp_cluster_support_diis_error_vector_queue_isempty( sp_sys ) == SP_SUPPORT_FALSE )
    {   sp_sys->diis_error_vector_queue_front = (sp_sys->diis_error_vector_queue_front+1)%sp_sys->diis_max_depth;
        sp_sys->diis_cur_depth--;
    }
    #ifdef DEBUG_SUPPORT
    const int front = sp_sys->diis_error_vector_queue_front;
    const int rear  = sp_sys->diis_error_vector_queue_rear;
    printf("%d\t%d\n",front,rear);
    #endif

    return;     // if queue is empty do nothing
}

double sp_cluster_support_diis_get_error_vector_gnorm( void* sp_sys_void )      // returns latest .. i.e., the one at rear
{   sp_cluster_system* sp_sys = (sp_cluster_system*)sp_sys_void;
    double ret;
    if( sp_cluster_support_diis_error_vector_queue_isempty( sp_sys ) == SP_SUPPORT_FALSE )
    {   ret = gsl_blas_dnrm2( sp_sys->diis_error_vector[ sp_sys->diis_error_vector_queue_rear ] );
        return ret/sp_sys->number_of_sp_ion;
    }
    return 0.;
}

void sp_cluster_support_diis_solve_least_square_problem( void* sp_sys_void )
{   sp_cluster_system* sp_sys = (sp_cluster_system*)sp_sys_void;

    // int count = front > rear ? (MAX - front + rear) : (rear - front);
    const int queue_capacity = sp_sys->diis_max_depth;
    const int front = sp_sys->diis_error_vector_queue_front;
    const int rear  = sp_sys->diis_error_vector_queue_rear;
    const int cur_size = front > rear ? (queue_capacity - front + rear) : (rear - front);       // number of elements in the queue (error_vector)
    
/*
gsl_vector*  diis_coefficient_vector;
gsl_permutation* diis_ws_p;          
gsl_matrix*  diis_error_matrix;      
gsl_matrix*  diis_error_matrix_inv;  
*/
    size_t len = cur_size + 1;

printf("\n#### IN LSP SOLVER / len : %zu\n",len);
printf("ws size : ErrorMat / %zu %zu\n",sp_sys->diis_error_matrix->size1,sp_sys->diis_error_matrix->size2);
printf("cursize/max_depth/cur_depth : %4d%4d%4d\n", cur_size, sp_sys->diis_max_depth - 1, sp_sys->diis_cur_depth );
    //if( !(cur_size == (sp_sys->diis_max_depth - 1)) )
    //if( !(cur_size == (sp_sys->diis_max_depth)) )
    //if( !(len == (sp_sys->diis_max_depth)) )
    if( !(len == sp_sys->diis_error_matrix->size1) )
    {  
printf("Resizing!!\n");
        // if current size is not same with actual max_depth in use ... i.e., max_depth - 1, since circular queue in use
        // have to resize the workspaces ... diis_ws_p (gsl_permutation*) / diis_error_matrix (gsl_matrix*) / diis_error_matrix_inv (gsl_matrix*) 

        gsl_vector_free(sp_sys->diis_coefficient_vector);
        gsl_vector_free(sp_sys->diis_least_square_condition);
        gsl_permutation_free(sp_sys->diis_ws_p);
        gsl_matrix_free(sp_sys->diis_error_matrix);
        gsl_matrix_free(sp_sys->diis_error_matrix_inv);

        sp_sys->diis_coefficient_vector = gsl_vector_calloc(len);
        sp_sys->diis_least_square_condition = gsl_vector_calloc(len);
        sp_sys->diis_ws_p = gsl_permutation_calloc(len);
        sp_sys->diis_error_matrix = gsl_matrix_calloc(len,len);
        sp_sys->diis_error_matrix_inv = gsl_matrix_calloc(len,len);      // len = cur_size (error_vector length) + 1 // 1 is to represent least square sense of coeff
    }

    int ii=0; int jj=0;     // tags to get actual matrix indices !
    double elem_tmp;
    for(int i=(front+1)%queue_capacity; i!=(rear+1)%queue_capacity ; i=(i+1)%queue_capacity)    // for loop circulate over queue elements
    {
        for(int j=(front+1)%queue_capacity; j!=(rear+1)%queue_capacity ; j=(j+1)%queue_capacity)
        {   
            #ifdef DEBUG_SUPPORT
            printf("i,j / tag_i,tag_j : %d\t%d / %d\t%d\n",i,j,ii,jj);
            //printf("size1 / size2 : %d\t%d\n", sp_sys->diis_error_matrix->size1, sp_sys->diis_error_matrix->size2);
            #endif

            gsl_blas_ddot(sp_sys->diis_error_vector[j],sp_sys->diis_error_vector[i],&elem_tmp);
            gsl_matrix_set(sp_sys->diis_error_matrix,ii,jj,elem_tmp);
            jj++;
        }
        jj=0;
        ii++;
    }
    printf("ERROR - 0 ?\n");
    ///// SET RESTS 
     for(int i=0;i<len;i++)
    {   gsl_matrix_set(sp_sys->diis_error_matrix,i,len-1,-1.);
        gsl_matrix_set(sp_sys->diis_error_matrix,len-1,i,-1.);        }
    printf("ERROR?\n");
    gsl_matrix_set(sp_sys->diis_error_matrix,len-1,len-1,0.);
    printf("ERROR?\n");

    ///// SET least Square Condition vector
    gsl_vector_set(sp_sys->diis_least_square_condition,len-1,-1.);

    //// calculate inverse
    int signum; // variable for LU_decomp

    #ifdef DEBUG_SUPPORT
    printf("ErrorMatrix\n");
    sp_cluster_support_matrix_view( sp_sys->diis_error_matrix );
    #endif

    /*  Say,
        Error matrix                    E
        Coefficient Vector              C
        RHS Least Square Condition      L

        Need to find 'C' vector

        Solve : EC = L, i.e., C = E_Inverse  L
    */

    gsl_linalg_LU_decomp(sp_sys->diis_error_matrix,sp_sys->diis_ws_p,&signum);
    gsl_linalg_LU_invert(sp_sys->diis_error_matrix,sp_sys->diis_ws_p,sp_sys->diis_error_matrix_inv);    // invert is saved in 'sp_sys->diis_error_matrix_inv'
    // Inverse Found

    // int gsl_blas_dgemv(CBLAS_TRANSPOSE_t TransA, double alpha, const gsl_matrix *A, const gsl_vector *x, double beta, gsl_vector *y)
    gsl_blas_dgemv(CblasNoTrans,1.,sp_sys->diis_error_matrix_inv,sp_sys->diis_least_square_condition,0.,sp_sys->diis_coefficient_vector);
    // Solve LeaseSquare Problem Answer saved in 'sp_sys->diis_coefficient_vector'

    #ifdef DEBUG_SUPPORT
    //sp_cluster_support_vector_view( sp_sys->diis_least_square_condition );
    printf("coefficient vector, ignore the last dummy lambda value\n");
    sp_cluster_support_vector_view( sp_sys->diis_coefficient_vector );
    //sp_cluster_support_matrix_view( sp_sys->diis_error_matrix );
    printf("InverseMatrix\n");
    sp_cluster_support_matrix_view( sp_sys->diis_error_matrix_inv );
    printf("cursize/max_depth/cur_depth : %d\t%d\t%d\n", cur_size, sp_sys->diis_max_depth - 1, sp_sys->diis_cur_depth );

    

    printf(" ---- FINALISE LEAST SQUARE SOLVER \n\n");
    #endif
    return;
}


void sp_cluster_support_diis_least_square_result_update( void* sp_sys_void )
{   sp_cluster_system* sp_sys = (sp_cluster_system*)sp_sys_void;
    int low_idx;
    int jj = 0;
    double cs, cx, cy, cz, norm_factor;

    const int queue_capacity = sp_sys->diis_max_depth;
    const int front = sp_sys->diis_error_vector_queue_front;
    const int rear  = sp_sys->diis_error_vector_queue_rear;

    for(int i=0;i<sp_sys->number_of_sp_ion;i++)
    {
        low_idx = sp_cluster_support_get_lowest_state(sp_sys->sp_ion[i].eigen_value);
 
        // results are in 'sp_sys->diis_coefficient_vector' ... data structure sp_sys->diis_coefficient_vector->size or 'stride'?
        // size can also be obtained by 'sp_sys->diis_cur_depth'

        jj = 0;
        cs = 0.; cx = 0.; cy = 0.; cz = 0.;

        for(int j=(front+1)%queue_capacity; j!=(rear+1)%queue_capacity ; j=(j+1)%queue_capacity)
        {   
            cs += gsl_vector_get( sp_sys->diis_coefficient_vector,jj) * gsl_vector_get(sp_sys->diis_prev_eigen_vector[j],i*4+0);
            cx += gsl_vector_get( sp_sys->diis_coefficient_vector,jj) * gsl_vector_get(sp_sys->diis_prev_eigen_vector[j],i*4+1);
            cy += gsl_vector_get( sp_sys->diis_coefficient_vector,jj) * gsl_vector_get(sp_sys->diis_prev_eigen_vector[j],i*4+2);
            cz += gsl_vector_get( sp_sys->diis_coefficient_vector,jj) * gsl_vector_get(sp_sys->diis_prev_eigen_vector[j],i*4+3);
            jj++;

            printf("index j / tag jj: %d\t%d\n",j,jj);
        }
        
#ifdef DEBUG_SUPPORT
printf("%dth\t%12.4lf%12.4lf%12.4lf%12.4lf%20.4lf\n",i+1,cs,cx,cy,cz,sqrt(cs*cs+cx*cx+cy*cy+cz*cz));
#endif     
        norm_factor = sqrt(cs*cs+cx*cx+cy*cy+cz*cz);  
        norm_factor = 1.;

        gsl_matrix_set(sp_sys->sp_ion[i].eigen_vector,0,low_idx,cs/norm_factor);
        gsl_matrix_set(sp_sys->sp_ion[i].eigen_vector,1,low_idx,cx/norm_factor);
        gsl_matrix_set(sp_sys->sp_ion[i].eigen_vector,2,low_idx,cy/norm_factor);
        gsl_matrix_set(sp_sys->sp_ion[i].eigen_vector,3,low_idx,cz/norm_factor);
        // set eigenvector ... updated!!!

        // finish up the rest ...
    }


    return;
}





















