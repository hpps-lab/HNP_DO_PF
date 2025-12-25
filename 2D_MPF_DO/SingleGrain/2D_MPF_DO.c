#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//=================================================-=================================================//
//*****************************************Sub Function*****************************************//
//=================================================-=================================================//

#define max2(a, b) ((a) > (b) ? (a) : (b))

//=======================================================================

void swap(double **f, double **fn)
{
    double *tmp = *f;
    *f = *fn;
    *fn = tmp;
}

//=======================================================================

void output_vtkdata
(
	char *filename,
    const int number_of_phase,
	const int nx, 
	const int ny,
	const double dx,
	const double dy,
	const double *pf
)
{
	char fname[256];
	
	snprintf(fname,sizeof(fname),"%s.vtk",filename);
	FILE *fp = fopen(fname,"w");
	
	fprintf(fp, "# vtk DataFile Version 3.0\n");
	fprintf(fp, "output\n");
	fprintf(fp, "ASCII\n");
	fprintf(fp, "DATASET STRUCTURED_GRID\n");
	fprintf(fp, "DIMENSIONS %6d %6d %6d\n",nx,ny,1);
	fprintf(fp, "POINTS %10d float\n",nx*ny*1);
	
	for (int j=0; j<ny; j++) {
	for (int i=0; i<nx; i++) {
	    fprintf(fp, "%12.4e %12.4e %3.1f\n",(double)i*dx, (double)j*dy, 0.);
	}
	}
   
	fprintf(fp, "POINT_DATA %10d\n",nx*ny*1);
   
	for (int m=0; m<number_of_phase; m++) {
	    fprintf(fp, "SCALARS Phase-field%02d float\n", m);
	    fprintf(fp, "LOOKUP_TABLE default\n");
	
	    for (int j=0; j<ny; j++) {
	    for (int i=0; i<nx; i++) {
	        const int ix = j*nx + i ;
	        fprintf(fp, "%6.3f\n",pf[m*nx*ny + ix]);
	    }
	    }
	}
	
	fprintf(fp, "SCALARS phase-field_square_sum float\n");
	fprintf(fp, "LOOKUP_TABLE default\n");
	
	for (int j=0; j<ny; j++) {
	for (int i=0; i<nx; i++) {
	    const int ix = j*nx + i ;

		double pf_sq_sum = 0.0;
		for (int m=0; m<number_of_phase; m++) pf_sq_sum += (0.5*(pf[m*nx*ny + ix] + 1.0))*(0.5*(pf[m*nx*ny + ix] + 1.0));

	    fprintf(fp, "%6.3f\n",pf_sq_sum);
	}
	}
	fclose(fp);
}

//=======================================================================

double x_axis_interface_position_exploration
(
	const int nx, 
	const int ny, 
	const int m,
	const int j,
	const double dx,
	const double *pf
)
{
    double interface_pos_x = 0.0;
    for (int i=0; i<nx-1; i++) {

        const int ix = j*nx + i;

        if(0.0 <= pf[m*nx*ny + ix] && pf[m*nx*ny + ix+1] < 0.0){

            interface_pos_x = i*dx + (pf[m*nx*ny + ix])/(pf[m*nx*ny + ix]-pf[m*nx*ny + ix+1])*dx;
            break;
        }
    }
    return interface_pos_x;
}

//=======================================================================

void Multi_Phase_field
(
    const int number_of_phase,
	const int  nx,
	const int  ny,
	const double pi,
	const double rdx,
	const double rdy,
	const double dt,
	const double pf_max,
	const double pf_min,
	const double *aaa_DO,
	const double *www_DO,
	const double *pmb_DO,
	const double *df_ij,
	const double *pf,
	      double *pff
)
{
	for(int j=0; j<ny; j++){
	for(int i=0; i<nx; i++){	

		const int ix = j*nx + i;

		double pfl[number_of_phase];
		for(int m0=0; m0<number_of_phase; m0++){

			const double pfi_ix = pf[m0*nx*ny+ix];
			const double pfi_im = (i!=0   ) ? pf[m0*nx*ny+ix-1 ] : pf[m0*nx*ny+ix+1 ];
			const double pfi_ip = (i!=nx-1) ? pf[m0*nx*ny+ix+1 ] : pf[m0*nx*ny+ix-1 ];
			const double pfi_jm = (j!=0   ) ? pf[m0*nx*ny+ix-nx] : pf[m0*nx*ny+ix+nx];
			const double pfi_jp = (j!=ny-1) ? pf[m0*nx*ny+ix+nx] : pf[m0*nx*ny+ix-nx];

			pfl[m0] = (pfi_ip + pfi_im - 2.0*pfi_ix)*rdx*rdx + (pfi_jp + pfi_jm - 2.0*pfi_ix)*rdy*rdy;
		}

		int negativephase_flag[number_of_phase];
        for(int m0=0; m0<number_of_phase; m0++) negativephase_flag[m0] = 0;

		for(int m0=0; m0<number_of_phase; m0++){ //itration for removing negativephase 

			int number_of_positivephase = 0;
        	for(int m0=0; m0<number_of_phase; m0++){
				if(negativephase_flag[m0]==0) number_of_positivephase++;
			}

			for(int m0=0; m0<number_of_phase; m0++){

				double       pffi_ix = pf_min;
				const double pfi_ix = pf[m0*nx*ny+ix];

				if(number_of_positivephase!=1 && negativephase_flag[m0]==0){

					double dp = 0.0;
					for(int m1=0; m1<number_of_phase; m1++){

						if(m0!=m1 && negativephase_flag[m1]==0){

							const double pfj_ix = pf[m1*nx*ny+ix];

							const double M_ij = pmb_DO[number_of_phase*m1 + m0];

							double       dp1 = 0.0;
							double       dp2 = 0.0;
							const double dp3 = -8.0/pi*sqrt((pfi_ix+1.0)*(pfj_ix+1.0))*df_ij[ number_of_phase*m1 + m0];

							for(int m2=0; m2<number_of_phase; m2++){

								if(negativephase_flag[m2]==0){

									const double pfk_ix = pf[m2*nx*ny+ix];

									const double w_ik = www_DO[number_of_phase*m2 + m0];
									const double w_jk = www_DO[number_of_phase*m2 + m1];
									const double a_ik = aaa_DO[number_of_phase*m2 + m0];
									const double a_jk = aaa_DO[number_of_phase*m2 + m1];

									dp1 += (w_ik-w_jk)*(pfk_ix+1.0);
									dp2 += 0.5*(a_ik*a_ik-a_jk*a_jk)*pfl[m2];
								}
							}//end m2 loop

							dp += M_ij*(dp1 + dp2 + dp3);
						}
					}//end m1 loop
					
					pffi_ix = pfi_ix -2.0/number_of_positivephase*dp*dt;

				}else{

					pffi_ix = (negativephase_flag[m0]==0) ? pf_max : pf_min;
				}
				
				pff[m0*nx*ny+ix] = pffi_ix;

			}//end m0 loop

			int number_of_positivephase_0 = 0;
			for(int m0=0; m0<number_of_phase; m0++){

				if(negativephase_flag[m0]==0 && pff[m0*nx*ny+ix] < pf[m0*nx*ny+ix] && pff[m0*nx*ny+ix] <= pf_min){

					negativephase_flag[m0] = -1;
				}

				if(negativephase_flag[m0]==0) number_of_positivephase_0++;
			}

			if(number_of_positivephase==number_of_positivephase_0){

				double psum = 0.0;
				number_of_positivephase_0 = 0;
				for(int m0=0; m0<number_of_phase; m0++){

					if(pff[m0*nx*ny+ix]>=pf_max){

						pff[m0*nx*ny+ix] = pf_max;

						for(int m1=0; m1<number_of_phase; m1++){

							if(m0!=m1) pff[m1*nx*ny+ix] = pf_min;
						}

						number_of_positivephase_0 = 1;
						break;
					}else if(pff[m0*nx*ny+ix]<=pf_min){

						pff[m0*nx*ny+ix] = pf_min;

					}else{

						number_of_positivephase_0++;
						psum += (pff[m0*nx*ny+ix]+1.0)*0.5;
					}
				
				}

				if(number_of_positivephase_0!=1){

					for(int m0=0; m0<number_of_phase; m0++){

						if(negativephase_flag[m0]==0){

							pff[m0*nx*ny+ix] = (pff[m0*nx*ny+ix]+1.0)/psum - 1.0;
						}
					}
				}
				
				break; 
			}
		} //for loop end (itration for removing negativephase)
	}
    }
}

//=================================================-=================================================//
//********************************************Main Function******************************************//
//=================================================-=================================================//

int main(int argc, char** argv)
{
    double *PF;
    double *PPF;
    double *GAM_ij, *MOB_ij, *DF_ij;
    double *AAA_DO, *WWW_DO, *PMB_DO;

   const int nend = 3000;
    const int nout =   150;
    const int nvel =   150;

	const int ncx = 256;
	const int ncy = 256;
	const int ngx = ncx+1;
	const int ngy = ncy+1;

//<<parameter setting>>//
	const double pi  = M_PI;
	const double HHH = 200.0e-6;    //Domain size x [m]
	const double LLL = HHH/ncy*ncx; //Domain size y [m]
	const double dx  = LLL/ncx;     //Grid spacing [m]
	const double dy  = HHH/ncy;     //Grid spacing [m]

    const double delta    = 4.0*dx;  //Interface width of DO model [m]
	const double mobility = 1.0e-12; //Reference boundary mobility[m^4/(Js)]
	const double gamma    = 1.0;     //Reference boundary energy [J/m^2]
	const double df       = 5.0e+5;  //Reference driving force [J/m^3]

    const double pf_max   = +1.0;                        //Max value of phi
    const double pf_min   = -1.0;                        //Min value of phi

    const int number_of_phase = 2;
    GAM_ij = (double *) malloc(number_of_phase*number_of_phase*sizeof(double)); if(GAM_ij == NULL){fprintf(stderr,"Malloc error\n"); exit(1);}
    MOB_ij = (double *) malloc(number_of_phase*number_of_phase*sizeof(double)); if(MOB_ij == NULL){fprintf(stderr,"Malloc error\n"); exit(1);}
    DF_ij  = (double *) malloc(number_of_phase*number_of_phase*sizeof(double)); if(DF_ij  == NULL){fprintf(stderr,"Malloc error\n"); exit(1);}

 	for(int m0=0; m0<number_of_phase; m0++){
	for(int m1=0; m1<number_of_phase; m1++){

        const int mij = number_of_phase*m1 + m0;

        if(m0!=m1){

            MOB_ij[mij] = mobility;
            GAM_ij[mij] = gamma;
            DF_ij[ mij] = (m0==0) ? df : ((m1==0) ? -df : 0.0);

        }else{

            MOB_ij[mij] = 0.0;
            GAM_ij[mij] = 0.0;
            DF_ij[ mij] = 0.0;

        }
    }
    }

    AAA_DO = (double *) malloc(number_of_phase*number_of_phase*sizeof(double)); if(AAA_DO == NULL){fprintf(stderr,"Malloc error\n"); exit(1);}
    WWW_DO = (double *) malloc(number_of_phase*number_of_phase*sizeof(double)); if(WWW_DO == NULL){fprintf(stderr,"Malloc error\n"); exit(1);}
    PMB_DO = (double *) malloc(number_of_phase*number_of_phase*sizeof(double)); if(PMB_DO == NULL){fprintf(stderr,"Malloc error\n"); exit(1);}

 	for(int m0=0; m0<number_of_phase; m0++){
	for(int m1=0; m1<number_of_phase; m1++){

        const int mij = number_of_phase*m1 + m0;

        if(m0!=m1){

			AAA_DO[mij] = 2.0*sqrt(2.0*delta*GAM_ij[mij])/pi;
			WWW_DO[mij] = 4.0*GAM_ij[mij]/delta;
			PMB_DO[mij] = pi*pi*MOB_ij[mij]/(8.0*delta);

		}else{

			AAA_DO[mij] = 0.0;
			WWW_DO[mij] = 0.0;
			PMB_DO[mij] = 0.0;
		}
	}
	}

    double mob_max = 0.;
    double gam_max = 0.;
    for(int mij=0; mij<number_of_phase*number_of_phase; mij++){
        if(mob_max < MOB_ij[mij]) mob_max = MOB_ij[mij];
        if(gam_max < GAM_ij[mij]) gam_max = GAM_ij[mij];
    }

    const double dt = dx*dx/(6.0*mob_max*gam_max);

	PF  = (double *) malloc(number_of_phase*ngx*ngy*sizeof(double)); if(PF  == NULL){fprintf(stderr,"Malloc error PF\n"); exit(1);}
	PPF = (double *) malloc(number_of_phase*ngx*ngy*sizeof(double)); if(PPF == NULL){fprintf(stderr,"Malloc error PPF\n"); exit(1);}

//<<Initial profile Setting>>//
	for(int j=0; j<ngy; j++){
	for(int i=0; i<ngx; i++){
	
		const int ix = j*ngx + i;
		
		double xx = ((double)i - ncx/2.0)*dx;
		double yy = ((double)j - ncy/2.0)*dy;
		double rr = sqrt(xx*xx + yy*yy) - ncx/16.0*dx;
		
		if(rr < -0.5*delta){
			PF[0*ngx*ngy + ix] = pf_max;
		}else if(0.5*delta < rr){
			PF[0*ngx*ngy + ix] = pf_min;
		}else{
			PF[0*ngx*ngy + ix] = -sin(pi*rr/delta);
		}
	
		if(0.5*delta < rr){
			PF[1*ngx*ngy + ix] = pf_max;
		}else if(rr < -0.5*delta){
			PF[1*ngx*ngy + ix] = pf_min;
		}else{
			PF[1*ngx*ngy + ix] = sin(pi*rr/delta);
		}
	
		const double pf0 = 1.0 - (PF[0*ngx*ngy + ix]+1.0)*0.5;
		for(int m=1; m<number_of_phase; m++){

			PF[m*ngx*ngy + ix] = (PF[m*ngx*ngy + ix]+1.0)*pf0-1.0;
		}
	}
	}

	//Initial data output
    if(1){
        FILE *fini = fopen("computational_condition.dat", "w");
        if(fini == NULL){
        	fprintf(stderr,"can't open computational_condition.dat\n");
        	exit(1);
        } 

        fprintf(fini, "<<<<< Setting Data >>>>> \n");
        fprintf(fini, "ngx                 [-] =  %10d\n", ngx);
        fprintf(fini, "ngy                 [-] =  %10d\n", ngy);

        fprintf(fini, "<<<<< Dimensional values >>>>> \n");
        fprintf(fini, "Domain size x (L)            [m]        =  %16.7e\n", ncx*dx);
        fprintf(fini, "Domain size y (H)            [m]        =  %16.7e\n", ncy*dy);
        fprintf(fini, "Interface thickness          [m]        =  %16.7e\n", delta);
        fprintf(fini, "Mobirity M_ij(i!=0, j!=0)    [m^4/(Js)] =  %16.7e\n", mobility);
        fprintf(fini, "Bound. Ene. g_0j, g_i0       [J/m^2]    =  %16.7e\n", gamma);
        fprintf(fini, "Lattice size                 [m]        =  %16.7e\n", dx);
        fprintf(fini, "Time increment               [s]        =  %16.7e\n", dt);

        fclose(fini);
    }

    //<<interface velocity along x-axies initialaize>>//
    double interface_pos_x1_at_y00h = x_axis_interface_position_exploration(ngx, ngy, 0,    0, dx, PF);
    double interface_pos_x0_at_y00h = interface_pos_x1_at_y00h;
    double interface_pos_x1_at_y05h = x_axis_interface_position_exploration(ngx, ngy, 0, ncy/2, dx, PF);
    double interface_pos_x0_at_y05h = interface_pos_x1_at_y05h;

//<<<<<<<<<<<<<<<<<<<<<<<<< main loop start >>>>>>>>>>>>>>>>>>>>>>>>>//

	for(int nstep=1; nstep<=nend ; nstep++){

        if(nstep%100==0) printf("nstep = %8d\n", nstep);

        //Finite difference computation of multi-phase-field equation
        Multi_Phase_field(number_of_phase,ngx,ngy, pi, 1.0/dx, 1.0/dy, dt, pf_max,pf_min,
            AAA_DO,WWW_DO,PMB_DO, DF_ij, PF, PPF);
		swap(&PF,  &PPF);

        //<<interface velocity along x-axies>>//
		if(nstep%nvel==0){
            if(1){
                interface_pos_x0_at_y00h = interface_pos_x1_at_y00h;
                interface_pos_x1_at_y00h = x_axis_interface_position_exploration(ngx, ngy, 0,    0, dx, PF);
                double interface_vel_x1_at_y00h = (interface_pos_x1_at_y00h-interface_pos_x0_at_y00h)/(nvel*dt);
                
		        FILE *fp_v = fopen("interface_velocity_x_at_y00h.dat","a");
		        fprintf(fp_v, "%6d %16.7e %16.7e %16.7e\n", nstep, (double)nstep*dt, interface_pos_x1_at_y00h, interface_vel_x1_at_y00h);
		        fclose(fp_v);
            }

            if(1){
                interface_pos_x0_at_y05h = interface_pos_x1_at_y05h;
                interface_pos_x1_at_y05h = x_axis_interface_position_exploration(ngx, ngy, 0, ncy/2, dx, PF);
                double interface_vel_x1_at_y05h = (interface_pos_x1_at_y05h-interface_pos_x0_at_y05h)/(nvel*dt);

		        FILE *fp_v = fopen("interface_velocity_x_at_y05h.dat","a");
		        fprintf(fp_v, "%6d %16.7e %16.7e %16.7e\n", nstep, (double)nstep*dt, interface_pos_x1_at_y05h, interface_vel_x1_at_y05h);
		        fclose(fp_v);
            }
        }

        if(nstep%nout==0 || nstep==1){
            char fvti[128];
            snprintf(fvti,sizeof(fvti),"paraview_%08d",nstep);	
            output_vtkdata(fvti,number_of_phase, ngx, ngy, dx, dy, PF);
        }

        
		if(max2(interface_pos_x1_at_y00h,interface_pos_x1_at_y05h) > 0.95*LLL){
            printf("grain 0 has reached the right end\n");
            break;
        }
    }
    
	return 0;
}

//=================================================-=================================================//
//************************************************END************************************************//
//=================================================-=================================================//
