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
	const double *pf,
	const double *ps
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
	fprintf(fp, "POINTS %10d double\n",nx*ny*1);
	
	for (int j=0; j<ny; j++) {
	for (int i=0; i<nx; i++) {
	    fprintf(fp, "%12.4e %12.4e %3.1f\n",(double)i*dx, (double)j*dy, 0.);
	}
	}
   
	fprintf(fp, "POINT_DATA %10d\n",nx*ny*1);
   
	for (int m=0; m<number_of_phase; m++) {
	    fprintf(fp, "SCALARS Phase-field%02d double\n", m);
	    fprintf(fp, "LOOKUP_TABLE default\n");
	
	    for (int j=0; j<ny; j++) {
	    for (int i=0; i<nx; i++) {
	        const int ix = j*nx + i ;
	        fprintf(fp, "%6.3f\n",pf[m*nx*ny + ix]);
	    }
	    }
	}
	
	for (int m=0; m<number_of_phase; m++) {
	    fprintf(fp, "SCALARS Nonlinear-phase-field%02d double\n", m);
	    fprintf(fp, "LOOKUP_TABLE default\n");
	
	    for (int j=0; j<ny; j++) {
	    for (int i=0; i<nx; i++) {
	        const int ix = j*nx + i ;
	        fprintf(fp, "%7.3f\n",ps[m*nx*ny + ix]);
	    }
	    }
	}
	
	fprintf(fp, "SCALARS phase-field_square_sum double\n");
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

int active_phase_check
(
    const int ix,
    const int nx,
    const int ny,
	const int number_of_phase, 
    const double pi,
	const double theta_th, 
	const double *ps,
          int *activephase_flag
)
{
	int number_of_activephase = 0;
	for(int m0=0; m0<number_of_phase; m0++){

		if(-0.5*pi+theta_th < ps[m0*nx*ny+ix]){

			number_of_activephase++;

			if(ps[m0*nx*ny+ix] < 0.5*pi-theta_th){
                
				activephase_flag[m0] = 0;

			}else{

				number_of_activephase = 1;
				activephase_flag[m0] = +1;
				for(int m1=0; m1<number_of_phase; m1++) if(m0!=m1) activephase_flag[m1] = -1;

				break;
			}
		}else{
			activephase_flag[m0] = -1;
		}
	}

    return number_of_activephase;
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
	const double ps_max,
	const double ps_min,
	const double theta_th,
	const double *aaa_DO,
	const double *www_DO,
	const double *pmb_DO,
	const double *aaa_DW,
	const double *www_DW,
	const double *pmb_DW,
	const double *df_ij,
	const double *pf,
	      double *pff,
	const double *ps,
	      double *pss
)
{
	for(int j=0; j<ny; j++){
	for(int i=0; i<nx; i++){	

		const int ix = j*nx + i;

        int number_of_activephase = 0;
		int activephase_flag[number_of_phase];
        number_of_activephase = active_phase_check(ix,nx,ny,number_of_phase, pi, theta_th, ps, activephase_flag);

		double psn[number_of_phase];
		double psl[number_of_phase];
		for(int m0=0; m0<number_of_phase; m0++){

			const double psi_ix = ps[m0*nx*ny+ix];
			const double psi_im = (i!=0   ) ? ps[m0*nx*ny+ix-1 ] : ps[m0*nx*ny+ix+1 ];
			const double psi_ip = (i!=nx-1) ? ps[m0*nx*ny+ix+1 ] : ps[m0*nx*ny+ix-1 ];
			const double psi_jm = (j!=0   ) ? ps[m0*nx*ny+ix-nx] : ps[m0*nx*ny+ix+nx];
			const double psi_jp = (j!=ny-1) ? ps[m0*nx*ny+ix+nx] : ps[m0*nx*ny+ix-nx];

			double psx = (psi_ip-psi_im)*0.5*rdx;
			double psy = (psi_jp-psi_jm)*0.5*rdy;

			psn[m0] = psx*psx + psy*psy;
			psl[m0] = (psi_ip + psi_im + - 2.0*psi_ix)*rdx*rdx + (psi_jp + psi_jm - 2.0*psi_ix)*rdy*rdy;
		}

		for(int m0=0; m0<number_of_phase; m0++){

			double       pssi_ix = ps_min;
			double       pffi_ix = pf_min;
			const double psi_ix = ps[m0*nx*ny+ix];
			const double pfi_ix = pf[m0*nx*ny+ix];

			if(number_of_activephase!=1 && activephase_flag[m0]==0){

				if(number_of_activephase==2){

                    for(int m1=0; m1<number_of_phase; m1++){

                        if(m0!=m1 && activephase_flag[m1]==0){

				            const double M_ij = pmb_DO[number_of_phase*m1 + m0];
				            const double w_ij = www_DO[number_of_phase*m1 + m0];
				            const double a_ij = aaa_DO[number_of_phase*m1 + m0];
				            const double dfij = df_ij[ number_of_phase*m1 + m0];

				            const double dp1 = a_ij*a_ij*psl[m0];
				            const double dp2 = (2.0*w_ij - a_ij*a_ij*psn[m0])*tan(psi_ix);
				            const double dp3 = 8.0/pi*dfij;

				            pssi_ix = psi_ix + M_ij*(dp1+dp2+dp3)*dt;
                            break;
                        }
                    }

				}else{

					double dp = 0.0;
					for(int m1=0; m1<number_of_phase; m1++){

						if(m0!=m1 && activephase_flag[m1]==0){

							const double pfj_ix = pf[m1*nx*ny+ix];

							const double M_ij = pmb_DO[number_of_phase*m1 + m0];

							double       dp1 = 0.0;
							double       dp2 = 0.0;
							const double dp3 = -8.0/pi*sqrt((pfi_ix+1.0)*(pfj_ix+1.0))*df_ij[ number_of_phase*m1 + m0];

							for(int m2=0; m2<number_of_phase; m2++){

								if(activephase_flag[m2]==0){

									const double pfk_ix = pf[m2*nx*ny+ix];

									const double w_ik = www_DO[number_of_phase*m2 + m0];
									const double w_jk = www_DO[number_of_phase*m2 + m1];
									const double a_ik = aaa_DO[number_of_phase*m2 + m0];
									const double a_jk = aaa_DO[number_of_phase*m2 + m1];

									dp1 += (w_ik-w_jk)*(pfk_ix+1.0);
									dp2 += 0.5*(a_ik*a_ik-a_jk*a_jk)*(psl[m2]*sqrt(1.0-pfk_ix*pfk_ix) - pfk_ix*psn[m2]);
								}
							}//end m2 loop

							dp += M_ij*(dp1 + dp2 + dp3);
						}
					}//end m1 loop
		
					pssi_ix = psi_ix -2.0/(number_of_activephase*sqrt(1.0-pfi_ix*pfi_ix))*dp*dt;
				}

			}else{

				int number_of_coexistencephase = 0;
				double ps_sum   = 0.;
				double M_sum    = 0.;
				double a_sq_sum = 0.;
				double w_sum    = 0.;
				double M_ave    = 0.;
				double a_sq_ave = 0.;
				double w_ave    = 0.;

				for(int m1=0; m1<number_of_phase; m1++){

					if(m1!=m0 && ps_min<ps[m1*nx*ny+ix]){

						number_of_coexistencephase++;
                        const double dp = tanh(ps[m1*nx*ny+ix]) - pf_min;
						ps_sum   += dp;
						M_sum    += pmb_DW[number_of_phase*m1 + m0]*dp;
						a_sq_sum += aaa_DW[number_of_phase*m1 + m0]*aaa_DW[number_of_phase*m1 + m0]*dp;
						w_sum    += www_DW[number_of_phase*m1 + m0]*dp;
					}
					if(m1!=m0){
						M_ave    += pmb_DW[number_of_phase*m1 + m0];
						a_sq_ave += aaa_DW[number_of_phase*m1 + m0]*aaa_DW[number_of_phase*m1 + m0];
						w_ave    += www_DW[number_of_phase*m1 + m0];
                    }
				}

				const double M_ij = (number_of_coexistencephase==0) ? M_ave/(number_of_phase-1)          : M_sum/ps_sum;
				const double a_ij = (number_of_coexistencephase==0) ? sqrt(a_sq_ave/(number_of_phase-1)) : sqrt(a_sq_sum/ps_sum);
				const double w_ij = (number_of_coexistencephase==0) ? w_ave/(number_of_phase-1)          : w_sum/ps_sum;

				const double dp1 = psl[m0]*a_ij*a_ij;
				const double dp2 = (w_ij-2.0*a_ij*a_ij*psn[m0])*tanh(psi_ix);

				pssi_ix = psi_ix + M_ij*(dp1 + dp2)*dt;
			}

			if(ps_max < pssi_ix){
				pssi_ix = ps_max;
			}else if(pssi_ix < ps_min){
				pssi_ix = ps_min;
			}

			if(pssi_ix < -0.5*pi){
				pffi_ix = pf_min;
			}else if(0.5*pi < pssi_ix){
				pffi_ix = pf_max;
			}else{
				pffi_ix = sin(pssi_ix);
			}

			pss[m0*nx*ny+ix] = pssi_ix;
			pff[m0*nx*ny+ix] = pffi_ix;
		
		}//end m0 loop
		
        number_of_activephase = active_phase_check(ix,nx,ny,number_of_phase, pi, theta_th, pss, activephase_flag);

		if(number_of_activephase!=1){

			double psum = 0.0;
			for(int m0=0; m0<number_of_phase; m0++){
		
				if(activephase_flag[m0]==0) psum += (pff[m0*nx*ny+ix]+1.0)*0.5;
			}

			for(int m0=0; m0<number_of_phase; m0++){

				if(activephase_flag[m0]==0){

					pff[m0*nx*ny+ix] = (pff[m0*nx*ny+ix]+1.0)/psum - 1.0;
					pss[m0*nx*ny+ix] = asin(pff[m0*nx*ny+ix]);

				}
			}
		}
	}
    }
}

//=================================================-=================================================//
//********************************************Main Function******************************************//
//=================================================-=================================================//

int main(int argc, char** argv)
{
    double *PF, *PS;
    double *PPF, *PPS;
    double *GAM_ij, *MOB_ij, *DF_ij;
    double *AAA_DO, *WWW_DO, *PMB_DO;
    double *AAA_DW, *WWW_DW, *PMB_DW;

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
    const double ps_max   = +(4.0*dx/delta*pi + 0.5*pi); //Max value of psi
    const double ps_min   = -(4.0*dx/delta*pi + 0.5*pi); //Min value of psi
    const double theta_th =  3./180.*pi;

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

    AAA_DW = (double *) malloc(number_of_phase*number_of_phase*sizeof(double)); if(AAA_DW == NULL){fprintf(stderr,"Malloc error\n"); exit(1);}
    WWW_DW = (double *) malloc(number_of_phase*number_of_phase*sizeof(double)); if(WWW_DW == NULL){fprintf(stderr,"Malloc error\n"); exit(1);}
    PMB_DW = (double *) malloc(number_of_phase*number_of_phase*sizeof(double)); if(PMB_DW == NULL){fprintf(stderr,"Malloc error\n"); exit(1);}


 	for(int m0=0; m0<number_of_phase; m0++){
	for(int m1=0; m1<number_of_phase; m1++){

        const int mij = number_of_phase*m1 + m0;

        if(m0!=m1){

			AAA_DO[mij] = 2.0*sqrt(2.0*delta*GAM_ij[mij])/pi;
			WWW_DO[mij] = 4.0*GAM_ij[mij]/delta;
			PMB_DO[mij] = pi*pi*MOB_ij[mij]/(8.0*delta);

			AAA_DW[mij] = sqrt(3.0*GAM_ij[mij]/(pi/delta));
			WWW_DW[mij] = 6.0*(pi/delta)*GAM_ij[mij];
			PMB_DW[mij] = (pi/delta)*MOB_ij[mij]/3.0;

		}else{

			AAA_DO[mij] = 0.0;
			WWW_DO[mij] = 0.0;
			PMB_DO[mij] = 0.0;

			AAA_DW[mij] = 0.0;
			WWW_DW[mij] = 0.0;
			PMB_DW[mij] = 0.0;
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
	PS  = (double *) malloc(number_of_phase*ngx*ngy*sizeof(double)); if(PS  == NULL){fprintf(stderr,"Malloc error PS\n"); exit(1);}
	PPF = (double *) malloc(number_of_phase*ngx*ngy*sizeof(double)); if(PPF == NULL){fprintf(stderr,"Malloc error PPF\n"); exit(1);}
	PPS = (double *) malloc(number_of_phase*ngx*ngy*sizeof(double)); if(PPS == NULL){fprintf(stderr,"Malloc error PPS\n"); exit(1);}

//<<Initial profile Setting>>//
	for(int j=0; j<ngy; j++){
	for(int i=0; i<ngx; i++){
	
		const int ix = j*ngx + i;
		
		double xx = ((double)i - ncx/2.0)*dx;
		double yy = ((double)j - ncy/2.0)*dy;
		double rr = sqrt(xx*xx + yy*yy) - ncx/16.0*dx;
		
		if(pi*rr/delta < ps_min){
			PS[0*ngx*ngy + ix] = ps_max;
		}else if(ps_max < pi*rr/delta){
			PS[0*ngx*ngy + ix] = ps_min;
		}else{
			PS[0*ngx*ngy + ix] = -pi*rr/delta;
		}
	
		if(ps_max < pi*rr/delta){
			PS[1*ngx*ngy + ix] = ps_max;
		}else if(pi*rr/delta < ps_min){
			PS[1*ngx*ngy + ix] = ps_min;
		}else{
			PS[1*ngx*ngy + ix] = +pi*rr/delta;
		}
		
		for(int m=0; m<number_of_phase; m++){
			if(PS[m*ngx*ngy + ix] < -0.5*pi){
				PF[m*ngx*ngy + ix] = pf_min;
			}else if(0.5*pi < PS[m*ngx*ngy + ix]){
				PF[m*ngx*ngy + ix] = pf_max;
			}else{
				PF[m*ngx*ngy + ix] = sin(PS[m*ngx*ngy + ix]);
			}
		}
	
		const double pf0 = 1.0 - (PF[0*ngx*ngy + ix]+1.0)*0.5;
		for(int m=1; m<number_of_phase; m++){

			if(-0.5*pi < PS[m*ngx*ngy + ix] && PS[m*ngx*ngy + ix] < 0.5*pi){

				PF[m*ngx*ngy + ix] = (PF[m*ngx*ngy + ix]+1.0)*pf0-1.0;
				PS[m*ngx*ngy + ix] = asin(PF[m*ngx*ngy + ix]);
			}
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
        fprintf(fini, "<<<<< Dimensionless values >>>>> \n");
        fprintf(fini, "Psi_max                   [pi]    =  %16.7e\n", ps_max/pi);
        fprintf(fini, "Theta_th                  [pi]    =  %16.7e\n", theta_th/pi);

        fclose(fini);
    }

    //<<interface velocity along x-axies initialaize>>//
    double interface_pos_x1_at_y00h = x_axis_interface_position_exploration(ngx, ngy, 0,    0, dx, PS);
    double interface_pos_x0_at_y00h = interface_pos_x1_at_y00h;
    double interface_pos_x1_at_y05h = x_axis_interface_position_exploration(ngx, ngy, 0, ncy/2, dx, PS);
    double interface_pos_x0_at_y05h = interface_pos_x1_at_y05h;

//<<<<<<<<<<<<<<<<<<<<<<<<< main loop start >>>>>>>>>>>>>>>>>>>>>>>>>//

	for(int nstep=1; nstep<=nend ; nstep++){

        if(nstep%100==0) printf("nstep = %8d\n", nstep);

        //Finite difference computation of multi-phase-field equation
        Multi_Phase_field(number_of_phase,ngx,ngy, pi, 1.0/dx, 1.0/dy, dt, pf_max,pf_min,ps_max,ps_min,theta_th,
            AAA_DO,WWW_DO,PMB_DO, AAA_DW,WWW_DW,PMB_DW, DF_ij, PF, PPF, PS, PPS);
		swap(&PF,  &PPF);
		swap(&PS,  &PPS);

        //<<interface velocity along x-axies>>//
		if(nstep%nvel==0){
            if(1){
                interface_pos_x0_at_y00h = interface_pos_x1_at_y00h;
                interface_pos_x1_at_y00h = x_axis_interface_position_exploration(ngx, ngy, 0,    0, dx, PS);
                double interface_vel_x1_at_y00h = (interface_pos_x1_at_y00h-interface_pos_x0_at_y00h)/(nvel*dt);
                
		        FILE *fp_v = fopen("interface_velocity_x_at_y00h.dat","a");
		        fprintf(fp_v, "%6d %16.7e %16.7e %16.7e\n", nstep, (double)nstep*dt, interface_pos_x1_at_y00h, interface_vel_x1_at_y00h);
		        fclose(fp_v);
            }

            if(1){
                interface_pos_x0_at_y05h = interface_pos_x1_at_y05h;
                interface_pos_x1_at_y05h = x_axis_interface_position_exploration(ngx, ngy, 0, ncy/2, dx, PS);
                double interface_vel_x1_at_y05h = (interface_pos_x1_at_y05h-interface_pos_x0_at_y05h)/(nvel*dt);

		        FILE *fp_v = fopen("interface_velocity_x_at_y05h.dat","a");
		        fprintf(fp_v, "%6d %16.7e %16.7e %16.7e\n", nstep, (double)nstep*dt, interface_pos_x1_at_y05h, interface_vel_x1_at_y05h);
		        fclose(fp_v);
            }
        }

        if(nstep%nout==0 || nstep==1){
            char fvti[128];
            snprintf(fvti,sizeof(fvti),"paraview_%08d",nstep);	
            output_vtkdata(fvti,number_of_phase, ngx, ngy, dx, dy, PF, PS);
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
