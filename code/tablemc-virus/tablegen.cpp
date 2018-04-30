#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h> //floating point limits, need FLT_MAX
#include "rotations.h"
#include "tables.h"
#include "mc.h"
#include "mt.h"
#include "util.h"
#ifndef __unix__
//Cannot do memory-mapped tables when not UNIX.
//Actually we could use MapViewOfFile on Windows, but we don't want to bother.
#define NO_MMAP_TABLES
#endif
#ifndef NO_MMAP_TABLES
#include <sys/types.h>
#include <sys/mman.h>
#include <err.h>
#include <fcntl.h>
#include <unistd.h>
#endif

//For smoothing.  Cells outside CUTOFF*rho do not count towards average.
#define MIN_WEIGHT 1E-6
#define TABLESIZE  131072
#define LMAX       1000

table::table(const char * fname, int part, int numparts)
{
    hdr.totalpoints=0;
    energy=NULL;
    generate_table(fname,part,numparts);
}

void table::read_table_header_info(const char * fname)
{
    /*char cwd[1024];
    if (getcwd(cwd, sizeof(cwd)) != NULL)
        printf("Current working dir: %s\n", cwd);*/
    int sub_dip,rdie,i,ir;
    char buffer[255];
    char * p;
    hdr.headersize=sizeof(hdr);
    hdr.version=VERSION;
    hdr.options=0; /*for now*/
    printf("Reading control file: %s\n",fname);
    FILE * f;
    f=fopen(fname,"r");
    if (f==NULL) {
        printf("FATAL ERROR: file %s is not found\n",fname);
        die();
    }
    fscanf(f,"%lg %lg %lg %lg %lg\n",&hdr.minr,&hdr.maxr,&hdr.dr,&hdr.dsph,&hdr.deuler);
    //printf("*** %g %g %g %g\n",hdr.maxx,hdr.dx,hdr.dphipsi,hdr.dtheta);
    hdr.dsph=hdr.dsph*DEG_TO_RAD;
    hdr.deuler=hdr.deuler*DEG_TO_RAD;
    //hdr.nr=ceil(hdr.maxr/hdr.dr)+1;
    hdr.nsphtheta=floor(M_PI/hdr.dsph);
    hdr.nsphphi=hdr.nsphtheta*2;
    hdr.ntheta=floor(M_PI/hdr.deuler);
    hdr.nphipsi=hdr.ntheta*2;
    //hdr.dr=hdr.maxr/(hdr.nr-1); /*ensure that resolutions are consistent with number of grid points*/
    hdr.dsph=M_PI/hdr.nsphtheta; // Not clear how to handle the case when
    //hdr.deuler=2.0*M_PI/hdr.nphipsi;
    hdr.deuler=M_PI/hdr.ntheta;
    hdr.decr=1; //not decimated.
    hdr.decsph=1;
    hdr.deceuler=1;
    //hdr.ntrans=hdr.nr * hdr.nsphtheta * hdr.nsphphi;
    hdr.norient=hdr.nphipsi * hdr.nphipsi * hdr.ntheta;
    //hdr.clashpoints=hdr.nsphtheta*hdr.nsphphi*hdr.nphipsi*hdr.ntheta*hdr.nphipsi;
    hdr.rfactor=1+hdr.dr/hdr.minr;
    hdr.logrfactor=log(hdr.rfactor);
    hdr.nr=(long int) (log(hdr.maxr/hdr.minr)/hdr.logrfactor+0.5)+1;
    hdr.ntrans=hdr.nr * hdr.nsphtheta * hdr.nsphphi;
    hdr.totalpoints=(long long int) hdr.ntrans * (long long int) hdr.norient;
    hdr.clashpoints=0;
    hdr.invalidpoints=0;
    //hdr.totalpoints=(long long int) hdr.ntrans * (long long int) hdr.norient;
    /*fscanf(f,"%lg %lg %lg %lg %lg\n",&hdr.go_hardcore,&hdr.go_cutoff,&hdr.go_delta,&hdr.native_energy,&hdr.nonnative_energy);*/
    fgets(buffer,sizeof(buffer),f);
    p=strtok(buffer," "); //put a null after the end of the filename
    strncpy(hdr.go_model_map_fname,buffer,sizeof(hdr.go_model_map_fname));
    trim_string(hdr.go_model_map_fname);
    p+=strlen(p)+1;//advanced past the null
    read_go_params(p,&hdr.go_params);
    /*fgets(hdr.go_model_map_fname,sizeof(hdr.go_model_map_fname),f);
    trim_string(hdr.go_model_map_fname);*/
    fscanf(f,"%lg %lg %lg %d %lg\n",&hdr.beta,&hdr.sph_alpha,&hdr.euler_alpha,&hdr.n_invalid_dev,&hdr.en_invalid_margin); //, &hdr.n_invalid_dev, &hdr.en_clash,&hdr.en_invalid_margin);
    hdr.sph_alpha=hdr.sph_alpha*DEG_TO_RAD;
    hdr.sph_alpha=hdr.sph_alpha/SPHERICAL_DIFFUSION_CONST;
    hdr.sph_alpha=hdr.sph_alpha*hdr.sph_alpha;
    hdr.euler_alpha=hdr.euler_alpha*DEG_TO_RAD;
    hdr.euler_alpha=hdr.euler_alpha/SO3_DIFFUSION_CONST;
    hdr.euler_alpha=hdr.euler_alpha*hdr.euler_alpha;
    if (hdr.beta<=0.0) hdr.beta=300.0; /*default 300 K for association free energy*/
    hdr.beta=1/(KBOLTZ*hdr.beta);
    hdr.en_smooth_limit=DUMMY_ENERGY;
    if (hdr.en_invalid_margin<=0.0) hdr.n_invalid_dev=0;
    if (hdr.n_invalid_dev<=0) hdr.en_invalid_margin=0.0;
    /*hdr.en_invalid_margin=0.0;
    hdr.n_invalid_dev=0;*/
    //hdr.en_clash=0.0;
    //if ((hdr.nsmoothr>0) || (hdr.nsmoothsph>0) || (hdr.nsmootheuler>0)) hdr.options = hdr.options | OP_BOLTZMANN_AVERAGED;
    if ((hdr.sph_alpha>0) || (hdr.euler_alpha>0)) hdr.options = hdr.options | OP_BOLTZMANN_AVERAGED;
    if ((hdr.options & OP_BOLTZMANN_AVERAGED) && (hdr.en_invalid_margin>0.0)) { //sanity check
        printf("You may not apply smoothing and mark invalid cells at the same time.\n");
        die();
    }
    //fscanf(f,"%d %ld %lg\n",&hdr.nboltzmannaverage,&hdr.seed,&hdr.boltzmannbeta)
    //->boltzmannbeta = 1/(KBOLTZ*hdr.boltzmannbeta);
    //if (sub_dip) hdr.options = hdr.options | OP_DIPOLE_SUBTRACTED;
    //if (rdie) hdr.options = hdr.options | OP_DIST_DEP_DIELECTRIC;
    if (sizeof(energy_t)==sizeof(float)) hdr.options = hdr.options | OP_SINGLE_PRECISION;
    /*if (hdr.nboltzmannaverage!=0) {
        hdr.options = hdr.options | OP_BOLTZMANN_AVERAGED;
        init_genrand(hdr.seed);
    }*/
    time(&hdr.created); /*current date/time*/
    fgets(hdr.ref_frag_fname,sizeof(hdr.ref_frag_fname),f);
    trim_string(hdr.ref_frag_fname);
    fgets(hdr.other_frag_fname,sizeof(hdr.other_frag_fname),f);
    trim_string(hdr.other_frag_fname);
    fclose(f);
    radial_table =(double *) checkalloc(hdr.nr,sizeof(double));//destroyed by destructor
    for (ir=0; ir<hdr.nr; ir++) radial_table[ir]=hdr.minr*pow(hdr.rfactor,ir);
    one_minr=1.0/hdr.minr;
    one_logrfactor=1.0/hdr.logrfactor;
    one_dsph=1.0/hdr.dsph;
    one_deuler=1.0/hdr.deuler;
}

/*void table::fill_clash_table(void)
{
    long int iphi,itheta,ipsi,isphtheta,isphphi;
    long long int index,points;
    double rmin,rmax,r,x[3];
    double sphtheta,sphphi,phi,theta,psi;
    int clash;
        hdr.minr=hdr.maxr;
        hdr.max_clash=0;
        oldcenter[0]=0.0;
        oldcenter[1]=0.0;
        oldcenter[2]=0.0;
        oldorient[0]=1.0;
        oldorient[1]=0.0;
        oldorient[2]=0.0;
        oldorient[3]=0.0;
        points=0;
        update_coords(0,oldcenter,oldorient,oldcoords);
        for (iphi=0; iphi<hdr.nphipsi; iphi++)
            for (itheta=0; itheta<hdr.ntheta; itheta++)
                for (ipsi=0; ipsi<hdr.nphipsi; ipsi++) {
                    phi=iphi*hdr.deuler;
                    theta=itheta*hdr.deuler;
                    psi=ipsi*hdr.deuler;
            for (isphtheta=0; isphtheta<hdr.nsphtheta; isphtheta++)
                for (isphphi=0; isphphi<hdr.nsphphi; isphphi++) {
                    rmin=0.0;
                    rmax=hdr.maxr;
                    sphtheta=isphtheta*hdr.dsph;
                    sphphi=isphphi*hdr.dsph;
                    while (fabs(rmax-rmin)>0.01) {
                        r=(rmax+rmin)/2.0;
                        sph_to_cart(r,sphtheta,sphphi,&x[0]);
                        oldcenter[3]=x[0];
                        oldcenter[4]=x[1];
                        oldcenter[5]=x[2];
                        euler_to_quat(phi,theta,psi,&oldorient[4]);
                        update_coords(1,oldcenter,oldorient,oldcoords);
                        clash=check_clash(hdr.clashfactor,frags[0].natom,&frags[0].types[0],&oldcoords[3*fragstart[0]],frags[1].natom,&frags[1].types[0],&oldcoords[3*fragstart[1]]);
                        if (clash) rmin=r; else rmax=r;
                    }
                    //minimum clash distance is between rmin and rmax
                    index=clash_table_index(isphtheta,isphphi,iphi,itheta,ipsi);
                    clash_table[index]=rmin;
                    if (rmin<hdr.minr) hdr.minr=rmin;
                    if (rmin>hdr.max_clash) hdr.max_clash=rmin;
                    points++;
                    if ((points%100000)==0) printf("%lld\n",points);
                }
        }
    printf("Minimum and maximum clash distances:         %.2f A- %.2f A\n",hdr.minr,hdr.max_clash);
    //We now need to set up the r factor and nr count.
    hdr.rfactor=1+hdr.dr/hdr.minr;
    hdr.logrfactor=log(hdr.rfactor);
    hdr.nr=(long int) (log(hdr.maxr/hdr.minr)/hdr.logrfactor+0.5)+1;
    hdr.ntrans=hdr.nr * hdr.nsphtheta * hdr.nsphphi;
    hdr.norient=hdr.nphipsi * hdr.nphipsi * hdr.ntheta;
    hdr.totalpoints=(long long int) hdr.ntrans * (long long int) hdr.norient;
}*/



double table::get_pair_energy(go_model_info * go_model, fragmenttype * frag1, fragmenttype * frag2, double * coords1, double * coords2, double r, double sphtheta, double sphphi, double phi, double theta, double psi)
{
    double x[3],q[4];
                        sph_to_cart(r,sphtheta,sphphi,&x[0]);


                        euler_to_quat(phi,theta,psi,&q[0]);
                        //update_coords(1,oldcenter,oldorient,oldcoords);
                        frag2->get_coords(&x[0],&q[0],coords2);
                            //for (i=0;i<frag2.natom;i++){
                            //    translated[3*i]=rotated[3*i]+x;
                            //    translated[3*i+1]=rotated[3*i+1]+y;
                            //    translated[3*i+2]=rotated[3*i+2]+z;
                            //}
                            //return interaction_energy(0,1,oldcenter,oldorient,oldcoords);
    //return ffield->exact_interaction_energy(FALSE,0.0,0.0,hdr.eps,(hdr.options & OP_DIST_DEP_DIELECTRIC)!=0,frag1,coords1,frag2,coords2);
    return go_model->energy(false,0,0,&hdr.go_params,frag1,coords1,frag2,coords2);
}

void table::fill_table(go_model_info * go_model, fragmenttype * frag1, fragmenttype * frag2)
{
    long int ir,isphtheta,isphphi,iphi,itheta,ipsi;
    long long int index,points,index2;
    int i,nclash,ndev;
    double r,sphtheta,sphphi,x[3],q[4],phi,theta,psi,en[13],endev,diff,enx,edip,enmin,rotmatrix[3][3];
    double * coords1;
    double * coords2;
    FILE * test_output;
    coords1=(double *) malloc(3*frag1->natom*sizeof(double));
    coords2=(double *) malloc(3*frag2->natom*sizeof(double));
    enmin=10000000.0;
    //ix=7; iy=17; iz=14; iphi=9; itheta=3; ipsi=11;
    //ix=6; iy=8; iz=6; iphi=6; itheta=1; ipsi=2;
    points=1;
    hdr.invalidpoints=0;
    hdr.clashpoints=0;
        /*Set fragment 0 to the standard position/orientation.*/
        x[0]=0.0;
        x[1]=0.0;
        x[2]=0.0;
        q[0]=1.0;
        q[1]=0.0;
        q[2]=0.0;
        q[3]=0.0;
        //update_coords(0,oldcenter,oldorient,oldcoords);
        frag1->get_coords(&x[0],&q[0],coords1);
        ir=23; isphtheta=1; isphphi=5; iphi=2; itheta=2; ipsi=3;
        for (iphi=0; iphi<hdr.nphipsi; iphi++)
            for (itheta=0; itheta<hdr.ntheta; itheta++)
                for (ipsi=0; ipsi<hdr.nphipsi; ipsi++) {
                    phi=iphi*hdr.deuler;
                    theta=(itheta+0.5)*hdr.deuler;
                    psi=ipsi*hdr.deuler;
                    //euler_to_matrix(phi,theta,psi,rotmatrix);
                    //print_matrix(rotmatrix);
                    //for (i=0;i<frag2.natom;i++) matmul(rotmatrix,&frag2.refgeom[3*i],&rotated[3*i]);

        for (isphtheta=0; isphtheta<hdr.nsphtheta; isphtheta++)
            for (isphphi=0; isphphi<hdr.nsphphi; isphphi++) {
                for (ir=0; ir<hdr.nr; ir++) {
                    index=calculate_index(ir,isphtheta,isphphi,iphi,itheta,ipsi);
                    if (index>=hdr.totalpoints) {
                        printf("Error in table generation.\n");
                        die();
                    }
                    if ((index<hdr.startindex) || (index>hdr.endindex)) continue; //out of range
                    r=radial_table[ir];
                    //printf("*** %d %.4f\n",ir,r);
                    sphtheta=(isphtheta+0.5)*hdr.dsph;
                    sphphi=isphphi*hdr.dsph;
                    en[0]=get_pair_energy(go_model,frag1,frag2,coords1,coords2,r,sphtheta,sphphi,phi,theta,psi);
                    /*if (hdr.en_invalid_margin>0) {
                        en[1]=get_pair_energy(go_model,frag1,frag2,coords1,coords2,r*sqrt(hdr.rfactor),sphtheta,sphphi,phi,theta,psi);
                        en[2]=get_pair_energy(go_model,frag1,frag2,coords1,coords2,r/sqrt(hdr.rfactor),sphtheta,sphphi,phi,theta,psi);
                        en[3]=get_pair_energy(go_model,frag1,frag2,coords1,coords2,r,sphtheta+0.5*hdr.dsph,sphphi,phi,theta,psi);
                        en[4]=get_pair_energy(go_model,frag1,frag2,coords1,coords2,r,sphtheta-0.5*hdr.dsph,sphphi,phi,theta,psi);
                        en[5]=get_pair_energy(go_model,frag1,frag2,coords1,coords2,r,sphtheta,sphphi+0.5*hdr.dsph,phi,theta,psi);
                        en[6]=get_pair_energy(go_model,frag1,frag2,coords1,coords2,r,sphtheta,sphphi-0.5*hdr.dsph,phi,theta,psi);
                        en[7]=get_pair_energy(go_model,frag1,frag2,coords1,coords2,r,sphtheta,sphphi,phi+0.5*hdr.deuler,theta,psi);
                        en[8]=get_pair_energy(go_model,frag1,frag2,coords1,coords2,r,sphtheta,sphphi,phi-0.5*hdr.deuler,theta,psi);
                        en[9]=get_pair_energy(go_model,frag1,frag2,coords1,coords2,r,sphtheta,sphphi,phi,theta+0.5*hdr.deuler,psi);
                        en[10]=get_pair_energy(go_model,frag1,frag2,coords1,coords2,r,sphtheta,sphphi,phi,theta-0.5*hdr.deuler,psi);
                        en[11]=get_pair_energy(go_model,frag1,frag2,coords1,coords2,r,sphtheta,sphphi,phi,theta,psi+0.5*hdr.deuler);
                        en[12]=get_pair_energy(go_model,frag1,frag2,coords1,coords2,r,sphtheta,sphphi,phi,theta,psi-0.5*hdr.deuler);
                        endev=0.0;
                        ndev=0;
                        for (i=1; i<=12; i++) {
                            //if (en[i]>hdr.en_clash) nclash++;
                            diff=fabs(en[i]-en[0]);
                          	     if (diff>hdr.en_invalid_margin) ndev++;
                        	}
                    }*/

                            /*if (nclash==13) {
                                energy[index]=DUMMY_ENERGY;
                                hdr.clashpoints++;
                            } else */
                        /*if ((hdr.en_invalid_margin>0.0) && (ndev>=hdr.n_invalid_dev)) {
                            energy[index]=INVALID_ENERGY; //aries too much
                            hdr.invalidpoints++;
                        } else*/ if (en[0]>FLT_MAX) {
                            energy[index]=DUMMY_ENERGY;
                        } else energy[index]=en[0];

                        //if (en[0]>DUMMY_ENERGY) en[0]=DUMMY_ENERGY;

				//}
                            //exact_interaction_energy(hdr.eps,frag1.natom,&frag1.types[0],&frag1.refgeom[0],frag2.natom,&frag2.types[0],translated);
                            //if (table[index]<enmin) enmin=table[index];
                            if (en[0]<enmin) enmin=en[0];
                            /*if (en<-1000.0) {
                                printf("index: %lld %lld\n",index,points);
                                test_output=fopen("fragments.pdb","w");
                                write_frame_pdb(test_output,0,oldcoords);
                                fclose(test_output);
                                printf("%.2f %.2f %.2f %.2f %.2f %.2f %g %g\n",r,sphtheta*RAD_TO_DEG,sphphi*RAD_TO_DEG,phi*RAD_TO_DEG,theta*RAD_TO_DEG,psi*RAD_TO_DEG,en,enmin);
                                die();
                            }*/
                            if ((points%1000)==0) {
                                printf("%lld points completed\n",points);
                                fflush(stdout);
                            }
                            points++;
                            //printf("%.2f %.2f %.2f %.2f %.2f %.2f %.2f\n",x,y,z,phi*RAD_TO_DEG,theta*RAD_TO_DEG,psi*RAD_TO_DEG,table[index]);
                            //printf("%d %d %d %d %d %d %.2f\n",ix,iy,iz,iphi,itheta,ipsi,table[index]);
                    }
            }
    }
    free(coords1);
    free(coords2);
    printf("index: %lld\n",index);
    /*test_output=fopen("fragments.pdb","w");
    write_pair_pdb2(test_output,frag1,frag2,coords1,coords2);
    fclose(test_output);*/
    printf("Minimum interaction energy: %g\n",enmin);
    printf("Clash points: %lld\n",hdr.clashpoints);
    printf("Invalid points: %lld (%.2f%%)\n",hdr.invalidpoints,100*((double)hdr.invalidpoints/(double)hdr.totalpoints));

}


void table::fake_fill_table(void)
{
	long long int index;
	long int iphi, itheta, ipsi, ir, isphtheta, isphphi;
	for (index=0; index<hdr.totalpoints; index++) energy[index]=DUMMY_ENERGY;
	/*for (ir=0; ir<hdr.nr; ir++)
	for (isphtheta=0; isphtheta<hdr.nsphtheta; isphtheta++)
	for (isphphi=0; isphphi<hdr.nsphphi; isphphi++) {
		index=calculate_index(ir, isphtheta, isphphi,0,0,0);
		energy[index]=0.0;
	}*/
	for (iphi=0; iphi<hdr.nphipsi; iphi++)
	for (itheta=0; itheta<hdr.ntheta; itheta++)
	for (ipsi=0; ipsi<hdr.nphipsi; ipsi++) {
		index=calculate_index(0,0,0,iphi,itheta,ipsi);
		energy[index]=0.0;
	}
}


//center -- indices of cell, displacement, indices of neighbor
//returns whether the cell exists or not.
//r, theta, phi, phi', theta', psi'
int table::find_neighbor(long int * indices, long int * disp, long int * neighbor)
{
     int i,swap,boundscheck;
     for (i=0; i<6; i++) neighbor[i]=indices[i]+disp[i];
     //off the end of the r coordinate?
     if ((neighbor[0]<0) || (neighbor[0]>=hdr.nr)) return FALSE;
     //if we are close to the north or south pole, but we move over it
     if (neighbor[1]<0) {
         neighbor[1]=-neighbor[1];
         neighbor[2]+=hdr.nsphphi/2;
     }
     //Go from nsphtheta, nsphtheta+1, etc. to nsphtheta-1, nsphtheta-2, etc.
     if (neighbor[1]>=hdr.nsphtheta) {
         neighbor[1]=(hdr.nsphtheta)-(neighbor[1]-(hdr.nsphtheta-1));
         neighbor[2]+=hdr.nsphphi/2;
     }
     //If we change the sign of theta', q_1 and q_2 change sign but q_0 and q_3 don't.
     //Therefore we must change the sign of (phi'-psi') while keeping (phi'+psi') constant.
     //Interchange phi' and psi'.
     if (neighbor[4]<0) {
         neighbor[4]=-neighbor[4];
         swap=neighbor[3];
         neighbor[3]=neighbor[5];
         neighbor[5]=swap;
     }
     if (neighbor[4]>=hdr.ntheta) {
         neighbor[4]=(hdr.ntheta)-(neighbor[4]-(hdr.ntheta-1));
         swap=neighbor[3];
         neighbor[3]=neighbor[5];
         neighbor[5]=swap;
     }
     /*reduce phi, phi', and psi' modulo */
     if (neighbor[2]<0) neighbor[2]+=hdr.nsphphi;
     if (neighbor[3]<0) neighbor[3]+=hdr.nphipsi;
     if (neighbor[5]<0) neighbor[5]+=hdr.nphipsi;
     neighbor[2]=neighbor[2]%hdr.nsphphi;
     neighbor[3]=neighbor[3]%hdr.nphipsi;
     neighbor[5]=neighbor[5]%hdr.nphipsi;
     //Check that everything is within bounds.  If not, print a warning and return FALSE.
     /*boundscheck=((neighbor[0]>=0) && (neighbor[0]<hdr.nr) && (neighbor[1]>=0) && (neighbor[1]<hdr.nsphtheta) && (neighbor[2]>=0) && (neighbor[2]<hdr.nsphphi));
     boundscheck=boundscheck && ((neighbor[3]>=0) && (neighbor[3]<hdr.nphipsi) && (neighbor[4]>=0) && (neighbor[4]<hdr.ntheta) && (neighbor[5]>=0) && (neighbor[5]<hdr.nphipsi));
     if (!boundscheck) {
        printf("BOUNDS CHECK FAILURE in smoothing\n");
        printf("indices: %ld %ld %ld %ld %ld %ld\n",indices[0],indices[1],indices[2],indices[3],indices[4],indices[5]);
        printf("displacements: %ld %ld %ld %ld %ld %ld\n",disp[0],disp[1],disp[2],disp[3],disp[4],disp[5]);
        printf("neighbor: %ld %ld %ld %ld %ld %ld\n",neighbor[0],neighbor[1],neighbor[2],neighbor[3],neighbor[4],neighbor[5]);
        return FALSE;
     }*/
     return TRUE;
}

/*void table::get_neighbors(int nsmooth, long int * ix, long int * dstart, long int * dstop)
{
    //Ensure circular region at the poles.
    for (i=0; i<6; i++) dstart=-nsmooth;
    for (i=0; i<6; i++) dstop=nsmooth;
    if (ix[1]==0) {
        dstart[1]=0;
        dstart[2]=0;
        dstop[2]=hdr.nsphphi;
    }
    if (ix[1]==(hdr.nsphtheta-1)) {
        dstop[1]=0;
        dstart[2]=0;
        dstop[2]=hdr.nsphphi;
    }
    if (ix[4]==0) {
        dstart[4]=0;
        dstart[3]=0;
        dstop[3]=hdr.nsphphi;
        dstart[5]=0; //Since there is only one DOF
        dstop[5]=0;
    }
    if (ix[4]==(hdr.nsphtheta-1)) {
        dstop[1]=0;
        dstart[2]=0;
        dstop[2]=hdr.nsphphi;
    }
}*/

/*void table::mark_invalid_cells(void)
{
    long int ir,isphtheta,isphphi,iphi,itheta,ipsi,idiff;
    long long int indices[13];
    long int isphphim1,isphphip1,iphim1,iphip1,ipsim1,ipsip1;
    double energies[13],diff,maxdiff;
    int i;
}*/
//For efficiency reasons this subroutine does not use the sph_to_cart, but uses its own formulas,
//to avoid repeatedly calling the trig functions.
void table::boltzmann_average_trans(void)
{
    //long int ix[6],d[6],dstart[6],dstop[6],neigh[6],idiff;
    long int ir,isphtheta,isphphi,iphi,itheta,ipsi,idiff;
    long int ineighr,ineighsphtheta,ineighsphphi,ineighphi,ineightheta,ineighpsi;
    int i,nvalid,exist,nneigh,nsmoothrplus,nsmoothrminus,nsmoothsph,bin;
    double r,sphtheta,sphphi,phi,theta,psi;
    double sum,sum2,en,volume,neighvolume,v1,betainv,weight,weightsum,x1[3],x2[3],cosgamma,sin2gamma,sigma2;
    long long int index,neighborindex,transindex,points;
    double * prob;
    double * costheta, *sintheta, *cosphi, *sinphi;
    double spherical_diffusion_kernel[TABLESIZE],work1[LMAX+1];
    printf("Creating spherical diffusion kernel...\n");
    create_spherical_diffusion_kernel(TABLESIZE,LMAX,hdr.sph_alpha,spherical_diffusion_kernel,work1);
    //double debugsum1, debugsum2;
    costheta=(double *) checkalloc(hdr.nsphtheta,sizeof(double));
    sintheta=(double *) checkalloc(hdr.nsphtheta,sizeof(double));
    for (i=0; i<hdr.nsphtheta; i++) {
        sphtheta=(i+0.5)*hdr.dsph;
        costheta[i]=cos(sphtheta);
        sintheta[i]=sin(sphtheta);
    }
    cosphi=(double *) checkalloc(hdr.nsphphi,sizeof(double));
    sinphi=(double *) checkalloc(hdr.nsphphi,sizeof(double));
    for (i=0; i<hdr.nsphphi; i++) {
        sphphi=i*hdr.dsph;
        cosphi[i]=cos(sphphi);
        sinphi[i]=sin(sphphi);
    }
    points=0;
    prob=(double *) malloc(hdr.totalpoints*sizeof(double));
    if (prob==NULL) {
        printf("Could not allocate space for Boltzmann factors.\n");
        die();
    }
    printf("Computing Boltzmann factors...\n");
    for (index=0; index<hdr.totalpoints; index++) {
        if (energy[index]>50.0) prob[index]=0.0; else prob[index]=exp(-hdr.beta*energy[index]);
    }
    printf("Applying boltzmann averaging in translational space (angular coordinates only).\n");
    betainv=1/hdr.beta;
    v1=(hdr.rfactor-1)*hdr.dsph*hdr.dsph*hdr.deuler*hdr.deuler*hdr.deuler;
    //sigma2=hdr.trans_alpha*hdr.trans_alpha;
    //debugsum1=0.0;
    //debugsum2=0.0;
    for (iphi=0; iphi<hdr.nphipsi; iphi++)
    for (itheta=0; itheta<hdr.ntheta; itheta++)
    for (ipsi=0; ipsi<hdr.nphipsi; ipsi++)
    for (ir=0; ir<hdr.nr; ir++)
    for (isphtheta=0; isphtheta<hdr.nsphtheta; isphtheta++)
    for (isphphi=0; isphphi<hdr.nsphphi; isphphi++) {
        //printf("%d %d %d %d %d %d\n",ir,isphtheta,isphphi,iphi,itheta,ipsi);
        index=calculate_index(ir,isphtheta,isphphi,iphi,itheta,ipsi);
        r=radial_table[ir];
        //sphtheta=(isphtheta+0.5)*hdr.dsph;
        //sphphi=isphphi*hdr.dsph;
        //theta=(itheta+0.5)*hdr.deuler;
        volume=r*r*r*sintheta[isphtheta]*v1; //drop factor of sin(theta)
        nneigh=0;
        sum=0.0;
        weightsum=0.0;
        //sph_to_cart(r,sphtheta,sphphi,x1);
        /*x1[0]=r*sintheta[isphtheta]*cosphi[isphphi];
        x1[1]=r*sintheta[isphtheta]*sinphi[isphphi];
        x1[2]=r*costheta[isphtheta];*/
        ineighphi=iphi;
        ineightheta=itheta;
        ineighpsi=ipsi;
        ineighr=ir;
        //for (ineighr=0; ineighr<hdr.nr; ineighr++)
        for (ineighsphtheta=0; ineighsphtheta<hdr.nsphtheta; ineighsphtheta++)
        for (ineighsphphi=0; ineighsphphi<hdr.nsphphi; ineighsphphi++) {
            //nneigh++;
            r=radial_table[ineighr];
            //sphtheta=(ineighsphtheta+0.5)*hdr.dsph;
            //sphphi=ineighsphphi*hdr.dsph;
            //sph_to_cart(r,sphtheta,sphphi,x2);
            idiff=abs(isphphi-ineighsphphi);
            if (idiff>(hdr.nsphphi/2)) idiff=hdr.nsphphi-idiff; //need to go short way around
            cosgamma=costheta[isphtheta]*costheta[ineighsphtheta]+sintheta[isphtheta]*sintheta[ineighsphtheta]*cosphi[idiff];
            bin=(int) (((cosgamma+1.0)/2.0)*TABLESIZE);
            if (bin<0) bin=0;
            if (bin>=TABLESIZE) bin=TABLESIZE-1;
            weight=spherical_diffusion_kernel[bin];
            if (weight<=MIN_WEIGHT) continue;
            neighborindex=calculate_index(ineighr,ineighsphtheta,ineighsphphi,ineighphi,ineightheta,ineighpsi);

            //if ((index==1) || (neighborindex==1)) printf("neighbors: %lld %lld\n",index,neighborindex);
            //theta=(ineightheta+0.5)*hdr.deuler;
            neighvolume=r*r*r*sintheta[ineighsphtheta]*v1;
            if (neighvolume==0.0) continue;
            //weight=gausstable[bin]; //exp(-dist/rho2);
            //weight=exp(-sin2gamma/sigma2);
            weightsum+=weight*neighvolume;
            sum+=neighvolume*weight*prob[neighborindex];
            //if ((index==1) || (neighborindex==1)) printf("neighbors: %lld %lld %.4f\n",index,neighborindex,weight);
#ifdef DEBUG
            if (index==15) {
                printf("***: %lld %lld %ld %ld %ld %ld %.4f %.4f %.4f %.4f %.4f %.4f %.4f\n",index,neighborindex,isphtheta,isphphi,ineighsphtheta,ineighsphphi,prob[index],prob[neighborindex],cosgamma,weight,sum,weightsum,sum/weightsum);
            }
#endif
        }
        //printf("totalweight: %lld %.10f %ld %ld %ld %ld %ld %ld$\n",index,weightsum,ir,isphtheta,isphphi,iphi,itheta,ipsi);
        //debugsum1+=sum;
        //w_i exp(-beta*F_i) = (1/n) sum_n w_j exp(-beta*F_j)
        //energy not changed for sum=0.
        sum=sum/(weightsum);
#ifdef DEBUG
        if ((sum>20.0) && (prob[index]<20.0)) printf("lowered: %lld\n",index);
#endif
        if (sum<0) {
            printf("Error in smoothing for indices %ld %ld %ld %ld %ld %ld\n",ir,isphtheta,isphphi,iphi,itheta,ipsi);
            die();
        } else if (sum>0) energy[index]=-betainv*log(sum);
        points++;
        //debugsum2+=sum*volume*weightsum;
        if ((points%10000)==0) {
                printf("%lld points averaged\n",points);
                fflush(stdout);
        }
    }
    //printf("debug sums: %.10f %.10f %.10f\n",debugsum1,debugsum2,volume);
    free(prob);
}

void table::boltzmann_average_orient(void)
{
    //long int ix[6],d[6],dstart[6],dstop[6],neigh[6];
    long int ir,isphtheta,isphphi,iphi,itheta,ipsi,ialpha,ibeta;
    long int ineighr,ineighsphtheta,ineighsphphi,ineighphi,ineightheta,ineighpsi;
    int i,nvalid,exist,nneigh,nsmoothrplus,nsmoothrminus,nsmoothsph,bin;
    double r,sphtheta,sphphi,phi,theta,psi;
    double sum,sum2,en,volume,neighvolume,v1,betainv,weight,weightsum,q1[4],q2[4],q1c[4],q3[4],q1a[4],cosrot;
    long long int index,neighborindex,transindex,points;
    double * prob;
    double * costheta, *sintheta, *cosphi, *sinphi;
    double so3_diffusion_kernel[TABLESIZE],work1[LMAX+1];
    printf("Creating SO(3) diffusion kernel...\n");
    create_so3_diffusion_kernel(TABLESIZE,LMAX,hdr.euler_alpha,so3_diffusion_kernel,work1);
    //double debugsum1, debugsum2;
    costheta=(double *) checkalloc(hdr.ntheta*4,sizeof(double));
    sintheta=(double *) checkalloc(hdr.ntheta*4,sizeof(double));
    for (i=0; i<hdr.ntheta*4; i++) {
        theta=i*hdr.deuler/4;
        costheta[i]=cos(theta);
        sintheta[i]=sin(theta);
    }
    cosphi=(double *) checkalloc(hdr.nphipsi*2,sizeof(double));
    sinphi=(double *) checkalloc(hdr.nphipsi*2,sizeof(double));
    for (i=0; i<hdr.nphipsi*2; i++) {
        phi=i*hdr.deuler/2;
        cosphi[i]=cos(phi);
        sinphi[i]=sin(phi);
    }
    points=0;
    prob=(double *) malloc(hdr.totalpoints*sizeof(double));
    if (prob==NULL) {
        printf("Could not allocate space for Boltzmann factors.\n");
        die();
    }
    printf("Computing Boltzmann factors...\n");
    for (index=0; index<hdr.totalpoints; index++) {
        if (energy[index]>50.0) prob[index]=0.0; else prob[index]=exp(-hdr.beta*energy[index]);
    }
    betainv=1/hdr.beta;
    v1=(hdr.rfactor-1)*hdr.dsph*hdr.dsph*hdr.deuler*hdr.deuler*hdr.deuler;
    //debugsum1=0.0;
    //debugsum2=0.0;
    for (iphi=0; iphi<hdr.nphipsi; iphi++)
    for (itheta=0; itheta<hdr.ntheta; itheta++)
    for (ipsi=0; ipsi<hdr.nphipsi; ipsi++)
    for (ir=0; ir<hdr.nr; ir++)
    for (isphtheta=0; isphtheta<hdr.nsphtheta; isphtheta++)
    for (isphphi=0; isphphi<hdr.nsphphi; isphphi++) {
        //printf("%ld %ld %ld %ld %ld %ld\n",ir,isphtheta,isphphi,iphi,itheta,ipsi);
        index=calculate_index(ir,isphtheta,isphphi,iphi,itheta,ipsi);
        //above the limit? do not smooth
        //if (energy[index]>hdr.en_smooth_limit) continue;
        r=radial_table[ir];
        //sphtheta=(isphtheta+0.5)*hdr.dsph;
        //sphphi=isphphi*hdr.dsph;
        //theta=(itheta+0.5)*hdr.deuler;
        //(itheta+0.5)*hdr.deuler = (4*itheta+2)*hdr.deuler/4
        volume=sintheta[4*itheta+2]*v1; //orientational part only
        nneigh=0;
        sum=0.0;
        weightsum=0.0;
        //euler_to_quat(phi,theta,psi,q1);
        ialpha=iphi-ipsi;
        if (ialpha<0) ialpha+=2*hdr.nphipsi;
        if (ialpha>=2*hdr.nphipsi) ialpha-=2*hdr.nphipsi;
        ibeta=iphi+ipsi;
        if (ibeta<0) ibeta+=2*hdr.nphipsi;
        if (ibeta>=2*hdr.nphipsi) ibeta-=2*hdr.nphipsi;
        //if ((ialpha<0) || (ialpha>=2*hdr.nphipsi) || (ibeta<0) || (ibeta>=2*hdr.nphipsi)) printf("***1\n");
        //(itheta+0.5)*hdr.deuler/2 = (2*itheta+1)*hdr.deuler/4
        q1[0]=cosphi[ibeta]*costheta[2*itheta+1];
        q1[1]=cosphi[ialpha]*sintheta[2*itheta+1];
        q1[2]=sinphi[ialpha]*sintheta[2*itheta+1];
        q1[3]=sinphi[ibeta]*costheta[2*itheta+1];
        //normalize_quat(q1);
#ifdef DEBUG
        phi=iphi*hdr.deuler;
        theta=(itheta+0.5)*hdr.deuler;
        psi=ipsi*hdr.deuler;
        euler_to_quat(phi,theta,psi,q1a);
        for (i=0; i<4; i++) if (fabs(q1[i]-q1a[i])>1e-6) {
            printf("q1: %.10f %.10f %.10f %.10f\n",q1[0],q1[1],q1[2],q1[3]);
            printf("q1a: %.10f %.10f %.10f %.10f\n",q1a[0],q1a[1],q1a[2],q1a[3]);
            die();
        }
#endif
        conjugate_quat(q1,q1c);
        ineighr=ir;
        ineighsphtheta=isphtheta;
        ineighsphphi=isphphi;
        for (ineighphi=0; ineighphi<hdr.nphipsi; ineighphi++)
        for (ineightheta=0; ineightheta<hdr.ntheta; ineightheta++)
        for (ineighpsi=0; ineighpsi<hdr.nphipsi; ineighpsi++) {
            //nneigh++;
            //if above limit, do not include in sum (maintains symmetry)
            //if (energy[neighborindex]>hdr.en_smooth_limit) continue;
            r=radial_table[ineighr];
            //sphtheta=(ineighsphtheta+0.5)*hdr.dsph;
            //sphphi=ineighsphphi*hdr.dsph;
            //sph_to_cart(r,sphtheta,sphphi,x2);
            neighvolume=sintheta[4*ineightheta+2]*v1;
            if (neighvolume==0.0) continue;
            ialpha=ineighphi-ineighpsi;
            if (ialpha<0) ialpha+=2*hdr.nphipsi;
            if (ialpha>=2*hdr.nphipsi) ialpha-=2*hdr.nphipsi;
            ibeta=ineighphi+ineighpsi;
            if (ibeta<0) ibeta+=2*hdr.nphipsi;
            if (ibeta>=2*hdr.nphipsi) ibeta-=2*hdr.nphipsi;
            if ((ialpha<0) || (ialpha>=2*hdr.nphipsi) || (ibeta<0) || (ibeta>=2*hdr.nphipsi)) printf("***2\n");
            q2[0]=cosphi[ibeta]*costheta[2*ineightheta+1];
            q2[1]=cosphi[ialpha]*sintheta[2*ineightheta+1];
            q2[2]=sinphi[ialpha]*sintheta[2*ineightheta+1];
            q2[3]=sinphi[ibeta]*costheta[2*ineightheta+1];
            //normalize_quat(q2); //unnecessary, normalizing q3 = q2 q1^-1 takes care of numerical errors
            //multiply_quat(q2,q1c,q3);  //We only need the real part of q3.
            q3[0]=q2[0]*q1c[0]-q2[1]*q1c[1]-q2[2]*q1c[2]-q2[3]*q1c[3];
            //normalize_quat(q3); //They should be normalized by construction, with only very small numerical errors.
            cosrot=2*q3[0]*q3[0]-1; //cosine of overall rotation angle
            //if ((cosrot<=-0.9999999) && (q3[1]<0)) continue; //Don't double count rotations by 180 degrees from the original orientaiton.
            bin=(int) (((cosrot+1.0)/2.0)*TABLESIZE);
            //bin=(int) (cosrot);
            if (bin<0) bin=0;
            if (bin>=TABLESIZE) bin=TABLESIZE-1;
            weight=so3_diffusion_kernel[bin];
            neighborindex=calculate_index(ineighr,ineighsphtheta,ineighsphphi,ineighphi,ineightheta,ineighpsi);
            if (weight<=MIN_WEIGHT) continue;
            //if ((index==1) || (neighborindex==1)) printf("neighbors: %lld %lld\n",index,neighborindex);
            //theta=(ineightheta+0.5)*hdr.deuler;

            //weight=gausstable[bin]; //exp(-dist/rho2);
#ifdef DEBUG
            //if ((index==20000) || (neighborindex==20000)) printf("***: %lld %lld %.10f %d %.10f\n",index,neighborindex,cosrot,bin,weight);
            /*if (cosrot<=-0.99999999) {
                printf("debugging pad\n");
            }*/
#endif
            weightsum+=weight*neighvolume;
            sum+=neighvolume*weight*prob[neighborindex];
            //if ((index==1) || (neighborindex==1)) printf("neighbors: %lld %lld %.4f\n",index,neighborindex,weight);
        }
        //printf("totalweight: %lld %.10f %ld %ld %ld %ld %ld %ld$\n",index,weightsum,ir,isphtheta,isphphi,iphi,itheta,ipsi);
        //debugsum1+=sum;
        //w_i exp(-beta*F_i) = (1/n) sum_n w_j exp(-beta*F_j)
        //energy not changed for sum=0.
        sum=sum/(weightsum);
        if (sum<0) {
            printf("Error in smoothing for indices %ld %ld %ld %ld %ld %ld\n",ir,isphtheta,isphphi,iphi,itheta,ipsi);
            die();
        } else if (sum>0) energy[index]=-betainv*log(sum);
        points++;
        //debugsum2+=sum*volume*weightsum;
        if ((points%10000)==0) {
                printf("%lld points averaged\n",points);
                fflush(stdout);
        }
    }
    //printf("debug sums: %.10f %.10f %.10f\n",debugsum1,debugsum2,volume);
    free(prob);
}

double table::assoc_free_energy(double beta)
{
    long int ir,isphtheta,isphphi,iphi,itheta,ipsi;
    int i;
    double r,sphtheta,sphphi,phi,theta,psi;
    double sum,sum2,en,volume,v1,betainv,totalvolume;
    long long int index;
    betainv=1/beta;
    v1=(hdr.rfactor-1)*hdr.dsph*hdr.dsph*hdr.deuler*hdr.deuler*hdr.deuler;
    //totalvolume=32.0*M_PI*M_PI*M_PI*hdr.maxr*hdr.maxr*hdr.maxr/3.0;
    sum=0.0;
    totalvolume=0.0;
    for (iphi=0; iphi<hdr.nphipsi; iphi++)
    for (itheta=0; itheta<hdr.ntheta; itheta++)
    for (ipsi=0; ipsi<hdr.nphipsi; ipsi++)
    for (ir=0; ir<hdr.nr; ir++)
    for (isphtheta=0; isphtheta<hdr.nsphtheta; isphtheta++)
    for (isphphi=0; isphphi<hdr.nsphphi; isphphi++) {
        //printf("%d %d %d %ld %ld %ld\n",ir,isphtheta,isphphi,iphi,itheta,ipsi);
        index=calculate_index(ir,isphtheta,isphphi,iphi,itheta,ipsi);
        r=hdr.minr*pow(hdr.rfactor,((double)ir +0.5));
        sphtheta=(isphtheta+0.5)*hdr.dsph;
        theta=(itheta+0.5)*hdr.deuler;
        volume=r*r*r*sin(sphtheta)*sin(theta)*v1;
        sum+=volume*exp(-beta*energy[index]);
      	totalvolume+=volume;
    }
    //Sum is the partition function Z for interaction between the two fragments.
    return -betainv*log(sum/totalvolume);
}



void table::generate_table(const char * control_file, int part, int numparts)
{
    long int ir,isphtheta,isphphi,iphi,itheta,ipsi;
    long long int index,points,clash_table_points,index2;
    int i;
    double r,sphtheta,sphphi,x[3],phi,theta,psi,en,enmin,rotmatrix[3][3];
    double * energy2;
    double sum,phi2,theta2,psi2;
    double * coords1;
    char pmffname[255];
    FILE * test_output;
    go_model_info * go_model;
    fragmenttype * frag1;
    fragmenttype * frag2;
    //fragment frag1, frag2;
    printf("Generating table -- part %d of %d.\n",part,numparts);
    read_table_header_info(control_file);
    go_model = new go_model_info();
    go_model->read_file(hdr.go_model_map_fname);
    go_model->set_parameters(&hdr.go_params);

    hdr.go_model_entries=go_model->nentries;
    hdr.go_model_native=go_model->nnative;
    //Temporarily we do not use fragment names when constructing tables.  Maybe we should record the fragment names in the header file?
    frag1=new fragmenttype(go_model->ifragname,hdr.ref_frag_fname);
    frag2=new fragmenttype(go_model->jfragname,hdr.other_frag_fname);
    //set up multi-parts for incomplete tables
    hdr.startindex=(long long int) (hdr.totalpoints*((double) part)/((double) numparts));
    hdr.endindex=(long long int) (hdr.totalpoints*((double) (part+1))/((double) numparts))-1;
    if ((hdr.startindex!=0) || (hdr.endindex!=(hdr.totalpoints-1))) hdr.options|=OP_INCOMPLETE;
    //Check to make sure the reference and "other" fragments match the first and second in the go model.
    //As a matter of policy we require that the first fragment in the go model be the reference fragment, and the second fragment type be the other fragment.
    //this is useless, we don't specify fragment names here
    //if ((strcasecmp(frag1->fragname,go_model->ifragname)!=0) || (strcasecmp(frag2->fragname,go_model->jfragname)!=0)) {
        //printf("Table generation fragment types
    /*clash_table = (double *) malloc(hdr.clashpoints*sizeof(double));
    if (clash_table == NULL) {
        printf("Could not allocate memory for clash table\n");
        die();
    }
    exact=TRUE;
    printf("Filling clash table with %lld points.\n",hdr.clashpoints);*/
    //fill_clash_table();
    print_header_info();
    if ((hdr.options & OP_INCOMPLETE) && ((hdr.sph_alpha>0) || (hdr.euler_alpha>0))) {
        printf("Cannot smooth an incomplete table.\n");
        die();
    }
    /*pbc=FALSE;
    cutoff2=100000000000000.0;*/
    energy = (energy_t *) malloc(hdr.totalpoints*sizeof(energy_t));
    if (energy == NULL) {
        printf("Could not allocate memory for table\n");
        die();
    }
    for (index=0; index<hdr.totalpoints; index++) energy[index]=NAN;
    fill_table(go_model,frag1,frag2);
    //Only compute free energy if no invalid points.
    //if (hdr.invalidpoints==0) {
    if (!(hdr.options & OP_INCOMPLETE)) printf("Association free energy at temperature %.2f K: %.6f kcal/mol\n",1/(KBOLTZ*hdr.beta),assoc_free_energy(hdr.beta));
    //snprintf(pmffname,sizeof(pmffname),"pmf-%s",
    //calculate_dist_pmf(hdr.beta,"pmf.dat"); //this is less than ideal
    //the only other place to put it is in write_table

    if (hdr.sph_alpha>0) {
        boltzmann_average_trans();
        printf("Association free energy at temperature %.2f K: %.6f kcal/mol\n",1/(KBOLTZ*hdr.beta),assoc_free_energy(hdr.beta));
    }
    if (hdr.euler_alpha>0) {
        boltzmann_average_orient();
        printf("Association free energy at temperature %.2f K: %.6f kcal/mol\n",1/(KBOLTZ*hdr.beta),assoc_free_energy(hdr.beta));
    }
    delete go_model;
    delete frag1;
    delete frag2;
    hdr.crc32_data=digital_crc32((unsigned char *) energy,hdr.totalpoints*sizeof(energy_t));
}

void table::do_smooth(double smooth_temp, double trans_scale, double orient_scale)
{
    int ir;
    energy_t * energy2;
    hdr.sph_alpha=trans_scale;
    hdr.sph_alpha=hdr.sph_alpha*DEG_TO_RAD;
    hdr.sph_alpha=hdr.sph_alpha/SPHERICAL_DIFFUSION_CONST;
    hdr.sph_alpha=hdr.sph_alpha*hdr.sph_alpha;
    hdr.euler_alpha=orient_scale;
    hdr.euler_alpha=hdr.euler_alpha*DEG_TO_RAD;
    hdr.euler_alpha=hdr.euler_alpha/SO3_DIFFUSION_CONST;
    hdr.euler_alpha=hdr.euler_alpha*hdr.euler_alpha;
    hdr.beta=smooth_temp;
    if (hdr.beta<=0.0) hdr.beta=300.0; /*default 300 K for association free energy*/
    hdr.beta=1/(KBOLTZ*hdr.beta);
    time(&hdr.created); /*current date/time*/
    if ((hdr.sph_alpha>0) || (hdr.euler_alpha>0)) hdr.options = hdr.options | OP_BOLTZMANN_AVERAGED;
    if ((hdr.options & OP_BOLTZMANN_AVERAGED) && (hdr.en_invalid_margin>0.0)) { //sanity check
        printf("You may not apply smoothing and mark invalid cells at the same time.\n");
        die();
    }
#ifndef NO_MMAP_TABLES
    //We need to get rid of the memory mapping on the table data so we can write to it without modifying the original table.
    energy2 = (energy_t *) checkalloc(hdr.totalpoints,sizeof(energy_t));
    memcpy(energy2, energy, hdr.totalpoints*sizeof(energy_t));
    munmap(map,mapsize);
    energy=energy2;
#endif
    radial_table =(double *) checkalloc(hdr.nr,sizeof(double));//destroyed by destructor
    for (ir=0; ir<hdr.nr; ir++) radial_table[ir]=hdr.minr*pow(hdr.rfactor,ir);
    printf("-------------------------------------------------------------------------\n");
    printf("Smoothing table.\n");
    print_header_info();
    printf("Association free energy at temperature %.2f K: %.6f kcal/mol\n",1/(KBOLTZ*hdr.beta),assoc_free_energy(hdr.beta));
    if (hdr.sph_alpha>0) {
        boltzmann_average_trans();
        printf("Association free energy at temperature %.2f K: %.6f kcal/mol\n",1/(KBOLTZ*hdr.beta),assoc_free_energy(hdr.beta));
    }
    if (hdr.euler_alpha>0) {
        boltzmann_average_orient();
        printf("Association free energy at temperature %.2f K: %.6f kcal/mol\n",1/(KBOLTZ*hdr.beta),assoc_free_energy(hdr.beta));
    }
    hdr.crc32_data=digital_crc32((unsigned char *) energy,hdr.totalpoints*sizeof(energy_t));
}


void table::write_dx(char * fname, double phi, double theta, double psi, double xmax, double dx)
{
    double  r, sphtheta, sphphi, x[3], en;
    int nx,i,j,k;
    long int points;
    FILE * output;
    output=fopen(fname,"w");
    nx=((int) ceil(2*xmax/dx))+1;
    dx=2*xmax/(nx-1);
    /*DX file header*/
    fprintf(output,"# Produced from table with phi, theta, phi = %.2f %.2f %.2f\n",phi,theta,psi);
    phi=phi*DEG_TO_RAD;
    theta=theta*DEG_TO_RAD;
    psi=psi*DEG_TO_RAD;
    fprintf(output,"object 1 class gridpositions counts %d %d %d\n",nx,nx,nx);
    fprintf(output,"origin %.4f %.4f %.4f\n",-xmax,-xmax,-xmax);
    fprintf(output,"delta %.4f %.4f %.4f\n",dx,0.0,0.0);
    fprintf(output,"delta %.4f %.4f %.4f\n",0.0,dx,0.0);
    fprintf(output,"delta %.4f %.4f %.4f\n",0.0,0.0,dx);
    fprintf(output,"object 2 class gridconnections counts %d %d %d\n",nx,nx,nx);
    fprintf(output,"object 3 class array type double rank 0 items :\n");
    printf("total %ld points\n",(long int) nx*nx*nx);
    points=1;
    for (i=0; i<nx; i++)
        for (j=0; j<nx; j++)
            for (k=0; k<nx; k++) {
                x[0]=-xmax+i*dx+1e-6; //cart_to_sph does not cover the case where x=0 exactly, for efficiency reasons
                x[1]=-xmax+j*dx;
                x[2]=-xmax+k*dx;
                r=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);

                if (r>hdr.maxr) en=0;
                else if (r<hdr.minr) en=DUMMY_ENERGY;
                else {
                    cart_to_sph(&x[0],r,&sphtheta,&sphphi);
                    en=get_energy(FALSE,r,sphtheta,sphphi,phi,theta,psi);
                }
                fprintf(output,"%20.10e",en);
                if ((points%3)==0) fprintf(output,"\n");
                if ((points%1000000)==0) {
                    printf("%ld points written\n",points);
                    fflush(output);
                }
                points++;
            }
    points--;
    if ((points%3)!=0) fprintf(output,"\n");
    fprintf(output,"attribute \"dep\" string \"positions\"\n");
    fprintf(output,"object \"regular positions regular connections\"\n");
    fprintf(output,"component \"positions\" value 1\n");
    fprintf(output,"component \"connections\" value 2\n");
    fprintf(output,"component \"data\" value 3\n");
    fclose(output);
}


void table::write_dx_exact(char * fname, double phi, double theta, double psi, double xmax, double dx)
{
    double  r, sphtheta, sphphi, x[3], q1[4], q2[4], en;
    int nx,i,j,k;
    long int points;
    go_model_info * go_model;
    fragmenttype * frag1;
    fragmenttype * frag2;
    double * coords1;
    double * coords2;
    FILE * output;
    go_model = new go_model_info();
    go_model->read_file(hdr.go_model_map_fname);
    go_model->set_parameters(&hdr.go_params);
    frag1=new fragmenttype("\0",hdr.ref_frag_fname);
    frag2=new fragmenttype("\0",hdr.other_frag_fname);
    coords1=(double *) malloc(3*frag1->natom*sizeof(double));
    coords2=(double *) malloc(3*frag2->natom*sizeof(double));
    output=fopen(fname,"w");
    nx=((int) ceil(2*xmax/dx))+1;
    dx=2*xmax/(nx-1);
    ///DX file header
    fprintf(output,"# Produced from table with phi, theta, phi = %.2f %.2f %.2f\n",phi,theta,psi);
    phi=phi*DEG_TO_RAD;
    theta=theta*DEG_TO_RAD;
    psi=psi*DEG_TO_RAD;
    fprintf(output,"object 1 class gridpositions counts %d %d %d\n",nx,nx,nx);
    fprintf(output,"origin %.4f %.4f %.4f\n",-xmax,-xmax,-xmax);
    fprintf(output,"delta %.4f %.4f %.4f\n",dx,0.0,0.0);
    fprintf(output,"delta %.4f %.4f %.4f\n",0.0,dx,0.0);
    fprintf(output,"delta %.4f %.4f %.4f\n",0.0,0.0,dx);
    fprintf(output,"object 2 class gridconnections counts %d %d %d\n",nx,nx,nx);
    fprintf(output,"object 3 class array type double rank 0 items :\n");
    printf("total %ld points\n",(long int) nx*nx*nx);
    points=1;
    q1[0]=1.0;
    q1[1]=0.0;
    q1[2]=0.0;
    q1[3]=0.0;
    x[0]=0.0;
    x[1]=0.0;
    x[2]=0.0;
    frag1->get_coords(&x[0],&q1[0],coords1);
    for (i=0; i<nx; i++)
        for (j=0; j<nx; j++)
            for (k=0; k<nx; k++) {
                x[0]=-xmax+i*dx+1e-6; //cart_to_sph does not cover the case where x=0 exactly, for efficiency reasons
                x[1]=-xmax+j*dx;
                x[2]=-xmax+k*dx;
                r=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
                cart_to_sph(&x[0],r,&sphtheta,&sphphi);
                en=get_pair_energy(go_model,frag1,frag2,coords1,coords2,r,sphtheta,sphphi,phi,theta,psi);
                fprintf(output,"%20.10e",en);
                if ((points%3)==0) fprintf(output,"\n");
                if ((points%1000)==0) {
                    printf("%ld points written\n",points);
                    fflush(output);
                }
                points++;
            }
    points--;
    if ((points%3)!=0) fprintf(output,"\n");
    fprintf(output,"attribute \"dep\" string \"positions\"\n");
    fprintf(output,"object \"regular positions regular connections\"\n");
    fprintf(output,"component \"positions\" value 1\n");
    fprintf(output,"component \"connections\" value 2\n");
    fprintf(output,"component \"data\" value 3\n");
    fclose(output);
    free(coords1);
    free(coords2);
    delete frag1;
    delete frag2;
}

//This visualizes in orientation space.  It is a 3-dimensional ball of radius 180 degrees. A point within the sphere represents
//a rotation along the axis connecting the point to the origin, by an angle equal to the distance of the point from the origin.
void table::write_dx_orient(char * fname, double r, double sphtheta, double sphphi, double dx)
{
    double phi,theta,psi, x[3], axis[3],q[4],angle,en, xmax;
    int nx,i,j,k,l;
    long int points;
    FILE * output;
    output=fopen(fname,"w");
    sphtheta=sphtheta*DEG_TO_RAD;
    sphphi=sphphi*DEG_TO_RAD;
    xmax=180.0;
    nx=((int) ceil(2*xmax/dx))+1;
    dx=2*xmax/(nx-1);
    /*DX file header*/
    fprintf(output,"# Orientational. Produced from table with r, theta, phi = %.2f %.2f %.2f\n",r,sphtheta,sphphi);
    fprintf(output,"object 1 class gridpositions counts %d %d %d\n",nx,nx,nx);
    fprintf(output,"origin %.4f %.4f %.4f\n",-xmax,-xmax,-xmax);
    fprintf(output,"delta %.4f %.4f %.4f\n",dx,0.0,0.0);
    fprintf(output,"delta %.4f %.4f %.4f\n",0.0,dx,0.0);
    fprintf(output,"delta %.4f %.4f %.4f\n",0.0,0.0,dx);
    fprintf(output,"object 2 class gridconnections counts %d %d %d\n",nx,nx,nx);
    fprintf(output,"object 3 class array type double rank 0 items :\n");
    printf("total %ld points\n",(long int) nx*nx*nx);
    points=1;
    for (i=0; i<nx; i++)
        for (j=0; j<nx; j++)
            for (k=0; k<nx; k++) {
                //the entries of x[] are in degrees
                x[0]=-xmax+i*dx+1e-6;
                x[1]=-xmax+j*dx;
                x[2]=-xmax+k*dx;
                angle=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
                if (angle>180.0) {
                    en=0;
                } else {
                    for (l=0; l<3; l++) axis[l]=x[l]/angle;
                    angle=angle*DEG_TO_RAD;
                    axisangle_to_quat(angle,axis,q);
                    quat_to_euler(q,&phi,&theta,&psi);
                    en=get_energy(FALSE,r,sphtheta,sphphi,phi,theta,psi);
                }
                fprintf(output,"%20.10e",en);
                if ((points%3)==0) fprintf(output,"\n");
                if ((points%100000)==0) {
                    printf("%ld points written\n",points);
                    fflush(output);
                }
                points++;
            }
    points--;
    if ((points%3)!=0) fprintf(output,"\n");
    fprintf(output,"attribute \"dep\" string \"positions\"\n");
    fprintf(output,"object \"regular positions regular connections\"\n");
    fprintf(output,"component \"positions\" value 1\n");
    fprintf(output,"component \"connections\" value 2\n");
    fprintf(output,"component \"data\" value 3\n");
    fclose(output);
}

//Calculate the pmf in terms of r directly from the table by integration.
void table::calculate_dist_pmf(double temp, char * fname)
{
    long int ir,isphtheta,isphphi,iphi,itheta,ipsi;
    int i;
    double r,sphtheta,sphphi,phi,theta,psi;
    double sum,sum2,en,volume,v1,betainv,beta,totalvolume;
    long long int index;
    FILE * output;
    output=fopen(fname,"w");
    if (output==NULL) {
        printf("Could not open file %s\n",fname);
        die();
    }
    betainv=KBOLTZ*temp;
    beta=1/betainv;
    v1=(hdr.rfactor-1)*hdr.dsph*hdr.dsph*hdr.deuler*hdr.deuler*hdr.deuler;
    for (ir=0; ir<hdr.nr; ir++) {
    //totalvolume=32.0*M_PI*M_PI*M_PI*hdr.maxr*hdr.maxr*hdr.maxr/3.0;
        sum=0.0;
        totalvolume=0.0;
        r=hdr.minr*pow(hdr.rfactor,((double)ir+0.5));
        for (iphi=0; iphi<hdr.nphipsi; iphi++)
        for (itheta=0; itheta<hdr.ntheta; itheta++)
        for (ipsi=0; ipsi<hdr.nphipsi; ipsi++)
        for (isphtheta=0; isphtheta<hdr.nsphtheta; isphtheta++)
        for (isphphi=0; isphphi<hdr.nsphphi; isphphi++) {
            //printf("%d %d %d %ld %ld %ld\n",ir,isphtheta,isphphi,iphi,itheta,ipsi);
            index=calculate_index(ir,isphtheta,isphphi,iphi,itheta,ipsi);
            sphtheta=(isphtheta+0.5)*hdr.dsph;
            theta=(itheta+0.5)*hdr.deuler;
            volume=r*r*r*sin(sphtheta)*sin(theta)*v1;
            sum+=volume*exp(-beta*energy[index]);
            totalvolume+=volume;
       }
       if ((totalvolume>0) && (sum>0)) fprintf(output,"%.4f %.6f\n",r,-betainv*log(sum/totalvolume));
    }
    fclose(output);
    printf("Association free energy at temperature %.2f K: %.6f kcal/mol\n",1/(KBOLTZ*beta),assoc_free_energy(beta));
    //Sum is the partition function Z for interaction between the two fragments.
    //return -betainv*log(sum/totalvolume);
}


//Special "decimation" constructor for constructing a smaller table from a larger one.
//decr, etc. -- decimation factors.
table::table(table * orig_table, int decr, int decsph, int deceuler)
{
    long long int this_index, orig_index;
    long int ir, isphtheta, isphphi, iphi, itheta, ipsi;
    printf("Original table:\n");
    orig_table->print_header_info();
    hdr=orig_table->hdr;
    hdr.options |= OP_DECIMATED;
    hdr.decr=decr;
    hdr.decsph=decsph;
    hdr.deceuler=deceuler;
    hdr.dr*=decr;
    hdr.dsph*=decsph;
    hdr.deuler*=deceuler;
    time(&hdr.created); /*current date/time*/
    hdr.nr/=decr;
    hdr.nsphtheta/=decsph;
    hdr.nsphphi/=decsph;
    hdr.nphipsi/=deceuler;
    hdr.ntheta/=deceuler;
    //the following are repeated from read_table_header_info
    hdr.norient=hdr.nphipsi * hdr.nphipsi * hdr.ntheta;
    //hdr.clashpoints=hdr.nsphtheta*hdr.nsphphi*hdr.nphipsi*hdr.ntheta*hdr.nphipsi;
    hdr.rfactor=1+hdr.dr/hdr.minr;
    hdr.logrfactor=log(hdr.rfactor);
    hdr.nr=(long int) (log(hdr.maxr/hdr.minr)/hdr.logrfactor+0.5)+1;
    hdr.ntrans=hdr.nr * hdr.nsphtheta * hdr.nsphphi;
    hdr.totalpoints=(long long int) hdr.ntrans * (long long int) hdr.norient;
    hdr.clashpoints=0;
    hdr.invalidpoints=0;
    //now allocate space for the new table and fill it with values from the old table
    printf("\n");
    printf("--------------------------------------------------------------------------------\n");
    printf("New table:\n");
    print_header_info();
    energy = (energy_t *) malloc(hdr.totalpoints*sizeof(energy_t));
    if (energy == NULL) {
        printf("Could not allocate memory for table\n");
        die();
    }
    for (iphi=0; iphi<hdr.nphipsi; iphi++)
    for (itheta=0; itheta<hdr.ntheta; itheta++)
    for (ipsi=0; ipsi<hdr.nphipsi; ipsi++)
    for (ir=0; ir<hdr.nr; ir++)
    for (isphtheta=0; isphtheta<hdr.nsphtheta; isphtheta++)
    for (isphphi=0; isphphi<hdr.nsphphi; isphphi++) {
        this_index=calculate_index(ir,isphtheta,isphphi,iphi,itheta,ipsi);
        orig_index=orig_table->calculate_index(ir*decr,isphtheta*decsph,isphphi*decsph,iphi*deceuler,itheta*deceuler,ipsi*deceuler);
        energy[this_index]=orig_table->energy[orig_index];
        if (energy[this_index]==INVALID_ENERGY) hdr.invalidpoints++;
    }
    hdr.crc32_data=digital_crc32((unsigned char *) energy,hdr.totalpoints*sizeof(energy_t));
}

//Combine tables.
table::table(int ntables, table * * tables)
{
    int itable;
    long long int index,start,end;
    hdr=tables[0]->hdr;
    hdr.options &= ~OP_INCOMPLETE;
    hdr.startindex=0;
    hdr.endindex=hdr.totalpoints-1;
    printf("\n");
    //printf("--------------------------------------------------------------------------------\n");
    printf("New table:\n");
    print_header_info();
    energy = (energy_t *) malloc(hdr.totalpoints*sizeof(energy_t));
    if (energy == NULL) {
        printf("Could not allocate memory for table\n");
        die();
    }
    for (index=0; index<hdr.totalpoints; index++) energy[index]=NAN;
    for (index=0; index<hdr.totalpoints; index++) {
        //Find the table that does not have a NAN in its entry.
        for (itable=0; itable<ntables; itable++) {
            start=tables[itable]->hdr.startindex;
            end=tables[itable]->hdr.endindex;
            if ((index>=start) && (index<=end)) energy[index]=tables[itable]->energy[index-start]; //account for shift when loading
        }
    }
    //check to see if there are any nan's left
    for (index=0; index<hdr.totalpoints; index++) if (isnan(energy[index])) {
        printf("Error: table still incomplete after combination.\n");
        printf("first empty index %lld\n",index);
        die();
    }
    hdr.crc32_data=digital_crc32((unsigned char *) energy,hdr.totalpoints*sizeof(energy_t));
}
