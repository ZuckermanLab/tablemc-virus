#ifdef EXCHANGE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "mc.h"
#include "util.h"
#include "mt.h"
#include "fragments.h"


//Assumes that we are in parallel mode. Also assumes (for now) replicas == nodes.
//need to take care of splitting the trajectory files.
void simulation::exchange_init(FILE * input) //, char * fragfmt)
{
    char tablefmt[255], rexfname[255], xyzfmt[255], quatfmt[255], restfmt[255];
    char * tablefmts;
    double * temps;
    double * lambdas;
    double t,l,total_table_size;
    int irep,iirep,i;
    myrep=mynod;
    nrep=numnod;
    //read all data.
    //tablefmts=(char *) checkalloc(nrep,sizeof(tablefmt));
    betas=(double *) checkalloc(nrep,sizeof(double));
    //lambdas=(double *) checkalloc(nrep,sizeof(double));
    //All nodes read the file.
    //fscanf(input,"%s %s %s\n",xyzfmt,quatfmt,restfmt);
    strncpy(xyzfmt,xyzfname,sizeof(xyzfmt));
    snprintf(xyzfname,sizeof(xyzfname),xyzfmt,myrep+1);
    trim_string(xyzfname);
    strncpy(quatfmt,quatfname,sizeof(quatfmt));
    snprintf(quatfname,sizeof(quatfname),quatfmt,myrep+1);
    trim_string(quatfname);
    strncpy(restfmt,endrestartfname,sizeof(restfmt));
    snprintf(endrestartfname,sizeof(quatfname),restfmt,myrep+1);
    trim_string(endrestartfname);
    fscanf(input,"%d %s\n",&exchfreq,rexfname);
    printf("Enabling replica exchange for replica %d of %d\n",myrep+1,nrep);
    printf("Exchange frequency:           %d\n",exchfreq);
    printf("Exchange long file name:      %s\n",rexfname);
    fflush(stdout);
    for (iirep=0; iirep<nrep; iirep++) {
        if (feof(input)) {
            printf("Error reading replica exchange section of input file.\n");
            die();
        }
        fscanf(input,"%d  %lg\n",&irep,&t);
	if (irep>nrep) {
            printf("Incorrect replica number in input file.\n");
            die();
        }
        irep--; //1-based replica number
        betas[irep]=1/(KBOLTZ*t);
        printf("Replica %d will use and temperature %.2f K\n",irep+1,1/(KBOLTZ*betas[irep]));
        fflush(stdout);
    }
    //broadcast all data
    /*MPI_Bcast(&exchfreq,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(rexfname,sizeof(rexfname),MPI_CHAR,0,MPI_COMM_WORLD);
    MPI_Bcast(tablefmts,nrep*sizeof(tablefmt),MPI_CHAR,0,MPI_COMM_WORLD);
    MPI_Bcast(betas,nrep,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(lambdas,nrep,MPI_DOUBLE,0,MPI_COMM_WORLD);*/
    //Open teh log file across all replicas.
    MPI_File_delete(rexfname,MPI_INFO_NULL);
    MPI_File_open(MPI_COMM_WORLD,rexfname,MPI_MODE_CREATE | MPI_MODE_WRONLY,MPI_INFO_NULL,&rexlog);
    //printf("lambda1 %.2f\n",lambdas[myrep]);
    /*tables_lambda=lambdas[myrep];
    //printf("lambda2 %.2f %.2f %d\n",tables_lambda,lambdas[myrep],myrep);
    if (tables_lambda<0.0) tables_lambda=0.0;
    if (tables_lambda>1.0) tables_lambda=1.0;*/
    //0 = fully exact, 1 = fully tabulated.
    /*use_std_tables=(tables_lambda>0.0);
    if ((tables_lambda>0.0) && (tables_lambda<1.0)) {
        printf("Will mix exact/tabulated energies in this simulation.  Fraction %.2f exact and %.2f tabulated.\n",(1.0-tables_lambda),tables_lambda);
    } else if (tables_lambda==0.0) {
        printf("Will calculate all nonbonded interactions exactly.\n");
    } else if (tables_lambda==1.0) {
        printf("Will use noncovalent tables in this simulation.\n");
    }
    use_cov_tables=false;
    printf("Will calculate all peptide covalent interactions exactly.\n");*/
    beta=betas[myrep];
    printf("Temperature set to %.2f K.\n",1/(KBOLTZ*beta));
    /*strncpy(tablefmt,&tablefmts[myrep*sizeof(tablefmt)],sizeof(tablefmt));
    trim_string(tablefmt);
    if (use_std_tables || enwrite) {
        tables = (table * *) checkalloc(nfragtypes*nfragtypes,sizeof(table *));
        for (i=0; i<nfragtypes*nfragtypes; i++) tables[i]=NULL;
        load_tables(tablefmt,fragfmt,tables);
    } //else tablefmt[0]='\0';
    total_table_size=0.0;
    if (use_std_tables || enwrite) for (i=0; i<nfragtypes*nfragtypes; i++) if (tables[i]!=NULL) total_table_size+=tables[i]->getsize();
    if (use_std_tables || use_cov_tables || enwrite) printf("Total table size:     %.2f MB\n",total_table_size);
    free(lambdas);
    free(tablefmts);*/
}

#define NMSG 5

void simulation::exchange(int icycle, double * cum_energy, double * center, double * orient, double * coords)
{
    int rep_swap,ifrag,partner,i; //rep_swap and rep_swap+1 will exchange
    double * exchcenter;
    double * exchorient;
    double * exchcoords;
    double eii,eij,ejj,eji,db,p,r;
    bool accept,new_nb_list;
    bool * dummy_moved;
    //MPI_Status status;
    //MPI_Request req;
    char buffer[255];
    //printf("node %d entering exchange\n",mynod);
    //Node 0 chooses the replica to swap and broadcasts.
    if (mynod==0) rep_swap=int((nrep-1)*genrand_real3());
    MPI_Bcast(&rep_swap,1,MPI_INT,0,MPI_COMM_WORLD);
    printf("Attempting to swap replicas %d and %d\n",rep_swap+1,rep_swap+2);
    if ((myrep!=rep_swap) && (myrep!=(rep_swap+1))) return; //we're done
    exchcenter=(double *) checkalloc(3*nfrag,sizeof(double));
    exchorient=(double *) checkalloc(4*nfrag,sizeof(double));
    exchcoords=(double *) checkalloc(3*natom,sizeof(double));
    dummy_moved=(bool *) checkalloc(nfrag,sizeof(bool));
    for (ifrag=0; ifrag<nfrag; ifrag++) dummy_moved[ifrag]=true;

    //Exchange energies, centers, and orientations.  Receive
    if (myrep==rep_swap) partner=rep_swap+1;
    if (myrep==(rep_swap+1)) partner=rep_swap;
    //printf("node %d partner %d\n",myrep,partner);
    //fflush(stdout);
    eii=*cum_energy;
    //printf("node %d about to excahnge with node %d\n",mynod,partner);
    //fflush(stdout);
    MPI_Sendrecv(&eii,1,MPI_DOUBLE,partner,NMSG*icycle+1,&ejj,1,MPI_DOUBLE,partner,NMSG*icycle+1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    //printf("node %d exchanged energies with node %d\n",mynod,partner);
    //fflush(stdout);
    MPI_Sendrecv(center,3*nfrag,MPI_DOUBLE,partner,NMSG*icycle+2,exchcenter,3*nfrag,MPI_DOUBLE,partner,NMSG*icycle+2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    //printf("node %d exchanged centers with node %d\n",mynod,partner);
    //fflush(stdout);
    MPI_Sendrecv(orient,4*nfrag,MPI_DOUBLE,partner,NMSG*icycle+3,exchorient,4*nfrag,MPI_DOUBLE,partner,NMSG*icycle+3,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    //printf("node %d exchanged orients with node %d\n",mynod,partner);
    //fflush(stdout);
    //compute cross energies
    if (exact) for (ifrag=0; ifrag<nfrag; ifrag++) update_coords(ifrag,exchcenter,exchorient,exchcoords);
         if (use_nb_list) {
#ifdef TIMERS
             switch_timer(TIMER_NB_LIST);
#endif

             new_nb_list=frag_nblist->check_nb_list(pbc,halfboxsize,boxsize,cutoff,nfrag,dummy_moved,exchcenter);
            if (new_nb_list) {
                //printf("New nonbond list at step %ld\n",istep);
                //if (pbc) recenter();
                frag_nblist->create_nonbond_list(pbc,halfboxsize,boxsize,nfrag,newcenter);
                //create_nb_atom_exact_list(exact,nb_list_per_frag,nb_list_count,nonbond_list,&nb_atom_list);
            }
#ifdef TIMERS
             switch_timer(TIMER_OTHER);
#endif
         }
    eji=total_energy(exchcenter,exchorient,exchcoords);
    //exchange these
    printf("Cross energy: %.4f kcal/mol\n",eji);
    //print_energies(stdout,FALSE,"Cross: ",icycle,energies,eji);
    //printf("node %d about to exchange energies with node %d\n",mynod,partner);
    //fflush(stdout);
    MPI_Sendrecv(&eji,1,MPI_DOUBLE,partner,NMSG*icycle+4,&eij,1,MPI_DOUBLE,partner,NMSG*icycle+4,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    //printf("node %d phase 1 complete\n",mynod);
    //fflush(stdout);
    //determine whether or not to swap -- detailed balance condition.  eii and eji on hamiltonian i, eij and ejj on hamiltonian j
    if (myrep==rep_swap) {
        db=-betas[myrep]*(eji-eii)-betas[partner]*(eij-ejj);
        //printf("%.6f %.6f %.6f %.6f %.6f %.6f %.6f\n",eii,ejj,eij,eji,betas[myrep],betas[partner],db);
        if (db>0) p=1; else p=exp(db);
        r=genrand_real3();
        accept=(r<p);
        snprintf(buffer,sizeof(buffer),"%d %d %d %.6f %.6f %.6f %.6f %.6f %c\n",icycle,rep_swap+1,rep_swap+2,eii,ejj,eij,eji,p,yesno(accept));
        MPI_File_write_shared(rexlog,buffer,strlen(buffer),MPI_CHAR,MPI_STATUS_IGNORE);
        MPI_Send(&accept,1,MPI_INT,partner,NMSG*icycle+5,MPI_COMM_WORLD);
    } else { //must be rep_swap+1
        MPI_Recv(&accept,1,MPI_INT,partner,NMSG*icycle+5,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }
    //printf("node %d phase 2 complete, accept=%c\n",mynod,yesno(accept));
    //fflush(stdout);
    if (accept) {
        //finish the swap -- since mcloop assumes that "old" and "new" coordinates are the same most of the time, we need to update both of them.
        //for (i=0; i<EN_TERMS; i++) cum_energies[i]=energies[i];
        *cum_energy=eji;
        for (ifrag=0; ifrag<nfrag; ifrag++) copy_frag(ifrag,exchcenter,exchorient,exchcoords,newcenter,neworient,newcoords);
        for (ifrag=0; ifrag<nfrag; ifrag++) copy_frag(ifrag,exchcenter,exchorient,exchcoords,oldcenter,oldorient,oldcoords);
         if (use_nb_list) {
#ifdef TIMERS
             switch_timer(TIMER_NB_LIST);
#endif

             new_nb_list=frag_nblist->check_nb_list(pbc,halfboxsize,boxsize,cutoff,nfrag,dummy_moved,newcenter);
            if (new_nb_list) {
                //printf("New nonbond list at step %ld\n",istep);
                //if (pbc) recenter();
                frag_nblist->create_nonbond_list(pbc,halfboxsize,boxsize,nfrag,newcenter);
                //create_nb_atom_exact_list(exact,nb_list_per_frag,nb_list_count,nonbond_list,&nb_atom_list);
            }
#ifdef TIMERS
             switch_timer(TIMER_OTHER);
#endif
         }
    }
    fflush(stdout);
    free(dummy_moved);
}

#endif
