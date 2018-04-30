#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#if defined(PARALLEL) || defined(EXCHANGE)
#include <mpi.h>
#endif
#include "mc.h"
#include "tables.h"
#include "mt.h"
#include "rotations.h"
#include "util.h"
#ifdef __unix__
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/sysinfo.h>
#endif

/*Update coordinates for fragment ifrag.*/
//consider making this a member function of fragment
void simulation::update_coords(int ifrag, double * centers, double * orients, double * coords)
{
    fragmenttype * frag;
    int iatomstart;
    frag = fragtypes[frags[ifrag].type];
    iatomstart=frags[ifrag].start;
    frag->get_coords(&centers[3*ifrag],&orients[4*ifrag],&coords[3*iatomstart]);
}


/*use xyz file format.*/



/*Copy fragment i's center, orientation, and coordinates from "1" to "2".*/
void simulation::copy_frag(int ifrag, double * center1, double * orient1, double * coords1, double * center2, double * orient2, double * coords2)
{
    int iatomstart,iatomend,iatom,k;
    iatomstart=frags[ifrag].start;
    iatomend=iatomstart+fragtypes[frags[ifrag].type]->natom-1;
    //printf("fragment %d start %d end %d\n",ifrag,iatomstart,iatomend);
    for (k=0; k<3; k++) center2[3*ifrag+k]=center1[3*ifrag+k];
    for (k=0; k<4; k++) orient2[4*ifrag+k]=orient1[4*ifrag+k];
    //for (k=3*iatomstart; k<=(3*iatomend+2); k++) coords2[k]=coords1[k];
    if (exact) for (iatom=iatomstart; iatom<=iatomend; iatom++)
        for (k=0; k<3; k++) coords2[3*iatom+k]=coords1[3*iatom+k];
    //update_coords(ifrag,center2,orient2,coords2);
}

double simulation::interaction_energy(int pbc, int ifrag, int jfrag, double * center, double * orient,double * coords)
{
    //Switches between exact_interaction_energy and table_interaction_energy.
    int reffrag,otherfrag,reftype,othertype,itype,jtype,iatomstart,jatomstart,inatom,jnatom,k;
    double en,enexact,entable,rij[3],r2,rdiff;
    table * tbl;
    go_model_info * go_model;
    //Master spherical cutoff ensures same for both table-based and exact simulations.
#ifdef TIMERS
    switch_timer(TIMER_CHECK_CUTOFF);
#endif
    //this cannot yet be replaced with pbc_distance2 because we need the displacement vector rij
    for (k=0; k<3; k++) {
        rij[k]=center[3*jfrag+k]-center[3*ifrag+k];
        if (pbc) {
            if (rij[k]>halfboxsize) rij[k]-=boxsize;
            if (rij[k]<-halfboxsize) rij[k]+=boxsize;
        }
    }
    r2=rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2];
    if (r2>cutoff2) {
#ifdef TIMERS
        switch_timer(TIMER_INT_OTHER);
#endif
        return 0.0;
    }
    entable=0;
    enexact=0;
    enevalcount++;
    itype=frags[ifrag].type;
    jtype=frags[jfrag].type;
    if ((!exact) || enwrite) {
#ifdef TIMERS
        switch_timer(TIMER_INT_PREP);
#endif
        //figure out which fragment will be the reference.
        //ensure interactions are always taken in a consistent order

        if (itype==jtype) { //it really doesn't matter which is the reference fragment
            tbl=tables[itype][jtype];
            if (tbl==NULL) {
                printf("Table for fragment types %s and %s not loaded.\n",fragtypes[itype]->fragname, fragtypes[jtype]->fragname);
                die();
            }
            if (ifrag<jfrag) {
                reffrag=ifrag;
                otherfrag=jfrag;
            } else if (ifrag>jfrag) {
                reffrag=jfrag;
                otherfrag=ifrag;
                for (k=0; k<3; k++) rij[k]=-rij[k];
            } else { //shouldn't happen
                printf("Interaction energy error.\n");
                die();
            }
        } else {
            //We need to consult the table to see which fragment should be the reference.
            //This guarantees consistency in the choice of which fragment is the reference between table generation and table use.
            tbl=tables[itype][jtype];
            if (tbl==NULL) tbl=tables[jtype][itype];
            if (tbl==NULL) {
                printf("Table for fragment types %s and %s not loaded.\n",fragtypes[itype]->fragname, fragtypes[jtype]->fragname);
                die();
            }
            if (itype==tbl->reffragtype) {
                reffrag=ifrag;
                otherfrag=jfrag;
            } else if (jtype==tbl->reffragtype) {
                reffrag=jfrag;
                otherfrag=ifrag;
                for (k=0; k<3; k++) rij[k]=-rij[k];
            } else { //shouldn't happen
                printf("Interaction energy error.\n");
                die();
            }
        }
        //double table_interaction_energy(double r, double * rij, int reftype, double * qref, int othertype, double * qother)
        //scalefactors should be symmetric
        entable=scalefactors[itype][jtype]*tbl->table_interaction_energy(enwrite,interp,r2,&rij[0],&orient[4*reffrag],&orient[4*otherfrag],&rdiff);
        en=entable;
        /*if (en_by_table!=NULL) {
            if (itype<=jtype) tableindex=itype*nfragtypes+jtype; else tableindex=jtype*nfragtypes+itype;
            en_by_table[tableindex]+=entable;
        }*/
        entablecount++;
    }
    if (exact || enwrite || (entable<=INVALID_ENERGY)) {
        inatom=fragtypes[itype]->natom;
        jnatom=fragtypes[jtype]->natom;
        iatomstart=frags[ifrag].start;
        jatomstart=frags[jfrag].start;
        //double exact_interaction_energy(double eps, int natom1, int * types1, double * coords1, int natom2, int * types2, double * coords2)
        //if (ifrag<jfrag) {
#ifdef TIMERS
        switch_timer(TIMER_INT_EXACT);
#endif
        //enexact=ffield->exact_interaction_energy(pbc,halfboxsize,boxsize,eps,rdie,fragtypes[itype],&coords[3*iatomstart],fragtypes[jtype],&coords[3*jatomstart]
        //match itype and jtype to corresponding for go model.  Go models are not symmetric as to fragment types.
        if (go_models[itype][jtype]!=NULL) {
            go_model=go_models[itype][jtype];
            enexact=go_model->energy(pbc,halfboxsize,boxsize,&params,fragtypes[itype],&coords[3*iatomstart],fragtypes[jtype],&coords[3*jatomstart]);
        } else if (go_models[jtype][itype]!=NULL) { //must switch the two fragments
            go_model=go_models[jtype][itype];
            enexact=go_model->energy(pbc,halfboxsize,boxsize,&params,fragtypes[jtype],&coords[3*jatomstart],fragtypes[itype],&coords[3*iatomstart]);
        } else {
            //something's wrong
            printf("Go model not loaded for fragment types %s and %s.\n",fragtypes[itype]->fragname,fragtypes[jtype]->fragname);
            die();
        }
        en=enexact;
        enexactcount++;
#ifdef TIMERS
        switch_timer(TIMER_INT_OTHER);
#endif
            //=exact_interaction_energy(eps,jnatom,&frags[jtype].types,&coords[3*jatomstart],inatom,&frags[itype].types,&coords[3*iatomstart]);
        /*} else {
            //en=exact_interaction_energy(eps,inatom,&frags[itype].types,&coords[3*iatomstart],jnatom,&frags[jtype].types,&coords[3*jatomstart]);
            en=exact_interaction_energy(eps,jnatom,&frags[jtype].types[0],&coords[3*jatomstart],inatom,&frags[itype].types[0],&coords[3*iatomstart]);
        }*/
    }
    /*if (enwrite) {
        fprintf(energy_output,"%d %d %.4f %.10f %.10f %.10f\n",ifrag,jfrag,rdiff,enexact,entable,en);
        if (fabs(enexact-entable)>5.0) write_pair_pdb(pairs_output,ifrag,jfrag,coords);
    }*/
    //if (exact) return enexact; else return entable;
#ifdef DEBUG
    printf("Interaction energy: %d %s %d %s %.16f\n",ifrag,fragtypes[frags[ifrag].type]->fragname,jfrag,fragtypes[frags[jfrag].type]->fragname,en);
#endif
    return en;
}

#ifdef UMBRELLA
double simulation::umbrella_energy(double * center)
{
    int k;
    double r,r2,dx[3];
    r=sqrt(pbc_distance2(pbc,halfboxsize,boxsize,&center[3*ifragumb],&center[3*jfragumb]));
    return 0.5*kumb*(r-rumb)*(r-rumb);
}
#endif

/*Change in energy upon moving fragment.*/
/*this assumes only one fragment has been moved.*/
//check nonbond list.
double simulation::moved_energy(int imoved, double * center, double * orient, double * coords)
{
    double en;
    int j,jfrag;
    //eold=0.0;
#ifdef TIMERS
    switch_timer(TIMER_INT_OTHER);
#endif
    en=0.0;
    if (use_nb_list) {
        for (j=0; j<frag_nblist->nb_list_count[imoved]; j++) {
            jfrag=frag_nblist->nonbond_list[frag_nblist->nb_list_per_frag*imoved+j];
            //printf("Interaction %d %d\n",imoved,jfrag);
                //fprintf(energy_output,"NEW\n");
            en+=interaction_energy(pbc,imoved,jfrag,center,orient,coords);
                //fprintf(energy_output,"OLD\n");
            //eold+=interaction_energy(pbc,imoved,jfrag,oldcenter,oldorient,oldcoords);
                //fprintf(energy_output,"moved_energy: imoved, j, enew, eold, de = %d %d %.4f %.4f %.4f\n",imoved,j,enew,eold,enew-eold);
        }
    } else {
        for (jfrag=0; jfrag<nfrag; jfrag++) if (jfrag!=imoved) en+=interaction_energy(pbc,imoved,jfrag,center,orient,coords);
    }
#ifdef TIMERS
    switch_timer(TIMER_OTHER);
#endif
#ifdef UMBRELLA
     en+=umbrella_energy(center);
#endif
     return en;
}

//this only works on "new" coordinates.
double simulation::total_energy(double * center, double * orient, double * coords)
{
    int i,j;
    double energy, inte;
#ifdef TIMERS
    switch_timer(TIMER_INT_OTHER);
#endif

    energy=0.0;
    for (i=0; i<nfrag; i++)
        for (j=(i+1); j<nfrag; j++)
            energy+=interaction_energy(pbc,i,j,center,orient,coords);

        /*{
            inte=interaction_energy(i,j,newcenter,neworient,newcoords);
            fprintf(energy_output,"total_energy: i, j, inte, energy = %d %d %.4f %.4f\n",i,j,inte,energy);
            energy+=inte;
        }*/
#ifdef TIMERS
    switch_timer(TIMER_OTHER);
#endif
#ifdef UMBRELLA
     energy+=umbrella_energy(newcenter);
#endif
    return energy;
}

/*Trial move. Either translational or orientational move.*/
void simulation::mcmove(int * imoved, int * movetype, double * center, double * orient, double * coords)
{
    int ifrag,k,move;
    double r,x[3],m,q[4],newq[4];
    r=genrand_real3();
    ifrag=floor(r*nfrag);
    *imoved=ifrag;
    r=genrand_real3();
    for (move=1; move<=NUM_MOVES; move++) {
        if (r<=cumprob[move]) break;
    }
    *movetype = move;
    switch(move) {
        case MOVE_TRANS:
        /*it's a translational move*/
        //printf("translational move\n");
            //use a spherical distribution, maybe mayke this gaussian in future?
            do {
               m=0;
               for (k=0; k<3; k++) {
                   x[k]=2.0*genrand_real3()-1.0;
                   m+=x[k]*x[k];
               }
            } while (m>=1.0);
            //x[k] is a random vector within the unit sphere;
            for (k=0; k<3; k++) center[3*ifrag+k]+=x[k]*dtrans;
            break;
        case MOVE_ORIENT:
        /*it's orientational*/
        //printf("orientational move\n");
            rand_small_quat(dorient,q);
            multiply_quat(&orient[4*ifrag],q,newq);
            normalize_quat(newq);
            for (k=0; k<4; k++) orient[4*ifrag+k]=newq[k];
            break;
	/*case MOVE_UNIF_ORIENT:
            rand_unif_quat(q);
            multiply_quat(&orient[4*ifrag],q,newq);
            normalize_quat(newq);
            for (k=0; k<4; k++) orient[4*ifrag+k]=newq[k];
            break;*/
        default: /*This should never happen.*/
            printf("Error in switch statement.\n");
            die();
    }
    //we only need to update the coordinates if doing an exact simulation
    if (exact) update_coords(ifrag,center,orient,coords);
    //if (pbc) recenter();
}


//Recenter fragments within the box.
void simulation::recenter(void)
{
    int ifrag,k,flag;
    for (ifrag=0; ifrag<nfrag; ifrag++)
    {
        flag=FALSE;
        for (k=0; k<3; k++) {
            if (oldcenter[3*ifrag+k]>halfboxsize) {
                oldcenter[3*ifrag+k]-=boxsize;
                flag=TRUE;
            }
            if (oldcenter[3*ifrag+k]<-halfboxsize) {
                oldcenter[3*ifrag+k]+=boxsize;
                flag=TRUE;
            }
            if (flag) {
                update_coords(ifrag,oldcenter,oldorient,oldcoords);
                copy_frag(ifrag,oldcenter,oldorient,oldcoords,newcenter,neworient,newcoords);
            }
        }
    }
}

/*The Monte Carlo Loop.*/
//void mcloop(FILE * xyzoutput,FILE * quatoutput)
void simulation::mcloop(void)
{
    long int istep,nacc[NUM_MOVES+1],natt[NUM_MOVES+1];
    int movedfrag,movetype,i,ifrag,new_nb_list;
    double cum_energy,fresh_energy,enew,eold,de,de2,r,p,accrate,mctime,elapsedtime,c[3],q[4],rmsd;
    clock_t starttime;
    time_t start,end;
#ifdef __unix__
    struct rusage usage;
    struct sysinfo si;
    struct rlimit as_limit;
#endif
    /*At the beginning "new" and "old" coordinates are the same.*/
    //new_nonbond_list=NULL;
    //use_nb_list=(listcutoff>cutoff);
    if (nsave_xyz>0) {
        xyzoutput=fopen(xyzfname,"wb"); //dcd file
        if (xyzoutput==NULL) {
            printf("Cannot open trajectory file %s.\n",xyzfname);
            die();
        }
        write_dcd_header(xyzoutput);
        printf("Will write coordinates to file %s\n",xyzfname);
    } else {
        printf("No coordinate trajectory file will be written.\n");
    }
    trim_string(quatfname);
    quatoutput=fopen(quatfname,"w");
    if (quatoutput==NULL) {
        printf("Cannot open trajectory file %s.\n",quatfname);
        die();
    }
    printf("Will write centers/quaternions to file %s\n",quatfname);
    if (use_nb_list) {
        frag_nblist=new fragment_nblist(nfrag,listcutoff);
        frag_nblist->create_nonbond_list(pbc,halfboxsize,boxsize,nfrag,newcenter);
    }
    fresh_energy=total_energy(newcenter,neworient,newcoords);
    cum_energy=fresh_energy;
    for (i=1; i<=NUM_MOVES; i++) {
        nacc[i]=0;
        natt[i]=0;
    }
#ifdef __unix__
    getrusage(RUSAGE_SELF,&usage);
    sysinfo(&si);
    getrlimit(RLIMIT_AS,&as_limit);
    printf("Total %.2f MB used.\n",((double) usage.ru_maxrss)/1024);
    printf("Total %.2f MB free out of %.2f MB available.\n",((double)(si.freeram+si.bufferram))/MB,((double)si.totalram)/MB);
    printf("Address space limit %.2f MB (soft), %.2f MB (hard).\n",((double) as_limit.rlim_cur)/MB, ((double)as_limit.rlim_max)/MB);
#endif
#ifdef TIMERS
    init_timers();
#endif
    time(&start);
    printf("Starting Monte Carlo at %s\n",ctime(&start));
    enexactcount=0;
    entablecount=0;
    enevalcount=0;
    printf("Step %ld: Cumulative energy = %.4f   Fresh energy = %.4f\n",0,cum_energy,fresh_energy);
    //die();
    starttime=clock();
    for (istep=1; istep<=nmcstep; istep++){
#ifdef TIMERS
         switch_timer(TIMER_MC_MOVE);
#endif
         mcmove(&movedfrag,&movetype,newcenter,neworient,newcoords);
         eold=moved_energy(movedfrag,oldcenter,oldorient,oldcoords);
         if (use_nb_list) {
#ifdef TIMERS
             switch_timer(TIMER_NB_LIST);
#endif

            new_nb_list=frag_nblist->check_nb_list(pbc,halfboxsize,boxsize,cutoff,nfrag,movedfrag,newcenter);
            if (new_nb_list) {
                //printf("New nonbond list at step %ld\n",istep);
                //if (pbc) recenter();
                frag_nblist->create_nonbond_list(pbc,halfboxsize,boxsize,nfrag,newcenter);
            }
#ifdef TIMERS
             switch_timer(TIMER_OTHER);
#endif
         } else { //still need to recenter if no list
            if (pbc) recenter();
         }
         enew=moved_energy(movedfrag,newcenter,neworient,newcoords);
         de=enew-eold;
         natt[movetype]++;
         //fresh_energy=total_energy(newcenter,neworient);
         //de2=fresh_energy-cum_energy;
         //printf("mcloop: step, de, de2 = %d %.4f %.4f\n",istep,de,de2);
         //This guards against floating point errors.
         if (de<0) p=1.0; else p=exp(-beta*de);
         r=genrand_real3();
         if (r<p) {
             //ACCEPTED
             cum_energy+=de;
             nacc[movetype]++;
             copy_frag(movedfrag,newcenter,neworient,newcoords,oldcenter,oldorient,oldcoords);
             if (de<=-0.5*DUMMY_ENERGY) cum_energy=total_energy(newcenter,neworient,newcoords); //This guards against numerical errors related to "declashing."
         } else {
             //REJECTED
             copy_frag(movedfrag,oldcenter,oldorient,oldcoords,newcenter,neworient,newcoords);
             //If we made a new nonbond list, we now need to get rid of it.  This is duplicative, but avoids the need to manage two nonbond lists.
             if (use_nb_list && new_nb_list) {
#ifdef TIMERS
                 switch_timer(TIMER_NB_LIST);
#endif
                 frag_nblist->create_nonbond_list(pbc,halfboxsize,boxsize,nfrag,oldcenter);
#ifdef TIMERS
                 switch_timer(TIMER_OTHER);
#endif
             }
         }
         //check_nb_list(istep);
         /*Now, "new" and "old" should be the same again, and on the MC trajectory.*/
         /*Handle saving and printing.*/
         if ((nsave_xyz>0) && ((istep%nsave_xyz)==0)) {
              //we don't update coordinates on the moves (for performance reasons), so we need to update them now
              for (ifrag=0; ifrag<nfrag; ifrag++) update_coords(ifrag,newcenter,neworient,newcoords);
              write_dcd_frame(xyzoutput,newcoords);
         }
         //if ((istep%nsave_xyz)==0) write_frame_pdb(xyzoutput,istep,newcoords);
         if ((istep%nsave_quat)==0) write_frame_quat(quatoutput,nprevstep+istep,newcenter,neworient);
         if ((istep%nprint)==0) {
             fresh_energy=total_energy(newcenter,neworient,newcoords);
             printf("Step %ld: Cumulative energy = %.4f   Fresh energy = %.4f Deviation %.6f\n",istep,cum_energy,fresh_energy,cum_energy-fresh_energy);
             for (i=1; i<=NUM_MOVES; i++) {
                if (natt[i]>0) accrate=((double) nacc[i]/(double) natt[i])*100.0; else accrate=0.0;
                printf("%s moves   Attempted %ld   Accepted %ld   Acceptance rate=%.2f%%\n",mc_move_names[i],natt[i],nacc[i],accrate);
             }
             printf("Exact/table/overall energy evaluations: %ld %ld %ld\n",enexactcount,entablecount,enevalcount);
             if (initcenter!=NULL) {
                 rmsd_fit(nfrag,mass,initcenter,newcenter,c,q,&rmsd);
                 printf("Fragment center RMSD: %.3f A\n",rmsd);
             }
             if (fabs(cum_energy-fresh_energy)>10.0) {
                 //print_energies(stdout,FALSE,"Cum. energy:",istep,cum_energies,cum_energy);
                 printf("Too much deviation between cumulative and fresh energies.\n");
                 if (strlen(endrestartfname)>0) write_restart(nprevstep+istep,endrestartfname);
                 fflush(stdout);
                 die();
             }
             fflush(stdout);
         }
#ifdef EXCHANGE
         if ((istep%exchfreq)==0) {
                exchange(istep/exchfreq,&cum_energy,newcenter,neworient,newcoords);
                //print_energies(stdout,FALSE,"Cum: ",istep,cum_energies,cum_energy);
                fresh_energy=total_energy(newcenter,neworient,newcoords);
                //print_energies(stdout,FALSE,"Energy:",istep,fresh_energies,fresh_energy);
        }
#endif
    } //end of main loop "for(istep=1; istep<=nmcstep; istep++)"
#if defined(PARALLEL) || defined(EXCHANGE)
    //printf("node %d before final barrier\n",mynod);
    //fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);
    //printf("node %d after final barrier\n",mynod);
    //fflush(stdout);
#endif
#ifdef TIMERS
    switch_timer(TIMER_NONE);
#endif
    mctime=((double)(clock()-starttime))/CLOCKS_PER_SEC;
    time(&end);
    elapsedtime=difftime(end,start);
    printf("Total time %.2f seconds.\n",mctime);
    printf("%.2f MC steps per second.\n",nmcstep/mctime);
    printf("Total elapsed time %.2f seconds.\n",elapsedtime);
#ifdef __unix__
    //print some statistics regarding memory usage
    getrusage(RUSAGE_SELF,&usage);
    printf("CPU time: %.2f sec user mode, %.2f sec system mode.\n",convtime(usage.ru_utime),convtime(usage.ru_stime));
    printf("Total %.2f MB used, %ld page faults.\n",((double)usage.ru_maxrss)/1024,usage.ru_majflt);
#endif
#ifdef TIMERS
    print_timers();
#endif
    if (strlen(endrestartfname)>0) write_restart(nprevstep+nmcstep,endrestartfname);
    if (nsave_xyz>0) fclose(xyzoutput);
    fclose(quatoutput);
#ifdef EXCHANGE
    //printf("node %d closing replica log\n",mynod);
    //fflush(stdout);
    MPI_File_close(&rexlog);
    //printf("node %d closed replica log.\n",mynod);
    //fflush(stdout);
#endif
    if (enwrite) {
        fclose(energy_output);
        fclose(pairs_output);
    }
}


