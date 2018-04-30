#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "rotations.h"
#include "util.h"
#include "mt.h"
#include "morph.h"


using namespace std;

#define MAX_FRAGMENT_NAME  32
void read_frame_quat(FILE * input, int nfrag, long long int * istep, int * fragtypes, double * center, double * orient)
{
    int ifrag,ifragx,k,itype;
    double x[3],q[4];
    char fragname[MAX_FRAGMENT_NAME];
    for (ifragx=0; ifragx<nfrag; ifragx++){
        fscanf(input,"%lld %d %s %lg %lg %lg %lg %lg %lg %lg\n",istep,&ifrag,fragname,&x[0],&x[1],&x[2],&q[0],&q[1],&q[2],&q[3]);
        ifrag--; //zero-based
        //ifrag=ifrag-1;
        fragtypes[ifrag]=-1;
        for (itype=0; itype<NFRAGTYPES; itype++) if (strcasecmp(fragname,fragnames[itype])==0) fragtypes[ifrag]=itype;
        if ((ifrag>=nfrag) || (fragtypes[ifrag]<0)) {
            printf("error reading trajectory file\n");
            die();
        }
        normalize_quat(q);
        for (k=0; k<3; k++) center[3*ifrag+k]=x[k];
        for (k=0; k<4; k++) orient[4*ifrag+k]=q[k];
        normalize_quat(&orient[4*ifrag]);
    }
}

void read_template(char * fname, int * ntmpfrag, int * * tmpfragtypes, double * * tmpcenter, double * * tmporient)
{
    int ifrag,i,itype,inewfrag,natomx,k;
    long int nprevstep;
    double x[3],q[4];
    char fragname[MAX_FRAGMENT_NAME];
    FILE * f;
    f=fopen(fname,"r");
    if (f==NULL) {
        printf("Could not open restart file %s\n",fname);
        die();
    }
    //printf("Reading template file %s.\n",fname);
    fscanf(f,"%d %d %ld\n",&natomx,ntmpfrag,&nprevstep);
    *tmpcenter=(double *) checkalloc(3*(*ntmpfrag),sizeof(double));
    *tmporient=(double *) checkalloc(4*(*ntmpfrag),sizeof(double));
    *tmpfragtypes=(int *) checkalloc(*ntmpfrag,sizeof(int));
    /*Read the fragment, center, and quaternion.*/
    /*if (nprevstep>0) {
        printf("Reading RNG state vector from file %s.\n",fname);
        read_rng_state(f);
        reseed=false;
    }*/
    while (!feof(f)) {
    //for (i=0;i<nfrag;i++){
        fscanf(f,"%d %s %lg %lg %lg %lg %lg %lg %lg\n",&ifrag,fragname,&x[0],&x[1],&x[2],&q[0],&q[1],&q[2],&q[3]);
        /*itype=frag_type_by_name(name);
        if (itype<0) {
            printf("Unrecognized type in restart file %s.\n",fname);
            die();
        }
        inewfrag=insert_fragment(itype);*/
        ifrag=ifrag-1;
        /*if (ifrag!=inewfrag) {
            printf("Error reading restart file %s.\n",fname);
            die();
        }*/
        (*tmpfragtypes)[ifrag]=-1;
        for (itype=0; itype<NFRAGTYPES; itype++) if (strcasecmp(fragname,fragnames[itype])==0) (*tmpfragtypes)[ifrag]=itype;
        normalize_quat(q);
        for (k=0; k<3; k++) (*tmpcenter)[3*ifrag+k]=x[k];
        for (k=0; k<4; k++) (*tmporient)[4*ifrag+k]=q[k];
        normalize_quat(&(*tmporient)[4*ifrag]);
        //update_coords(ifrag,oldcenter,oldorient,oldcoords);
        //copy_frag(ifrag,oldcenter,oldorient,oldcoords,newcenter,neworient,newcoords);
    }
    fclose(f);
}

//arguments: nfrag boxsize template traj output
int main(int argc, char * argv[])
{
    int nfrag,ntemplatefrag;
    double boxsize,halfboxsize,cutoff,anglecutoff,distratio;
    bool pbc;
    FILE * traj;
    FILE * templatefile;
    FILE * output;
    FILE * pdb_output;
    double * center;
    double * orient;
    int * fragtypes;
    double * templatecenter;
    double * templateorient;
    int * tmpfragtypes;
    int * mapping;
    long long int istep,iframe;
    double rmsd;
    clusters * clus;
    char buffer[255];
    time_t now,start;
    double seconds;
    int ifrag,i,j,k;
    nfrag=atoi(argv[1]);
    ntemplatefrag=120;
    boxsize=atof(argv[2]);
    //cutoff=atof(argv[3]);
    distratio=atof(argv[3]);
    anglecutoff=atof(argv[4])*DEG_TO_RAD;
    pbc=(boxsize>0);
    halfboxsize=0.5*boxsize;
    /*templatefile=fopen(argv[3],"r");
    if (templatefile==NULL) {
        printf("Could not open file %s\n",argv[3]);
        die();
    }*/
    //trim_string(argv[4]);
    traj=fopen(argv[6],"r");
    if (traj==NULL) {
        printf("Could not open file %s\n",argv[6]);
        die();
    }
    output=fopen(argv[7],"w");
    if (output==NULL) {
        printf("Could not open file %s\n",argv[7]);
        die();
    }
    pdb_output=NULL;
    if (argc>6) pdb_output=fopen(argv[8],"w");
    center=(double *) checkalloc(3*nfrag,sizeof(double));
    orient=(double *) checkalloc(4*nfrag,sizeof(double));
    fragtypes=(int *) checkalloc(nfrag,sizeof(int));
    //templatecenter=(double *) checkalloc(3*ntemplatefrag,sizeof(double));
    //templateorient=(double *) checkalloc(4*ntemplatefrag,sizeof(double));
    mapping=(int *) checkalloc(nfrag,sizeof(int));
    for (ifrag=0; ifrag<nfrag; ifrag++) mapping[ifrag]=ifrag; //default mapping
    read_template(argv[5],&ntemplatefrag,&tmpfragtypes,&templatecenter,&templateorient);
    /*if (ntemplatefrag!=nfrag) {
        printf("trajectory does not match template, nfrag = %d, ntemplatefrag = %d\n",nfrag,ntemplatefrag);
        die();
    }*/
    //fclose(templatefile);
    time(&now);
    strncpy(buffer,ctime(&now),sizeof(buffer));
    printf("Clique-based analysis starting at %s\n",buffer);
    start=now; //start of current frame;
    iframe=1;
    while (!feof(traj)) {
        read_frame_quat(traj,nfrag,&istep,fragtypes,center,orient);

        //separate half the fragments by translating them along the z-axis by 20 A for testing purposes
        //for (ifrag=nfrag/2; ifrag<nfrag; ifrag++) center[3*ifrag+2]+=20;
        //printf("Step %lld\n",istep);
        clus=new clusters(pbc,boxsize,halfboxsize,nfrag,fragtypes,center,orient,ntemplatefrag,tmpfragtypes,templatecenter,templateorient,distratio,anglecutoff);
        //fprintf(output,"Step %lld\n",istep);
        clus->report(iframe,output);
        if (pdb_output!=NULL) clus->write_pdb_frame(iframe,center,pdb_output);
        fflush(output);
        //output information... or have "clusters" do it
        delete clus;
        if ((iframe%1)==0) {
            time(&now);
            strncpy(buffer,ctime(&now),sizeof(buffer));
            seconds=difftime(now,start);
            printf("frame %lld completed at %s, took %.2f seconds\n",iframe,buffer,seconds);
            start=now;
        }
        iframe++;
        //mapped_pbc_rmsd_fit(pbc,boxsize,halfboxsize,nfrag,center,ntemplatefrag,templatecenter,mapping,&trans,&rmsd);
        //fprintf(output,"%lld %.3f\n",istep,rmsd);
        //break;
    }
    fclose(traj);
    fclose(output);
    if (pdb_output!=NULL) fclose(pdb_output);
    free(center);
    free(orient);
    free(mapping);
    return 0;
}
