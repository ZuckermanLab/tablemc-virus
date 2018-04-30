#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include "rotations.h"
#include "util.h"
#include "mt.h"
#include "morph.h"

//maximum number of iterations
#define MAX_ITER       10
#define MAX_TOTAL_ITER 5 //per fragment
//"pre iterations" -- use an infinite cutoff, do not change clusters, do not declare convergence
#define PRE_ITER   3
using namespace std;

//this is needed for the i/o routines
void pdb_coordinate(double val, char * str)
{
     if (fabs(val)>1000) sprintf(str,"%8.2f",val); else sprintf(str,"%8.3f",val);
}

//Compute rmsd between two point sets of unequal size using the mapping
//mapping[icoord1] is the entry in coords2 that corresponds to coord1.
//The mapping may have entries of "-1" which mean that the corresponding coords in
/*void mapped_pbc_rmsd_fit(bool pbc, double boxsize, double halfboxsize, int ncoord1, double * coords1, int ncoord2, double * coords2, int * mapping, rigid_trans * trans, double * rmsd)
{
    int nactualcoord,icoord1,icoord2,iactualcoord,k,first_mapped_frag;
    double * mappedcoords1;
    double * mappedcoords2;
    double * weights;
    double dx[3];
    nactualcoord=0;
    for (icoord1=0; icoord1<ncoord1; icoord1++) if (mapping[icoord1]>=0) nactualcoord++;
    mappedcoords1=(double *) checkalloc(3*nactualcoord,sizeof(double));
    mappedcoords2=(double *) checkalloc(3*nactualcoord,sizeof(double));
    weights=(double *) checkalloc(nactualcoord,sizeof(double));
    iactualcoord=0;
    for (icoord1=0; icoord1<ncoord1; icoord1++) if (mapping[icoord1]>=0) {
        for (k=0; k<3; k++) mappedcoords1[3*iactualcoord+k]=coords1[3*icoord1+k];
        icoord2=mapping[icoord1];
        for (k=0; k<3; k++) mappedcoords2[3*iactualcoord+k]=coords2[3*icoord2+k];
        if (pbc) for (k=0; k<3; k++) { //Ensure that all fragments are the same periodic cell as the first fragment.
            dx[k]=mappedcoords1[3*iactualcoord+k]-mappedcoords1[k];
            if (dx[k]>halfboxsize) mappedcoords1[3*iactualcoord+k]-=boxsize;
            if (dx[k]<-halfboxsize) mappedcoords1[3*iactualcoord+k]+=boxsize;
        }
        weights[iactualcoord]=1.0;
        iactualcoord++;
    }
    rmsd_fit(nactualcoord,weights,mappedcoords1,mappedcoords2,&trans->disp[0],&trans->rot[0],rmsd);
    free(mappedcoords1);
    free(mappedcoords2);
    free(weights);
}*/
//Find clusters (could be partially formed capsids)
clusters::clusters(bool pbc, double boxsize, double halfboxsize, int _nfrag, int * fragtypes, double * center, double * orient,
    int _ntmpfrag, int * tmpfragtypes, double * tmpcenter, double * tmporient, graph * tmpcontacts, double cutoff, double anglecutoff)
{
    int ifrag,iclus,inewfrag,jnewfrag,iter,total_iter,nassigned;
    bool done,conv,preiter;
    double thecutoff;
    FILE * debug_output;
    //double rmsdconv,rmsddiff;
    int * old_cluster_id;
    int * old_mapping;
    int max_total_iter;
    done=false;
    nclusters=0;
    nfrag=_nfrag;
    ntmpfrag=_ntmpfrag;
    info=NULL;
    mindist=NULL;
    max_total_iter=MAX_TOTAL_ITER*nfrag;
    //init_cluster_info(&info[nclusters-1]);
    cluster_id=(int *) checkalloc(nfrag,sizeof(int));
    mapping=(int *) checkalloc(nfrag,sizeof(int));
    old_cluster_id=(int *) checkalloc(nfrag,sizeof(int));
    old_mapping=(int *) checkalloc(nfrag,sizeof(int));
    for (ifrag=0; ifrag<nfrag; ifrag++) {
        cluster_id[ifrag]=-1; //not assigned yet
        mapping[ifrag]=-1;
        old_cluster_id[ifrag]=-1;
        old_mapping[ifrag]=-1;
    }
    //set up the first cluster
    //guess_clusters(pbc,boxsize,halfboxsize,center,cutoff); //heuristic on cutoff for distance matrix
    add_new_cluster(pbc,boxsize,halfboxsize,0,fragtypes,center,orient,tmpfragtypes,tmpcenter,tmporient);
    //insert_cluster(fragtypes,tmpfragtypes);
    //find_transformations(pbc,boxsize,halfboxsize,center,tmpcenter);
    //assignment algorithm:
#ifdef DEBUG
    debug_output=fopen("debug.pdb","w");
#endif
    done=false;
    total_iter=1;
    while (!done && (total_iter<max_total_iter)) {
        //Alternately reassign the fragments to clusters and update the transformations
        //until the mapping converges.
        conv=false;
        iter=1;
        while ((!conv) && (iter<MAX_ITER)) {
            /*preiter=(iter<=PRE_ITER);
            if (preiter) thecutoff=1e10; else*/ thecutoff=cutoff;
#ifdef DEBUG
            fprintf(debug_output,"REMARK frame %d\n",total_iter);
#endif
            reassign(pbc,boxsize,halfboxsize,fragtypes,center,orient,tmpfragtypes,tmpcenter,tmporient,tmpcontacts,thecutoff,anglecutoff,debug_output);
#ifdef DEBUG
            fprintf(debug_output,"END\n");
            fflush(debug_output);
#endif
            find_transformations(pbc,boxsize,halfboxsize,center,orient,tmpcenter,tmporient);
            //this sorts the clusters in descending order by size, relabeling both cluster_id and old_cluster_id.
            sort_clusters(old_cluster_id);
            //remove empty clusters (which have been placed at the end of the list)
            while (info[nclusters-1].size==0) nclusters--;
            //check to see if we still have at least one fragment assigned
            nassigned=0;
            for (ifrag=0; ifrag<nfrag; ifrag++) if (cluster_id[ifrag]>=0) nassigned++;
            if (nassigned==0) {
                printf("error: lost all assignments at iter=%d total_iter=%d\n",iter,total_iter);
                die();
            }
            //check convergence
            //if (!preiter) {
                conv=true;
                for (ifrag=0; ifrag<nfrag; ifrag++) if ((cluster_id[ifrag]!=old_cluster_id[ifrag]) || (mapping[ifrag]!=old_mapping[ifrag])) conv=false;
                //update old ids and mappings
                for (ifrag=0; ifrag<nfrag; ifrag++) {
                    old_cluster_id[ifrag]=cluster_id[ifrag];
                    old_mapping[ifrag]=mapping[ifrag];
                }
            //} else conv=false;
            iter++;
            total_iter++;
            if ((total_iter%10)==0) printf("%d iterations, total assigned = %d, clusters = %d, max cluster = %d\n",
                                            total_iter,nassigned,nclusters,info[0].size);
        }
        //Check to see if all fragments have been assigned a cluster.  If not, create a new cluster and add the first unassigned fragment to it.
        inewfrag=0;
        while ((inewfrag<nfrag) && (cluster_id[inewfrag]>=0)) inewfrag++;
        done=(inewfrag>=nfrag); //all fragments have been assigned stably
        if (!done) add_new_cluster(pbc,boxsize,halfboxsize,inewfrag,fragtypes,center,orient,tmpfragtypes,tmpcenter,tmporient);
        //check to see if they are all assigned
        /*done=true;
        for (inewfrag=0; inewfrag<nfrag; inewfrag++) done=done && (cluster_id[inewfrag]>=0);
        if (done) break;
        //if not done, insert a new cluster
        insert_cluster(fragtypes,tmpfragtypes);
        find_transformations(pbc,boxsize,halfboxsize,center,tmpcenter);*/
        /*if (!done) {
            //see if there is a second new fragment.  If so we can try to assign them both.
            jnewfrag=inewfrag+1;
            while ((jnewfrag<nfrag) && (cluster_id[jnewfrag]>=0)) jnewfrag++;
            if (jnewfrag<nfrag) add_new_cluster(inewfrag,jnewfrag,center,tmpcenter); else add_new_cluster(inewfrag,center,tmpcenter);
            //find_transformations(pbc,boxsize,halfboxsize,center,tmpcenter);
        }*/
    }
    if (!done) {
        printf("Warning: failed to converge and assign all fragments. %d fragments assigned.\n",nassigned);
        printf("Asssigning all remaining fragments to their own clusters.");
        while (!done) {
            inewfrag=0;
            while ((inewfrag<nfrag) && (cluster_id[inewfrag]>=0)) inewfrag++;
            done=(inewfrag>=nfrag); //all fragments have been assigned stably
            if (!done) add_new_cluster(pbc,boxsize,halfboxsize,inewfrag,fragtypes,center,orient,tmpfragtypes,tmpcenter,tmporient);
        }
    }
    free(old_cluster_id);
    free(old_mapping);
    sort_clusters(NULL);
    get_mindist(pbc,halfboxsize,boxsize,center,orient,tmpcenter,tmporient);
#ifdef DEBUG
    fclose(debug_output);
#endif
}

clusters::~clusters()
{
    int iclus;
    if (info!=NULL) {
        /*for (iclus=find_transformations(pbc,boxsize,halfboxsize,center,tmpcenter);; iclus<nclusters; iclus++) {
            free(info[iclus].fragments);
            free(info[iclus].tmpfragments);
            }*/
        free(info);
    }
    free(cluster_id);
    free(mapping);
    if (mindist!=NULL) free(mindist);
}


//force a single fragment in cluster to fit perfectly
void clusters::force_fit_fragment(bool pbc, double boxsize, double halfboxsize, int ifrag, double * center, double * orient, int itmpfrag, double * tmpcenter, double * tmporient, int iclus)
{
    double q[4],v[3];
    int k;
    //set up rotation to bring one into the other (=tmporient[itmpfrag]^-1 * orient[ifrag])
    conjugate_quat(&tmporient[4*itmpfrag],&q[0]);
    multiply_quat(&q[0],&orient[4*ifrag],&info[iclus].rot[0]);
    //set up displacement to make fragment "ifrag" match template fragment itmpfrag after rotation
    rotate_vector_by_quat(&info[iclus].rot[0],&tmpcenter[3*itmpfrag],&v[0]);
    if (pbc) for (k=0; k<3; k++) {
        info[iclus].disp[k]=center[3*ifrag+k]-v[k];
        if (info[iclus].disp[k]>halfboxsize) info[iclus].disp[k]-=boxsize;
        if (info[iclus].disp[k]<-halfboxsize) info[iclus].disp[k]+=boxsize;
    }

}

//same thing, but for two fragments (finds the rotation about the interfragment axis that makes ifrag's orientation coincide as closely as possible with itmpfrag's
//this relies on an initial rotation from rmsd_fit
void clusters::force_fit_fragment(bool pbc, double boxsize, double halfboxsize, int ifrag, int jfrag, double * center, double * orient, int itmpfrag, int jtmpfrag, double * tmpcenter, double * tmporient, int iclus)
{
    double q1[4],q2[4],q3[4],q4[4],v[3],axis[3],com[3],tmpcom[3],r2,dot;
    int k;
    //start by finding unit vector along axis between two fragments
    r2=0.0;
    for (k=0; k<3; k++) {
        axis[k]=center[3*jfrag+k]-center[3*ifrag+k];
        if (pbc) {
            if (axis[k]>halfboxsize) axis[k]-=boxsize;
            if (axis[k]<-halfboxsize) axis[k]+=boxsize;
        }
        r2+=axis[k]*axis[k];
    }
    r2=sqrt(r2);
    for (k=0; k<3; k++) axis[k]/=r2;
    //conjugate_quat(&tmporient[4*itmpfrag],&q[0]);
    //This finds the net additional rotation needed to make ifrag and itmpfrag coincide.
    //q3 = (template orientation * orig. rotation)^-1 * (frame orientation)
    multiply_quat(&tmporient[4*itmpfrag],&info[iclus].rot[0],&q1[0]);
    conjugate_quat(&q1[0],&q2[0]);
    multiply_quat(&q2[0],&orient[4*ifrag],&q3[0]);
    //Now we "project" this rotation onto the axis by projecting its imaginary component onto the axis and renormalizing.
    dot=0.0;
    for (k=0; k<3; k++) dot+=q3[k+1]*axis[k]; //=Im(q3).axis
    for (k=0; k<3; k++) q3[k+1]=dot*axis[k]; //=(Im(q3).axis) axis
    normalize_quat(&q3[0]);
    //now include this additional rotation in the net rotation (q4)
    multiply_quat(&info[iclus].rot[0],&q3[0],&q4[0]);
    for (k=0; k<4; k++) info[iclus].rot[k]=q4[k];
    //find the COM of both the template and frame cluster.
    for (k=0; k<3; k++) {
        com[k]=(center[3*ifrag+k]+center[3*jfrag+k])/2;
        tmpcom[k]=(tmpcenter[3*itmpfrag+k]+tmpcenter[3*jtmpfrag+k])/2;
    }
    //set up displacement to make them match after rotation
    rotate_vector_by_quat(&info[iclus].rot[0],&tmpcom[0],&v[0]);
    if (pbc) for (k=0; k<3; k++) {
        info[iclus].disp[k]=com[k]-v[k];
        if (info[iclus].disp[k]>halfboxsize) info[iclus].disp[k]-=boxsize;
        if (info[iclus].disp[k]<-halfboxsize) info[iclus].disp[k]+=boxsize;
    }
    //we don't need to do this, called from find_transformations
    /*cluster_id[ifrag]=iclus;
    mapping[ifrag]=itmpfrag;
    cluster_id[jfrag]=iclus;
    mapping[jfrag]=jtmpfrag;*/
}
//Guess the transformation for the next cluster. Returns whether all the clusters have been assigned
//Create a new cluster for each unassigned fragment.
// bool clusters::add_new_clustr(double * center, double * tmpcenter)
void clusters::add_new_cluster(bool pbc, double boxsize, double halfboxsize, int ifrag, int * fragtypes, double * center, double * orient, int * tmpfragtypes, double * tmpcenter, double * tmporient)
{
    int k,inewclus,itmpfrag;
    double q[4],v[3];
    //For each fragment that hasn't been assigned to a cluster.
    //find the first template fragment of matching type
    itmpfrag=0;
    while (fragtypes[ifrag]!=tmpfragtypes[itmpfrag]) itmpfrag++;
        //set up a new cluster
    info=(cluster_info *) checkrealloc(info,nclusters+1,sizeof(cluster_info));
    nclusters++;
    inewclus=nclusters-1;
    info[inewclus].size=1; //will be updated by find_transformations
    info[inewclus].rmsd=0; //will also be updated
    force_fit_fragment(pbc,boxsize,halfboxsize,ifrag,center,orient,itmpfrag,tmpcenter,tmporient,inewclus);
    cluster_id[ifrag]=inewclus;
    mapping[ifrag]=itmpfrag;
}

void clusters::insert_cluster(int * fragtypes, int * tmpfragtypes)
{
    int fragsadded,ifrag,itype,itmpfrag[NFRAGTYPES],inewclus,k;
    bool done;
    info=(cluster_info *) checkrealloc(info,nclusters+1,sizeof(cluster_info));
    nclusters++;
    inewclus=nclusters-1;
    info[inewclus].size=0; //will be updated by find_transformations
    info[inewclus].rmsd=0; //will also be updated
    info[inewclus].rot[0]=1.0; //start with unit quaternion
    for (k=1; k<4; k++) info[inewclus].rot[k]=0.0;
    //assign the next ntmpfrag unassigned fragments to the new cluster.
    //find the first fragment in the template of each type
    for (itype=0; itype<NFRAGTYPES; itype++) {
        itmpfrag[itype]=0;
        while (tmpfragtypes[itmpfrag[itype]]!=itype) itmpfrag[itype]++;
    }
    ifrag=0;
    done=false;
    while (!done) {
        //locate the next unassigned fragment
        while ((ifrag<nfrag) && (cluster_id[ifrag]>=0)) ifrag++;
        if (ifrag>=nfrag) break; //we're done
        itype=fragtypes[ifrag];
        //itmpfrag[itype] is the next fragment from the template of type itype to be assigned.
        if (itmpfrag[itype]>=ntmpfrag) {ifrag++; continue;};  //out of template fragments of this type, can't assign this fragment
        cluster_id[ifrag]=inewclus;
        mapping[ifrag]=itmpfrag[itype]; //the next fragment of the same type to be assigned
        ifrag++; //start searching with the next fragment
        //find the next template fragment of this type
        itmpfrag[itype]++;
        while ((itmpfrag[itype]<ntmpfrag) && (tmpfragtypes[itmpfrag[itype]]!=itype)) itmpfrag[itype]++;
        //check to see if we are out of all fragment types
        done=true;
        for (itype=0; itype<NFRAGTYPES; itype++) done=done && (itmpfrag[itype]>=ntmpfrag);
    }
    //need to finish by calling find-transformations
}
//void clusters::guess_clusters(
//same, but add two new fragments to the cluster.  Does not fill in correct transformation, need to apply find_transformations afterward
/*void clusters::add_new_cluster(int ifrag, int jfrag, double * center, double * tmpcenter)
{
    int k,inewclus;
    //For each fragment that hasn't been assigned to a cluster.
        //set up a new cluster
    info=(cluster_info *) checkrealloc(info,nclusters+1,sizeof(cluster_info));
    nclusters++;
    inewclus=nclusters-1;
    info[inewclus].size=2; //will be updated by find_transformations
    info[inewclus].rmsd=0; //will also be updated
    info[inewclus].rot[0]=1.0; //start with unit quaternion
    for (k=1; k<4; k++) info[inewclus].rot[k]=0.0;
    //set up displacement to make fragment "ifrag" match template fragment 0
    for (k=0; k<3; k++) info[inewclus].disp[k]=center[3*ifrag+k]-tmpcenter[k];
    //temporarily assign the second new fragment to the second fragment in the template.  This may not be correct, we rely on reassign to help us.
    cluster_id[ifrag]=inewclus;
    mapping[ifrag]=0;
    cluster_id[jfrag]=inewclus;
    mapping[jfrag]=1;
}*/

//For each cluster, RMSD fit given the current assignment and mapping.
void clusters::find_transformations(bool pbc, double boxsize, double halfboxsize, double * center, double * orient, double * tmpcenter, double * tmporient)
{
    double * mapped_from_template;
    double * mapped_from_frame;
    double * weights;
    //int * backmap; //gives the fragment for each "actual"
    int iclus, ifrag, itmpfrag,i,k;
    double dx[3],q[4];
    //make sure "info" is properly allocated
    info=(cluster_info *) checkrealloc(info,nclusters,sizeof(cluster_info));
    mapped_from_template=(double *) checkalloc(3*ntmpfrag,sizeof(double));
    mapped_from_frame=(double *) checkalloc(3*ntmpfrag,sizeof(double));
    weights=(double *) checkalloc(ntmpfrag,sizeof(double));
    //backmap=(int *) checkalloc(ntmpfrag,sizeof(int));
    for (iclus=0; iclus<nclusters; iclus++) {
        //clear out the mappedcoords array, put a nan in so that bad stuff can be detected
        for (itmpfrag=0; itmpfrag<3*ntmpfrag; itmpfrag++) mapped_from_template[itmpfrag]=NAN;
        for (itmpfrag=0; itmpfrag<3*ntmpfrag; itmpfrag++) mapped_from_frame[itmpfrag]=NAN;
        //pair off coordinates in the frame and in the template
        i=0;
        for (ifrag=0; ifrag<nfrag; ifrag++) if ((cluster_id[ifrag]==iclus) && (mapping[ifrag]>=0)) { //must have been assigned to a template fragment wtihin the cluster.
            //backmap[i]=ifrag;
            //info[iclus].fragments[i]=ifrag;
            //info[iclus].tmpfragments[i]=mapping[ifrag];
            for (k=0; k<3; k++) mapped_from_template[3*i+k]=tmpcenter[3*mapping[ifrag]+k]; //from the frame
            for (k=0; k<3; k++) mapped_from_frame[3*i+k]=center[3*ifrag+k]; //from the template
            if (pbc) for (k=0; k<3; k++) { //Ensure that all fragments are the same periodic cell as the first fragment in the cluster
                dx[k]=mapped_from_frame[3*i+k]-mapped_from_frame[k];
                if (dx[k]>halfboxsize) mapped_from_frame[3*i+k]-=boxsize;
                if (dx[k]<-halfboxsize) mapped_from_frame[3*i+k]+=boxsize;
            }
            weights[i]=1.0;
            if (i<2) { //see note on these fields in morph.h
                info[iclus].frag[i]=ifrag;
                info[iclus].tmpfrag[i]=mapping[ifrag];
            }
            i++;
        }
        info[iclus].size=i;
        //Done pairing off coordinates, perform the RMSD fit.  Reference is first.
        if (info[iclus].size>1) {
            rmsd_fit(info[iclus].size,weights,mapped_from_template,mapped_from_frame,&info[iclus].disp[0],&info[iclus].rot[0],&info[iclus].rmsd);
        } else if (info[iclus].size==1) {
            //need to fix rotation
            //ifrag=info[iclus].frag[0];
            //itmpfrag=info[iclus].tmpfrag[0];
            force_fit_fragment(pbc,boxsize,halfboxsize,info[iclus].frag[0],center,orient,info[iclus].tmpfrag[0],tmpcenter,tmporient,iclus);
        } /*else { //zero size, ok
            printf("error\n");
            die();
        }*/
        if (info[iclus].size==2) {
            //handle the tricky case where two fragments need to be "force-fitted"
            //this assumes an RMSD fitting has been performed to get an initial rotation
            //probably should change the name of this routine
            force_fit_fragment(pbc,boxsize,halfboxsize,info[iclus].frag[0],info[iclus].frag[1],center,orient,
                info[iclus].tmpfrag[0],info[iclus].tmpfrag[1],tmpcenter,tmporient,iclus);
        }
        //need to do something about two fragments, buut this is more complicated
        //else printf("empty cluster %d\n",iclus);
    }
    free(mapped_from_template);
    free(mapped_from_frame);
    free(weights);
    //free(backmap);
}


struct entry {
    int frag;
    int tmpfrag;
    int clus;
    double dist2;
};

int compare_entries(const void * e1, const void * e2)
{
    double d1 = ((entry *) e1)->dist2;
    double d2 = ((entry *) e2)->dist2;
    if (d1<d2) return -1;
    if (d1>d2) return 1;
    return 0;
}
#define MAX_CLUSTER 10
void clusters::reassign(bool pbc, double boxsize, double halfboxsize, int * fragtypes, double * center, double * orient,
            int * tmpfragtypes, double * tmpcenter, double * tmporient, graph * tmpcontacts, double cutoff, double anglecutoff, FILE * debug_output)
{
    double * transtemplates;
    bool * assigned;
    double com[3],v[3],dx[3],rotmatrix[3][3],cutoff2,dist2,cacutoff,q1[4],cangle;
    int iclus, ifrag,itmpfrag,index,closest,i,k,ncomp,icomp;
    const char * pdbatomfmt = "ATOM  %5d %4s %3s %c%4d    %8s%8s%8s%6.2f%6.2f      %-4d\n";
    char buffer[255],xval[9],yval[9],zval[9];
    vector<entry> entries;
    entry e;
    subset tmpfragments_in_clus;
    long int * comps;
    entries.reserve(nfrag);
    transtemplates=(double *) checkalloc(3*nclusters*ntmpfrag,sizeof(double));
    assigned=(bool *) checkalloc(nclusters*ntmpfrag,sizeof(bool));
    cutoff2=cutoff*cutoff;
    cacutoff=cos(anglecutoff/2); //acutoff = angular cutoff for deviation of
    //unassign all fragments if changing clusters
    for (ifrag=0; ifrag<nfrag; ifrag++) {
        cluster_id[ifrag]=-1;
        mapping[ifrag]=-1;
    }
    //Prepare a set of coordinates containing a copy of the template for each cluster,
    //transformed through the transformation.  In effect this is a "model" of the frame.
    for (iclus=0; iclus<nclusters; iclus++) {
        //transform the template through the transformation
        quat_to_matrix(info[iclus].rot,&rotmatrix[0][0]);
        for (itmpfrag=0; itmpfrag<ntmpfrag; itmpfrag++) {
            index=iclus*ntmpfrag+itmpfrag;
            matmul(&rotmatrix[0][0],&tmpcenter[3*itmpfrag],v);
            for (k=0; k<3; k++) transtemplates[3*index+k]=v[k]+info[iclus].disp[k];
        }
    }
#ifdef DEBUG
    //write this to a pdb file.
    for (index=0; index<MAX_CLUSTER*ntmpfrag; index++)  {
        //Find the closest template fragment to fragment ifrag.
        iclus=index/ntmpfrag;
        itmpfrag=index%ntmpfrag;
        for (k=0; k<3; k++) v[k]=0.0;
        if (iclus<nclusters) for (k=0; k<3; k++) v[k]=transtemplates[3*index+k];
        pdb_coordinate(v[0],xval);
        pdb_coordinate(v[1],yval);
        pdb_coordinate(v[2],zval);
        snprintf(buffer,sizeof(buffer),pdbatomfmt,index+1,"X","XXX",'X',itmpfrag+1,xval,yval,zval,1.0,0.0,iclus+1);
        fputs(buffer,debug_output);
    }
#endif
    //for each fragment, find the atom in the "model" that is closest.  If closer than the cutoff, assign it.
    for (index=0; index<nclusters*ntmpfrag; index++)  {
        //Find the closest template fragment to fragment ifrag.
        iclus=index/ntmpfrag;
        itmpfrag=index%ntmpfrag;
        for (ifrag=0; ifrag<nfrag; ifrag++) if (fragtypes[ifrag]==tmpfragtypes[itmpfrag]) {
            dist2=pbc_distance2(pbc,halfboxsize,boxsize,&center[3*ifrag],&transtemplates[3*index]);
            /*for (k=0; k<3; k++) {
                dx[k]=center[3*ifrag+k]-transtemplates[3*index+k];
                m=1000;
                if (pbc) {
                    if (dx[k]>halfboxsize) dx[k]-=boxsize;
                    if (dx[k]<-halfboxsize) dx[k]+=boxsize;
                }
                dist2+=dx[k]*dx[k];
            }*/
            if (dist2<cutoff2) {
                //rotate template fragment through rotation to bring it into fragme;
                multiply_quat(&tmporient[4*itmpfrag],&info[iclus].rot[0],&q1[0]);
                normalize_quat(&q1[0]);
                // angular difference to acutal orientation
                cangle=dist(&q1[0],&orient[4*ifrag]);
                if (cangle>cacutoff) {
                    e.frag=ifrag;
                    e.tmpfrag=itmpfrag;
                    e.clus=iclus;
                    e.dist2=dist2;
                    entries.push_back(e);
                }
            }
        }
    }
    qsort(&entries[0],entries.size(),sizeof(entry),compare_entries);
    for (index=0; index<nclusters*ntmpfrag; index++) assigned[index]=false;
    for (i=0; i<entries.size(); i++) {
        index=entries[i].clus*ntmpfrag+entries[i].tmpfrag;
        ifrag=entries[i].frag;
        if ((cluster_id[ifrag]<0) && !assigned[index]) {
            cluster_id[ifrag]=entries[i].clus;
            mapping[ifrag]=entries[i].tmpfrag;
            assigned[index]=true;
        }
    }
    //For each cluster, see if the template fragments that are mapped as part of that cluster
    //constitute a connected subgraph of tmpcontacts (the contact graph for the template).
    //If not, divide the subgraph into connected components, and break up the cluster accordingly.
    tmpfragments_in_clus.init(ntmpfrag);
    comps = (long int *) checkalloc(ntmpfrag,sizeof(long int));
    iclus=0;
    //for (iclus=0; iclus<nclusters; iclus++) {
    //use a while loop here because the number of clusters may be increasing as we go along.
    while (iclus<nclusters) {
        tmpfragments_in_clus.clear();
        for (ifrag=0; ifrag<nfrag; ifrag++) {
            if (cluster_id[ifrag]==iclus) tmpfragments_in_clus+=mapping[ifrag];
        }
        ncomp=tmpcontacts->connected_components(tmpfragments_in_clus,comps);
        //if this "cluster" has more than one connected component...
        if (ncomp>1) {
            //for each fragment, if it is in the correct cluster...
            for (ifrag=0; ifrag<nfrag; ifrag++) if (cluster_id[ifrag]==iclus) {
                icomp=comps[mapping[ifrag]];
                //reassign clusters as follows: if icomp is zero (the first component), leave alone
                //otherwise it belongs to one of (ncomp-1) new clusters, numbered (nclusters) thru (nclusters + ncomp - 2)
                if (icomp>0) cluster_id[ifrag]=nclusters+icomp-1;
            }
            //we just added (ncomp-1) additional clusters.
            nclusters+=ncomp-1;
        }
        iclus++;
    }
    //since the number of clusters may have changed, need to reallocate "info"
    info=(cluster_info *) checkrealloc(info,nclusters,sizeof(cluster_info));
    free(transtemplates);
    free(assigned);
    free(comps);
}


//Compare two pointers to cluster info, based on size.
//We want to sort in reverse order, so return -1 if size1>size2.
int compare_cluster_info(const void * info1, const void * info2)
{
    int size1 = ((cluster_info * ) info1)->size;
    int size2 = ((cluster_info * ) info2)->size;
    if (size1>size2) return -1;
    if (size1<size2) return 1;
    return 0; //must be equal
}

void clusters::sort_clusters(int * old_cluster_id)
{
    int iclus,ifrag;
    int * orig_to_new;
    for (iclus=0; iclus<nclusters; iclus++) info[iclus].orig_index=iclus;
    qsort(info,nclusters,sizeof(cluster_info),compare_cluster_info);
    orig_to_new=(int *) checkalloc(nclusters,sizeof(int));
    //change cluster_id's to "new" indices
    for (iclus=0; iclus<nclusters; iclus++) orig_to_new[info[iclus].orig_index]=iclus;
    for (ifrag=0; ifrag<nfrag; ifrag++) {
        if (cluster_id[ifrag]>=0) cluster_id[ifrag]=orig_to_new[cluster_id[ifrag]];
        if ((old_cluster_id!=NULL) && (old_cluster_id[ifrag]>=0)) old_cluster_id[ifrag]=orig_to_new[old_cluster_id[ifrag]];
    }
    free(orig_to_new);
}

//assume they are already sorted.
void clusters::report(long int iframe,FILE * output)
{
    int iclus,ifrag;
    //cluster_info * * sorted_clusters;
    //sorted_clusters=(cluster_info * *) checkalloc(nclusters,sizeof(cluster_info *));
    //for (iclus=0; iclus<nclusters; iclus++) sorted_clusters[iclus]=&info[iclus];
    //sort the pointers to avoid disturbing their order in the actual "info"
    //qsort(sorted_clusters,nclusters,sizeof(cluster_info *),compare_cluster_info);
    //find the min distance between
    //Can we get away with hiding the empty clusters?
    for (iclus=0; iclus<nclusters; iclus++) if (info[iclus].size>0) fprintf(output,"cluster %ld %d %d %.3f\n",iframe,iclus+1,info[iclus].size,info[iclus].rmsd);
    if (mindist!=NULL) for (ifrag=0; ifrag<nfrag; ifrag++) fprintf(output,"fragment %ld %d %d %.3f %.3f\n",iframe,ifrag,cluster_id[ifrag]+1,mindist[ifrag],angle[ifrag]*RAD_TO_DEG);
}

void clusters::identify_edge_frags(int pbc, double halfboxsize, double boxsize, double * tmpcenter, bool * is_edge)
{
    bool * largest_template;
    int ifrag,jtmpfrag,itmpfrag,count;
    int * closest;
    double * minffdist;
    double d;
    for (itmpfrag=0; itmpfrag<ntmpfrag; itmpfrag++) is_edge[itmpfrag]=false;
    //This identifies which fragments in the template correspond to the largest cluster.
    largest_template=(bool *) checkalloc(ntmpfrag,sizeof(bool));
    for (itmpfrag=0; itmpfrag<ntmpfrag; itmpfrag++) largest_template[itmpfrag]=false;
    for (ifrag=0; ifrag<nfrag; ifrag++) if (cluster_id[ifrag]==0) largest_template[mapping[ifrag]]=true;
    count=0;
    for (itmpfrag=0; itmpfrag<ntmpfrag; itmpfrag++) if (largest_template[itmpfrag]) count++;
    //now, for each fragment not in the template identify the closest fragment in the template
    closest=(int *) checkalloc(ntmpfrag,sizeof(int));
    minffdist=(double *) checkalloc(ntmpfrag,sizeof(double));
    for (itmpfrag=0; itmpfrag<ntmpfrag; itmpfrag++) {
        minffdist[itmpfrag]=1e20;
        closest[itmpfrag]=-1;
    }
    for (jtmpfrag=0; jtmpfrag<ntmpfrag; jtmpfrag++) if (largest_template[jtmpfrag])
        for (itmpfrag=0; itmpfrag<ntmpfrag; itmpfrag++) if (!largest_template[itmpfrag]) {
            d=pbc_distance2(pbc,halfboxsize,boxsize,&tmpcenter[3*itmpfrag],&tmpcenter[3*jtmpfrag]);
            if (d<minffdist[itmpfrag]) {
                minffdist[itmpfrag]=d;
                closest[itmpfrag]=jtmpfrag;
                //printf("mindist: %d %d %.4f\n",itmpfrag,jtmpfrag,sqrt(d));
            }
        }
    //mark all closest fragments so identified as being edge fragments.
    for (itmpfrag=0; itmpfrag<ntmpfrag; itmpfrag++) if (closest[itmpfrag]>=0) is_edge[closest[itmpfrag]]=true;
    //for (itmpfrag=0; itmpfrag<ntmpfrag; itmpfrag++) if (closest[itmpfrag]>=0) printf("mindist: %d %d %.4f\n",itmpfrag,closest[itmpfrag],sqrt(minffdist[itmpfrag]));
    printf("edge fragments: ");
    for (itmpfrag=0; itmpfrag<ntmpfrag; itmpfrag++) if (is_edge[itmpfrag]) printf("%d ",itmpfrag+1);
    printf("\n");
    free(minffdist);
    free(closest);
}

void clusters::get_mindist(int pbc, double halfboxsize, double boxsize, double * center, double * orient, double * tmpcenter, double * tmporient)
{
    double d;
    int iclus,ifrag,jfrag,itmpfrag,k;
    double * transtemplate_center;
    double * transtemplate_orient;
    bool * occupied;
    double rotmatrix[3][3],v[3],q[4];

    transtemplate_center=(double *) checkalloc(3*ntmpfrag,sizeof(double));
    transtemplate_orient=(double *) checkalloc(4*ntmpfrag,sizeof(double));
    occupied=(bool *) checkalloc(ntmpfrag,sizeof(bool));
    quat_to_matrix(info[0].rot,&rotmatrix[0][0]);
    for (itmpfrag=0; itmpfrag<ntmpfrag; itmpfrag++) {
        matmul(&rotmatrix[0][0],&tmpcenter[3*itmpfrag],v);
        for (k=0; k<3; k++) transtemplate_center[3*itmpfrag+k]=v[k]+info[0].disp[k];
        multiply_quat(&tmporient[4*itmpfrag],&info[0].rot[0],&transtemplate_orient[4*itmpfrag]);
        occupied[itmpfrag]=false;
    }
    //omit from the distance analysis all fragments "occupied" in maximum size cluster"
    for (ifrag=0; ifrag<nfrag; ifrag++) if (cluster_id[ifrag]==0) occupied[mapping[ifrag]]=true;
    //edge_fragments=(bool *) checkalloc(ntmpfrag,sizeof(bool));
    //identify_edge_frags(pbc,halfboxsize,boxsize,tmpcenter,edge_fragments);
    //for (iclus=0; iclus<nclusters; iclus++) info[iclus].mindist=1e20;
    mindist=(double *) checkalloc(nfrag,sizeof(double));
    angle=(double *) checkalloc(nfrag,sizeof(double));
    for (ifrag=0; ifrag<nfrag; ifrag++) {
        mindist[ifrag]=1e20;
        angle[ifrag]=M_PI;
    }
    for (ifrag=0; ifrag<nfrag; ifrag++) //if (cluster_id[ifrag]!=0)
        for (itmpfrag=0; itmpfrag<ntmpfrag; itmpfrag++) if (!occupied[itmpfrag]) { //not largest cluster
            d=pbc_distance2(pbc,halfboxsize,boxsize,&center[3*ifrag],&transtemplate_center[3*itmpfrag]);
            //if (d<info[cluster_id[jfrag]].mindist) info[cluster_id[jfrag]].mindist=d;
            if (d<mindist[ifrag]) {
                mindist[ifrag]=d;
                angle[ifrag]=dist(&orient[4*ifrag],&transtemplate_orient[4*itmpfrag]);
            }

    }
    //for (iclus=0; iclus<nclusters; iclus++) info[iclus].mindist=sqrt(info[iclus].mindist);
    for (ifrag=0; ifrag<nfrag; ifrag++) {
        mindist[ifrag]=sqrt(mindist[ifrag]);
        angle[ifrag]=2*acos(angle[ifrag]);
    }
    free(transtemplate_center);
    free(transtemplate_orient);
}

//PDB frame contains a pseudoatom for each fragment, colored by its cluster id.
void clusters::write_pdb_frame(long long int frame, double * center, FILE * output)
{
    int ifrag;
    const char * pdbatomfmt = "ATOM  %5d %4s %3s %c%4d    %8s%8s%8s%6.2f%6.2f      %-4d\n";;
    char buffer[255],xval[9],yval[9],zval[9],chain;
    fprintf(output,"REMARK frame %lld\n",frame);
    //fprintf(output,"MODEL     %4d\n",frame
    for (ifrag=0; ifrag<nfrag; ifrag++) {
        pdb_coordinate(center[3*ifrag],xval);
        pdb_coordinate(center[3*ifrag+1],yval);
        pdb_coordinate(center[3*ifrag+2],zval);
        if (cluster_id[ifrag]>=26) chain='Z'; else chain=(char)(cluster_id[ifrag]+'A');
        snprintf(buffer,sizeof(buffer),pdbatomfmt,ifrag+1,"X","XXX",chain,
            ifrag+1,xval,yval,zval,1.0,0.0,cluster_id[ifrag]);
        fputs(buffer,output);
    }
    fprintf(output,"END\n");
}

void clusters::write_xyz_frame(long long int frame, double * center, FILE * output)
{
    int ifrag;
    fprintf(output,"%d\n",nfrag);
    fprintf(output,"frame %lld\n",frame);
    //fprintf(output,"MODEL     %4d\n",frame
    for (ifrag=0; ifrag<nfrag; ifrag++) fprintf(output,"%c%d %.5f %.5f %.5f\n",(char)(cluster_id[ifrag]+'A'),
        ifrag+1,center[3*ifrag],center[3*ifrag+1],center[3*ifrag+2]);
    //fprintf(output,"END\n");
}
//Guess the clusters by constructing a
void clusters::guess_clusters(bool pbc, double boxsize, double halfboxsize, double * center, double cutoff)
{
    int ifrag,jfrag,k,inewclus;
    double dx[3],dist2,cutoff2;
    bool * close; //a matrix of pairs of fragments that are close, to be used as an incidence matrix
    //init_cluster_id=(int *) checkalloc(nfrag,sizeof(int));
    //init_mapping=(int *) checkalloc(nfrag,sizeof(int));
    /*for (ifrag=0; ifrag<nfrag; ifrag++) {
        init_cluster_id[ifrag]=-1;
        init_mapping[ifrag]=-1;
    }*/
    close=(bool *) checkalloc(nfrag*nfrag,sizeof(bool));
    cutoff2=cutoff*cutoff;
    for (ifrag=0; ifrag<nfrag*nfrag; ifrag++) close[ifrag]=false;
    for (ifrag=0; ifrag<nfrag; ifrag++) {
        close[ifrag*nfrag+ifrag]=true;
        for (jfrag=ifrag+1; jfrag<nfrag; jfrag++) {
            dist2=0.0;
            for (k=0; k<3; k++) {
                dx[k]=center[3*jfrag+k]-center[3*ifrag+k];
                if (pbc) {
                    if (dx[k]>halfboxsize) dx[k]-=boxsize;
                    if (dx[k]<-halfboxsize) dx[k]+=boxsize;
                }
                dist2+=dx[k]*dx[k];
            }
            close[ifrag*nfrag+jfrag]=(dist2<cutoff2);
            close[jfrag*nfrag+ifrag]=(dist2<cutoff2);
        }
    }
    nclusters=0;
    while (true) {
        //find first unassigned fragment
        ifrag=0;
        while ((ifrag<nfrag) && (cluster_id[ifrag]>=0)) ifrag++;
        if (ifrag>=nfrag) break; //no unassigned fragments, we're done
        nclusters++;
        inewclus=nclusters-1;
        info=(cluster_info *) checkrealloc(info,nclusters,sizeof(cluster_info));
        cluster_id[ifrag]=inewclus;
        mapping[ifrag]=0;
        info[inewclus].size=1;
        //root call for recursion
        identify_component(nfrag,ntmpfrag,ifrag,inewclus,close);
    }
    free(close);

}

//recursive depth-first search for connected components, operates on init_cluster_id, puts no more than ntmpfrag components in (by using the info[].size elements)
void clusters::identify_component(int nfrag, int ntmpfrag, int ifrag, int iclus, bool * close)
{
    int jfrag;
    for (jfrag=0; jfrag<nfrag; jfrag++) if ((cluster_id[jfrag]<0)&& close[ifrag*nfrag+jfrag]) {
        if (info[iclus].size>=ntmpfrag) return; //can't add any more
        cluster_id[jfrag]=iclus;
        mapping[jfrag]=info[iclus].size;
        info[iclus].size++;
        //recursive call
        identify_component(nfrag,ntmpfrag,jfrag,iclus,close);
    }
}

