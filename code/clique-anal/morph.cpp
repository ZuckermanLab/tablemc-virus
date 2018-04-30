#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include "rotations.h"
#include "util.h"
#include "mt.h"
#include "morph.h"
#include "graph.h"

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

struct pair_info {
    int ifrag, jfrag;
    double dist;
    double rel_orient[4]; //quaternion giving relative orientation
};

void get_pairs(bool pbc, double boxsize, double halfboxsize, int nfrag, double * center, double * orient, vector<pair_info> * pairs)
{
    int ifrag, jfrag;
    //double distcutoff2 = distcutoff*distcutoff;
    pair_info pair;
    double d;
    double q1[4];
    pairs->clear();
    //this is slightly inefficient, we could make half the trips through the loop (for jfrag=ifrag+1;...) and add two pairs each time
    for (ifrag=0; ifrag<nfrag; ifrag++) for (jfrag=0; jfrag<nfrag; jfrag++) if (ifrag!=jfrag) {
        d=pbc_distance2(pbc,boxsize,halfboxsize,&center[3*ifrag],&center[3*jfrag]);
        //if (d<distcutoff2) {
            pair.ifrag=ifrag;
            pair.jfrag=jfrag;
            pair.dist=sqrt(d);
            //compute q_j q_i^-1
            conjugate_quat(&orient[4*ifrag],&q1[0]);
            multiply_quat(&orient[4*jfrag],&q1[0],&pair.rel_orient[0]);
            normalize_quat(&pair.rel_orient[0]);
            pairs->push_back(pair);
        //}
    }
}

inline long int vertex_index(int nfrag, int ntmpfrag, int ifrag, int itmpfrag)
{
    return ntmpfrag*ifrag+itmpfrag;
}

inline void vertex_to_fragments(long int ivertex, int nfrag, int ntmpfrag, int * ifrag, int * itmpfrag)
{
    *ifrag=ivertex/ntmpfrag; //integer divide
    *itmpfrag=ivertex%ntmpfrag;
}

clusters::clusters(bool pbc, double boxsize, double halfboxsize, int _nfrag, int * fragtypes, double * center, double * orient,
    int _ntmpfrag, int * tmpfragtypes, double * tmpcenter, double * tmporient, double distratio, double anglecutoff)
{
    int ifrag,iclus,itmpfrag,inewfrag;
    FILE * debug_output;
    bool done,conv,preiter;
    bool * template_mapped;
    double cosacutoff = cos(anglecutoff/2);
    //double rmsdconv,rmsddiff;
    vector<pair_info> frame_pairs;
    vector<pair_info> tmp_pairs;
    vector<clique_info> cliques;
    graph * g;
    long int iframepair,itmppair;
    long int i,ivertex,jvertex,iclique,edges;
    double a,r;
    nclusters=0;
    nfrag=_nfrag;
    ntmpfrag=_ntmpfrag;
    info=NULL;
    mindist=NULL;
    //init_cluster_info(&info[nclusters-1]);
    cluster_id=(int *) checkalloc(nfrag,sizeof(int));
    mapping=(int *) checkalloc(nfrag,sizeof(int));
    for (ifrag=0; ifrag<nfrag; ifrag++) {
        cluster_id[ifrag]=-1; //not assigned yet
        mapping[ifrag]=-1;
    }
    get_pairs(pbc,boxsize,halfboxsize,nfrag,center,orient,&frame_pairs);
    get_pairs(pbc,boxsize,halfboxsize,ntmpfrag,tmpcenter,tmporient,&tmp_pairs);
    printf("Comparing %ld fragment pairs in frame to %ld fragment pairs in template.\n",frame_pairs.size(),tmp_pairs.size());
    //Each vertex corresponds to a fragment in the frame, and a fragment in the tmeplate
    g=new graph(nfrag*ntmpfrag);
    //For each pair in the frame and for each pair in the template
    edges=0;
    for (iframepair=0; iframepair<frame_pairs.size(); iframepair++)
        for (itmppair=0; itmppair<tmp_pairs.size(); itmppair++)
            //check matching types and distances
            if ((fragtypes[frame_pairs[iframepair].ifrag]==tmpfragtypes[tmp_pairs[itmppair].ifrag]) &&
                (fragtypes[frame_pairs[iframepair].jfrag]==tmpfragtypes[tmp_pairs[itmppair].jfrag])) {
                    //check the distance ratio
                    r=frame_pairs[iframepair].dist/tmp_pairs[itmppair].dist;
                    if (r<1) r=1/r;
                    if (r<distratio) {
                        //check the relative orientations are within the angle cutoff
                        a=dist(&frame_pairs[iframepair].rel_orient[0],&tmp_pairs[itmppair].rel_orient[0]);
                        if (a>cosacutoff) {
                            //compute vertices and add the edge
                            ivertex=vertex_index(nfrag,ntmpfrag,frame_pairs[iframepair].ifrag,tmp_pairs[itmppair].ifrag);
                            jvertex=vertex_index(nfrag,ntmpfrag,frame_pairs[iframepair].jfrag,tmp_pairs[itmppair].jfrag);
                            if (ivertex!=jvertex) g->add_edge(ivertex,jvertex);
                            edges++;
                            if (edges%10000==0) printf("%ld edges added\n",edges);
                        }
                    }

                }
    printf("total %ld edges in graph\n",edges);
    //detect the cliques
#ifdef DEBUG
    debug_output=fopen("graph.txt","w");
    g->output(debug_output);
    fclose(debug_output);
#endif
    g->detect_cliques(&cliques);
    //they are sorted in order from biggest to smallest
#ifdef DEBUG
    debug_output=fopen("cliques.txt","w");
#endif
    nclusters=cliques.size();
    template_mapped=(bool *) checkalloc(ntmpfrag,sizeof(bool));
    info=(cluster_info *) checkalloc(nclusters,sizeof(cluster_info));
    for (iclique=0; iclique<cliques.size(); iclique++) {
        //no template fragments have been mapped yet for this clique
        for (itmpfrag=0; itmpfrag<ntmpfrag; itmpfrag++) template_mapped[itmpfrag]=false;
        for (i=0; i<cliques[iclique].count; i++) {
            ivertex=cliques[iclique].vertices[i];
            vertex_to_fragments(ivertex,nfrag,ntmpfrag,&ifrag,&itmpfrag);
#ifdef DEBUG
            fprintf(debug_output,"clique %d vertex %d ifrag %d itmpfrag %d\n",iclique,ivertex,ifrag,itmpfrag);
#endif
            //do not assign a fragment from the frame or templatemore than once
            if ((cluster_id[ifrag]<0) && !template_mapped[itmpfrag]) {
                cluster_id[ifrag]=iclique;
                mapping[ifrag]=itmpfrag;
                template_mapped[itmpfrag]=true;
            }
        }
    }
#ifdef DEBUG
    fclose(debug_output);
#endif
    free(template_mapped);
    //This is causing an error at the moment. Need to fix.
    delete g;
    //assign all remaining clusters
    done=false;
    while (!done) {
        inewfrag=0;
        while ((inewfrag<nfrag) && (cluster_id[inewfrag]>=0)) inewfrag++;
        done=(inewfrag>=nfrag); //all fragments have been assigned stably
        if (!done) add_new_cluster(pbc,boxsize,halfboxsize,inewfrag,fragtypes,center,orient,tmpfragtypes,tmpcenter,tmporient);
    }
    //this finds the transformations and the size of each cluster;
    find_transformations(pbc,boxsize,halfboxsize,center,orient,tmpcenter,tmporient);
    //sort clusters in order by size from largest to smallest
    sort_clusters(NULL);
    //remove empty clusters (which have been placed at the end of the list)
    while (info[nclusters-1].size==0) nclusters--;
    get_mindist(pbc,halfboxsize,boxsize,center,orient,tmpcenter,tmporient);
}

clusters::~clusters()
{
    int iclus;
    if (info!=NULL) {
        /*for (iclus=0find_transformations(pbc,boxsize,halfboxsize,center,tmpcenter);; iclus<nclusters; iclus++) {
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

//For each cluster, RMSD fit given the current assignment and mapping.
void clusters::find_transformations(bool pbc, double boxsize, double halfboxsize, double * center, double * orient, double * tmpcenter, double * tmporient)
{
    double * mapped_from_template;
    double * mapped_from_frame;
    double * weights;
    int * backmap; //gives the fragment for each "actual"
    int iclus, ifrag, itmpfrag,i,k;
    double dx[3],q[4];
    mapped_from_template=(double *) checkalloc(3*ntmpfrag,sizeof(double));
    mapped_from_frame=(double *) checkalloc(3*ntmpfrag,sizeof(double));
    weights=(double *) checkalloc(ntmpfrag,sizeof(double));
    backmap=(int *) checkalloc(ntmpfrag,sizeof(int));
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
    free(backmap);
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


