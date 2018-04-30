#include "nblist.h"
#include "util.h"
#include "fragments.h"
#include <vector>
#include <math.h>
#include <stdlib.h>

#define NBFACTOR   10.0

/*nblist::nblist(int top->nfrag)
{
    nb_list_count=(int *) malloc(top->nfrag*sizeof(int):
    nonbond_list=NULL;
}*/
//Creates/updates nonbond list based on fragment positions. "By cubes" algorithm.
//We need to maintain two nonbond lists, one for the "new" configuration and one for the "old" configuration.
//void nblist::update_list(double listcutoff, int pbc, double boxsize, int top->nfrag, double * center)
grid::grid()
{
    boxlist=NULL;
    boxcount=NULL;
    boxneighcount=NULL;
    boxneighbors=NULL;
    //create_nonbond_list will allocate the first time through.
    //new_nonbond_list=NULL;
    nbox[0]=0; //new boxes first time
    nbox[1]=0;
    nbox[2]=0;
    nboxtot=0;
    max_per_box=0;
    volume=0.0;
    density=0.0;
}

grid::~grid()
{
    free(boxlist);
    free(boxcount);
    free(boxneighbors);
    free(boxneighcount);
}

void grid::create_grid(bool pbc, double listcutoff, int ncoords, double * coords)
{
    double xmin[3],xmax[3],dx[3];
    int noldbox[3];
    int ibox[3],jbox[3],jjbox[3];
    int icoord,k,box,box2;
    int bad,new_boxes;
    //Start by finding a bounding box for the system.
    for (k=0; k<3; k++) {
        xmin[k]=1000000;
        xmax[k]=-1000000;
    }
    for (icoord=0; icoord<ncoords; icoord++)
        for (k=0; k<3; k++) {
            if (coords[3*icoord+k]<xmin[k]) xmin[k]=coords[3*icoord+k];
            if (coords[3*icoord+k]>xmax[k]) xmax[k]=coords[3*icoord+k];
        }
    volume=1.0;
    nboxtot=1;
    new_boxes=FALSE;
    for (k=0; k<3; k++) {
        //ensure finite distance in each direction (fix bug with ace-gly-nme
        if ((xmax[k]-xmin[k])<0.1) {
            xmax[k]+=0.05;
            xmin[k]-=0.05;
        }
        volume*=(xmax[k]-xmin[k]);
        noldbox[k]=nbox[k];
        nbox[k]=((int) ((xmax[k]-xmin[k])/listcutoff));
        if (nbox[k]<1) nbox[k]=1;
        new_boxes=new_boxes || (nbox[k]!=noldbox[k]);
        dx[k]=(xmax[k]-xmin[k])/nbox[k];//this makes teh boxes a little bigger than listcutoff
        nboxtot*=nbox[k];
    }
    density=((double) ncoords)/volume;
    max_per_box=(int) (density*dx[0]*dx[1]*dx[2]*NBFACTOR)+1;
    //We only need to do the following if the number of boxes in any direction has changed.
    //They only depend on the geometrical relationship between the boxes, not the coordinates.
    if (new_boxes) {
#ifdef DEBUG
        printf("New boxes.\n");
#endif
    	boxlist = (int *) checkrealloc(boxlist,max_per_box*nboxtot,sizeof(int));
    	boxcount = (int *) checkrealloc(boxcount,nboxtot,sizeof(int));
    	boxneighcount = (int *) checkrealloc(boxneighcount,nboxtot,sizeof(int));
    	boxneighbors = (int *) checkrealloc(boxneighbors,27*nboxtot,sizeof(int));
    //}
        for (box=0; box<nboxtot; box++) boxneighcount[box]=0;
        for (ibox[0]=0; ibox[0]<nbox[0]; ibox[0]++)
        for (ibox[1]=0; ibox[1]<nbox[1]; ibox[1]++)
        for (ibox[2]=0; ibox[2]<nbox[2]; ibox[2]++) {
            box=ibox[2]*nbox[0]*nbox[1]+ibox[1]*nbox[0]+ibox[0];
            for (jbox[0]=ibox[0]-1; jbox[0]<=ibox[0]+1; jbox[0]++)
            for (jbox[1]=ibox[1]-1; jbox[1]<=ibox[1]+1; jbox[1]++)
            for (jbox[2]=ibox[2]-1; jbox[2]<=ibox[2]+1; jbox[2]++) {
                bad=FALSE;
                for (k=0; k<3; k++) {
                    jjbox[k]=jbox[k];
                    if (jjbox[k]<0) {
                        //If there are only two boxes in this direction, we need only two iterations.
                        if (pbc && (nbox[k]>=3)) jjbox[k]+=nbox[k]; else bad=TRUE;
                    }
                    if (jjbox[k]>=nbox[k]) {
                        if (pbc && (nbox[k]>=3)) jjbox[k]-=nbox[k]; else bad=TRUE;
                    }
                }
                if (bad) continue; //nonexistent cell.
                box2=jjbox[2]*nbox[0]*nbox[1]+jjbox[1]*nbox[0]+jjbox[0];
                if (box<=box2) {
                    boxneighbors[27*box+boxneighcount[box]]=box2;
                    boxneighcount[box]++;
                    if (box!=box2) {
                        boxneighbors[27*box2+boxneighcount[box2]]=box;
                        boxneighcount[box2]++;
                    }
                }
            }
        }
    }
    for (box=0; box<nboxtot; box++) boxcount[box]=0;
    //Make a list of fragments in each box.
    for (icoord=0; icoord<ncoords; icoord++) {
        for (k=0; k<3; k++) {
            ibox[k]=((int) ((coords[3*icoord+k]-xmin[k])/dx[k]));
            if (ibox[k]<0) ibox[k]=0;
            if (ibox[k]>=nbox[k]) ibox[k]=nbox[k]-1;
        }
        box=ibox[2]*nbox[0]*nbox[1]+ibox[1]*nbox[0]+ibox[0];
        boxlist[box*max_per_box+boxcount[box]]=icoord;
        boxcount[box]++;
        if (boxcount[box]>=max_per_box) {
            printf("Too many fragments in box, increase NBFACTOR\n");
            die();
        }
    }
}

fragment_nblist::fragment_nblist(int nfrag, double listcut)
{
    grd=new grid();
    nonbond_list=NULL;
    last_nb_list_center=(double *) checkalloc(3*nfrag,sizeof(double));
    nb_list_per_frag=0;
    nb_list_size=0;
    listcutoff=listcut;
    nb_list_count=(int*) checkalloc(nfrag,sizeof(int)); //for the first time
}

fragment_nblist::~fragment_nblist()
{
    delete grd;
    free(nb_list_count);
    if (nonbond_list!=NULL) free(nonbond_list);
    if (last_nb_list_center!=NULL) free(last_nb_list_center);
}

void fragment_nblist::create_nonbond_list(bool pbc, double halfboxsize, double boxsize, int nfrag, double * center)
{
    int box,box2,ifrag,jfrag,i,j,k,m;
    double listcutoff2,rij[3],r2;
    grd->create_grid(pbc,listcutoff,nfrag,center);
    //Allocate space for the nonbond list, if needed.
    nb_list_per_frag=((int) (grd->density*(4/3)*M_PI*listcutoff*listcutoff*listcutoff*NBFACTOR))+1;
    nb_list_size=nb_list_per_frag*nfrag;
    nonbond_list = (int*) checkrealloc(nonbond_list,nb_list_size,sizeof(int));
    if (nonbond_list==NULL) {
        printf("Could not allocate memory for nonbond list.\n");
        die();
    }
    for (ifrag=0; ifrag<nfrag; ifrag++) nb_list_count[ifrag]=0;
    listcutoff2=listcutoff*listcutoff;
    for (box=0; box<grd->nboxtot; box++)
        for (m=0; m<grd->boxneighcount[box]; m++) {
            box2=grd->boxneighbors[27*box+m];
            for (i=0; i<grd->boxcount[box]; i++)
                for (j=0; j<grd->boxcount[box2]; j++) {
                    ifrag=grd->boxlist[box*grd->max_per_box+i];
                    jfrag=grd->boxlist[box2*grd->max_per_box+j];
                    if (ifrag>=jfrag) continue;
                    for (k=0; k<3; k++) {
                        rij[k]=center[3*jfrag+k]-center[3*ifrag+k];
                        if (pbc) {
                            if (rij[k]>halfboxsize) rij[k]-=boxsize;
                            if (rij[k]<-halfboxsize) rij[k]+=boxsize;
                        }
                    }
                    r2=rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2];
                    //only count each fragment pair once.
                    if (r2<listcutoff2) {
                        //Add to the nonbond list, twice, so that a nonbond list is available for each fragment.
                        nonbond_list[nb_list_per_frag*ifrag+nb_list_count[ifrag]]=jfrag;
                        nb_list_count[ifrag]++;
                        nonbond_list[nb_list_per_frag*jfrag+nb_list_count[jfrag]]=ifrag;
                        nb_list_count[jfrag]++;
                        if ((nb_list_count[ifrag]>nb_list_per_frag) || (nb_list_count[jfrag]>nb_list_per_frag)) {

                            printf("Nonbond list too small, increase NBFACTOR\n");
                            die();
                        }
                    }

                }
        }
    //Record last positions at which the nonbond list was done.
    for (i=0; i<nfrag; i++) last_nb_list_center[i]=center[i];

}

//Check to see if any moved fragment has moved more than (listcutoff-cutoff)/2.
//int nblist::check_list(int istep, double cutoff, double listcutoff, int pbc, double boxsize, int top->nfrag, double * center);
bool fragment_nblist::check_nb_list(bool pbc, double halfboxsize, double boxsize, double cutoff, int nfrag, bool * moved, double * center)
{
    bool redo_nb_list;
    int ifrag,k;
    double rij[3],r2,margin,margin2;
    margin=(listcutoff-cutoff)/2.0;
    margin2=margin*margin;
    redo_nb_list=FALSE;
    if (last_nb_list_center==NULL) {
        last_nb_list_center = (double *) checkrealloc(last_nb_list_center,3*nfrag,sizeof(double));
        return true;
    }
    for (ifrag=0; ifrag<nfrag; ifrag++) if (moved[ifrag]) {
        for (k=0; k<3; k++) {
            rij[k]=center[3*ifrag+k]-last_nb_list_center[3*ifrag+k];
            if (pbc) {
                if (rij[k]>halfboxsize) rij[k]-=boxsize;
                if (rij[k]<-halfboxsize) rij[k]+=boxsize;
            }
        }
        r2=rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2];
        if (r2>margin2) {
            redo_nb_list=TRUE;
            break;
        }
    }
    return redo_nb_list;
}

bool fragment_nblist::check_nb_list(bool pbc, double halfboxsize, double boxsize, double cutoff, int nfrag, int imoved, double * center)
{
    bool redo_nb_list;
    int ifrag,k;
    double rij[3],r2,margin,margin2;
    margin=(listcutoff-cutoff)/2.0;
    margin2=margin*margin;
    redo_nb_list=FALSE;
    if (last_nb_list_center==NULL) {
        last_nb_list_center = (double *) checkrealloc(last_nb_list_center,3*nfrag,sizeof(double));
        return true;
    }

        for (k=0; k<3; k++) {
            rij[k]=center[3*imoved+k]-last_nb_list_center[3*imoved+k];
            if (pbc) {
                if (rij[k]>halfboxsize) rij[k]-=boxsize;
                if (rij[k]<-halfboxsize) rij[k]+=boxsize;
            }
        }
        r2=rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2];
        if (r2>margin2) {
            redo_nb_list=true;
        }
    return redo_nb_list;
}



