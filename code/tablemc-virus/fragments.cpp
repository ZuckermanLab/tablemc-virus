#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h>
#include "fragments.h"
#include "rotations.h"
#include "util.h"
//i/o function for reading a pdb line
void read_pdb_line(const char * buf, atominfo * atom, double * coords, double * mass)
{
    int i;
    char el;
    double weight;
    sscanf(buf+12,"%4s",atom->name);
    atom->name[4]='\0';
    sscanf(buf+22,"%3d",&atom->res);
    sscanf(buf+17,"%3s",atom->resname);
    sscanf(buf+60,"%lf",&weight);
    atom->resname[3]='\0';
    atom->chain=buf[21];
    atom->type=0;
    atom->is_site=(weight!=0);
    sscanf(buf+30,"%lf%lf%lf",&coords[0],&coords[1],&coords[2]);
    i=0;
    el=' ';
    while (isspace(el) && (i<4)) {
        el=atom->name[i];
        i++;
    }
    *mass = atom_mass(el);
}

fragmenttype::fragmenttype(const char * name, const char * fname)
{
    FILE * f;
    double xx,yy,zz,q,inertia[3][3],axes[3][3],cross[3],det,temp[3];
    char buf[255];
    char aname[4];
    char el;
    int i,j,k,junk2,minres,maxres;
    printf("Reading coordinates for fragment type %s from file %s.\n",name,fname);
    //printf("*** %d\n",strlen(fname));
    strncpy(fragname,name,sizeof(fragname));
    strncpy(fragfname,fname,sizeof(fragfname));
    f=fopen(fname,"r");
    if (f==NULL) {
        printf("Could not open fragment file %s\n",fname);
        die();
    }
    atoms=NULL;
    refgeom=NULL;
    mass=NULL;
    totalmass=0;
    natom=0;
    nsite=0;
    xx=0;
    yy=0;
    zz=0;
    startres=100000000;
    endres=0;
    /*it is a pdb file*/
    while (!feof(f)) {
        fgets(buf,sizeof(buf),f);
        if (feof(f)) break;
        if (strncasecmp("END",buf,3)==0) break;
        if (strncasecmp("ATOM  ",buf,6)!=0) continue;
        atoms=(atominfo *) checkrealloc(atoms,natom+1,sizeof(atominfo));
        refgeom=(double *)checkrealloc(refgeom,3*(natom+1),sizeof(double));
        mass=(double *)checkrealloc(mass,natom+1,sizeof(double));
        read_pdb_line(buf,&atoms[natom],&refgeom[3*natom],&mass[natom]);
        atoms[natom].fragment=-1;
        if (atoms[natom].res<startres) startres=atoms[natom].res;
        if (atoms[natom].res>endres) endres=atoms[natom].res;
        xx+=refgeom[3*natom];
        yy+=refgeom[3*natom+1];
        zz+=refgeom[3*natom+2];
        totalmass+=mass[natom];
        if (atoms[natom].is_site) nsite++;
        natom++;
    }
    fclose(f);
    /*place barycenter at origin*/
    xx/= natom;
    yy/= natom;
    zz/= natom;
    for (i=0; i<natom; i++) {
        refgeom[3*i]-=xx;
        refgeom[3*i+1]-=yy;
        refgeom[3*i+2]-=zz;
    }
    //compute moment of inertia and align with axes
    for (i=0; i<3; i++) for (j=0; j<3; j++) inertia[i][j]=0.0;
    for (i=0; i<natom; i++) {
        xx=refgeom[3*i];
        yy=refgeom[3*i+1];
        zz=refgeom[3*i+2];
        inertia[0][0]+=mass[i]*(yy*yy+zz*zz);
        inertia[0][1]+=-mass[i]*xx*yy;
        inertia[0][2]+=-mass[i]*xx*zz;
        inertia[1][1]+=mass[i]*(xx*xx+zz*zz);
        inertia[1][2]+=-mass[i]*yy*zz;
        inertia[2][2]+=mass[i]*(xx*xx+yy*yy);
    }
    inertia[1][0]=inertia[0][1];
    inertia[2][0]=inertia[0][2];
    inertia[2][1]=inertia[1][2];
    jacobi(3,&inertia[0][0],&axes[0][0]);
    //the axes are {axes[0][0], axes[1][0], axes[2][0]}, etc.
    //make sure it's a rotation matrix (det axes = +1)
    cross[0]=axes[1][1]*axes[2][2]-axes[2][1]*axes[1][2];
    cross[1]=axes[2][1]*axes[0][2]-axes[0][1]*axes[2][2];
    cross[2]=axes[0][1]*axes[1][2]-axes[1][1]*axes[0][2];
    det=axes[0][0]*cross[0]+axes[1][0]*cross[1]+axes[2][0]*cross[2];
    if (det<0) for (i=0; i<3; i++) axes[i][2]*=-1.0;
    //this actually multiplies by the transpose of axes,
    for (i=0; i<natom; i++) {
        matmul(&axes[0][0],&refgeom[3*i],&temp[0]);
        for (k=0; k<3; k++) refgeom[3*i+k]=temp[k];
    }

    /*qtot=0.0;
    for (k=0; k<3; k++) dipole[k]=0.0;
    for (i=0; i<natom; i++) {
        q=ffield->chargeParams[types[i]];
        qtot+=q;
        for (k=0; k<3; k++) dipole[k]+=q*(refgeom[3*i+k]);
    }
    dipolemag=0.0;
    for (k=0; k<3; k++) dipolemag+=dipole[k]*dipole[k];
    dipolemag=sqrt(dipolemag);
    printf("Total charge for fragment %s: %g\n",fname,qtot);
    printf("Dipole moment for fragment %s: %.2f %.2f %.2f\n",fname,dipole[0],dipole[1],dipole[2]);
    printf("Magnitude of dipole:           %.2f\n",dipolemag);*/
    printf("Fragment type %s read. Total %d atoms and %d sites.\n",fragname,natom,nsite);
}

fragmenttype::~fragmenttype()
{
    free(atoms);
    free(refgeom);
    free(mass);
}

void fragmenttype::get_coords(double * center, double * orient, double * coords)
{
    double x[3],xx[3],rotmatrix[3][3];
    int i,k;
    quat_to_matrix(&orient[0],&rotmatrix[0][0]);
    for (i=0; i<natom; i++) {
        /*copy reference coordinates from fragment atom i */
        for (k=0; k<3; k++) x[k]=refgeom[3*i+k];
        matmul(&rotmatrix[0][0],x,xx);
        coords[3*i]=xx[0]+center[0];
        coords[3*i+1]=xx[1]+center[1];
        coords[3*i+2]=xx[2]+center[2];
    }
}


//returns the RMSD
void fragmenttype::fit_fragment(double * coords, double * center, double * orient, double * rmsd)
{
    fit_fragment(coords,center,orient,mass,rmsd);
}

void fragmenttype::fit_fragment(double * coords, double * center, double * orient, double * weights, double * rmsd)
{
    //I think it will be ok if only 1 or 2 atoms.  Will check this later.
    //Do the fit using only heavy atoms.
    int i;
    rmsd_fit(natom,weights,refgeom,coords,center,orient,rmsd);
}

int get_fragtype_by_name(char * name, int nfragtypes, fragmenttype * * fragtypes)
{
    int itype;
    for (itype=0; itype<nfragtypes; itype++)
        if (strncasecmp(name,fragtypes[itype]->fragname,sizeof(fragtypes[itype]->fragname))==0) return itype;
    return -1; //not found
}

//same, but for fragment file names.  Used by the table constructor
//changed to "strstr" to permit being in another directory
//changed further to compare only the base names
int get_fragtype_by_fname(char * fname, int nfragtypes, fragmenttype * * fragtypes)
{
    int itype;
    char fnamecopy[255];
    char fragfnamecopy[255];
    for (itype=0; itype<nfragtypes; itype++) {
        strncpy(fnamecopy,fname,sizeof(fnamecopy));
        strncpy(fragfnamecopy,fragtypes[itype]->fragfname,sizeof(fragfnamecopy));
        //if (strncasecmp(fname,fragtypes[itype]->fragfname,sizeof(fragtypes[itype]->fragfname))==0) return itype;
        if (strcmp(basename(fnamecopy),basename(fragfnamecopy))==0) return itype;
    }
    return -1; //not found
}

void write_pair_pdb2(FILE * output, fragmenttype * frag1, fragmenttype * frag2, double * coords1, double * coords2)
{
    char buf[32];
    int i,iatom,ii;
    const char * pdbcrysfmt = "CRYST1%9.3%9.3%9.3%7.2%7.2%7.2\n";
    const char * pdbatomfmt = "ATOM  %6d%4s F%02d %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n";
    //fprintf(output,"%d\n",natom);
    //fprintf(output,"step %ld\n",istep);
    //for (ifrag=0; ifrag<nfrag; ifrag++) {
    ii=1;
        //frag=frags[fragtypes[ifrag]];
        for (i=0; i<frag1->natom; i++) {
            //iatom=fragstart[ifrag]+i;
            //snprintf(buf,sizeof(buf),"%s_%d",frag->names[i],ifrag+1);
            fprintf(output,pdbatomfmt,ii,frag1->atoms[i].name,0,1,coords1[3*i],coords1[3*i+1],coords1[3*i+2],1.0,0.0);
            ii++;
        }
    //}
        //frag=frags[fragtypes[jfrag]];
        for (i=0; i<frag2->natom; i++) {
            //iatom=fragstart[jfrag]+i;
            //snprintf(buf,sizeof(buf),"%s_%d",frag->names[i],jfrag+1);
            fprintf(output,pdbatomfmt,ii,frag2->atoms[i].name,1,2,coords2[3*i],coords2[3*i+1],coords2[3*i+2],1.0,0.0);
            ii++;
        }

    fprintf(output,"END\n");
    fflush(output);
}

