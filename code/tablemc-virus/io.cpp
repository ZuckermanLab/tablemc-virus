#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "fragments.h"
#include "mc.h"
#include "rotations.h"
#include "mt.h"
#include "util.h"

//DCD related stuff.
#define AKMA               4.88882129E-02
#define CHARMM_VERSION     34
//Writing the random number generator state vector
const char * charmm_signature = "CORD";
//Taken from the ICNTRL array in subroutine WRITCV of dynio.src in CHARMM.
struct dcd_header {
    char signature[4]; //="CORD"
    unsigned int nframes;
    unsigned int begin; //Number of steps first frame.
    unsigned int skip; //Number of steps interval between frames.
    unsigned int nstep; //Total number of steps.
    unsigned int icntrl6to8[3];
    unsigned int ndegf; //Number of degrees of freedom.
    unsigned int nfixed; //Number of fixed atoms.
    float timestep_akma; //Timestep in AKMA units.
    unsigned int crystal; //Whether or not crystal data is written (side lengths for PBC
    unsigned int dim4; //whether or not 4-dimensinal dynamics
    unsigned int qcg; //whether or not charges written (CHEQ)
    unsigned int icntrl14to19[6];
    unsigned int version; //CHARMM version.
};

struct dcd_titles {
    unsigned int ntitle;
    char title[80];
};

//The offsets are one less than indicated on the PDB standard, because
//C arrays are zero-based.
const char * origin = "Produced by Tabulated Monte Carlo code";

//Fortran writes the record length as a 32-bit integer before and after writing each section of data.
//For DCDs we must emulate this behavior exactly.
void fortran_fwrite(const void * ptr, size_t size, size_t count, FILE * stream)
{
    unsigned int nbytes;
    nbytes=size*count;
    fwrite(&nbytes,sizeof(nbytes),1,stream);
    fwrite(ptr,size,count,stream);
    fwrite(&nbytes,sizeof(nbytes),1,stream);
}


void simulation::write_dcd_header(FILE * dcdfile)
{
    int i;
    dcd_header hdr;
    dcd_titles titles;
    strncpy(hdr.signature,charmm_signature,sizeof(hdr.signature));
    hdr.nframes = (nmcstep/nsave_xyz);
    hdr.begin = nprevstep + nsave_xyz; //I think this is correct for multiple files, but not absolutely sure
    hdr.skip = nsave_xyz;
    hdr.nstep = nmcstep;
    for (i=0; i<3; i++) hdr.icntrl6to8[i]=0;
    hdr.ndegf = 6*nfrag-6; //not necessarily exact degrees of freedom
    hdr.nfixed = 0;
    hdr.timestep_akma = 0.001f/AKMA; //Fake time step of 1 fs per MC step.
    hdr.crystal = pbc; //maybe should set equal to "pbc" and write crystal data?
    hdr.dim4 = 0;
    hdr.qcg = 0;
    for (i=0; i<6; i++) hdr.icntrl14to19[i]=0;
    hdr.version = CHARMM_VERSION;
    titles.ntitle = 1;
    for (i=0; i<sizeof(titles.title); i++) titles.title[i]=' ';
    strncpy(titles.title,origin,strlen(origin));
    fortran_fwrite(&hdr,sizeof(hdr),1,dcdfile);
    fortran_fwrite(&titles,sizeof(titles),1,dcdfile);
    fortran_fwrite(&natom,sizeof(natom),1,dcdfile);
    fflush(dcdfile);
}

void simulation::write_dcd_frame(FILE * dcdfile, double * coords)
{
    int i;
    float * xwrite;
    float * ywrite;
    float * zwrite;
    double xtlabc[6];
    xwrite=(float *) checkalloc(natom,sizeof(float));
    ywrite=(float *) checkalloc(natom,sizeof(float));
    zwrite=(float *) checkalloc(natom,sizeof(float));
    if (pbc) {
        for (i=0; i<6; i++) xtlabc[i]=0.0;
        //first, third, and sixth element are the diagonal elements of the matrix
        xtlabc[0]=boxsize;
        xtlabc[2]=boxsize;
        xtlabc[5]=boxsize;
        fortran_fwrite(xtlabc,sizeof(double),6,dcdfile);
    }
    for (i=0; i<natom; i++) {
        xwrite[i]=coords[3*i];
        ywrite[i]=coords[3*i+1];
        zwrite[i]=coords[3*i+2];
    }
    fortran_fwrite(xwrite,sizeof(float),natom,dcdfile);
    fortran_fwrite(ywrite,sizeof(float),natom,dcdfile);
    fortran_fwrite(zwrite,sizeof(float),natom,dcdfile);
    fflush(dcdfile);
    free(xwrite);
    free(ywrite);
    free(zwrite);
    //fflush(dcdfile);
}


void simulation::expand_trajectory(char * datfname, char * dcdfname, long int nframes, int freq, double _boxsize)
{
    FILE * input;
    FILE * output;
    double x[3],q[4];
    long int istep,iframe;
    int ifrag,ifragx,k;
    input=fopen(datfname,"r");
    if (input==NULL) {
        printf("Could not open input file %s.\n",datfname);
        die();
    }
    output=fopen(dcdfname,"wb");
    if (output==NULL) {
        printf("Could not open output file %s.\n",dcdfname);
        die();
    }
    printf("Expanding trajectory every %d frames from %s to %s.\n",freq,datfname,dcdfname);
    //set up parameters needed for DCD header -- read first line to d
    nprevstep=0;
    boxsize=_boxsize;
    halfboxsize=0.5*boxsize;
    pbc=(boxsize>0);
    fscanf(input,"%ld %d %lf %lf %lf %lf %lf %lf %lf\n",&istep,&ifragx,&x[0],&x[1],&x[2],&q[0],&q[1],&q[2],&q[3]);
    nsave_xyz=istep*freq; //assume the first step count  is the # of steps between frames
    nmcstep=istep*nframes;
    write_dcd_header(output);
    rewind(input);
    iframe=0;
    while (!feof(input)) {
        //use the "new" coordinates to preserve any "old" coordinates.
        read_frame_quat(input,&istep,newcenter,neworient);
        if ((istep%nsave_xyz)==0) { //every nth frame
            for (ifrag=0; ifrag<nfrag; ifrag++) update_coords(ifrag,newcenter,neworient,newcoords);
            write_dcd_frame(output,newcoords);
            printf("Frame %d step %ld written\n",iframe,istep);
        }
        iframe++;
        if (iframe>nframes) break;
    }
    //ensure consistent state after trajectory reading, probably doesn't matter
    for (ifrag=0; ifrag<nfrag; ifrag++) copy_frag(ifrag,oldcenter,oldorient,oldcoords,newcenter,neworient,newcoords);
    fclose(input);
    fclose(output);
}



//allows rereading "center-quaternion" files to obtain their energies.  More exact than reading a pdb file.
//this assumes each frame is written as a block of lines equal to the number of fragments
void simulation::read_frame_quat(FILE * input, long int * istep, double * center, double * orient)
{
    int ifrag,ifragx,k;
    double x[3],q[4];
    char fragname[MAX_FRAGMENT_NAME];
    for (ifragx=0; ifragx<nfrag; ifragx++){
        fscanf(input,"%ld %d %s %lg %lg %lg %lg %lg %lg %lg\n",istep,&ifrag,fragname,&x[0],&x[1],&x[2],&q[0],&q[1],&q[2],&q[3]);
        ifrag=ifrag-1;
        if ((ifrag>=nfrag) || (strncasecmp(fragname,fragtypes[frags[ifrag].type]->fragname,MAX_FRAGMENT_NAME)!=0)) {
            printf("error reading trajectory file\n");
            die();
        }
        normalize_quat(q);
        for (k=0; k<3; k++) center[3*ifrag+k]=x[k];
        for (k=0; k<4; k++) orient[4*ifrag+k]=q[k];
        normalize_quat(&orient[4*ifrag]);
    }
}

void simulation::write_frame_xyz(FILE * output, long int istep, double * coords)
{
    char buf[32];
    int ifrag,i,iatom;
    fragmenttype * frag;
    fprintf(output,"%d\n",natom);
    fprintf(output,"step %ld\n",istep);
    for (ifrag=0; ifrag<nfrag; ifrag++) {
        frag=fragtypes[frags[ifrag].type];
        for (i=0; i<frag->natom; i++) {
            iatom=frags[ifrag].start+i;
            snprintf(buf,sizeof(buf),"%s_%d",frag->atoms[i].name,ifrag+1);
            fprintf(output, "%s %.10f %.10f %.10f\n",buf,coords[3*iatom],coords[3*iatom+1],coords[3*iatom+2]);
        }
    }
    fflush(output);
}

void pdb_coordinate(double val, char * str)
{
     //*str=' ';
     if ((val>10000) || (val<-1000)) sprintf(str,"%8.1f",val); 
     else if ((val>1000) || (val<-100)) sprintf(str,"%8.2f",val); 
     else sprintf(str,"%8.3f",val);
     //if (fabs(val)>10000) sprintf(str+1,"%7.1f",val); else if (fabs(val)>1000) sprintf(str+1,"%7.2f",val); else sprintf(str+1,"%7.3f",val);
}

void atom_number(int iatom, char * str)
{
     if (iatom<=99999) sprintf(str,"%5d",iatom); else strncpy(str,"*****",5);
}

void simulation::write_frame_pdb(FILE * output, long int istep, double * coords)
{
    char buf[81],xval[9],yval[9],zval[9],atomnum[6];//room for terminating null
    int ifrag,i,iatom;
    fragmenttype * frag;
    const char * pdbcrysfmt = "CRYST1%9.3%9.3%9.3%7.2%7.2%7.2\n";
    //const char * pdbatomfmt = "ATOM  %5d %4s %3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f      %-4d\n";
    const char * pdbatomfmt = "ATOM  %5s %4s %3s %c%4d    %8s%8s%8s%6.2f%6.2f       %-4d\n";
    /*fprintf(output,"%d\n",natom);
    fprintf(output,"step %ld\n",istep);*/
    for (ifrag=0; ifrag<nfrag; ifrag++) {
        frag=fragtypes[frags[ifrag].type];
        for (i=0; i<frag->natom; i++) {
            iatom=frags[ifrag].start+i;
            //snprintf(buf,sizeof(buf),"%s_%d",frag->names[i],ifrag+1);
            pdb_coordinate(coords[3*iatom],xval);
            pdb_coordinate(coords[3*iatom+1],yval);
            pdb_coordinate(coords[3*iatom+2],zval);
            atom_number(iatom+1,atomnum);
            snprintf(buf,sizeof(buf),pdbatomfmt,atomnum,frag->atoms[i].name,frag->atoms[i].resname,fragtypes[frags[ifrag].type]->fragname[0],
                        frag->atoms[i].res,xval,yval,zval,1.0,0.0,ifrag+1);
            fputs(buf,output);
        }
    }
    fprintf(output,"END\n");
    fflush(output);
}
/*void simulation::write_pair_pdb(FILE * output, int ifrag, int jfrag, double * coords)
{
    char buf[32];
    int i,iatom,ii;
    fragmenttype * frag;
    const char * pdbcrysfmt = "CRYST1%9.3%9.3%9.3%7.2%7.2%7.2\n";
    const char * pdbatomfmt = "ATOM  %6d%4s %3s %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n";
    //fprintf(output,"%d\n",natom);
    //fprintf(output,"step %ld\n",istep);
    //for (ifrag=0; ifrag<nfrag; ifrag++) {
    ii=1;
        frag=frags[fragtypes[ifrag]];
        for (i=0; i<frag->natom; i++) {
            iatom=fragstart[ifrag]+i;
            //snprintf(buf,sizeof(buf),"%s_%d",frag->names[i],ifrag+1);
            fprintf(output,pdbatomfmt,ii,frag->atoms[i].name,fragtypes[ifrag],1,coords[3*iatom],coords[3*iatom+1],coords[3*iatom+2],1.0,0.0);
            ii++;
        }
    //}
        frag=frags[fragtypes[jfrag]];
        for (i=0; i<frag->natom; i++) {
            iatom=fragstart[jfrag]+i;
            //snprintf(buf,sizeof(buf),"%s_%d",frag->names[i],jfrag+1);
            fprintf(output,pdbatomfmt,ii,frag->atoms[i].name,fragtypes[jfrag],2,coords[3*iatom],coords[3*iatom+1],coords[3*iatom+2],1.0,0.0);
            ii++;
        }

    fprintf(output,"END\n");
    fflush(output);
}*/



void simulation::write_frame_quat(FILE * output, long int istep, double * center, double * orient)
{
    int ifrag;
    for (ifrag=0; ifrag<nfrag; ifrag++) fprintf(output,"%ld %d %s %.10f %.10f %.10f %.10f %.10f %.10f %.10f\n",
        istep,ifrag+1,fragtypes[frags[ifrag].type]->fragname,center[3*ifrag],center[3*ifrag+1],center[3*ifrag+2],orient[4*ifrag],orient[4*ifrag+1],orient[4*ifrag+2],orient[4*ifrag+3]);
    fflush(output);
}


//The MTRIXn records in the PDB file give the "crystallographic asymmetric unit", not the "biological assembly".
//The "biological assembly" is given by REMARK 350  BIOMT records.
#define PDB_MATRIX "REMARK 350   BIOMT" //to character 18 in the file
#define PDB_ATOM   "ATOM"
#define PDB_END    "END"
//we may need a special "initatoms" array to give atom info for the PDB file, compared with the system as constructed
void simulation::read_pdb_file(char * fname)
{
    FILE * input;
    char buf[255],buf2[255],buf3[255];
    double matrix[3][3];
    initatoms=NULL;
    trans=NULL;
    initcoords=NULL;
    mass=NULL;
    ninitatom=0;
    ntrans=0;
    input = fopen(fname,"r");
    if (input==NULL) {
        printf("Could not open file %s\n",fname);
        die();
    }
    printf("Reading coordinates and transformations from PDB file %s.\n",fname);
    while (!feof(input)) {
        fgets(buf,sizeof(buf),input);
        if (strncasecmp(buf,PDB_MATRIX,18)==0) {
            //It is a transformation.  Read two more lines.
            fgets(buf2,sizeof(buf2),input);
            fgets(buf3,sizeof(buf3),input);
            if ((strncasecmp(buf2,PDB_MATRIX,18)!=0) || (strncasecmp(buf2,PDB_MATRIX,18)!=0)) {
                printf("Error reading transformation.\n");
                die();
            }
            trans=(transinfo *) checkrealloc(trans,ntrans+1,sizeof(transinfo));
            //read the rotation matrix, convert to a quaternion and store
            sscanf(buf+24,"%lf%lf%lf",&matrix[0][0],&matrix[0][1],&matrix[0][2]);
            sscanf(buf2+24,"%lf%lf%lf",&matrix[1][0],&matrix[1][1],&matrix[1][2]);
            sscanf(buf3+24,"%lf%lf%lf",&matrix[2][0],&matrix[2][1],&matrix[2][2]);
            matrix_to_quat(matrix,&trans[ntrans].q[0]);
            //read the trnaslation
            sscanf(buf+59,"%lf",&trans[ntrans].center[0]);
            sscanf(buf2+59,"%lf",&trans[ntrans].center[1]);
            sscanf(buf3+59,"%lf",&trans[ntrans].center[2]);
            ntrans++;
        } else if (strncasecmp(buf,PDB_ATOM,4)==0) { //it's an atom
            initatoms=(atominfo *) checkrealloc(initatoms,ninitatom+1,sizeof(atominfo));
            initcoords=(double *) checkrealloc(initcoords,3*(ninitatom+1),sizeof(double));
            mass=(double *) checkrealloc(mass,ninitatom+1,sizeof(double));
            read_pdb_line(buf,&initatoms[ninitatom],&initcoords[3*ninitatom],&mass[ninitatom]);
            ninitatom++;
        } else if (strncasecmp(buf,PDB_END,3)==0) {
            break; //end of file;
        } else continue; //something else
    }
    printf("Total %d atoms and %d transformations in initial structure.\n",ninitatom,ntrans);
}

void simulation::read_restart(char * fname)
{
    int ifrag,i,itype,inewfrag,natomx,nfragx,k;
    double x[3],q[4];
    char name[MAX_FRAGMENT_NAME];
    FILE * f;
    f=fopen(fname,"r");
    if (f==NULL) {
        printf("Could not open restart file %s\n",fname);
        die();
    }
    printf("Reading restart file %s.\n",fname);
    fscanf(f,"%d %d %ld\n",&natomx,&nfragx,&nprevstep);
    /*Read the fragment, center, and quaternion.*/
    /*if (nprevstep>0) {
        printf("Reading RNG state vector from file %s.\n",fname);
        read_rng_state(f);
        reseed=false;
    }*/
    while (!feof(f)) {
    //for (i=0;i<nfrag;i++){
        fscanf(f,"%d %s %lg %lg %lg %lg %lg %lg %lg\n",&ifrag,name,&x[0],&x[1],&x[2],&q[0],&q[1],&q[2],&q[3]);
        itype=frag_type_by_name(name);
        if (itype<0) {
            printf("Unrecognized type in restart file %s.\n",fname);
            die();
        }
        inewfrag=insert_fragment(itype);
        ifrag=ifrag-1;
        if (ifrag!=inewfrag) {
            printf("Error reading restart file %s.\n",fname);
            die();
        }
        normalize_quat(q);
        for (k=0; k<3; k++) oldcenter[3*ifrag+k]=x[k];
        for (k=0; k<4; k++) oldorient[4*ifrag+k]=q[k];
        normalize_quat(&oldorient[4*ifrag]);
        //update_coords(ifrag,oldcenter,oldorient,oldcoords);
        //copy_frag(ifrag,oldcenter,oldorient,oldcoords,newcenter,neworient,newcoords);
    }
    fclose(f);
}

void simulation::write_restart(long int istep, char * fname)
{
    FILE * output;
    int ifrag;
    output=fopen(fname,"w");
    if (output==NULL) {
        printf("Failed to open restart file %s.\n",fname);
        die();
    }
    printf("Writing restart file %s.\n",fname);
    fprintf(output,"%d %d %ld\n",natom,nfrag,istep);
    //if (istep>0) write_rng_state(output);
    //changed from "new" to "old" coordinates so that we can write restart files during setup phase.
    for (ifrag=0; ifrag<nfrag; ifrag++) fprintf(output,"%d %s %.10f %.10f %.10f %.10f %.10f %.10f %.10f\n",
        ifrag+1,fragtypes[frags[ifrag].type]->fragname,oldcenter[3*ifrag],oldcenter[3*ifrag+1],oldcenter[3*ifrag+2],
            oldorient[4*ifrag],oldorient[4*ifrag+1],oldorient[4*ifrag+2],oldorient[4*ifrag+3]);
    fclose(output);
}
