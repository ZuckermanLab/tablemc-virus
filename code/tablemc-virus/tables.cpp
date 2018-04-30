#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "tables.h"
#include <math.h>
#include "fragments.h"
#include "rotations.h"
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


table::table(const char * fname, int new_table)//, int newtable)
{
    FILE * f;
    memset(&hdr,0,sizeof(hdr));
    energy=NULL;
    radial_table=NULL;
#ifndef NO_MMAP_TABLES
    map=NULL;
#endif
    hdr.totalpoints=0;
    if (new_table) {
        generate_table(fname,0,1);
    } else {
        f=fopen(fname,"rb");
        if (f==NULL) { //still couldn't open?
            printf("Could not open table file %s.\n",fname);
            die();
        }
        alloc_read_table(f,fname);
        fclose(f);
    }
    one_minr=1.0/hdr.minr;
    one_logrfactor=1.0/hdr.logrfactor;
    one_dsph=1.0/hdr.dsph;
    one_deuler=1.0/hdr.deuler;
}



table::~table(void)
{
#ifndef NO_MMAP_TABLES
    //get rid of the map, do not free energy in the usual manner
    if (map!=NULL) {
        munmap(map,mapsize);
    } else {
#endif
    free(energy);
#ifndef NO_MMAP_TABLES
    }
#endif
    if (radial_table!=NULL) free(radial_table);
    //free(clash_table);
}

//returns table size in megabytes
double table::getsize(void)
{
    return ((double) hdr.totalpoints*sizeof(energy_t))/((double) 1024*1024);
}

bool table::verify_checksum(void)
{
    unsigned int current_crc32;
    current_crc32=digital_crc32((unsigned char *) energy,hdr.totalpoints*sizeof(energy_t));
    return (current_crc32==hdr.crc32_data);
}


void table::print_header_info(void)
{
    char buffer[255];
    double mem;
    printf("Minimum and maximum radius:                       %.2f-%.2f A\n",hdr.minr,hdr.maxr);
    //printf("Radial, angular, orientational resolutions:  %.2f A %.2f deg %.2f deg\n",hdr.dr,RAD_TO_DEG*hdr.dsph,RAD_TO_DEG*hdr.deuler);
    printf("Radial resolution:                                %.2f A\n",hdr.dr);
    printf("Translational angular resolution:                 %.2f degrees\n",hdr.dsph*RAD_TO_DEG);
    printf("Orientational resolution:                         %.2f degrees\n",hdr.deuler*RAD_TO_DEG);
    printf("Number of radial points:                          %ld\n",hdr.nr);
    printf("Number of angular theta/phi points:               %ld %ld\n",hdr.nsphtheta,hdr.nsphphi);
    printf("Number of Euler theta/phi-psi points:             %ld %ld\n",hdr.ntheta,hdr.nphipsi);
    printf("Total translational points:                       %ld\n",hdr.ntrans);
    printf("Total orientational points:                       %ld\n",hdr.norient);
    printf("Total points:                                     %lld\n",hdr.totalpoints);
    printf("Table size:                                       %.2f MB\n",getsize());
    //printf("Minimum and maximum clash distances:         %.2f A- %.2f A\n",hdr.min_clash,hdr.max_clash);
    //printf("Clash factor:                                %.2f\n",hdr.clashfactor);
    printf("Radial factor:                                    %.2f\n",hdr.rfactor);
    //printf("Energy clash criterion:                           %.2f kcal/mol\n",hdr.en_clash);
    if (hdr.en_invalid_margin>0.0) {
        printf("Clash points:                                     %lld\n",hdr.clashpoints);
        printf("Invalid points:                                   %lld (%.2f%%)\n",hdr.invalidpoints,100*((double)hdr.invalidpoints/(double)hdr.totalpoints));
        printf("Energy invalidity margin:                         %.2f kcal/mol\n",hdr.en_invalid_margin);
        printf("Number of deviating points for invalidity:        %d\n",hdr.n_invalid_dev);
    } else {
        printf("No invalid cells.\n");
    }
    printf("Reference fragment file:                          %s\n",hdr.ref_frag_fname);
    printf("Other fragment file:                              %s\n",hdr.other_frag_fname);
    printf("Contact map file:                                 %s\n",hdr.go_model_map_fname);
    printf("Go model total contacts:                          %d\n",hdr.go_model_entries);
    printf("Go model native contacts:                         %d\n",hdr.go_model_native);
    print_go_params(hdr.go_params);
    /*printf("Go model hardcore and cutoff distances:           %.2f %.2f A\n",hdr.go_hardcore,hdr.go_cutoff);
    printf("Go model delta:                                   %.2f A\n",hdr.go_delta);
    printf("Go model native and nonnative energies:           %.2f %.2f kcal/mol\n",hdr.native_energy,hdr.nonnative_energy);*/
    strncpy(buffer,ctime(&hdr.created),sizeof(buffer));

    if ((hdr.options & OP_BOLTZMANN_AVERAGED)!=0) {
        printf("Boltzmann averaging was used in the construction of this table.\n");
        //printf("Number of cells averaged (r, anguler, Euler):     %d %d %d\n",hdr.nsmoothr,hdr.nsmoothsph,hdr.nsmootheuler);
        printf("Angular scale for angular averaging:              %.2f degrees\n",SPHERICAL_DIFFUSION_CONST*RAD_TO_DEG*sqrt(hdr.sph_alpha));
        printf("Angular scale for orientational averaging:        %.2f degrees\n",SO3_DIFFUSION_CONST*RAD_TO_DEG*sqrt(hdr.euler_alpha));
        printf("Temperature for Boltzmann averaging:              %.2f K\n",1/(KBOLTZ*hdr.beta));
        //printf("Energy limit for smoothing:                       %.2f kcal/mol\n",hdr.en_smooth_limit);
        printf("Energy limit for smoothing disabled this version.\n");
        //printf("Seed used:                                   %ld\n",hdr.seed);
    }
    /*if ((hdr.options & OP_DIST_DEP_DIELECTRIC)!=0) {
        printf("This table used a distance dependent dielectric.\n");
    }*/
    if ((hdr.options & OP_SINGLE_PRECISION)!=0) {
        printf("This table contains single precision numbers.\n");
    }
    if ((hdr.options & OP_DECIMATED)!=0) {
        printf("This table has been decimated from a larger table.\n");
    }
    if ((hdr.options & OP_INCOMPLETE)!=0) {
        printf("This table is incomplete.\n");
        printf("Starting and ending indices:                      %lld %lld\n",hdr.startindex,hdr.endindex);
    }
    printf("This table file created on:                       %s\n",buffer);
    printf("CRC-32 of table data:                             %.8lx\n",hdr.crc32_data);
    printf("Created by version:                               %.1f\n",hdr.version);
    fflush(stdout);
}



//each index zero-based
/*inline long long int table::clash_table_index(long int isphtheta, long int isphphi, long int iphi, long int itheta, long int ipsi)
{
    long int itrans,iorient;
    //itrans=hdr.nx*hdr.nx*ix+hdr.nx*iy+iz;
    itrans=(hdr.nsphphi*isphtheta)+isphphi;
    iorient=hdr.nphipsi*hdr.ntheta*iphi+hdr.nphipsi*itheta+ipsi;
    //printf("*** %ld %ld\n",itrans,iorient);
    return hdr.nsphtheta*hdr.nsphphi*((long long int) iorient)+itrans;
}*/

//Format a file name and open a table file for reading.
//Also sets the "reffragtype" and "otherfragtype" public members.
table::table(const char * fmt, const char * fragname1, const char * fragname2, int nfragtype, fragmenttype * * fragtypes)
{

    FILE * f;
    char fname[255],fragfname1[255],fragfname2[255];
    energy=NULL;
    radial_table=NULL;
#ifndef NO_MMAP_TABLES
    map=NULL;
#endif
    hdr.totalpoints=0;
    //First compose the file name (try both ways) and see which one exists and can be opened.
    snprintf(fname,sizeof(fname),fmt,fragname1,fragname2);
    strlower(fname);
    trim_string(fname);
    f=fopen(fname,"rb");
    if (f==NULL) { //try the other way
        snprintf(fname,sizeof(fname),fmt,fragname2,fragname1);
        strlower(fname);
        trim_string(fname);
        f=fopen(fname,"rb");
        if (f==NULL) { //still couldn't open?
            printf("Could not open table file %s.\n",fname);
            die();
        }
    }
    alloc_read_table(f,fname);
    fclose(f);
    //ref_frag_fname is the "reference" fragment.  Determine the fragment type number.
    /*snprintf(fragfname1,sizeof(fragfname1),fragfmt,fragname1);
    snprintf(fragfname2,sizeof(fragfname2),fragfmt,fragname2);
    reffragtype=-1; otherfragtype=-1;
    if (strcasecmp(fragfname1,hdr.ref_frag_fname)==0) {
        reffragtype=top->fragtypebyname(fragname1);
    } else if (strcasecmp(fragfname2,hdr.ref_frag_fname)==0) {
        reffragtype=top->fragtypebyname(fragname2);
    }
    if (strcasecmp(fragfname1,hdr.other_frag_fname)==0) {
        otherfragtype=top->fragtypebyname(fragname1);
    } else if (strcasecmp(fragfname2,hdr.other_frag_fname)==0) {
        otherfragtype=top->fragtypebyname(fragname2);
    }*/
    reffragtype=get_fragtype_by_fname(hdr.ref_frag_fname,nfragtype,fragtypes);
    otherfragtype=get_fragtype_by_fname(hdr.other_frag_fname,nfragtype,fragtypes);
    if ((reffragtype<0) || (otherfragtype<0)) {
        printf("Table fragment types in header do not match.\n");
        die();
    }
    one_minr=1.0/hdr.minr;
    one_logrfactor=1.0/hdr.logrfactor;
    one_dsph=1.0/hdr.dsph;
    one_deuler=1.0/hdr.deuler;

}

//in the event we are reading an incomplete table, the indices will be shifted by "startindex".  We account for this when combining tables.
void table::alloc_read_table(FILE * f, const char * fname)
{
    int ir;
    long long int actual_totalpoints;
#ifndef NO_MMAP_TABLES
    size_t filesize, pagesize;
#endif
    printf("Loading table file %s:\n",fname);
    //printf("*** %d\n",sizeof(*hdr));
    fread(&hdr,sizeof(hdr),1,f);
    if (ferror(f)) {
        printf("Error reading table file %s\n",fname);
        die();
    }
    if ((hdr.headersize!=sizeof(hdr)) || (hdr.version < VERSION)) {
        printf("File format error in table file %s\n",fname);
        die();
    }
    if ((sizeof(float)==sizeof(energy_t)) && !(hdr.options & OP_SINGLE_PRECISION)) {
        printf("This table uses double precision numbers, but program is compiled for single precision tables.\n");
        die();
    }
    if ((sizeof(double)==sizeof(energy_t)) && !(hdr.options & OP_SINGLE_PRECISION)) {
        printf("This table uses single precision numbers, but program is compiled for double precision tables.\n");
        die();
    }
    print_header_info();
    if (hdr.options & OP_INCOMPLETE) actual_totalpoints=hdr.endindex-hdr.startindex+1; else actual_totalpoints=hdr.totalpoints;
#ifndef NO_MMAP_TABLES

    //Determine the size of the whole file, and round up to the next page.
    filesize=hdr.headersize+actual_totalpoints*sizeof(energy_t);
    pagesize=sysconf(_SC_PAGE_SIZE);
    mapsize=filesize/pagesize;
    if ((filesize%pagesize)>0) mapsize++;
    mapsize*=pagesize;
    //map the file --   void * mmap(void *start, size_t length, int prot , int flags, int fd, off_t offset);
    //read only. MAP_POPULATE causes the entire table to be brought into memory.
    map=mmap(NULL, mapsize, PROT_READ, MAP_SHARED | MAP_POPULATE , fileno(f), 0);
    if (map==MAP_FAILED) {
        printf("Memory map failed.\n");
        die();
    }
    //map now points to the beginning of the table file.  Must advance to get to the beginning of the data.
    //this is ugly, since we cannot perform pointer arithmetic on void pointers
    energy=(energy_t *) ((char *) map + hdr.headersize);
#else
    /*clash_table = (double *) malloc(hdr.clashpoints*sizeof(double));
    if (clash_table == NULL) {
        printf("Could not allocate memory for clash table in file %s\n",fname);
        die();
    }*/
    energy = (energy_t *) malloc(actual_totalpoints*sizeof(energy_t)); //generate
    if (energy == NULL) {
        printf("Could not allocate memory for table in file %s\n",fname);
        die();
    }
    if (ferror(f)) {
        printf("Error reading table file %s\n",fname);
        die();
    }
    //fread(clash_table,sizeof(double),hdr.clashpoints,f);
    fread(energy,sizeof(energy_t),actual_totalpoints,f);
#endif
//Do not verify the CRC by default.  It is a time-consuming operation and not necessary on every WE segment.
#ifdef CRC
    if (!verify_checksum()) {
        printf("Checksum failure for table %s\n",fname);
        die();
    }
#else
    printf("Checksum not verified.\n");
#endif
    //fclose(f);
}

void table::write_table(const char * fname)
{
    FILE * f;
    char pmffname[255];
    f=fopen(fname,"wb");
    if (f==NULL) {
        printf("Could not open table file %s for writing\n",fname);
        die();
    }
    fwrite(&hdr,sizeof(hdr),1,f);
    //fwrite(clash_table,sizeof(double),hdr.clashpoints,f);
    if (hdr.options & OP_INCOMPLETE) {
        fwrite(energy+hdr.startindex,sizeof(energy_t),hdr.endindex-hdr.startindex+1,f);
    } else {
        fwrite(energy,sizeof(energy_t),hdr.totalpoints,f);
    }
    fclose(f);
    printf("Table written to file %s.\n",fname);
}


/*double table::get_clash_radius(double sphtheta, double sphphi, double phi, double theta, double psi)
{
    long int isphtheta,isphphi,iphi,itheta,ipsi;
    long long index;
    isphtheta=(long int) (sphtheta/hdr.dsph + 0.5);
    isphphi=((long int) (sphphi/hdr.dsph + 0.5)) % hdr.nsphphi;
    itheta=(long int) (theta/hdr.deuler + 0.5);
    iphi=((long int) (phi/hdr.deuler+0.5)) % hdr.nphipsi;
    ipsi=((long int) (psi/hdr.deuler+0.5)) % hdr.nphipsi;
    index=clash_table_index(isphtheta,isphphi,iphi,itheta,ipsi);
    return clash_table[index];
}

double table::interpolate_clash_radius(double sphtheta, double sphphi, double phi, double theta, double psi)
{
    long int ind[5],ix2[5];
    long long int index;
    int i,j,bit;
    double frac[6],sum,prod,e;
    const int mask[5] = {0x2, 0x1, 0x10, 0x8, 0x4};
    //if (r>hdr.maxr) return 0.0;
    frac[0]=sphtheta/hdr.dsph;
    frac[1]=sphphi/hdr.dsph;
    frac[2]=phi/hdr.deuler;
    frac[3]=theta/hdr.deuler;
    frac[4]=psi/hdr.deuler;

    for (i=0; i<5; i++) {
        ind[i]=(long int) frac[i];
        frac[i]=frac[i]-ind[i];
    }
    if (ind[0]>=(hdr.nsphtheta-1)) {
        ind[0]=hdr.nsphtheta-2;
        frac[0]=1.0;
    }
    if (ind[3]>=(hdr.ntheta-1)) {
        ind[3]=hdr.ntheta-2;
        frac[3]=1.0;
    }
    sum=0.0;
    for (i=0; i<32; i++) {
        prod=1.0;
        for (j=0; j<5; j++) { //so that low order bits in i control low order indices, optimizing memory access
            if ((i & mask[j])!=0) bit=1; else bit=0;
            ix2[j]=ind[j]+bit;
            prod *= (bit ? frac[j] : (1-frac[j]));
        }
        //phi and psi are periodic, so we reduce their indices modulo
        ix2[1]=ix2[1]%hdr.nsphphi;
        ix2[2]=ix2[2]%hdr.nphipsi;
        ix2[4]=ix2[4]%hdr.nphipsi;
        index=clash_table_index(ix2[0],ix2[1],ix2[2],ix2[3],ix2[4]);
        e=clash_table[index];
        sum+=prod*e;
        //printf("*** %d %ld %ld %ld %ld %ld %ld %lld %.4f %.4f %.4f\n",i,ix2[0],ix2[1],ix2[2],ix2[3],ix2[4],ix2[5],index,prod,e,sum);
    }
    return sum;
}*/
//Primarily to get a separate timing for the table lookup itself.
double table::get_energy(int enwrite, double r, double sphtheta, double sphphi, double phi, double theta, double psi)
{
    long int ir,isphtheta,isphphi,iphi,itheta,ipsi;
    long int irmin,irmax,irmid;
    long long int index;
    int i,boundscheck;
    double u;
    energy_t * p;
    //More efficient to multiply than divide.
    //if (r<hdr.min_clash) return DUMMY_ENERGY;
    isphtheta=(long int) (sphtheta*one_dsph);
    if (isphtheta>=hdr.nsphtheta) isphtheta=hdr.nsphtheta-1;
    isphphi=((long int) (sphphi*one_dsph + 0.5));
    if (isphphi>=hdr.nsphphi) isphphi-=hdr.nsphphi;
    itheta=(long int) (theta*one_deuler);
    if (itheta>=hdr.ntheta) itheta=hdr.ntheta-1;
    iphi=((long int) (phi*one_deuler+0.5));
    if (iphi>=hdr.nphipsi) iphi-=hdr.nphipsi;
    ipsi=((long int) (psi*one_deuler+0.5));
    if (ipsi>=hdr.nphipsi) ipsi-=hdr.nphipsi;
    /*index=calculate_index(hdr,0,isphtheta,isphphi,iphi,itheta,ipsi);
    rclash=table[index];
    if (r<*rclash) return DUMMY_ENERGY;*/
    //if (enwrite) {
        /*index=calculate_index(hdr,0,isphtheta,isphphi,iphi,itheta,ipsi);
        rclash=table[index];*/
    //    fprintf(energy_output,"r, rclash: %.4f %.4f\n",r,*rclash);
    //printf("indices: %d %d %d %d %d %d\n",ir,isphtheta,isphphi,iphi,itheta,ipsi);
    //printf(energy_output,"indices: %d %d %d %d %d %d\n",ir,isphtheta,isphphi,iphi,itheta,ipsi);
    /*Note that table has coarse resolution when actual clash distance is much higher than min_clash.*/
    u=log(r*one_minr)*one_logrfactor;
    ir=(long int) (u+0.5);
    //if (ir>=hdr.nr) ir=hdr.nr-1;
    index=calculate_index(ir,isphtheta,isphphi,iphi,itheta,ipsi);
    //if (table[index]<=INVALID_ENERGY) printf("hit invalid energy cell.\n");
#ifdef TIMERS
    switch_timer(TIMER_INT_LOOKUP);
#endif
    return energy[index];
	//return 0.0;
}

//This ordering of the math operations seems to be optimal in terms of performance.

double table::table_interaction_energy(int enwrite, int interp, double r2, double * rij,  double * qref,  double * qother, double * rdiff)
{
    double r,disp[3],qij[4],e,qconj[4],rotmatrix[3][3],rclash;
    int typelo, typehi;
    double phi,theta,psi,sphtheta, sphphi,tmp,edip;
    long int ix[3],iphi,itheta,ipsi,i;
    long long int index;
    r=sqrt(r2);
    if (r<hdr.minr) return DUMMY_ENERGY; /* Prevent collisions and protect entry 0 for use as a clash distance*/
    conjugate_quat(qref,&qconj[0]);
#ifdef TIMERS
    switch_timer(TIMER_INT_ORIENT);
#endif
    multiply_quat(qother,qconj,&qij[0]);
    if (qij[0]<0) { //Otherwise, we can on very rare occasions get numerical errors from quat_to_euler in the last decimal place, shifting the table entry.
	qij[0]=-qij[0];
	qij[1]=-qij[1];
	qij[2]=-qij[2];
	qij[3]=-qij[3];
    }
#ifdef TIMERS
    quat_to_euler(qij,&phi,&theta,&psi);
    switch_timer(TIMER_INT_TRANS);
    rotate_vector_by_quat(qconj,rij,disp);
#else
    //other order faster?
    rotate_vector_by_quat(qconj,rij,disp);
    quat_to_euler(qij,&phi,&theta,&psi);
#endif
    cart_to_sph(disp,r,&sphtheta,&sphphi);
#ifdef TIMERS
    switch_timer(TIMER_INT_INDEX);
#endif
    e=get_energy(enwrite,r,sphtheta,sphphi,phi,theta,psi);
#ifdef TIMERS
    switch_timer(TIMER_INT_OTHER);
#endif
    return e;
}


void table::get_random_cell(double * r, double * sphtheta, double * sphphi, double * phi, double * theta, double * psi)
{
    long int ir,isphtheta,isphphi,iphi,itheta,ipsi;
            /*Choose a random position/orientation FROM THE TABLE and set fragment 1 to it.*/
        /*Make sure we stay within the cutoff, otherwise both energies are 0 and the test is pointless.*/
        ir=(long int) (genrand_real3()*(hdr.nr));

        isphtheta=(long int) (genrand_real3()*hdr.nsphtheta);
        isphphi=(long int) (genrand_real3()*hdr.nsphphi);
        //r=ir*hdr->dr;
        *r=hdr.minr*pow(hdr.rfactor,ir);
        *sphtheta=(isphtheta+0.5)*hdr.dsph;
        *sphphi=isphphi*hdr.dsph;

        iphi=(long int) (genrand_real3()*(hdr.nphipsi));
        itheta=(long int) (genrand_real3()*(hdr.ntheta));
        ipsi=(long int) (genrand_real3()*(hdr.nphipsi));
        *phi=iphi*hdr.deuler;
        *theta=(itheta+0.5)*hdr.deuler;
        *psi=ipsi*hdr.deuler;
        printf("indices: %ld %ld %ld %ld %ld %ld\n",ir,isphtheta,isphphi,iphi,itheta,ipsi);
}

void table::volume_test(int ntest)
{
    double q[4],phi,theta,psi,volume,expected,ratio,chisq,z,sumvolume;
    long int iphi,itheta,ipsi,iorient,itest,i,ndegf,ntest2;
    int * counts;
    counts = (int *) checkalloc(hdr.norient,sizeof(int));
    for (i=0; i<hdr.norient; i++) counts[i]=0;
    for (itest=1; itest<=ntest; itest++) {
        rand_unif_quat(q);
        quat_to_euler(q,&phi,&theta,&psi);
        itheta=(long int) (theta/hdr.deuler);
        if (theta>=M_PI) itheta=hdr.ntheta-1;
        iphi=((long int) (phi/hdr.deuler+0.5));
        if (iphi>=hdr.nphipsi) iphi-=hdr.nphipsi;
        ipsi=((long int) (psi/hdr.deuler+0.5));
        if (ipsi>=hdr.nphipsi) ipsi-=hdr.nphipsi;
        iorient=(hdr.ntheta*iphi+itheta)*hdr.nphipsi+ipsi;
        counts[iorient]++;
    }
    chisq=0.0;
    ntest2=0;
    sumvolume=0.0;
    for (iphi=0; iphi<hdr.nphipsi; iphi++)
    for (itheta=0; itheta<hdr.ntheta; itheta++)
    for (ipsi=0; ipsi<hdr.nphipsi; ipsi++) {

        //phi=iphi*hdr.deuler;
        theta=(itheta+0.5)*hdr.deuler;
        //psi=ipsi*hdr.deuler;
        volume=sin(theta)*hdr.deuler*hdr.deuler*hdr.deuler;
        sumvolume+=volume;
    }
    for (iphi=0; iphi<hdr.nphipsi; iphi++)
    for (itheta=0; itheta<hdr.ntheta; itheta++)
    for (ipsi=0; ipsi<hdr.nphipsi; ipsi++) {
        iorient=(hdr.ntheta*iphi+itheta)*hdr.nphipsi+ipsi;
        //phi=iphi*hdr.deuler;
        theta=(itheta+0.5)*hdr.deuler;
        //psi=ipsi*hdr.deuler;
        volume=sin(theta)*hdr.deuler*hdr.deuler*hdr.deuler;
        //p=(((double) counts[iorient])/((double)ntest));
        expected=ntest*(volume/sumvolume);
        printf("iphi, itheta, ipsi = %ld %ld %ld, expected %.2f, actual %d\n",iphi,itheta,ipsi,expected,counts[iorient]);
        chisq+=((counts[iorient]-expected)*(counts[iorient]-expected))/expected;
        ntest2+=counts[iorient];
    }
    ndegf=hdr.norient-1;
    z=(chisq-ndegf)/sqrt(2*(double) ndegf); //normal approximzation to chi square distribution
    printf("chi square, degrees of freedom, z = %.4f %ld %.4f\n",chisq,ndegf,z);
}


