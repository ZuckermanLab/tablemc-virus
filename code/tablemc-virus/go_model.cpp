#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "go_model.h"
#include "util.h"
#include "fragments.h"
#include "rotations.h" //for pbc_distance2

//#define ALPHA_CARBON "CA"
#define CLASH_ENERGY 1e20

//i/o support routines
//format:  go_hardcore go_cutoff m n native_energy nonnative_energy
//m and n must be even and positive, with m>n
void read_go_params(char * line, go_model_params * params)
{
    sscanf(line,"%lg %lg %lg %d %d %lg %lg",&params->hardcore,&params->cutoff,&params->subcutoff,&params->m,&params->n,&params->native_energy,&params->nonnative_energy);
    if ((params->m%2!=0) || (params->m<params->n) || (params->n<0) || (params->n%2!=0)) {
        printf("Invalid set of Go parameters.\n");
        die();
    }
    params->hardcore2=4*params->hardcore*params->hardcore; //hardcore distance = 2 * hardcore radius
    params->cutoff2=params->cutoff*params->cutoff;
    params->ratio=((double) params->m)/(params->n);
    params->scaled_native_en=-params->native_energy/(1-params->ratio);
    params->scaled_nonnative_en=-params->nonnative_energy/(1-params->ratio);
    params->hn=params->n/2;
    params->hm=params->m/2;
    params->rsubcutoff=1/params->subcutoff;
}

//this can be called from the setup or from table::print_header_info
void print_go_params(go_model_params params)
{
    printf("Go model hardcore and cutoff distances:           %.2f %.2f A\n",params.hardcore,params.cutoff);
    printf("Go model subcutoff (units of native distance):    %.2f\n",params.subcutoff);
    //printf("Go model delta:                                   %.2f A\n",hdr.go_delta);
    printf("Go model exponents:                               %d %d\n",params.m,params.n);
    printf("Go model native and nonnative energies:           %.2f %.2f kcal/mol\n",params.native_energy,params.nonnative_energy);
    if ((params.m%2!=0) || (params.m<params.n) || (params.n<0) || (params.n%2!=0)) {
        printf("Invalid set of Go parameters.\n");
        die();
    }
}
//for passing to "qsort" and "bsearch"
int compare_entries(const void * _entry1, const void * _entry2)
{
    go_model_entry * entry1 = (go_model_entry *) _entry1;
    go_model_entry * entry2 = (go_model_entry *) _entry2;
    if (entry1->ires<entry2->ires) {
        return -1;
    } else if (entry1->ires>entry2->ires) {
        return 1;
    } else if (entry1->jres<entry2->jres) {
        return -1;
    } else if (entry1->jres>entry2->jres) {
        return 1;
    } else return 0;
}


go_model_info::go_model_info()
{
    nentries=0;
    nnative=0;
    //entries=NULL;
    distances=NULL;
    ifragname[0]='\0';
    jfragname[0]='\0';
}



go_model_info::~go_model_info()
{
    //if (entries!=NULL) free(entries);
    if (distances!=NULL) free(distances);
}

/*int go_model_info::find_entry(int ires, int jres) {
    go_model_entry entry;
    go_model_entry * p;
    //int low, high, mid,comp;
    entry.ires=ires;
    entry.jres=jres;
    if (entries==NULL) return -1; //no entries? not there
    p=(go_model_entry *) bsearch(&entry,entries,nentries,sizeof(go_model_entry),compare_entries);
    if (p==NULL) return -1; else return p-entries; //index
}*/

//utility routine to convert fragment atoms to residues and find the appropriate position in the distance table
int go_model_info::get_index(fragmenttype * ifragtype, int ifragatom, fragmenttype * jfragtype, int jfragatom)
{
    int ires,jres,swap,index;
    ires=ifragtype->atoms[ifragatom].res;
    jres=jfragtype->atoms[jfragatom].res;
    /*if (jres>ires) {//ensure symmetry by having only one entry
        swap=ires;
        ires=jres;
        jres=swap;
    }*/
    index=(inres*(jres-jfragtype->startres))+(ires-ifragtype->startres);
    return index;
}
//to prepare lists of atoms for making Go models

//Construct the contact map for a Go model between fragments itype and jtype. (fills in only the ifragatom, jfragatom, and native_distance in the entry)
//nuniq -- number of "unique" fragments (not symmetry related) of type itype.
//we only consider the first nuniq fragments of type itype when constructing the go model, on the assumption that all the others are symemtry rlated.
//We could use the fragment centers as an aid to avoid calculating distances between fragments too far away to matter.
void go_model_info::create_contact_map(int itype, fragmenttype * ifragtype, int jtype, fragmenttype * jfragtype,  int nfrag,  fragmentinfo * fraginfo, double * nativecoords)
{
    int iatom,jatom,ifragatom,jfragatom,ifrag,jfrag,iuniq,ientry,ientry2,k,counter,inatom,jnatom,ires,jres,swap,index,index2;
    int * fraglist;
    double distance2;
    nentries=0;
    inres=ifragtype->endres-ifragtype->startres+1;
    jnres=jfragtype->endres-jfragtype->startres+1;
    distances=(double *) checkalloc(inres*jnres,sizeof(double));
    for (index=0; index<inres*jnres; index++) distances[index]=1e20;
    strncpy(ifragname,ifragtype->fragname,sizeof(ifragname));
    strncpy(jfragname,jfragtype->fragname,sizeof(jfragname));
    inatom=ifragtype->natom;
    jnatom=jfragtype->natom;
    counter=0;
    /*iuniq=0;
    fraglist=(int *) checkalloc(nuniq,sizeof(int));
    for (ifrag=0; ifrag<nfrag; ifrag++) if (fraginfo[ifrag].type==itype) {
        if (iuniq>=nuniq) break;
        fraglist[iuniq]=ifrag;
        iuniq++;
    }*/
    //entries=checkalloc(ifragtype->natom*jfragtype->natom,sizeof(go_model_entry));
    /*for (iuniq=0; iuniq<nuniq; iuniq++) {
        ifrag=fraglist[iuniq];*/
    for (ifrag=0; ifrag<nfrag; ifrag++) if (fraginfo[ifrag].type==itype)
        for (jfrag=0; jfrag<nfrag; jfrag++) if ((fraginfo[jfrag].type==jtype)&&(ifrag!=jfrag)) {
            printf("Doing fragments %d and %d\n",ifrag,jfrag);
            for (ifragatom=0; ifragatom<inatom; ifragatom++) if (ifragtype->atoms[ifragatom].is_site)
                for (jfragatom=0; jfragatom<jnatom; jfragatom++) if (jfragtype->atoms[jfragatom].is_site) {
                    //Both atoms are in fragments of the correct type, and are alpha carbons.
                    //Determine native distance.
                    iatom=fraginfo[ifrag].start+ifragatom;
                    jatom=fraginfo[jfrag].start+jfragatom;

                    distance2=0;
                    for (k=0; k<3; k++) distance2+=(nativecoords[3*iatom+k]-nativecoords[3*jatom+k])*(nativecoords[3*iatom+k]-nativecoords[3*jatom+k]);
                    //if ((((ifragatom==24) && (jfragatom==127)) || ((ifragatom==127) && (jfragatom==24))) && (distance2<2000)) printf("***: %d %d %d %d %.16f\n",ifrag,ifragatom,jfrag,jfragatom,sqrt(distance2));
                    counter++;
                    if (counter%10000==0) printf("Calculated %d distances\n",counter);
                    //See if the entry is already in the model.
                    index=get_index(ifragtype,ifragatom,jfragtype,jfragatom);
                    if (distance2<distances[index]) distances[index]=distance2;
                }
         }

    //symmetrize if symmetric -- avoid numerical errors
    if (strcasecmp(ifragname,jfragname)==0)
        for (ires=1; ires<=inres; ires++) for (jres=1; jres<=jnres; jres++) {
            index=inres*(jres-1)+(ires-1);
            index2=inres*(ires-1)+(jres-1);
            distance2=distances[index];
            if (distances[index2]<distance2) distance2=distances[index2];
            distances[index]=distance2;
            distances[index2]=distance2;
        }
    /*for (ientry=0; ientry<nentries; ientry++) {
        ires=entries[ientry].ires;
        jres=entries[ientry].jres;
        distance2=entries[ientry].native_distance2;
        ientry2=find_entry(jres,ires);
        if (entries[ientry2].native_distance2<distance2) distance2=entries[ientry2].native_distance2;
        entries[ientry].native_distance2=distance2;
        entries[ientry2].native_distance2=distance2;
    }*/
    //free(fraglist);
    //find sigma_i, and sigma_j, the nearest
}


//sets up parameters, identifies fragment types
void go_model_info::set_parameters(go_model_params * params)
{
    int ientry,ifragatom,jfragatom,ires,jres,index;
    //double go_hardcore2=4*go_hardcore*go_hardcore; //(go_hardcore + go_hardcore)^2
    //double go_cutoff2=go_cutoff*go_cutoff;
    /*double * radii1;
    double * radii2;*/
    double aux;
    nnative=0;
    nentries=inres*jnres;
    printf("Setting up Go model with hard-core radius %g A, cutoff %g A \n",params->hardcore,params->cutoff);
    for (ires=1; ires<=inres; ires++) for (jres=1; jres<=jnres; jres++) {
        index=inres*(jres-1)+(ires-1);
        if (distances[index]<params->cutoff2) nnative++;
    }
    /*radii1=(double *) checkalloc(inatom,sizeof(double));
    radii2=(double *) checkalloc(jnatom,sizeof(double));
    for (ifragatom=0; ifragatom<inatom; ifragatom++) radii1[ifragatom]=1000000;
    for (jfragatom=0; jfragatom<jnatom; jfragatom++) radii1[jfragatom]=1000000;
    for (ientry=0; ientry<nentries; ientry++) {
        if (entries[ientry].native_distance2<radii1[entries[ientry].ifragatom]) radii1[entries[ientry].ifragatom]=entries[ientry].native_distance2;
        if (entries[ientry].native_distance2<radii1[entries[ientry].jfragatom]) radii1[entries[ientry].jfragatom]=entries[ientry].native_distance2;
    }
    //each "radius" is half the minimum distance to the nearest other alpha carbon
    for (ifragatom=0; ifragatom<inatom; ifragatom++) radii1[ifragatom]=sqrt(radii1[ifragatom])/2;
    for (jfragatom=0; jfragatom<jnatom; jfragatom++) radii2[jfragatom]=sqrt(radii2[ifragatom])/2;*/
    /*for (ientry=0; ientry<nentries; ientry++) {
        //it's a native entry if distance less than cutoff
        entries[ientry].native=(entries[ientry].native_distance2<params->cutoff2);
        if (entries[ientry].native) nnative++;*/
        /*if (entries[ientry].native) {
            //native interaction, see eq. 7
            entries[ientry].low_distance2=entries[ientry].native_distance2*(1-delta)*(1-delta);
            entries[ientry].hi_distance2=entries[ientry].native_distance2*(1+delta)*(1+delta);
        } else {
            //non-native interaction, see eq. 5
            //aux=(radii1[entries[ientry].ifragatom]+radii2[entries[ientry].jfragatom])/2;
            //aux=aux*(1-delta);
            //entries[ientry].low_distance2=aux*aux;
            entries[ientry].low_distance2=go_hardcore2;
            entries[ientry].hi_distance2=go_cutoff2;
        }*/
    //}
    //identify itype and jtype

    //free(radii1);
    //free(radii2);
}

//this needs the fragmenttype objects in order to look up residue numbers.
double go_model_info::energy(int pbc, double halfboxsize, double boxsize, go_model_params * params,  fragmenttype * ifragtype, double * icoords, fragmenttype * jfragtype, double * jcoords)
{
    double en,entot,dx[3],r2;
    int ientry,k,ifragatom,jfragatom,iexp,ires,jres;
    double a,am,an,native_distance2;
    entot=0.0;
    for (ifragatom=0; ifragatom<ifragtype->natom; ifragatom++) if (ifragtype->atoms[ifragatom].is_site)
        for (jfragatom=0; jfragatom<jfragtype->natom; jfragatom++) if (jfragtype->atoms[jfragatom].is_site) {
            r2=pbc_distance2(pbc,halfboxsize,boxsize,&icoords[3*ifragatom],&jcoords[3*jfragatom]);
            if (r2<=0) {
                printf("go_model_info::energy: error, ifragatom=%d jfragatom=%d r2=%.2f\n",ifragatom,jfragatom,r2);
                die();
            }
            /*if (r2<entries[ientry].low_distance2) return CLASH_ENERGY;
            if (r2<entries[ientry].hi_distance2) { //intermediate range in equations
                if (entries[ientry].native) en+=native_energy; else en+=nonnative_energy;
            }*/
            native_distance2=distances[get_index(ifragtype,ifragatom,jfragtype,jfragatom)];
            if (native_distance2<params->cutoff2) a=native_distance2/r2; else a=params->hardcore2/r2;
            if (a<params->rsubcutoff)  continue; //<1% of native energy
            an=1;
            for (iexp=1; iexp<=params->hn; iexp++) an*=a;
            //an now equals (r/r_nat)^n or (r/r_hc)^n
            am=an;
            if (params->m==2*params->n) am=an*an; else for (; iexp<=params->hm; iexp++) am*=a;
            //am now equals (r/r_nat)^m or (r/r_hc)^m
            if (native_distance2<params->cutoff2) {
                en=params->scaled_native_en*(am-params->ratio*an);
            } else {
                en=params->scaled_nonnative_en*am;
            }
#ifdef DEBUG
            //if (en!=0) printf("Go model: %d %d %.16f %c %.16f %.16f\n",ifragatom,jfragatom,r2,yesno(entries[ientry].native),a,en);
#endif
            entot+=en;
    }
    return entot;
}

int go_model_info::count_native_contacts(int pbc, double halfboxsize, double boxsize, double cutoff, double ratio, fragmenttype * ifragtype, double * icoords, fragmenttype * jfragtype, double * jcoords)
{
    double dx[3],r2,cutoff2,ratio2;
    int ientry,k,ifragatom,jfragatom,iexp;
    double a,am,an,native_distance2;
    int ncontact;
    ncontact=0;
    cutoff2=cutoff*cutoff;
    ratio2=ratio*ratio;

    for (ifragatom=0; ifragatom<ifragtype->natom; ifragatom++)
        for (jfragatom=0; jfragatom<jfragtype->natom; jfragatom++) {
            native_distance2=distances[get_index(ifragtype,ifragatom,jfragtype,jfragatom)];
            if (native_distance2<=cutoff2) {
                r2=pbc_distance2(pbc,halfboxsize,boxsize,&icoords[3*ifragatom],&jcoords[3*jfragatom]);
                if (r2<=0) {
                    printf("go_model_info::energy: error, ifragatom=%d jfragatom=%d r2=%.2f\n",ifragatom,jfragatom,r2);
                    die();
                }
                a=r2/native_distance2;
                if (a<=ratio2) ncontact++; /*else {
                    printf("nonnative contact: %d %d %.4f %.4f %.4f\n",ifragatom,jfragatom,sqrt(r2),sqrt(native_distance2),sqrt(a));
                }*/
            }
    }
    return ncontact;
}

//this restores to the same state as after create_contact_map.  need to call set_parameters to finish the go model before use.
void go_model_info::read_file(char * fname)
{
    FILE * input;
    int ires,jres,index;
    double r;
    input=fopen(fname,"r");
    if (input==NULL) {
        printf("Could not open file %s.\n",fname);
        die();
    }
    printf("Reading Go model information from file %s.\n",fname);
    fscanf(input,"%s %s %d %d\n",ifragname,jfragname,&inres,&jnres);
    distances=(double *) checkalloc(inres*jnres,sizeof(double));
    for (index=0; index<inres*jnres; index++) distances[index]=1.0e20;
    //fscanf(input,"%d\n",&nentries);
    while (!feof(input)) {
        fscanf(input,"%d %d %lg\n",&ires,&jres,&r);
        index=inres*(jres-1)+(ires-1);
        distances[index]=r*r;
    }
    fclose(input);
}

void go_model_info::write_file(char * fname)
{
    FILE * output;
    int ires,jres,index;
    output=fopen(fname,"w");
    if (output==NULL) {
        printf("Could not open file %s.\n",fname);
        die();
    }
    fprintf(output,"%s %s %d %d\n",ifragname,jfragname,inres,jnres);
    //fprintf(output,"%d\n",nentries);
    for (ires=1; ires<=inres; ires++) for (jres=1; jres<=jnres; jres++) {
        index=inres*(jres-1)+(ires-1);
        fprintf(output,"%d %d %.16f\n",ires,jres,sqrt(distances[index]));
    }
    fclose(output);
    printf("Go model information written to file %s.\n",fname);
}

