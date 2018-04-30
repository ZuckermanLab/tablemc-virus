#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include "tables.h"
#include "mc.h"
#include "mt.h"
#include <math.h>
#include "rotations.h"
#include "go_model.h"
#include "util.h"
#include "graph.h"

/*Input file format:
Number of steps, Save frequency, Print frequency, seed, temperature
Exact, PBC, box size
Probability of translations, Max. trans size
Probability of rotations, Max. orientation size
Number of fragments
For each fragment: Number, How many, Fragment file
For each pair of fragments i<j: i, j, table file*/
/*If exact, last stection reads: eps, param file*/
#if defined(PARALLEL) || defined(EXCHANGE)
    simulation::simulation(char* fname, int _mynod, int _numnod)
#else
    simulation::simulation(char * fname)
#endif
{
    FILE * f;
    unsigned long seed;
    double ptot,mctemp;
    int fragcounts[MAX_FRAGMENT_TYPES];
    char fname2[255],tablefmt[255],buffer[255],go_model_fname[255];
    char * p;
    go_model_info * go_model;
    int itype,ifrag,i,jtype,k;
    time_t now;
#if defined(PARALLEL) || defined(EXCHANGE)
    mynod=_mynod;
    numnod=_numnod;
#endif
    frag_nblist=NULL;
    for (itype=0; itype<MAX_FRAGMENT_TYPES; itype++) for (jtype=0; jtype<MAX_FRAGMENT_TYPES; jtype++) go_models[itype][jtype]=NULL;
    for (itype=0; itype<MAX_FRAGMENT_TYPES; itype++) for (jtype=0; jtype<MAX_FRAGMENT_TYPES; jtype++) tables[itype][jtype]=NULL;
    nfragtypes=0;
    nfrag=0;
    natom=0;
    ntrans=0;
    frags=NULL;
    reseed=true;
    oldcenter=NULL;
    oldorient=NULL;
    oldcoords=NULL;
    newcenter=NULL;
    neworient=NULL;
    newcoords=NULL;
    initcenter=NULL;
    for (itype=0; itype<MAX_FRAGMENT_TYPES; itype++) fragtypes[itype]=NULL;
    exact=true; //so that coordinates are actually updated and copied during initialization
    printf("Reading control file: %s\n",fname);
    f=fopen(fname,"r");
    if (f==NULL) {
        printf("FATAL ERROR: file %s is not found\n",fname);
        die();
    }
    process_commands(f);
    initcenter=(double *) checkalloc(3*nfrag,sizeof(double));
    if (mass!=NULL) free(mass);
    mass=(double *) checkalloc(nfrag,sizeof(double));
    for (ifrag=0; ifrag<nfrag; ifrag++) {
        update_coords(ifrag,oldcenter,oldorient,oldcoords);
        copy_frag(ifrag,oldcenter,oldorient,oldcoords,newcenter,neworient,newcoords);
        for (k=0; k<3; k++) initcenter[3*ifrag+k]=oldcenter[3*ifrag+k];
        mass[ifrag]=fragtypes[frags[ifrag].type]->totalmass;
    }
    fgets(fname2,sizeof(fname2),f);
    if (feof(f)) {
        printf("No simulation info. Quitting.\n");
        time(&now);    
        strncpy(buffer,ctime(&now),sizeof(buffer)); 
        printf("Finished at %s\n",buffer);      
        fflush(stdout);
        die();
    }
    //otherwise, we are doing a simulation. Initialize all state arrays starting from oldcenter/oldcoords.
    printf("Doing simulation.\n");
    sscanf(fname2,"%s %s %s\n",xyzfname,quatfname,endrestartfname);
    trim_string(xyzfname);
    /*trim_string(startrestartfname);
    trim_string(endrestartfname);*/
#ifdef EXCHANGE
    fscanf(f,"%ld %ld %ld %ld %ld %d\n",&nmcstep,&nsave_quat, &nsave_xyz,&nprint,&seed, &enwrite);
#else// EXCHANGE
    fscanf(f,"%ld %ld %ld %ld %ld %d %lg\n",&nmcstep,&nsave_quat, &nsave_xyz,&nprint,&seed, &enwrite, &mctemp);
#endif

    //if (!reseed) { //otherwise, the restart file contained the RNG state.
    //For WE, we do need to reseed the random number generator even if a restart file is read, to assure independence of trajectories.
        if (seed==0) {
/*            seed=time(NULL);
#ifdef UNIX
              seed^=getpid(); //To ensure uniqueness even if many are started at same time.
#endif*/
            seed=read_random_seed();
            printf("Initialized seed based on /dev/urandom.  Seed used is %lu\n",seed);
        } else {
            printf("Seed specified in file.  Seed used is %lu\n",seed);
        }
        init_genrand(seed);
    //}
    fscanf(f,"%d %d %lg %lg %lg\n",&exact,&pbc,&boxsize,&cutoff,&listcutoff);
    interp=INTERP_NONE;
    halfboxsize=boxsize*0.5;
    cutoff2=cutoff*cutoff;
#ifdef UMBRELLA
    fscanf(f,"%d %d %lg %lg\n",&ifragumb,&jfragumb,&kumb,&rumb);
    printf("Umbrella potential fragments:         %d %d\n",ifragumb,jfragumb);
    printf("Umbrella force constant:              %.2f kcal/mol-A^2\n",kumb);
    printf("Umbrella distance:                    %.2f A\n",rumb);
    ifragumb--; jfragumb--; //zero-based
#endif
    fscanf(f,"%lg %lg\n",&prob[MOVE_TRANS],&dtrans);
    fscanf(f,"%lg %lg\n",&prob[MOVE_ORIENT],&dorient);
    //fscanf(f,"%lg\n",&prob[MOVE_UNIF_ORIENT]);
    //prob[MOVE_UNIF_ORIENT]=0;
    ptot=0.0;
    for(i=1;i<=NUM_MOVES;i++)ptot+=prob[i];
    for(i=1;i<=NUM_MOVES;i++)prob[i]/=ptot;
    cumprob[1]=prob[1];
    for(i=2;i<=NUM_MOVES;i++)cumprob[i]=cumprob[i-1]+prob[i];
    printf("Number of monte carlo steps:      %ld\n",nmcstep);
    printf("Save frequencies (quat, coords):  %ld %ld\n",nsave_quat,nsave_xyz);
#ifndef EXCHANGE
    printf("Temperature:                      %.2f K\n",mctemp);
    beta=1/(KBOLTZ*mctemp);
#endif
    if (pbc) {
        printf("PBC is on. Box size =             %.2f A\n",boxsize);
    } else {
        printf("PBC is off.\n");
    }
    /*if (rmargin>0.0) {
        printf("Will calculate energy exactly when it is within %.2f A of clash radius in table.\n",rmargin);
    }*/
    /*switch (interp) {
        case INTERP_NONE:
            printf("No interpolation. Energy will be taken from the closest entry in the table.\n");
            break;
        case INTERP_3D:
            printf("Energy will be interpolated in all three translational directions.\n");
            break;
        case INTERP_GRAD_6D:
            printf("Energy will be approximated based on numerical gradients in all six dimensions.\n");
            break;
        case INTERP_6D:
            printf("Energy will be interpolated in six dimensions.\n");
            break;
        default:
            printf("Bad interpolation type.\n");
            die();
    }*/
    use_nb_list=(listcutoff>cutoff);
    if (use_nb_list) {
        printf("Nonbond list will be used.\n");
        printf("Spherical/list cutoff         %.2f %.2f\n",cutoff,listcutoff);
    } else {
        printf("Nonbond list will not be used.\n");
        printf("Spherical cutoff              %.2f\n",cutoff);
    }
    //printf("Dielectric constant:         %.2f\n",eps);
    printf("Translational moves:  Maximum size %.2f A       Fraction %.2f%%\n",dtrans,prob[MOVE_TRANS]*100);
    printf("Orientational moves: Maximum size %.2f degrees Fraction %.2f%%\n",dorient,prob[MOVE_ORIENT]*100);
    //printf("Uniform orientational moves:                   Fraction %.2f%%\n",prob[MOVE_UNIF_ORIENT]*100);
    dorient=dorient/180.0; //onvert from degrees to [0,1]
    /*fscanf(f,"%d\n",&nfragtypes);
    printf("There will be %d fragment types.\n",nfragtypes);
    for (itype=0; itype<nfragtypes; itype++) {
        fscanf(f,"%d %d %s\n",&i,&fragcounts[itype],fname2);
        //printf("*** %d %d %s\n",i,fragcounts[itype],fname2);
        if ((i-1)!=itype) {
            printf("You must specify fragment types in order.\n");
            die();
        }
        trim_string(fname2);
        printf("Reading fragment %d from file %s\n",itype+1,fname2);
        //read_fragment(fname2,&frags[itype]);
        frags[itype]=new fragmenttype("\0",fname2,ffield);
        printf("Fragment file read. Fragment contains %d atoms.\n",frags[itype]->natom);
        printf("There will be %d fragments of type %d.\n",fragcounts[itype],itype+1);
    }
    if (nfragtypes<=0) {
        printf("No fragment types.\n");
        die();
    }*/
    //need to load tables!
    if (!exact || enwrite) {
        //load tables
        fgets(tablefmt,sizeof(tablefmt),f);
        load_tables(tablefmt);
        read_scale_factors(f);
    } else {
        //read parameters for the go potential
        //fscanf(f,"%s %lg %lg %lg %lg %lg\n",go_model_fname,&go_hardcore,&go_cutoff,&delta,&native_energy,&nonnative_energy);
        printf("Exact simulation.\n");
        //to prevent errors, force all go potentials to use same parameters
        fgets(buffer,sizeof(buffer),f);
        read_go_params(buffer,&params);
        print_go_params(params);
        read_go_models(f,&params);

    }
#ifdef EXCHANGE
    exchange_init(f);
#endif // EXCHANGE
    fclose(f);
}

void simulation::read_go_models(FILE * f, go_model_params * params)
{
    int itype, jtype;
    go_model_info * go_model;
    char buffer[255],go_model_fname[255];
    while (true) { //Read go model descriptions from the rest of the file.
            fgets(buffer,sizeof(buffer),f);
            if (feof(f)) break;
            if (strncasecmp(buffer,"END",3)==0) break;
            /*p=strtok(buffer," "); //put a null after the end of the filename
            strncpy(go_model_fname,buffer,sizeof(go_model_fname));
            trim_string(go_model_fname);
            p+=strlen(p)+1;//advanced past the null
            read_go_params(p,&params);*/
            strncpy(go_model_fname,buffer,sizeof(go_model_fname));
            trim_string(go_model_fname);
            go_model=new go_model_info();
            go_model->read_file(go_model_fname);
            if (params!=NULL) go_model->set_parameters(params);
            printf("Go model for fragments %s and %s:\n",go_model->ifragname,go_model->jfragname);
            printf("Contact map filename:                   %s\n",go_model_fname);
            printf("Total and native Go model contacts:     %d %d\n",go_model->nentries, go_model->nnative);
            //print_go_params(params);
            itype=get_fragtype_by_name(go_model->ifragname,nfragtypes,fragtypes);
            jtype=get_fragtype_by_name(go_model->jfragname,nfragtypes,fragtypes);
            go_models[itype][jtype]=go_model;
        }
}

void simulation::read_scale_factors(FILE * f) {
    int itype,jtype;
    char ifragname[MAX_FRAGMENT_NAME],jfragname[MAX_FRAGMENT_NAME],buffer[255];
    double sf;
    while (!feof(f)) {
            fgets(buffer,sizeof(buffer),f);
            if (strncasecmp(buffer,"END",3)==0) break;
            if (buffer[0]=='#') continue;
            sscanf(buffer,"%s %s %lg\n",ifragname,jfragname,&sf);
            itype=frag_type_by_name(ifragname);
            jtype=frag_type_by_name(jfragname);
            if ((itype<0) || (jtype<0)) {
                printf("Invalid type name.\n");
                die();
            }
            scalefactors[itype][jtype]=sf;
            scalefactors[jtype][itype]=sf;//just for good measure
    }
    for (itype=0; itype<nfragtypes; itype++)
        for (jtype=0; jtype<nfragtypes; jtype++)
            if (tables[itype][jtype]!=NULL) printf("Will scale table entries for fragments %s and %s by factor %.4f.\n",
                fragtypes[itype]->fragname,fragtypes[jtype]->fragname,scalefactors[itype][jtype]);
}

void simulation::print_summary_info(void) {
    printf("Total %d atoms and %d fragments in system so far.\n",natom,nfrag);
}


void simulation::process_commands(FILE * input) {
    char command[255],command2[255],fname[255],fmt[255],fragname[MAX_FRAGMENT_NAME],fragname2[MAX_FRAGMENT_NAME];
    char * token;
    const char * delim = " \t\n";
    char chain;
    int i,ifrag,itype,jtype,nuniq,startres,endres,count,itrans,norigfrag,freq;
    double ratio,frac,exclude_radius;
    long int nframes;
    bool rand_orientations;
    FILE * output;
    go_model_info * go_model;
    while (true) {
        fgets(command,sizeof(command),input);
        if (feof(input)) break; //end of file
        printf("Processing command: %s\n",command);
        //for (i=0; i<strlen(command); i++) command[i]=toupper(command[i]);
        strncpy(command2,command,sizeof(command2));
        token=strtok(command2,delim);
        if (token==NULL) continue; //blank line
        else if (*token=='#') continue; //comment
        else if (strncmp("END",token,3)==0) return; //end of commands
        else if (strncmp("READ",token,6)==0) { //read a PDB file
            token=strtok(NULL,delim);
            strncpy(fmt,token,sizeof(fmt));
            token=strtok(NULL,delim);
            strncpy(fname,token,sizeof(fname));
            if (strncmp("PDB",fmt,3)==0) {
                read_pdb_file(fname);
            } else if (strncmp("REST",fmt,4)==0) {
#ifdef EXCHANGE
                if (strstr(fname,"%d")!=NULL) { //fill in replica number
                    strncpy(fmt,fname,sizeof(fmt));
                    snprintf(fname,sizeof(fname),fmt,myrep+1);
                }
#endif // EXCHANGE
                read_restart(fname);
                print_summary_info();
            } else {
                printf("Unrecognized format.\n");
                die();
            }
        } else if (strncmp("WRITE",token,5)==0) { //Fill in coordinates and write a pdb file.
            token=strtok(NULL,delim);
            strncpy(fmt,token,sizeof(fmt));
            token=strtok(NULL,delim);
            strncpy(fname,token,sizeof(fname));
            if (strncmp("PDB",fmt,3)==0) {
                for (ifrag=0; ifrag<nfrag; ifrag++) update_coords(ifrag,oldcenter,oldorient,oldcoords);
                output=fopen(fname,"w");
                write_frame_pdb(output,0,oldcoords);
                fclose(output);
            } else if (strncmp("REST",fmt,4)==0) {
#ifdef EXCHANGE
                if (strstr(fname,"%d")!=NULL) { //fill in replica number
                    strncpy(fmt,fname,sizeof(fmt));
                    snprintf(fname,sizeof(fname),fmt,myrep+1);
                }
#endif // EXCHANGE
                write_restart(0,fname);
            } else {
                printf("Unrecognized format.\n");
                die();
            }
        } else if (strncmp("FRAG",token,4)==0) { //load a fragment type
            token=strtok(NULL,delim);
            strncpy(fragname,token,sizeof(fragname));
            token=strtok(NULL,delim);
            strncpy(fname,token,sizeof(fname));
            if (nfragtypes>MAX_FRAGMENT_TYPES) {
                printf("too many fragment types, increase MAX_FRAGMENT_TYPES\n");
                die();
            }
            fragtypes[nfragtypes]=new fragmenttype(fragname,fname);
            nfragtypes++;
        } else if (strncmp("INSERT",token,6)==0) { //insert fragments, do not initialize coordinates
            //while (true) { //allow multiple fragments.  Maybe we want the synatax of this to be a fragment name and then a count.
            //token=strtok(NULL,delim);
                //if (token==NULL) break;
            //strncpy(fragname,token,sizeof(fragname));
            token+=strlen(token)+1; //move past the token, access remainder of string
            sscanf(token,"%s %d",fragname,&count);
            itype=frag_type_by_name(fragname);
            if (itype<0) {
                printf("Unrecognized fragment type %s.\n",fragname);
                die();
            }
            for (i=1; i<=count; i++) insert_fragment(itype);
            print_summary_info();
            //}
        } else if (strncmp("FIT",token,3)==0) { //insert a fitted fragment
            token+=strlen(token)+1; //move past the token, access remainder of string
            sscanf(token,"%s %d %d",fragname,&startres,&endres);
            itype=frag_type_by_name(fragname);
            if (itype<0) {
                printf("Unrecognized fragment type %s.\n",fragname);
                die();
            }
            ifrag=insert_fitted_fragment(itype,startres,endres);
            print_summary_info();
        } else if (strncmp("BUILD",token,5)==0) { //Apply all the transformations to build the complete capsid
            norigfrag=nfrag;
            token=strtok(NULL,delim);
            if ((token==NULL) || (strncmp("ALL",token,3)==0)){ //apply them all
                apply_transformations();
            } else { //numbers specified
                while (token!=NULL) {
                    itrans=atoi(token)-1; //zero-based
                    apply_one_transformation(norigfrag,itrans);
                    token=strtok(NULL,delim);
                }
            }
            print_summary_info();
        } else if (strncmp("IMAGES",token,6)==0) {
            token=strtok(NULL,delim);
            strncpy(fname,token,sizeof(fname));
            write_charmm_images(fname);
        } else if (strncmp("LATTICE",token,7)==0) { //Set up fragments in a lattice.
             token=strtok(NULL,delim);
             if (strncmp("RANDOM",token,6)==0) { //Whether or not to randomize orientations
                 rand_orientations=true;
             } else if (strncmp("FIXED",token,5)==0) {
                 rand_orientations=false;
             } else {
                 printf("Need to specify RANDOM or FIXED orientations for lattice.\n");
                 die();
             }
             token+=strlen(token)+1;
             sscanf(token,"%lg %lg",&boxsize,&exclude_radius); //read the box size
             halfboxsize=0.5*boxsize;
             initialize_system_lattice(rand_orientations,exclude_radius);
        } else if (strncmp("CONTACTMAP",token,10)==0) { //generate contact map for the structure
            token+=strlen(token)+1;
            sscanf(token,"%s %s %s",fragname,fragname2,fname);
            itype=frag_type_by_name(fragname);
            jtype=frag_type_by_name(fragname2);
            for (ifrag=0; ifrag<nfrag; ifrag++) update_coords(ifrag,oldcenter,oldorient,oldcoords);
            go_model=new go_model_info();
            go_model->create_contact_map(itype,fragtypes[itype],jtype,fragtypes[jtype],nfrag,frags,oldcoords);
            go_model->write_file(fname);
            //testing only
            /*go_model->set_parameters(1.7,8.0,0.2);
            for (ifrag=0; ifrag<nfrag; ifrag++) copy_frag(ifrag,oldcenter,oldorient,oldcoords,newcenter,neworient,newcoords);
            cutoff2=100000000.0;
            exact=true;
            enwrite=false;
            native_energy=-1.0;
            nonnative_energy=0.3;
            printf("Total energy: %.4f\n",total_energy());*/
            delete go_model;
        } else if (strncmp("EXPAND",token,6)==0) {
            token+=strlen(token)+1;
            sscanf(token,"%s %ld %d %lg %s",quatfname,&nframes,&freq,&boxsize,xyzfname);
            trim_string(xyzfname);
            trim_string(quatfname);
            expand_trajectory(quatfname,xyzfname,nframes,freq,boxsize);
        } else if (strncmp("RXNCOORD",token,8)==0) {
            token+=strlen(token)+1;
            xyzfname[0]='\0';
            sscanf(token,"%s %ld %d %lg %lg %lg %lg %s %s",quatfname,&nframes,&freq,&boxsize,&cutoff,&ratio,&frac,fname,xyzfname);
            //further lines: contact map information
            //trim_string(xyzfname);
            trim_string(quatfname);
            read_go_models(input,NULL);
            reaction_coordinates_traj(quatfname,fname,xyzfname,nframes,freq,boxsize,cutoff,ratio,frac);
        } else {
            printf("Unrecognized command.\n");
            die();
        }
    }
}

int simulation::frag_type_by_name(char * name)
{
    return get_fragtype_by_name(name,nfragtypes,fragtypes);
}
int simulation::find_init_atom(char chain, int res, char aname[4])
{
    int iatom,ifound;
    ifound=-1;
    for (iatom=0; iatom<ninitatom; iatom++) {
        if ((initatoms[iatom].chain==chain) && (initatoms[iatom].res==res) && (strncasecmp(initatoms[iatom].name,aname,4)==0)) {
            ifound=iatom;
            break;
        }
    }
    return ifound;
}

//returns index of new fragment
int simulation::insert_fragment(int itype)
{
    fragmenttype * fragtype;
    int inewfrag;
    fragtype=fragtypes[itype];
    frags=(fragmentinfo *) checkrealloc(frags,nfrag+1,sizeof(fragmentinfo));
    inewfrag=nfrag;
    frags[inewfrag].type=itype;
    frags[inewfrag].start=natom;
    frags[inewfrag].end=natom+fragtype->natom-1;
    //frags[inewfrag].chain[0]='\0';
    frags[inewfrag].startres=0;
    frags[inewfrag].endres=0;
    nfrag++;
    natom+=fragtype->natom;
    oldcenter=(double *) checkrealloc(oldcenter,3*nfrag,sizeof(double));
    oldorient=(double *) checkrealloc(oldorient,4*nfrag,sizeof(double));
    oldcoords=(double *) checkrealloc(oldcoords,3*natom,sizeof(double));
    newcenter=(double *) checkrealloc(newcenter,3*nfrag,sizeof(double));
    neworient=(double *) checkrealloc(neworient,4*nfrag,sizeof(double));
    newcoords=(double *) checkrealloc(newcoords,3*natom,sizeof(double));
    return inewfrag;
}


//Insert a fragment, fitting it to the atoms indicated by chain, startres, endres
//Do multiple chains according to a correspondence.  We assume we use the same residue numbers in each chain.
//This matches chains in the fragment to the corresponding chains with the same letter in the PDB file.
int simulation::insert_fitted_fragment(int itype,  int startres, int endres)
{
    int ifrag,ifragatom,ires,iactualres,iactualatom,ichain,k;
    fragmenttype * fragtype;
    char * p;
    double * fitcoords;
    double rmsd;
    char chain;
    //ifrag=nfrag; //will be the index of the new fragment
    fragtype=fragtypes[itype];
    ifrag=insert_fragment(itype); //insert the fragment
    frags[ifrag].startres=startres;
    frags[ifrag].endres=endres;
    fitcoords=NULL;
    fitcoords=(double *) checkalloc(3*fragtype->natom,sizeof(double));
    //For each atom in the fragment, locate a corresponding atom in the initial coordinates with the same name, the specified chain,
    //and a corresponding residue number (starting from "startres") and transfer the coordinates to the fitcoords array.
    for (ifragatom=0; ifragatom<fragtype->natom; ifragatom++) {
        chain=fragtype->atoms[ifragatom].chain;
        //if (strchr(frags[ifrag].chains,chain)==NULL) strncat(frags[ifrag].chain*/
        iactualres=fragtype->atoms[ifragatom].res-fragtype->startres+startres;
        if (iactualres>endres) printf("warning: too few residues specified\n"); //Should we die here?
        iactualatom=find_init_atom(chain,iactualres,fragtype->atoms[ifragatom].name);
        for (k=0; k<3; k++) fitcoords[3*ifragatom+k]=initcoords[3*iactualatom+k];
    }
    fragtype->fit_fragment(fitcoords,&oldcenter[3*ifrag],&oldorient[4*ifrag],&rmsd);
    //iactualres is hopefully the last residue (we really should use max(iactualres) over the last loop
    printf("Fitted fragment of type %s to residues %d-%d from chain %c of PDB file.  Mass-weighted RMSD = %.3f A\n",fragtype->fragname,startres,iactualres,chain,rmsd);
    free(fitcoords);
    return ifrag;
}

//Apply all the transformations from the PDB file.
//We apply all the transformations to all the fragments, regardless of any "APPLY THE FOLLOWING TO CHAINS:" instructions.
//We may want to fix this later, but
void simulation::apply_one_transformation(int norigfrag, int itrans)
{
    int ioldfrag,inewfrag,k;
    double v1[3],v2[3],q[4];
    if (trans[itrans].q[0]>0.99999999) return; //don't apply an identity transformation
    printf("Applying transformation %d from PDB file.\n",itrans+1);
    for (ioldfrag=0; ioldfrag<norigfrag; ioldfrag++) {
        inewfrag=insert_fragment(frags[ioldfrag].type); //Add another copy of this fragment. The copy is now fragment nfrag-1
        //frags[inewfrag].chain=frags[ioldfrag].chain;
        frags[inewfrag].startres=frags[ioldfrag].startres;
        frags[inewfrag].endres=frags[ioldfrag].endres;
        for (k=0; k<3; k++) v1[k]=oldcenter[3*ioldfrag+k]-trans[itrans].center[k]; //displacement of orig. fragment center from trans. center
        rotate_vector_by_quat(&trans[itrans].q[0],v1,v2);
        for (k=0; k<3; k++) oldcenter[3*inewfrag+k]=v2[k]+trans[itrans].center[k]; //rotated displacement, add to trans. center
        //The new fragment's orientation is the old fragment's orientation, multiplied by the transformation's quaternion.
        multiply_quat(&oldorient[4*ioldfrag],&trans[itrans].q[0],&oldorient[4*inewfrag]);
    }
}

void simulation::apply_transformations(void)
{
    int norigfrag,itrans;
    norigfrag=nfrag; //the number of fragments in the original system.
    //Skip the identity transformation.
    printf("Applying all %d transformations from PDB file.\n",ntrans);
    for (itrans=0; itrans<ntrans; itrans++) apply_one_transformation(norigfrag,itrans);
}

//Write all transformations in charmm image file format.
void simulation::write_charmm_images(char * fname)
{
    double angle,axis[3];
    int itrans,counter,k;
    FILE * output;
    output=fopen(fname,"w");
    if (output==NULL) {
        printf("Could not open file %s\n",fname);
        die();
    }
    printf("Writing CHARMM images to file %s.\n",fname);
    counter=1;
    fprintf(output,"* written by tablemc\n");
    fprintf(output,"*\n");
    for (itrans=0; itrans<ntrans; itrans++) {
        fprintf(output,"IMAGE %d\n",itrans+1);
        quat_to_axisangle(&trans[itrans].q[0],&angle,&axis[0]);
        if (angle<1e-8) {
            //charmm requires a normalized vector
            axis[0]=0.0;
            axis[1]=0.0;
            axis[2]=1.0;
        }
        fprintf(output,"ROTATE %.4f %.4f %.4f %.4f\n",axis[0],axis[1],axis[2],angle*RAD_TO_DEG);
        //translations not supported -- for hep b they are all zero
    }
    fprintf(output,"END\n");
    fclose(output);
}

void simulation::load_tables(const char * fmt)
{
    int num_tables, count, frags_in_use, ifrag,jfrag,itype,jtype,i;
    double total_size;
    FILE * f;
    bool * need_table;
    char fname[255];
    need_table=(bool *) malloc(nfragtypes*nfragtypes*sizeof(bool));
    num_tables=0;
    for (i=0; i<nfragtypes*nfragtypes; i++) need_table[i]=false;
    for (ifrag=0; ifrag<nfrag; ifrag++)
        for (jfrag=0; jfrag<nfrag; jfrag++)
            if ((ifrag!=jfrag) /*&& (!closefragments[ifrag*nfrag+jfrag])*/) {
                itype=frags[ifrag].type;
                jtype=frags[jfrag].type;
                need_table[itype*nfragtypes+jtype]=true;
            }
    for (itype=0; itype<nfragtypes; itype++)
        for (jtype=itype; jtype<nfragtypes; jtype++)
            if (need_table[itype*nfragtypes+jtype]) num_tables++;
    //total_size=0.0;
    /*frags_in_use=0;
    for (ifrag=0; ifrag<nfragtypes; ifrag++) if (fragtypes[ifrag]->n_used>0) frags_in_use++;
    //num_tables=frags_in_use*(frags_in_use+1)/2;
    printf("Fragment types in use: %d\n",frags_in_use);*/
    printf("Need to load a total of %d tables.\n",num_tables);
    count=0;
    for (itype=0; itype<nfragtypes; itype++)
        for (jtype=itype; jtype<nfragtypes; jtype++) {
            scalefactors[itype][jtype]=1.0;
            if (need_table[itype*nfragtypes+jtype]) {
                //Maybe this code for finding the file name should be moved to table's constructor.
                count++;
                printf("\n");
                printf("--------------------------------------------------------------------------------\n");
                printf("Loading interaction table file %d of %d.\n",count,num_tables);
                //ifrag is the lesser fragment type.  Load the table! (Constructor will try both possible names.
                tables[itype][jtype]=new table(fmt,fragtypes[itype]->fragname,fragtypes[jtype]->fragname,nfragtypes,fragtypes);
                //tables[jtype*nfragtypes+itype]=tables[itype*nfragtypes+jtype];
                //total_size+=tables[itype*nfragtypes+jtype]->getsize();
                fflush(stdout);
            }
        }
    free(need_table);
    printf("Total interaction tables loaded: %d\n",count);

    //printf("Total standard table size:     %.2f MB\n",total_size);
}

/*void simulation::initialize_system_random(void)
{
    //double q[4];
    printf("Initializing to random configuration.\n");
    int ifrag,jfrag,reffrag,otherfrag,reftype,othertype,k,done,clash;
    double rij[3],r;
    for (ifrag=0; ifrag<nfrag; ifrag++) {
        oldcenter[3*ifrag]=boxsize*(genrand_real3()-0.5);
        oldcenter[3*ifrag+1]=boxsize*(genrand_real3()-0.5);
        oldcenter[3*ifrag+2]=boxsize*(genrand_real3()-0.5);
        rand_unif_quat(&oldorient[4*ifrag]);
        update_coords(ifrag,oldcenter,oldorient,oldcoords);
        //copy_frag(ifrag,oldcenter,oldorient,oldcoords,newcenter,neworient,newcoords);
    }
    done=FALSE;
    while (!done) {
        done=TRUE;
        for (ifrag=0; ifrag<nfrag; ifrag++)
            for (jfrag=(ifrag+1); jfrag<nfrag; jfrag++) {
                    for (k=0; k<3; k++) {
                        rij[k]=oldcenter[3*jfrag+k]-oldcenter[3*ifrag+k];
                        if (pbc) {
                            if (rij[k]>halfboxsize) rij[k]-=boxsize;
                            if (rij[k]<-halfboxsize) rij[k]+=boxsize;
                        }
                    }
                    r=sqrt(rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2]);
                    reffrag=ifrag;
                    otherfrag=jfrag;
                    if (fragtypes[jfrag]<fragtypes[ifrag]) {
                        reffrag=jfrag;
                        otherfrag=ifrag;
                    }
                    /*reftype=fragtypes[reffrag];
                    othertype=fragtypes[otherfrag];
                    clash=TRUE;
                    while (clash) {
                        clash=check_clash(headers[reftype][othertype].clashfactor,frags[reftype].natom,&frags[reftype].types[0],&oldcoords[3*fragstart[reffrag]],
                            frags[othertype].natom,&frags[othertype].types[0],&oldcoords[3*fragstart[otherfrag]]);
                    //if (r<2.0*headers[fragtypes[reffrag]][fragtypes[otherfrag]].dr) {
                        if (r<3.0) {
                            oldcenter[3*jfrag]=boxsize*(genrand_real3()-0.5);
                            oldcenter[3*jfrag+1]=boxsize*(genrand_real3()-0.5);
                            oldcenter[3*jfrag+2]=boxsize*(genrand_real3()-0.5);
                            rand_unif_quat(&oldorient[4*jfrag]);
                            update_coords(jfrag,oldcenter,oldorient,oldcoords);
                            //printf("updated: %d %d\n",
                            done=FALSE;
                        }
                    //}
            }
    }
}*/

//place fragments on a lattice.
//rand_orientations -- whether or not to randomize orientations
void simulation::initialize_system_lattice(bool rand_orientations, double exclude_radius)
{
    int n,i,j,k,ifrag,idx;
    bool * map;
    double spacing,volume,radius2,x[3];
    if (rand_orientations) init_genrand(read_random_seed());
    //n=(int) (ceil(pow((double)nfrag,1.0/3.0)));
    //spacing=boxsize/(double)n;
    radius2=exclude_radius*exclude_radius;
    volume=(boxsize*boxsize*boxsize-(4.0/3.0)*M_PI*exclude_radius*exclude_radius*exclude_radius);
    volume=volume/((double) nfrag);
    spacing=pow(volume,1.0/3.0);
    n=ceil(boxsize/spacing);
    spacing*=0.95; //this ensures that fragments do not touch across the periodic boundary cell walls.
    printf("Initializing to lattice with spacing %.2f A\n",spacing);
    ifrag=0;
    //initialize to random positions on the lattice to ensure an equal mixture of different types of fragments
    //if (rand_orientations) {
        map=(bool*) checkalloc(n*n*n,sizeof(bool));
        for (i=0; i<n*n*n; i++) map[i]=false;
        //exclude the sphere of radius exclude_radius around the origin
        for (i=0; i<n; i++)
            for (j=0; j<n; j++) 
               for (k=0; k<n; k++) {
                   x[0]=spacing*i-halfboxsize;
                   x[1]=spacing*j-halfboxsize;
                   x[2]=spacing*k-halfboxsize;
                   if ((x[0]*x[0]+x[1]*x[1]+x[2]*x[2])<radius2) {
                      idx=n*n*i+n*j+k;
                      map[idx]=true;
                   }
        }
        for (ifrag=0; ifrag<nfrag; ifrag++) {
            do {
               i=(int) (genrand_real3()*n);
               j=(int) (genrand_real3()*n);
               k=(int) (genrand_real3()*n);
               idx=n*n*i+n*j+k;
            } while (map[idx]);
            map[idx]=true;
            oldcenter[3*ifrag]=spacing*i-halfboxsize;         
            oldcenter[3*ifrag+1]=spacing*j-halfboxsize;
            oldcenter[3*ifrag+2]=spacing*k-halfboxsize;
            if (rand_orientations) rand_unif_quat(&oldorient[4*ifrag]);
            else {
                oldorient[4*ifrag]=1.0;
                for (k=1; k<4; k++) oldorient[4*ifrag+k]=0.0;
            }
        }
    /*} else {
        for (i=0; i<n; i++)
            for (j=0; j<n; j++)
                for (k=0; k<n; k++) {
                    oldcenter[3*ifrag]=spacing*i-halfboxsize;
                    oldcenter[3*ifrag+1]=spacing*j-halfboxsize;
                    oldcenter[3*ifrag+2]=spacing*k-halfboxsize;
                    oldorient[4*ifrag]=1.0;
                    oldorient[4*ifrag+1]=0.0;
                    oldorient[4*ifrag+2]=0.0;
                    oldorient[4*ifrag+3]=0.0;
                    //update_coords(ifrag,oldcenter,oldorient,oldcoords);
                    //copy_frag(ifrag,oldcenter,oldorient,oldcoords,newcenter,neworient,newcoords);
                ifrag++;
                if (ifrag>=nfrag) return;
            }
    }*/
}

/*int simulation::count_native_contacts(int pbc, double halfboxsize, double boxsize, double cutoff, double ratio, double * coords, FILE * verbose_output)*/
void simulation::reaction_coordinates(int pbc, double halfboxsize, double boxsize, double cutoff, double ratio, double frac, double * center, double * coords, int * ncomponents, int * max_component, double * _mindist, long int iframe, FILE * verbose_output)
{
    int ncontacts, nthis, ifrag, jfrag, itype, jtype, inatom, jnatom, iatomstart, jatomstart,ncomp,icomp,maxcomp,id_maxcomp;
    double thisfrac,mindist,dist2;
    go_model_info * go_model;
    int * ids;
    int * counts;
    graph * g;
    go_model_params go_params;
    go_params.cutoff=cutoff;
    go_params.cutoff2=cutoff*cutoff;
    go_params.hardcore=0;
    go_params.hardcore2=0;
    //determine number of native contacts
    for (itype=0; itype<nfragtypes; itype++) for (jtype=0; jtype<nfragtypes; jtype++) if (go_models[itype][jtype]!=NULL) go_models[itype][jtype]->set_parameters(&go_params);
    //ncontacts=0;
    g=new graph(nfrag);
    for (ifrag=0; ifrag<nfrag; ifrag++)
        for (jfrag=ifrag+1; jfrag<nfrag; jfrag++) {
            //printf("%d %d\n",ifrag,jfrag);
            itype=frags[ifrag].type;
            jtype=frags[jfrag].type;
            inatom=fragtypes[itype]->natom;
            jnatom=fragtypes[jtype]->natom;
            iatomstart=frags[ifrag].start;
            jatomstart=frags[jfrag].start;
            if (go_models[itype][jtype]!=NULL) {
                go_model=go_models[itype][jtype];
                nthis=go_model->count_native_contacts(pbc,halfboxsize,boxsize,cutoff,ratio,fragtypes[itype],&coords[3*iatomstart],fragtypes[jtype],&coords[3*jatomstart]);
            } else if (go_models[jtype][itype]!=NULL) { //must switch the two fragments
                go_model=go_models[jtype][itype];
                nthis=go_model->count_native_contacts(pbc,halfboxsize,boxsize,cutoff,ratio,fragtypes[jtype],&coords[3*jatomstart],fragtypes[itype],&coords[3*iatomstart]);
            } else {
                //something's wrong
                printf("Go model not loaded for fragment types %s and %s.\n",fragtypes[itype]->fragname,fragtypes[jtype]->fragname);
                die();
            }
            /*ncontacts+=nthis;
            if (verbose_output!=NULL) fprintf(verbose_output,"%d %d %d\n",ifrag,jfrag,nthis);*/
            //fraction of native contacts
            thisfrac=((double) nthis)/((double) go_model->nnative);
            /*if (nthis>0) printf("reaction coordinate: %d %d %s %s %d %d %.4f\n",ifrag,jfrag,
                    fragtypes[itype]->fragname,fragtypes[jtype]->fragname,nthis,go_model->nnative,thisfrac);*/
            if (thisfrac>=frac) g->add_edge(ifrag,jfrag);

    }
    //determine connected components and find the largest connected component
    ids=(int *) checkalloc(nfrag,sizeof(int));
    ncomp=g->connected_components(ids);
    counts=(int *) checkalloc(ncomp,sizeof(int));
    for (icomp=0; icomp<ncomp; icomp++) counts[icomp]=0;
    for (ifrag=0; ifrag<nfrag; ifrag++) counts[ids[ifrag]]++;
    maxcomp=0;
    id_maxcomp=-1;
    for (icomp=0; icomp<ncomp; icomp++) if (counts[icomp]>maxcomp) {
        maxcomp=counts[icomp];
        id_maxcomp=icomp;
    }
    //find minimum distance between fragments in the largest cluster and fragments not in the largest cluster.
    if (ncomp>=2) {
        mindist=1e20;
        for (ifrag=0; ifrag<nfrag; ifrag++) if (ids[ifrag]==id_maxcomp)
            for (jfrag=0; jfrag<nfrag; jfrag++) if (ids[jfrag]!=id_maxcomp) {
                dist2=pbc_distance2(pbc,halfboxsize,boxsize,&center[3*ifrag],&center[3*jfrag]);
                if (dist2<mindist) mindist=dist2;
            }
    } else mindist=0;
    *ncomponents=ncomp;
    *max_component=maxcomp;
    *_mindist=sqrt(mindist);
    if (verbose_output!=NULL) {
        for (icomp=0; icomp<ncomp; icomp++) fprintf(verbose_output,"%ld %d %d\n",iframe,icomp,counts[icomp]);
        fflush(verbose_output);
    }
    //return ncontacts;
}


void simulation::reaction_coordinates_traj(char * datfname, char * outfname, char * voutname, long int nframes, int freq, double _boxsize, double cutoff, double ratio, double frac)
{
    FILE * input;
    FILE * output;
    FILE * verbose_output;
    double x[3],q[4],mindist;
    long int istep,iframe;
    int ifrag,jfrag,ifragx,k,ncontacts,numcomp,maxcomp;
    boxsize=_boxsize;
    halfboxsize=0.5*boxsize;
    pbc=(boxsize>0);
    nsave_xyz=freq;
    input=NULL;
    output=NULL;
    verbose_output=NULL;
    if (strncasecmp(datfname,"RESTART",7)==0) {
        //use coordinates previously set up from restart file.
        for (ifrag=0; ifrag<nfrag; ifrag++) update_coords(ifrag,oldcenter,oldorient,oldcoords);
        //ncontacts=count_native_contacts(pbc,halfboxsize,boxsize,cutoff,ratio,oldcoords,NULL);
        reaction_coordinates(pbc,halfboxsize,boxsize,cutoff,ratio,frac,oldcenter,oldcoords,&numcomp,&maxcomp,&mindist,0,NULL);
        printf("Number and maximum size of components: %d %d %.4f\n",numcomp,maxcomp,mindist);
        return;
    } else {
        //open trajectory file.
        input=fopen(datfname,"r");
        if (input==NULL) {
            printf("Could not open input file %s.\n",datfname);
            die();
        }
        output=fopen(outfname,"w");
        if (output==NULL) {
            printf("Could not open output file %s.\n",outfname);
            die();
        }

        printf("Calculating native contacts for trajectory every %d frames from %s to %s.\n",freq,datfname,outfname);
        if ((voutname!=NULL) && (strlen(voutname)>0)) {
            verbose_output=fopen(voutname,"w");
            printf("Verbose output privided in %s.\n",voutname);
        }
    }
    iframe=0;
    while (!feof(input)) {
        //use the "new" coordinates to preserve any "old" coordinates.
        read_frame_quat(input,&istep,newcenter,neworient);
        iframe++;
        if ((iframe%nsave_xyz)==0) { //every nth frame
            for (ifrag=0; ifrag<nfrag; ifrag++) update_coords(ifrag,newcenter,neworient,newcoords);
            //ncontacts=count_native_contacts(pbc,halfboxsize,boxsize,cutoff,ratio,newcoords,verbose_output);
            reaction_coordinates(pbc,halfboxsize,boxsize,cutoff,ratio,frac,newcenter,newcoords,&numcomp,&maxcomp,&mindist,iframe,verbose_output);
            fprintf(output,"%ld %d %d %.4f\n",istep,numcomp,maxcomp,mindist);
            printf("Frame %d step %ld written\n",iframe,istep);
            fflush(output);
            fflush(verbose_output);
        }
        if (iframe>nframes) break;
    }
    //ensure consistent state after trajectory reading, probably doesn't matter
    for (ifrag=0; ifrag<nfrag; ifrag++) copy_frag(ifrag,oldcenter,oldorient,oldcoords,newcenter,neworient,newcoords);
    fclose(input);
    fclose(output);
    if (verbose_output!=NULL) fclose(verbose_output);
}

simulation::~simulation()
{
    int itype,jtype;
    //delete ffield;
    for (itype=0; itype<nfragtypes; itype++)
        for (jtype=0; jtype<nfragtypes; jtype++)
            if (go_models[itype][jtype]!=NULL) delete go_models[itype][jtype];
    if (frags!=NULL) free(frags);
    if (initatoms!=NULL) free(initatoms);
    if (mass!=NULL) free(mass);
    if (initcoords!=NULL) free(initcoords);
    if (oldcenter!=NULL) free(oldcenter);
    if (oldorient!=NULL) free(oldorient);
    if (oldcoords!=NULL) free(oldcoords);
    if (newcenter!=NULL) free(newcenter);
    if (neworient!=NULL) free(neworient);
    if (newcoords!=NULL) free(newcoords);
    if (initcenter!=NULL) free(initcenter);
    if (frag_nblist!=NULL) delete frag_nblist;
    for (itype=0; itype<nfragtypes; itype++) if (fragtypes[itype]!=NULL) delete fragtypes[itype];
    for (itype=0; itype<nfragtypes; itype++)
        for (jtype=(itype+1); jtype<nfragtypes; jtype++)
            if (tables[itype][jtype]!=NULL) delete tables[itype][jtype];
}


void simulation::comparison_test(void)
{
    int ntest,itest,k,ok;
    long int ir,isphtheta,isphphi,iphi,itheta,ipsi;
    double r,sphtheta,sphphi,x[3],phi,theta,psi,xx[3],disp[3],q[4],newq[4],enexact,entable,rotmatrix[3][3],err;
    go_model_info * go_model;
    //table_header * hdr;
    FILE * test_output;
    char fname[255];
    init_genrand(57);
    ntest=1000;
    nfrag=2;
    pbc=FALSE;
    enwrite=FALSE;
    go_model=new go_model_info();
    go_model->read_file("data/contact-map-all-sites-cd-cd");
    read_go_params("1.7 10 2.5 12 6 1.0 1.0",&params);
    go_model->set_parameters(&params);
    go_models[0][0]=go_model;
    //hdr=&tables[fragtypes[0]][fragtypes[1]]->hdr;
    for (itest=1; itest<=ntest; itest++) {
        printf("test %d\n",itest);
        //Set fragment 0 to the standard position/orientation.
        oldcenter[0]=0.0;
        oldcenter[1]=0.0;
        oldcenter[2]=0.0;
        oldorient[0]=1.0;
        oldorient[1]=0.0;
        oldorient[2]=0.0;
        oldorient[3]=0.0;
        update_coords(0,oldcenter,oldorient,oldcoords);
        tables[frags[0].type][frags[1].type]->get_random_cell(&r,&sphtheta,&sphphi,&phi,&theta,&psi);
        sph_to_cart(r,sphtheta,sphphi,&x[0]);
        printf("values: %.4f %.4f %.4f %.4f %.4f %.4f\n",r,sphtheta*RAD_TO_DEG,sphphi*RAD_TO_DEG,phi*RAD_TO_DEG,theta*RAD_TO_DEG,psi*RAD_TO_DEG);
        oldcenter[3]=x[0];
        oldcenter[4]=x[1];
        oldcenter[5]=x[2];
        euler_to_quat(phi,theta,psi,&oldorient[4]);
        update_coords(1,oldcenter,oldorient,oldcoords);
        //Rotate both fragments ABOUT THE ORIGIN by a random rotation.
        //xx[0]=(2.0*genrand_real3()-1.0)*boxsize;
        //xx[1]=(2.0*genrand_real3()-1.0)*boxsize;
        //xx[2]=(2.0*genrand_real3()-1.0)*boxsize;
        /*rand_unif_quat(&q[0]);
        quat_to_matrix(q,&rotmatrix[0][0]);
        //matmul(&rotmatrix[0][0],xx,disp);
        //oldcenter[0]+=disp[0];
        //oldcenter[1]+=disp[1];
        //oldcenter[2]+=disp[2];
        multiply_quat(&oldorient[0],q,newq);
        for (k=0; k<4; k++) oldorient[k]=newq[k];
        update_coords(0,oldcenter,oldorient,oldcoords);
        //Rotate fragment 1's center about the origin.
        matmul(&rotmatrix[0][0],&oldcenter[3],xx);
        oldcenter[3]=xx[0];
        oldcenter[4]=xx[1];
        oldcenter[5]=xx[2];
        multiply_quat(&oldorient[4],q,newq);
        //multiply_quat(q,&oldorient[4],newq);
        for (k=0; k<4; k++) oldorient[4+k]=newq[k];
        update_coords(1,oldcenter,oldorient,oldcoords);
        //Translate both fragments at random.
        xx[0]=(2.0*genrand_real3()-1.0)*boxsize;
        xx[1]=(2.0*genrand_real3()-1.0)*boxsize;
        xx[2]=(2.0*genrand_real3()-1.0)*boxsize;
        oldcenter[0]+=xx[0];
        oldcenter[1]+=xx[1];
        oldcenter[2]+=xx[2];
        oldcenter[3]+=xx[0];
        oldcenter[4]+=xx[1];
        oldcenter[5]+=xx[2];*/
        //Calculate both the exact and table-based energy.
        exact=TRUE;
        enexact=interaction_energy(FALSE,0,1,oldcenter,oldorient,oldcoords);
        exact=FALSE;
        entable=interaction_energy(FALSE,0,1,oldcenter,oldorient,oldcoords);
        snprintf(fname,sizeof(fname),"test/test-%d.pdb",itest);
        test_output=fopen(fname,"w");
        write_frame_pdb(test_output,0,oldcoords);
        fclose(test_output);
        printf("test %d energies: %.8f %.8f\n",itest,entable,enexact);
        err=(entable-enexact);
        if (enexact!=0) err=err/enexact;
        ok=(((entable>=DUMMY_ENERGY) && (enexact>50)) || (err<1.0e-6)) && (interp==INTERP_NONE);
        ok=ok || ((((entable>=DUMMY_ENERGY) && (enexact>50)) || (err<0.1)) && (interp!=INTERP_NONE));
        if (!ok) {
            printf("Test failed.\n");
            die();
        }
    }
}



