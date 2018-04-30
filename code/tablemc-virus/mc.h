#ifndef MC_H_INCLUDED
#define MC_H_INCLUDED

#include <stdio.h>
#if defined(PARALLEL) || defined(EXCHANGE)
#include <mpi.h>
#endif
#include "tables.h"
#include "fragments.h"
#include "nblist.h"
#include "go_model.h"




#define MOVE_TRANS 1
#define MOVE_ORIENT 2
//#define MOVE_UNIF_ORIENT 3
#define NUM_MOVES 2






#define INTERP_NONE 0
struct table;

static const char * mc_move_names[NUM_MOVES+1] = {"","Translational","Orientational"}; //,"Uniform orientational"};
/*File names: coordinate output, quaternion output, starting restart file (if needed), ending restart file*/
struct transinfo {
    double center[3]; //the "t" vector from the MTRIXn records.  I believe this is the center about which transformations are to take place.
    double q[4]; //rotation
};

class simulation {
private:
    table * tables[MAX_FRAGMENT_TYPES][MAX_FRAGMENT_TYPES];
    double scalefactors[MAX_FRAGMENT_TYPES][MAX_FRAGMENT_TYPES];
    int exact, pbc, interp, enwrite;
    double cutoff, cutoff2, boxsize, halfboxsize;
    //double go_hardcore,go_cutoff, delta, native_energy, nonnative_energy;
    go_model_params params;
    //char go_model_fname[255];
    go_model_info * go_models[MAX_FRAGMENT_TYPES][MAX_FRAGMENT_TYPES];
    long int entablecount, enexactcount, enevalcount;
    FILE * energy_output;
    FILE * pairs_output;
    char xyzfname[255],quatfname[255],startrestartfname[255],endrestartfname[255];
    FILE * xyzoutput;
    FILE * quatoutput;
    int nfragtypes;
    /*Reference geometries for the fragments -- indexed by type*/
    fragmenttype * fragtypes[MAX_FRAGMENT_TYPES];
    /*number of actual fragments, total number of atoms*/
    int nfrag,natom,ninitatom;
    /*which type of fragments each fragment is*/
    fragmentinfo * frags;
    int ntrans; //number of transformations needed to complete structure
    transinfo * trans; //the transformations
    /*current center and orientation of each fragment, and current coordinates*/
    double * initcoords;
    double * mass;
    atominfo * initatoms; //information about each atom

    double * newcenter;
    double * neworient;
    double * newcoords;
    double * oldcenter;
    double * oldorient;
    double * oldcoords;
    double * initcenter;
    /*other data*/
    double beta; /*Monte Carlo temperature, expressed as beta*/

    long int nmcstep, nprevstep, nsave_quat,nsave_xyz, nprint; /*number of MC steps, frequency of saving conformations, frequency of printing*/
    double dtrans,dorient; /*maximum translational/orientational displacement, probabilities of each*/
    double prob[NUM_MOVES+1],cumprob[NUM_MOVES+1];

    //int startoption;
    bool reseed; //whether or not to reseed the randon number generator.
    //Nonbond list stuff.
    /*int nb_list_per_frag,nb_list_size;
    int * nonbond_list;
    int * nb_list_count;
    double last_nb_list_center[3*MAX_FRAGMENTS];
    int nbox[3],nboxtot;
    int * boxlist;
    int * boxcount;
    int * boxneighbors;
    int * boxneighcount;*/
    double listcutoff;
    bool use_nb_list;
    fragment_nblist * frag_nblist;
#ifdef UMBRELLA
    //umbrella potential
    double kumb, rumb;
    int ifragumb, jfragumb;
    double umbrella_energy(double * center);
#endif
#if defined(PARALLEL) || defined(EXCHANGE)
    int mynod, numnod;
#endif
#ifdef EXCHANGE
    int myrep, nrep, exchfreq;
    double * betas;
    MPI_File rexlog;
    void exchange_init(FILE * input);
    void exchange(int icycle, double * cum_energy, double * center, double * orient, double * coords);
#endif
    void recenter(void);
    void mcmove(int * imoved, int * movetype, double * center, double * orient, double * coords);
    void update_coords(int ifrag, double * centers, double * orients, double * coords);
    void write_frame_xyz(FILE * output, long int istep, double * coords);
    void write_frame_pdb(FILE * output, long int istep, double * coords);
    void write_pair_pdb(FILE * output, int ifrag, int jfrag, double * coords);
    void read_frame_quat(FILE * input, long int * istep, double * center, double * orient);
    void write_frame_quat(FILE * output, long int istep, double * center, double * orient);
    void write_dcd_header(FILE * dcdfile);
    void write_dcd_frame(FILE * dcdfile, double * coords);
    void copy_frag(int ifrag, double * center1, double * orient1, double * coords1, double * center2, double * orient2, double * coords2);
    double interaction_energy(int pbc, int ifrag, int jfrag, double * center, double * orient,double * coords);
    double moved_energy(int imoved, double * center, double * orient, double * coords);
    double total_energy(double * center, double * orient, double * coords);
    void print_summary_info(void);
    int frag_type_by_name(char * name);
    int find_init_atom(char chain, int res, char aname[4]);
    void process_commands(FILE * input);
    int insert_fragment(int itype);
    int insert_fitted_fragment(int itype,  int startres, int endres);
    void apply_one_transformation(int norigfrag, int itrans);
    void apply_transformations(void);
    void write_charmm_images(char * fname);
    void expand_trajectory(char * datfname, char * dcdfname, long int nframes, int freq, double _boxsize);
    //int count_native_contacts(int pbc, double halfboxsize, double boxsize, double cutoff, double ratio, double * coords, FILE * verbose_output);
    void reaction_coordinates(int pbc, double halfboxsize, double boxsize, double cutoff, double ratio, double frac, double * center, double * coords, int * ncomponents, int * max_component, double * _mindist, long int iframe, FILE * verbose_output);
    void read_go_models(FILE * f, go_model_params * params);
    void read_scale_factors(FILE * f);
    void reaction_coordinates_traj(char * datfname, char * outfile, char * voutname, long int nframes, int freq, double _boxsize, double cutoff, double ratio, double frac);
    /*void initialize_system_from_pdb_file(char * fname);
    void initialize_system_random(void);*/
    void initialize_system_lattice(bool rand_orientations, double exclude_radius);
    void load_tables(const char * fmt);
    void read_pdb_file(char * fname);
    void read_restart(char * fname);
    void write_restart(long int istep, char * fname);
    void create_nonbond_list(double * center);
    int check_nb_list(int imoved, double * center);
public:
    void comparison_test(void);
#if defined(PARALLEL) || defined(EXCHANGE)
    simulation(char* fname, int _mynod, int _numnod);
#else
    simulation(char * fname);
#endif
    ~simulation();
    void mcloop(void);
};



#endif // MC_H_INCLUDED
