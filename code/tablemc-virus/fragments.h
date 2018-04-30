#ifndef FRAGMENTS_H_INCLUDED
#define FRAGMENTS_H_INCLUDED
#include <stdio.h>

#define MAX_ATOMS_PER_FRAGMENT 20
#define MAX_FRAGMENT_TYPES 10
#define MAX_FRAGMENTS 1000
#define MAX_ATOMS MAX_FRAGMENTS*MAX_ATOMS_PER_FRAGMENT
#define MAX_FRAGMENT_NAME 32

//topology info concerning atom
struct atominfo {
    int res;
    char resname[4];
    char name[5];
    int type;
    char chain;
    int fragment; //which fragment -- used only for the atominfo array in simulation
    bool is_site;
};

//needed for simulation, contact map.
struct fragmentinfo {
    int type; //index into fragtypes array
    int start; //starting atom in system
    int end; //ending atom, maybe not needed
    //char chains[10]; //original fitted chain
    int startres, endres; //original starting and ending residue
};


class fragmenttype {
private:
    //double refgeom[3*MAX_ATOMS_PER_FRAGMENT]; /*Reference geometry*/
    double * refgeom; //3*natom
    double * mass;
    //double qtot,dipolemag; /*total charge, magnitude of dipole moment*/
    //double dipole[3]; /*dipole moment in reference geometry*/
public:
    char fragname[MAX_FRAGMENT_NAME];
    char fragfname[255];
    int natom; /*Number of atoms in the fragment*/
    int nsite; //number of go sites
    double totalmass; //total mass
    int startres, endres; //Starting and ending residues as given in the fragment PDB file
    //char names[MAX_ATOMS_PER_FRAGMENT][4]; /*names of all the atoms */
    //int types[MAX_ATOMS_PER_FRAGMENT]; /*Types of atoms*/
    atominfo * atoms;
    fragmenttype(const char * name, const char * fname);
    ~fragmenttype();
    void get_coords(double * center, double * orient, double * coords);
    void fit_fragment(double * coords, double * center, double * orient, double * rmsd);
    void fit_fragment(double * coords, double * center, double * orient, double * weights, double * rmsd);
};

void read_pdb_line(const char * buf, atominfo * atom, double * coords, double * mass);
int get_fragtype_by_name(char * name, int nfragtypes, fragmenttype * * fragtypes);
int get_fragtype_by_fname(char * fname, int nfragtypes, fragmenttype * * fragtypes);
void write_pair_pdb2(FILE * output, fragmenttype * frag1, fragmenttype * frag2, double * coords1, double * coords2);
#endif // FRAGMENTS_H_INCLUDED
