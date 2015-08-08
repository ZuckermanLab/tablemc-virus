#ifndef MORPH_H_INCLUDED
#define MORPH_H_INCLUDED

#include <stdio.h>
//fragment types.  We hard code these for now.
#define NFRAGTYPES  2
//static const int nfragtypes = 2;

static const char * fragnames[2] = {"ab","cd"};

struct cluster_info {
    int size; //size of this cluster
    double rmsd; //rmsd of frame atoms to template atoms for this cluster
    //double oldrmsd; //on previous iteraton;
    //int * fragments; //unordered list of fragments in the cluster.
    //int * tmpfragments; //corresponding template fragments
    double disp[3]; //this is the transformation needed to take the template to an approximation of the cluster.
    double rot[4];
    int orig_index; //original index, keep track during sorting
    //In case there are 2 or fewer fragments in this cluster, we need to keep track of their pairings her
    //so that we can redetermine the net rotation fo
    int frag[2];
    int tmpfrag[2];
    //double mindist; //min. distance from the largest cluster to this one
};

struct clusters {
    int nclusters, nfrag, ntmpfrag;
    int * cluster_id; //Which cluster is the given fragment in?
    int * mapping; //Which fragment in the template corresponds to the given fragment?
    int * init_cluster_id; //initial guess;
    int * init_mapping;
    double * mindist;
    double * angle;
    cluster_info * info;
    //bool done;
    clusters(bool pbc, double boxsize, double halfboxsize, int _nfrag, int * fragtypes, double * center, double * orient,
        int _ntmpfrag, int * tmpfragtypes, double * tmpcenter, double * tmporient, double distratio, double anglecutoff, double distcutoff);
    void report(long int iframe, FILE * output);
    void write_pdb_frame(long long int istep, double * center, FILE * output);
    void write_xyz_frame(long long int istep, double * center, FILE * output);
    void get_mindist(int pbc, double halfboxsize, double boxsize, double * center, double * orient, double * tmpcenter, double * tmporient);
    ~clusters();
private:
    void find_transformations(bool pbc, double boxsize, double halfboxsize, double * center, double * orient, double * tmpcenter, double * tmporient);
    void force_fit_fragment(bool pbc, double boxsize, double halfboxsize, int ifrag, double * center, double * orient, int itmpfrag, double * tmpcenter, double * tmporient, int iclus);
    void force_fit_fragment(bool pbc, double boxsize, double halfboxsize, int ifrag, int jfrag, double * center, double * orient, int itmpfrag, int jtmpfrag, double * tmpcenter, double * tmporient, int iclus);
    void add_new_cluster(bool pbc, double boxsize, double halfboxsize, int ifrag, int * fragtypes, double * center, double * orient, int * tmpfragtypes, double * tmpcenter, double * tmporient);
    void add_new_cluster(int ifrag, int jfrag, double * center, double * tmpcenter);
    void insert_cluster(int * fragtypes, int * tmpfragtypes);
    void sort_clusters(int * old_cluster_id);
    void identify_edge_frags(int pbc, double halfboxsize, double boxsize, double * tmpcenter, bool * is_edge);

    //find_clusters(bool pbc, double boxsize, double halfboxsize, int nfrag, double * center, int ntmpfrag, double * tmpcenter);
};

#endif // MORPH_H_INCLUDED
