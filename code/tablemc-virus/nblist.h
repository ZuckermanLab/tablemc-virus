#ifndef NBLIST_H_INCLUDED
#define NBLIST_H_INCLUDED

struct pair {
    int i;
    int j;
};

struct grid {
    int nbox[3],nboxtot,max_per_box;
    int * boxlist;
    int * boxcount;
    int * boxneighbors;
    int * boxneighcount;
    double volume, density;
    grid();
    ~grid();
    void create_grid(bool pbc, double listcutoff, int ncoords, double * coords);
//    void create_nonbond_list(std::list<
};

struct fragment_nblist {
    grid * grd;
    double listcutoff;
    int nb_list_per_frag,nb_list_size;
    int * nonbond_list;
    int * nb_list_count;
    //nonbond list stuff
    double * last_nb_list_center;
    fragment_nblist(int nfrag, double listcut);
    ~fragment_nblist();
    void create_nonbond_list(bool pbc, double halfboxsize, double boxsize, int nfrag, double * center);
    bool check_nb_list(bool pbc, double halfboxsize, double boxsize, double cutoff, int nfrag, bool * moved, double * center);
    bool check_nb_list(bool pbc, double halfboxsize, double boxsize, double cutoff, int nfrag, int imoved, double * center);
};



#endif // NBLIST_H_INCLUDED
