#ifndef GRAPH_H_INCLUDED
#define GRAPH_H_INCLUDED

#include <vector>
#include <string.h>
#include "util.h"

using namespace std;

struct clique_info {
    long int count;
    long int * vertices;
};

class graph {
    char * incidence; //bitmap for incidence matrix
    vector<long int> * connections; //for each vertex, stores a list of connected vertices
    inline long int index(long int i, long int j)
    {
        if (i<j) return ((j*(j-1))/2)+i;
        else if (i==j) return -1; //no entry
        else if (i>j) return ((i*(i-1))/2)+j;
    }
    void id_connected_component(long int ivertex, long int component, long int * ids);
    void detect_cliques(int call_level, vector<clique_info> * cliques, char * r, char * p, char * x);
public:
    long int nvertices;
    graph(long int _nvertices);
    graph(graph * orig, long int del);
    ~graph();
    void add_edge(long int i, long int j);
    void remove_edge(long int i, long int j);
    bool is_edge(long int i, long int j);
    void output(FILE * output);
    long int degree(long int i);
    long int connected_components(long int * ids);
    void degeneracy_ordering(vector<long int> * order, long int * degen);
    void detect_cliques(vector<clique_info> * cliques);
};

//utility routines
inline bool is_set(char * bits, long int idx)
{
    long int b=idx>>3;
    long int bit=idx&7;
    if ((bits[b])&(1<<bit)) return true;
    return false;
}

inline void set_bit(char * bits, long int idx)
{
    long int b=idx>>3;
    long int bit=idx&7;
    if ((bits[b])&(1<<bit)) return; //ensures we don't double add edge to connections
    bits[b]|=(1<<bit);
}

inline void clear_bit(char * bits, long int idx)
{
    long int b=idx>>3;
    long int bit=idx&7;
    if (!((bits[b])&(1<<bit))) return; //ensures we don't remove an edge
    bits[b]&=~(1<<bit);
}

//determine if set1 and set2 are both empty
//this actually tests the extra bits on the end also, which should be clear
inline bool is_empty(long int n, char * set1, char * set2)
{
    long int nn=(n>>3)+1;
    long int i;
    for (i=0; i<nn; i++) if ((set1[i]|set2[i])!=0) return false;
    return true;
}

inline char * new_empty(long int n)
{
    long int nbytes=(n>>3)+1;
    char * p = (char *) checkalloc(nbytes,1);
    memset(p,0,nbytes);
    return p;
}

inline void clear_all(long int n, char * bits)
{
    long int nbytes=(n>>3)+1;
    memset(bits,0,nbytes);
}

inline char * dup_bits(long int n, char * bits)
{
    long int nbytes=(n>>3)+1;
    char * p = (char *) checkalloc(nbytes,1);
    memcpy(p,bits,nbytes);
    return p;
}

inline void dup_bits(long int n, char * bits, char * bits2)
{
    long int nbytes=(n>>3)+1;
    memcpy(bits2,bits,nbytes);
}


inline void print_set(char * name, long int n, char * bits)
{
    long int i;
    printf("%s = ",name);
    for (i=0; i<n; i++) if (is_set(bits,i)) printf("%d,",i);
    printf("\n");
}


#endif // GRAPH_H_INCLUDED
