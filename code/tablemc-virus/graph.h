#ifndef GRAPH_H_INCLUDED
#define GRAPH_H_INCLUDED

#include <vector>
using namespace std;
class graph {
    char * incidence; //bitmap for incidence matrix
    vector<int> * connections; //for each vertex, stores a list of connected vertices
    inline int index(int i, int j)
    {
        if (i<j) return ((j*(j-1))/2)+i;
        else if (i==j) return -1; //no entry
        else if (i>j) return ((i*(i-1))/2)+j;
    }
    void id_connected_component(int ivertex, int component, int * ids);
public:
    int nvertices;
    graph(int _nvertices);
    ~graph();
    void add_edge(int i, int j);
    void remove_edge(int i, int j);
    bool is_edge(int i, int j);
    void test_index(void);
    int degree(int i);
    int connected_components(int * ids);
};


#endif // GRAPH_H_INCLUDED
