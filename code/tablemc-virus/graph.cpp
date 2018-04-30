#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include "util.h"
#include "graph.h"

graph::graph(int _nvertices)
{
    int nbytes,i;
    nvertices=_nvertices;
    nbytes = (nvertices*(nvertices-1))/2;
    nbytes = nbytes/8+1;
    incidence = (char *) checkalloc(nbytes,1);
    memset(incidence,0,nbytes);
    connections = new vector<int>[nvertices];
    //for (i=0; i<nvertices; i++) checkalloc[i]=new vector_
}

graph::~graph()
{
    delete connections;//maybe this should be delete[] ?
    free(incidence);
}

int compare_ints(void * a, void * b)
{
    int ia = *((int *) a);
    int ib = *((int *) b);
    if (ia<ib) return -1;
    if (ia>ib) return 1;
    return 0;
}


bool graph::is_edge(int i, int j)
{
    int idx = index(i,j);
    int b=idx>>3;
    int bit=idx&7;
    if ((incidence[b])&(1<<bit)) return true;
    return false;
}

void graph::add_edge(int i, int j)
{
    int idx = index(i,j);
    int b=idx>>3;
    int bit=idx&7;
    if ((incidence[b])&(1<<bit)) return; //ensures we don't double add edge to connections
    incidence[b]|=(1<<bit);
    connections[i].push_back(j);
    connections[j].push_back(i);
    //qsort(&connections[i][0],
}

void graph::remove_edge(int i, int j)
{
    int idx = index(i,j);
    int b=idx>>3;
    int bit=idx&7;
    if (!((incidence[b])&(1<<bit))) return; //ensures we don't remove an edge
    incidence[b]&=~(1<<bit);
    vector<int>::iterator pos,found;
    for (pos=connections[i].begin(); pos!=connections[i].end(); pos++) if (*pos==j) found=pos;
    connections[i].erase(found); //this woun't work in the buggy state where the connections doesn't have the edge
    for (pos=connections[j].begin(); pos!=connections[j].end(); pos++) if (*pos==i) found=pos;
    connections[j].erase(found);
}

void graph::test_index(void)
{
    int i,j;
    //vector<int>::iterator pos;
    for (j=0; j<nvertices; j++) for (i=0; i<j; i++) {
        printf("%d %d %d %c\n",i,j,index(i,j),yesno(is_edge(i,j)));
    }
    for (i=0; i<nvertices; i++) {
        printf("Vertex %d:",i);
        //for (pos=connections[i].begin(); pos!=connections[i].end(); pos++) printf(" %d",*pos);
        for (j=0; j<connections[i].size(); j++) printf(" %d",connections[i][j]);
        printf("\n");
    }
}

int graph::degree(int i)
{
    return connections[i].size();
}

//depth first search
void graph::id_connected_component(int ivertex, int component, int * ids)
{
    int j,jvertex;
    ids[ivertex]=component;
    for (j=0; j<connections[ivertex].size(); j++) {
        jvertex=connections[ivertex][j];
        if (ids[jvertex]<0) id_connected_component(jvertex,component,ids);
    }
}

int graph::connected_components(int * ids)
{
    int ivertex,component;
    for (ivertex=0; ivertex<nvertices; ivertex++) ids[ivertex]=-1;
    ivertex=0;
    component=0;
    while (ivertex<nvertices) {
        id_connected_component(ivertex,component,ids);
        while ((ivertex<nvertices) && (ids[ivertex]>=0)) ivertex++;
        component++;
    }
    return component;
}

