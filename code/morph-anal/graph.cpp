#include <cstdio>
#include <cstdlib>
#include <vector>
#include "util.h"
#include "graph.h"


graph::graph(long int _nvertices)
{
    long int incsize,i;
    nvertices=_nvertices;
    incsize = (nvertices*(nvertices-1))/2;
    //incidence = new_empty(incsize);
    incidence.init(incsize);
    connections = new vector<long int>[nvertices];
    //for (i=0; i<nvertices; i++) checkalloc[i]=new vector_
}

graph::~graph()
{
    if (connections!=NULL) delete [] connections;//maybe this should be delete[] ?
    //if (incidence!=NULL) free(incidence);
}


long int compare_ints(void * a, void * b)
{
    long int ia = *((long int *) a);
    long int ib = *((long int *) b);
    if (ia<ib) return -1;
    if (ia>ib) return 1;
    return 0;
}


bool graph::is_edge(long int i, long int j)
{
    if ((i>=nvertices) || (j>=nvertices) || (i==j)) {
        printf("graph: vertex out of bounds\n");
        die();
    }
    return incidence[index(i,j)];
}

void graph::add_edge(long int i, long int j)
{
    if ((i>=nvertices) || (j>=nvertices) || (i==j)) {
        printf("graph: vertex out of bounds\n");
        die();
    }
    long int idx = index(i,j);
    incidence+=idx;
    connections[i].push_back(j);
    connections[j].push_back(i);
    //qsort(&connections[i][0],
}

void graph::remove_edge(long int i, long int j)
{
    if ((i>=nvertices) || (j>=nvertices) || (i==j)) {
        printf("graph: vertex out of bounds\n");
        die();
    }
    long int idx = index(i,j);
    if (!incidence[idx]) return; //don't remove edge twice
    incidence-=idx;
    vector<long int>::iterator pos,found;
    for (pos=connections[i].begin(); pos!=connections[i].end(); pos++) if (*pos==j) found=pos;
    connections[i].erase(found); //this woun't work in the buggy state where the connections doesn't have the edge
    for (pos=connections[j].begin(); pos!=connections[j].end(); pos++) if (*pos==i) found=pos;
    connections[j].erase(found);
}

void graph::output(FILE * output)
{
    long int i,j;
    //vector<long int>::iterator pos;
    /*for (j=0; j<nvertices; j++) for (i=0; i<j; i++) {
        printf("%d %ld %ld %c\n",i,j,index(i,j),yesno(is_edge(i,j)));
    }*/
    for (i=0; i<nvertices; i++) {
        fprintf(output,"Vertex %ld:",i);
        //for (pos=connections[i].begin(); pos!=connections[i].end(); pos++) printf(" %ld",*pos);
        for (j=0; j<connections[i].size(); j++) fprintf(output," %ld",connections[i][j]);
        fprintf(output,"\n");
    }
}

long int graph::degree(long int i)
{
    return connections[i].size();
}

void graph::get_incident_vertices(long int ivertex, subset& vertices)
{
    long int j;
    vertices.init(nvertices);
    for (j=0; j<connections[ivertex].size(); j++) vertices+=connections[ivertex][j];
}

//the following two subroutines identify connected components in the subgraph of the current graph
//induced by the subset "vertices"
//depth first search
void graph::id_connected_component(subset& vertices, long int ivertex, long int component, long int * ids)
{
    long int j,jvertex;
    if (!vertices[ivertex]) return;
    ids[ivertex]=component;
    for (j=0; j<connections[ivertex].size(); j++) {
        jvertex=connections[ivertex][j];
        if ((ids[jvertex]<0) && (vertices[jvertex])) id_connected_component(vertices,jvertex,component,ids);
    }
}
//connected components within the subset vertices
long int graph::connected_components(subset& vertices, long int * ids)
{
    long int ivertex,icomp,component;
    for (ivertex=0; ivertex<nvertices; ivertex++) ids[ivertex]=-1;
    ivertex=0;
    component=0;
    while (ivertex<nvertices) {
	//need to skip over any vertices not part of the subset
	while ((ivertex<nvertices) && ((ids[ivertex]>=0) || !vertices[ivertex])) ivertex++;
        id_connected_component(vertices,ivertex,component,ids);
        while ((ivertex<nvertices) && ((ids[ivertex]>=0) || !vertices[ivertex])) ivertex++;
        component++;
    }
    //compoennt is the number of components
    //sort in order by size
    return component;
}
//comparison routine for qsort to order by size, the parameters are really pointers to clique_info
//sort in order from largest to smallest
int compare_cliques(const void * clique1, const void * clique2)
{
    int count1=((clique_info *) clique1)->count;
    int count2=((clique_info *) clique2)->count;
    if (count1>count2) return -1;
    else if (count1<count2) return 1;
    else return 0;
}

//they aren't really cliques, but we reuse the structure
long int graph::connected_components(subset& vertices,vector<clique_info> * comps)
{
    long int * ids;
    long int ncomp,icomp,ivertex;
    ids=(long int *) checkalloc(nvertices,sizeof(long int));
    ncomp=connected_components(vertices,ids);
    comps->resize(ncomp);
    for (icomp=0; icomp<ncomp; icomp++) {
        (*comps)[icomp].count=0;
        (*comps)[icomp].vertices=NULL;
    }
    for (ivertex=0; ivertex<nvertices; ivertex++) if (vertices[ivertex]) {
        icomp=ids[ivertex];
        (*comps)[icomp].vertices=(long int *) checkrealloc((*comps)[icomp].vertices,(*comps)[icomp].count+1,sizeof(long int));
        (*comps)[icomp].vertices[(*comps)[icomp].count]=ivertex;
        (*comps)[icomp].count++;

    }
    qsort(&(*comps)[0],comps->size(),sizeof(clique_info),compare_cliques);
    return ncomp;
}
//Create a graph obtained by deleting one vertex and all edges involving it.. The "this" clause means
graph::graph(graph * orig, long int del)
{
    long int i, j, ii, jj,incsize;
    vector<long int>::iterator pos;
    //this duplicates code from "graph::graph" above, but can't call one constructor from another
    nvertices=orig->nvertices-1;
    incsize = (nvertices*(nvertices-1))/2;
    incidence.init(incsize);
    connections = new vector<long int>[nvertices];
    //now we have an empty graph. (i,j) are vertices in the original graph; (ii,jj) are vertices in the new graph.
    for (i=0; i<orig->nvertices; i++) if (i!=del) {
        if (i>del) ii=i-1; else ii=i;
        for (pos=orig->connections[i].begin(); pos!=orig->connections[i].end(); pos++) {
            j=*pos;
            if (j!=del) {
                if (j>del) jj=j-1; else jj=j;
                add_edge(ii,jj);
            }
        }
    }
}
void graph::read_from_file(char * fname)
{
    FILE * input;
    long int i,j,jj,_nvertices;
    input=fopen(fname,"r");
    if (input==NULL) {
        printf("could not open file %s\n",fname);
        die();
    }
    fscanf(input,"%ld\n",&_nvertices);
    if (nvertices!=_nvertices) {
        printf("error: wrong number of vertices\n");
        return;
    }
    incidence.clear();
    for (i=0; i<nvertices; i++) connections[i].clear();
    while (!feof(input)) {
        fscanf(input,"%ld %ld\n",&i,&j);
        add_edge(i,j);

    }
    fclose(input);
}



void graph::write_to_file(char * fname)
{
    FILE * output;
    long int i,j,jj;
    output=fopen(fname,"w");
    if (output==NULL) {
        printf("could not open file %s\n",fname);
        die();
    }
    fprintf(output,"%ld\n",nvertices);
    for (i=0; i<nvertices; i++) {
        for (jj=0; jj<connections[i].size(); jj++) {
            j=connections[i][jj];
            fprintf(output,"%ld %ld\n",i,jj);
        }
    }
    fclose(output);
}

/*void to_clique_info(long int n, subset& s, clique_info * info)
{
    int i;
    info->count=0;
    info->vertices=NULL;
    for (i=0; i<n; i++) if (s[i]) {
        info->vertices=(long int *) checkrealloc(info->vertices,info->count+1,sizeof(long int));
        info->vertices[info->count]=i;
        info->count++;
    }
}


//create a degeneracy ordering for the graph
void graph::degeneracy_ordering(vector<long int> * order, long int * degen)
{
    bool done;
    long int i;
    long int * d;
    long int max_d,k,v,w;
    vector<long int>::iterator pos;
    char * * vsets;
    char * added;
    added=new_empty(nvertices);
    d=(long int *) checkalloc(nvertices,sizeof(long int));
    max_d=0;
    //set up list of degrees and compute max degree
    for (i=0; i<nvertices; i++) d[i]=degree(i);
    max_d=0;
    for (i=0; i<nvertices; i++) if (d[i]>max_d) max_d=d[i];
    //this is a list of sets: each set contains the vertices with a given degree
    vsets=(subset *) checkalloc(max_d+1,sizeof(subset));
    for (i=0; i<=max_d; i++) vsets[i].init(nvertices);
    for (i=0; i<nvertices; i++) set_bit(vsets[d[i]],i);

    order->clear();
    k=0;
    for (;;) {
        //scan vsets to find one that is nonempty
        i=0;
        while ((i<=max_d) && (is_empty(nvertices,vsets[i],vsets[i]))) i++;
        if (i>max_d) break; //they're all empty, we're done
        //set k = max(k,i)
        if (i>k) k=i;
        //scan vsets[i] for the first vertex in it
        v=0;
        while ((v<nvertices) && (!is_set(vsets[i],v))) v++;
        //it should be that d[v]=i
        if ((v>=nvertices) || (d[v]!=i)) {
            printf("error");
            die();
        }
        //add v to the list and remove from vsets[i]
        order->push_back(v);
        clear_bit(vsets[i],v);
        set_bit(added,v);
        //for each w in the neighbors of v, subtract one from d[w] and move it from one list to another
        for (pos=connections[v].begin(); pos!=connections[v].end(); pos++) {
            w=*pos;
            if (!is_set(added,w)) {
                clear_bit(vsets[d[w]],w);
                if (d[w]>0) {
                    d[w]--;
                    set_bit(vsets[d[w]],w);
                }
            }
        }
    }
    *degen=k;
    //clean up
    free(added);
    for (i=0; i<max_d; i++) free(vsets[i]);
    free(vsets);
    free(d);
}

//The Bron-Kerbosch clique detection algorithm.  Find all maximal cliques that include all vertices in r, some of those in p, and none of those in x.
//r, p, and x are subsets of vertices, represented as bit strings.
void graph::detect_cliques(int call_level, vector<clique_info> * cliques, subset& r, subset& p, subset& x)
{
    long int i,pivot,deg,pivotdeg,nbytes,v;
    graph * newgraph;
    char * newp;
    char * newx;
    char * newr;
    char * neighbors;
    clique_info info;
    vector<long int>::iterator pos,pos2;
    if (is_empty(nvertices,p,x)) {
        //r is a maximal clique, identify in ids and return
        to_clique_info(nvertices,r,&info);
        cliques->push_back(info);
        //do not free newr, will be managed by caller
        if (cliques->size()%1000==0) {
            qsort(&(*cliques)[0],cliques->size(),sizeof(clique_info),compare_cliques);
            printf("%d cliques identified, maximum size %ld\n",cliques->size(),(*cliques)[0].count);
        }
        return;
    }
    //otherwise, choose as a pivot the vertex in union(p,x) of highest degree;
    pivotdeg=-1;
    for (i=0; i<nvertices; i++) if (is_set(p,i) || is_set(x,i)) {
        deg=degree(i);
        if (deg>pivotdeg) {
            pivotdeg=deg;
            pivot=i;
        }
    }
    neighbors=new_empty(nvertices);
    newp=new_empty(nvertices);
    newr=new_empty(nvertices);
    newx=new_empty(nvertices);
    //for each vertex in p but not neighbors(v)
    for (v=0; v<nvertices; v++) if (is_set(p,v) && ((pivot==v) || !is_edge(pivot,v))) {
        //create new sets newp, newx = intersect(p,neighbors(v)) or intersect(x,neighbors(v))
        if (call_level<=0) printf("call level %ld, vertex = %ld\n",call_level,v);
        clear_all(nvertices,neighbors);
        for (pos=connections[v].begin(); pos!=connections[v].end(); pos++) {
            i=*pos; //i is a neighbor of v
            set_bit(neighbors,i);
        }
        dup_bits(nvertices,r,newr);
        dup_bits(nvertices,p,newp);
        dup_bits(nvertices,x,newx);
        set_bit(newr,v); //now newr = union(r,{v})
        nbytes=(nvertices>>3)+1;
        for (i=0; i<nbytes; i++) {
            newp[i]=p[i]&neighbors[i];
            newx[i]=x[i]&neighbors[i];
        }
        //now newp = intersect(p,neighbors) and newx=intersect(x,neighbors)
        //test output

        //recursive call
        detect_cliques(call_level+1,cliques,newr,newp,newx);
        //get rid of newp, newr, newx, neighbors
        //remove v from p and add it to x
        clear_bit(p,v);
        set_bit(x,v);
    }
    free(neighbors);
    free(newp);
    free(newr);
    free(newx);
}
//root call
void graph::detect_cliques(vector<clique_info> * cliques)
{
    long int nbytes;
    long int i,j;
    char * p;
    char * r;
    char * x;
    char * newr;
    char * newp;
    char * newx;
    char * neighbors;
    long int v,degen;
    vector<long int> order; //degeneracy ordering
    vector<long int>::iterator pos;
    p=new_empty(nvertices);
    r=new_empty(nvertices);
    x=new_empty(nvertices);
    neighbors=new_empty(nvertices);
    newp=new_empty(nvertices);
    newr=new_empty(nvertices);
    newx=new_empty(nvertices);
    cliques->clear();
    //for (i=0; i<nvertices; i++) ids[i]=-1;
    //p needs to have every vertex set
    for (i=0; i<nvertices; i++) set_bit(p,i);
    //detect_cliques(0,cliques,r,p,x);
    degeneracy_ordering(&order,&degen);
    printf("Graph has degeneracy %ld\n",degen);
    //work with degeneracy ordering to reduce number of recursive calls (see wiki article)
    //this resembles the main loop in the recursive version of detect_cliques above
    for (j=order.size()-1; j>=0; j--) {
        v=order[j];
        if (j%100==0) printf("call level 0, index = %d, vertex = %ld\n",j,v);
        clear_all(nvertices,neighbors);
        for (pos=connections[v].begin(); pos!=connections[v].end(); pos++) {
            i=*pos; //i is a neighbor of v
            set_bit(neighbors,i);
        }
        dup_bits(nvertices,r,newr);
        dup_bits(nvertices,p,newp);
        dup_bits(nvertices,x,newx);
        set_bit(newr,v); //now newr = union(r,{v})
        nbytes=(nvertices>>3)+1;
        for (i=0; i<nbytes; i++) {
            newp[i]=p[i]&neighbors[i];
            newx[i]=x[i]&neighbors[i];
        }
        detect_cliques(1,cliques,newr,newp,newx);
        //remove v from p and add it to x
        clear_bit(p,v);
        set_bit(x,v);
    }
    qsort(&(*cliques)[0],cliques->size(),sizeof(clique_info),compare_cliques);
    printf("total %ld cliques identified\n",cliques->size());
    free(p);
    free(r);
    free(x);
    free(newp);
    free(newr);
    free(newx);
    free(neighbors);
}*/
