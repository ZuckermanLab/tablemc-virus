#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "string.h"
#include "mc.h"
#include "tables.h"
#include "rotations.h"
#include "mt.h"
#include "util.h"

#ifdef __unix__
#include <fenv.h>
#else
#include <float.h>
//#ifndef _EM_OVERFLOW
//#define _EM_OVERFLOW EM_OVERFLOW
//#endif
#endif

//Command line syntax: tablemc input_file  OR tablemc generate input_file table_file]
int main(int argc, char * argv[])
{
    int imoved,part,numparts,ntables,itable;
    double en,en2,de;
    char fname[255],buffer[255];
    time_t now;
    table * newtable;
    table * * tables; //for combination
    simulation * sim;
    fragmenttype * frag;
#ifdef __unix__
    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#else
//    _controlfp(0, _EM_ZERODIVIDE | _EM_INVALID | _EM_OVERFLOW);
#endif
    fill_trig_tables();
    /*double q[4],qq[4],r[3][3];
    init_genrand(57);
    for(;;) {
    rand_unif_quat(q);
    quat_to_matrix(q,&r[0][0]);
    matrix_to_quat(r,qq);
    }*/
    //sim = new simulation(argv[1]);
    //delete sim;
    //printf("%ld\n",symmetric_round(3.5));
    //printf("%ld\n",symmetric_round(-3.5));
    //printf("%ld\n",-7 % 5);
    //printf("%d %s %s %s %s\n",argc,argv[1],argv[2],argv[3],argv[4]);
    //printf("%.8f\n",acos(0.5));
    if (argc<2) {
        printf("Syntax: tablemc input_file\n");
        die();
    }
#if defined(PARALLEL) || defined(EXCHANGE)
    parallel_init(&mynod,&numnod,argv[argc-1]);
#endif
    time(&now);
    strncpy(buffer,ctime(&now),sizeof(buffer));
    printf("Tabulated Monte Carlo starting at %s\n",buffer);
    printf("Executable name: %s\n",argv[0]);
    fflush(stdout);
    if (strncasecmp(argv[1],"generate",8)==0) {
        if (argc<4) {
           printf("Syntax: tablemc generate input_file table_file\n");
           die();
        }
        part=0;
        numparts=1;
        if (argc==6) {
            part=atoi(argv[4]);
            numparts=atoi(argv[5]);
        }
        newtable=new table(argv[2],part,numparts);
        newtable->write_table(argv[3]);
        delete newtable;
#if !defined(PARALLEL) && !defined(EXCHANGE)
    } else if (strncasecmp(argv[1],"run",3)==0) {
        sim=new simulation(argv[2]);
        sim->mcloop();
        delete sim;
//Table functions other than generation not available in parallel mode. (Generation only provided in case I want to parallelize it.)
    } else if (strncasecmp(argv[1],"test",3)==0) {
        sim=new simulation(argv[2]);
        sim->comparison_test();
        delete sim;
    } else if (strncasecmp(argv[1],"dx",2)==0) {
        newtable=new table(argv[2],FALSE);
        newtable->write_dx(argv[8],atof(argv[3]),atof(argv[4]),atof(argv[5]),atof(argv[6]),atof(argv[7]));
        delete newtable;
    } else if (strncasecmp(argv[1],"smooth",6)==0) {
        newtable=new table(argv[2],FALSE);
        newtable->do_smooth(atof(argv[3]),atof(argv[4]),atof(argv[5]));
        newtable->write_table(argv[6]);
        delete newtable;
    } else if (strncasecmp(argv[1],"pmf",3)==0) {
        newtable=new table(argv[2],FALSE);
        newtable->calculate_dist_pmf(atof(argv[3]),argv[4]);
        delete newtable;
    } else if (strncasecmp(argv[1],"combine",7)==0) {
        //combine tables
        //arguments: tablemc combine new-table table1,...
        ntables=argc-3;
        tables=(table * *) checkalloc(ntables,sizeof(table *));
        printf("Combining tables.\n");
        for (itable=0; itable<ntables; itable++) {
            printf("------------------------------------------------\n");
            printf("Loading table %d of %d\n",itable+1,ntables);
            tables[itable]=new table(argv[itable+3],false);
        }
        printf("------------------------------------------------\n");
        newtable=new table(ntables,tables);
        newtable->write_table(argv[2]);
        for (itable=0; itable<ntables; itable++) delete tables[itable];
        free(tables);
#endif
    } else {
        //Default to running a simulation.
#if defined(PARALLEL) || defined(EXCHANGE)
        sim=new simulation(argv[1],mynod,numnod);
#else
        sim=new simulation(argv[1]);
#endif
        sim->mcloop();
        delete sim;
    }
    time(&now);
    strncpy(buffer,ctime(&now),sizeof(buffer));
    printf("Finished at %s\n",buffer);
    fflush(stdout);
#if defined(PARALLEL) || defined(EXCHANGE)
    parallel_finish();
#endif
   //comparison_test();
    //symmetry_test();
    //initialize_system();
    //read_restart("restart.txt");
    //mcloop();
    return 0;
}
