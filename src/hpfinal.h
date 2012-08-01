/*
 * final.h
 *
 *  Created on: Mar 25, 2012
 *      Author: hieppham
 */

#ifndef FINAL_H_
#define FINAL_H_

#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <map>
#include <math.h>
#include <algorithm>

#include <ga/ga.h>
#include <ga/std_stream.h>

#include "Customer.h"
#include "RouteInfo.h"
#include "HGAGenome.h"
#include "HPGradProjectGA.h"

#define max(a,b) (a>b)?a:b;
#define min(a,b) (a<b)?a:b;

const double NEGINF = -1000000;
const double POSINF = 1000000;
const unsigned int VERYLAST = 1000;

inline double distance(Customer* c1, Customer* c2){
    double dx = c1->x - c2->x;
    double dy = c1->y - c2->y;
    double _dist = double(sqrt(dx*dx + dy*dy));
    return _dist;
}

class HPGV{
public:
    static unsigned int nCus;    // number of customers
    static unsigned int mVeh;    // number of vehicles
    static unsigned int tDay;    // number of days in planning horizontal
    static unsigned int numRoute;   // pre-calculate number of routes

    static double maxDuration;  // maximum duration of a route
    static double maxLoad;      // maximum load of a vehicle

    // penalty parameters
    static double hPenalty;
    static double avgQ;
    static double avgD;
    static double avgW;

    // store best cost
    static double bestFeasibleCost;

    static unsigned int nPop;
    static unsigned int nKeep;

    // Solomon I1 heuristic
    static double I1Lambda;
    static double I1Mu;
    static double I1Alpha;

    // Solomon Time-Oriented, Nearest-Neighbor heuristic
    static double SGamma2;
    static double SGamma3;
    // Tabu search parameters
    static long utsIter;
    static double utsLambda;

    // RVNS parameters
    static unsigned int rvnsIter;
    static double rvnsPRev;
    static double rvnsTmax;

    static vector<Customer> gArrC; // global variable - array of visitors
    static vector<unsigned int> multipattern;   // IDs of customers that have more than 1 pattern
    static unsigned int numOfMP;    // size of multipattern array

    static Customer* gDepot;
    static Customer* gDynamic;
    static double gDynamicStart;

    static vector<vector<double> > gDistance;

    static int genCounter;

};

inline double metric(
        Customer* c1,
        Customer* c2,
        double firstBegin,           // start time of first customer
        vector<vector<double> >& gDistance
        ){
    double tmp = firstBegin + c1->d + gDistance[c1->id][c2->id];
    double secondBegin = (tmp > c2->e) ? tmp : c2->e;
    double T = secondBegin - (firstBegin + c1->d);
    tmp = (double)c2->l - tmp;
    return ((1 - HPGV::SGamma2 - HPGV::SGamma3) * gDistance[c1->id][c2->id] + HPGV::SGamma2 * T + HPGV::SGamma3 * tmp);
}

inline int numberOfSetBits(int i)
{
    i = i - ((i >> 1) & 0x55555555);
    i = (i & 0x33333333) + ((i >> 2) & 0x33333333);
    return (((i + (i >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
}

template <class C> void FreeClear( C & cntr ) {
    for ( typename C::iterator it = cntr.begin();
              it != cntr.end(); ++it ) {
        delete * it;
    }
    cntr.clear();
}

#if !defined(GALIB_USE_AUTO_INST)
#include <ga/GAList.C>
#include <ga/GAListGenome.C>
GALIB_INSTANTIATION_PREFIX GAList<int>;
GALIB_INSTANTIATION_PREFIX GAListGenome<int>;
#endif

#endif /* FINAL_H_ */
