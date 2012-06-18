//============================================================================
// Name        : hpfinal.cc
// Author      : Hiep Pham
// Version     : 1.0
// Copyright   : (C) 2012 by _HP
// Description : Main function of HP Final Project in C++, Ansi-style
//============================================================================

#include "hgafinal.h"
#include "cpu_time.c"
using namespace std;

// global variables
unsigned int gmVehicles = 1; // number of vehicles
unsigned int gnCustomers = 1; // number of customers
unsigned int gtDays = 1; // number of days (PVRP)
unsigned int numOfDepots = 1;
int typeOfVRP;
double maxDuration; // maximum duration of a route
double maxLoad; // maximum load of a vehicle
vector<Customer> gArrC; // global variable - array of visitors

HGAGenome bestSol(0);

Customer* gDepot;
Customer* gDynamic;
double gDynamicStart;
vector<Route > gRoute;
RouteData gRouteInfo;
vector<vector<double> > gDistance;
double sig1 = 0, sig2 = 0;   // coefficients for calculating metrics

// penalty parameters
double hPenalty = 0;
double aveQ = 0;
double aveD = 0;
double aveW = 0;

double bestFeasibleCost = 0;

unsigned int nPop = 10;
unsigned int nKeep = 4;

/*
 * ====================================
 */
int getInputData(fstream &ifs, char* filein) {
    string line;
    int tempIter;
    // temporary variables for input data
    int tid, tf, td, tq, ta, tl, te, tc;
    float tx, ty;

    ifs.open(filein, ios::in);
    if (ifs.is_open()) {
        if (ifs.good()) {
            ifs >> typeOfVRP >> gmVehicles >> gnCustomers >> gtDays >> numOfDepots;
        }
        for (unsigned int i = 0; i < gtDays; i++) {
            ifs >> maxDuration >> maxLoad;
        }
        if (maxDuration == 0){
            maxDuration = 10000;
        }
        ifs >> tid >> tx >> ty >> td >> tq >> tf >> ta;
        gDepot = new Customer(tid, tx, ty, td, tq, tf, ta);
        gDepot->a = gDepot->e = gDepot->pattern = 0;
        gDepot->comb.clear();
        ifs >> te >> tl;
        gDepot->setTime(te, tl);

        while (ifs.good() && gArrC.size() < gnCustomers) {
            tempIter = 0;
            ifs >> tid >> tx >> ty >> td >> tq >> tf >> ta;
            Customer* gCustomer = new Customer(tid, tx, ty, td, tq, tf, ta);
            gArrC.push_back(*gCustomer);
            gArrC[tid - 1].comb.clear();
            for (tempIter = 0; tempIter < ta; tempIter++) {
                ifs >> tc;
                gArrC[tid - 1].comb.push_back(tc);
            }
            gArrC[tid - 1].comb.resize(ta);
            ifs >> te >> tl;
            gArrC[tid - 1].setTime(te, tl);
            delete gCustomer;
        }
        ifs.close();
    }else{
        cerr << "Can not open \"" << filein << "\" for input.\n";
        exit(1);
    }
    return 0;
}

int writeOutputData(fstream &ofs, char* fileout, char* inputFileName, double totalTime){
    unsigned int iDay = 0;
    unsigned int iVeh = 0;
    unsigned int vod = 0;
    int hours = (int)(totalTime/3600);
    double remainSecOne = totalTime - 3600*hours;
    int minInHour = (int)(remainSecOne/60);
    remainSecOne -= (60*minInHour);
    int minAll = (int)(totalTime/60);
    double remainSecTwo = totalTime - (60*minAll);
    ofs.open(fileout, ios::out|ios::app);
    if (ofs.is_open()) {
        while (ofs.good()) {
            ofs << "Result of HGA for \"" << inputFileName << "\" (" << totalTime << " sec)\n";
            ofs << "(" << hours << " hr " << minInHour << " min " << remainSecOne << " sec" << ")";
            ofs << "(" << minAll << " min " << remainSecTwo << " sec" << ")\n";
            if (!bestSol.feasible){
                ofs << "* ";
            }
            ofs << bestSol.durationCost << endl;
                HGAGenome::globalCopyRoute(bestSol);
                for (iDay = 0; iDay < gtDays; iDay++){
                    for (iVeh = 0; iVeh < gmVehicles; iVeh++){
                        vod = iDay*gmVehicles + iVeh;
                        HGAGenome::delayDeparture(vod, gRoute[vod]);
                        Route::iterator uIter, endIter;
                        for (uIter = gRoute[vod].begin(), endIter = gRoute[vod].end(); uIter != endIter; ++uIter){
                            if (uIter == gRoute[vod].begin()){
                                ofs << (iDay+1) << "  " << (iVeh + 1) << "\t";
                                ofs << gRouteInfo[vod]->duration << "\t" << gRouteInfo[vod]->load << "\t";
                                ofs << "0(" << gRouteInfo[vod]->timeLeaveDepot << ")  ";
                                ofs << (*uIter)->cus->id << "[" << (*uIter)->cus->pattern << "]" << "(" << (*uIter)->timeStartService << ")  ";
                            }else{
                                ofs << (*uIter)->cus->id << "[" << (*uIter)->cus->pattern << "]" << "(" << (*uIter)->timeStartService << ")  ";
                            }
                        }
                        if (gRoute[vod].size() > 0){
                            uIter = --gRoute[vod].end();
                            ofs << "0 (" << (*uIter)->timeDeparture + gDistance[(*uIter)->cus->id][0] << ")\n";
                        }
                    }
                }
                ofs << "====================================" << endl;
            break;
        }
        ofs.close();
    }else{
        cerr << "Can not open \"" << fileout << "\" for output.\n";
        exit(1);
    }

    return 0;
}

/*
 * Here we cache distances between each customer and the others
 * ====================================
 */
void cacheDistances(vector<vector<double> >& gDistance, vector<Customer>& gArrC, unsigned int n){
    unsigned int i = 0;
    unsigned int j = 0;
    gDistance.resize(n+1);
    for (i = 0; i < n+1; i++){
        gDistance[i].resize(n+1);
        for (j = 0; j < n+1; j++){
            if (i == j){
                gDistance[i][j] = 0;
            }else{
                if (i < j){
                    if (i == 0){
                        gDistance[i][j] = distance(gDepot, &gArrC[j-1]);
                    }else{
                        gDistance[i][j] = distance(&gArrC[i-1], &gArrC[j-1]);
                    }
                }else{
                    gDistance[i][j] = gDistance[j][i];
                }
            }
        }
    }
}

int main(int argc, char** argv) {
    fstream ifs;
    fstream ofs;
    if (argc < 3) {
        cerr << "Invalid arguments! Please type names of input, settings and output files." << endl;
        cerr << "Usage:\n hgafinal [input] [settings] [output]" << endl;
        exit(1);
    } else {
        // read input file
        getInputData(ifs, argv[1]);
        cacheDistances(gDistance, gArrC, gnCustomers);
        HGAGenome genome(0);
        HPGradProjectGA ga(genome);
        HPGradScaling scaling;
        HPGradSelector select;

        ga.parameters(argv[2]); // read parameters from settings file
        nKeep = (int)(ga.populationSize() - ga.nReplacement());
        nPop = 2*((int)ga.nReplacement()/2);
        if (nPop > 2){
            nPop -= 2;
        }
        ga.minimize();          // minimize objective function
        ga.scaling(scaling);
        ga.selector(select);

        cout << "\nEvolving..." << endl;
        cpu_time();
        ga.evolve();
        double totalTime = cpu_time();

        if (bestFeasibleCost == 0){
            bestSol = (HGAGenome&) (ga.statistics().bestIndividual());
        }
        cout << "\nHGA finished! Total time: " << totalTime << "(sec)\n";
        writeOutputData(ofs, argv[3], argv[1], totalTime);

        exit(1);
    }
    return 0;
}
