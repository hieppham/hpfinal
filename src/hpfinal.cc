//============================================================================
// Name        : hpfinal.cc
// Author      : Hiep Pham
// Version     : 1.0
// Copyright   : (C) 2012 by _HP
// Description : Main function of HP Final Project in C++, Ansi-style
//============================================================================

#include "hpfinal.h"
#include "cpu_time.c"
using namespace std;

// global variables
unsigned int numOfDepots = 1;
int typeOfVRP;
// initialize static variables
unsigned int HPGV::nCus = 1;
unsigned int HPGV::mVeh = 1;
unsigned int HPGV::tDay = 0;
unsigned int HPGV::numRoute = 0;
unsigned int HPGV::nPop = 60;
unsigned int HPGV::nKeep = 40;

double HPGV::maxDuration = 0;
double HPGV::maxLoad = 0;

double HPGV::hPenalty = 0;
double HPGV::avgQ = 0;
double HPGV::avgD = 0;
double HPGV::avgW = 0;

double HPGV::bestFeasibleCost = 0;

double HPGV::I1Lambda = 0;
double HPGV::I1Mu = 0;
double HPGV::I1Alpha = 0;

double HPGV::SGamma2 = 0;
double HPGV::SGamma3 = 0;

long HPGV::utsIter = 80;
double HPGV::utsLambda = 0.015;

unsigned int HPGV::rvnsIter = 100;
double HPGV::rvnsPRev = 0.1;
double HPGV::rvnsTmax = 600;
HGAGenome bestSol(0);

vector<Customer> HPGV::gArrC(0);
vector<unsigned int> HPGV::multipattern(0);
unsigned int HPGV::numOfMP = 0;
Customer* HPGV::gDepot;
Customer* HPGV::gDynamic;
double HPGV::gDynamicStart = 0;

vector<vector<double> > HPGV::gDistance(1, vector<double>(0));

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
            ifs >> typeOfVRP >> HPGV::mVeh >> HPGV::nCus >> HPGV::tDay >> numOfDepots;
        }
        HPGV::numRoute = HPGV::tDay * HPGV::mVeh;
        for (unsigned int i = 0; i < HPGV::tDay; i++) {
            ifs >> HPGV::maxDuration >> HPGV::maxLoad;
        }
        if (HPGV::maxDuration == 0){
            HPGV::maxLoad = 10000;
        }
        ifs >> tid >> tx >> ty >> td >> tq >> tf >> ta;
        HPGV::gDepot = new Customer(tid, tx, ty, td, tq, tf, ta);
        HPGV::gDepot->a = HPGV::gDepot->e = HPGV::gDepot->pattern = 0;
        HPGV::gDepot->comb.clear();
        ifs >> te >> tl;
        HPGV::gDepot->setTime(te, tl);

        while (ifs.good() && HPGV::gArrC.size() < HPGV::nCus) {
            tempIter = 0;
            ifs >> tid >> tx >> ty >> td >> tq >> tf >> ta;
            Customer* gCustomer = new Customer(tid, tx, ty, td, tq, tf, ta);
            HPGV::gArrC.push_back(*gCustomer);
            HPGV::gArrC[tid - 1].comb.clear();
            if (ta > 1){
                HPGV::multipattern.push_back(tid);
            }
            for (tempIter = 0; tempIter < ta; tempIter++) {
                ifs >> tc;
                HPGV::gArrC[tid - 1].comb.push_back(tc);
            }
            HPGV::gArrC[tid - 1].comb.resize(ta);
            ifs >> te >> tl;
            HPGV::gArrC[tid - 1].setTime(te, tl);
            delete gCustomer;
        }
        HPGV::numOfMP = HPGV::multipattern.size();
        ifs.close();
    }else{
        cerr << "Can not open \"" << filein << "\" for input.\n";
        exit(1);
    }
    return 0;
}

int getParams(fstream &ifs, char* filein) {
    ifs.open(filein, ios::in);
    if (ifs.is_open()) {
        try {
            ifs >> HPGV::I1Lambda >> HPGV::I1Mu >> HPGV::I1Alpha;
            ifs >> HPGV::utsIter >> HPGV::utsLambda;
            ifs >> HPGV::rvnsIter >> HPGV::rvnsPRev >> HPGV::rvnsTmax;
            ifs >> HPGV::SGamma2 >> HPGV::SGamma3;
            ifs.close();
        } catch (...) {
            cerr << "Can not read parameters \"" << filein << "\" from file.\n";
            exit(1);
        }

    }else{
        cerr << "Can not open \"" << filein << "\" for getting parameters.\n";
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
            if (!bestSol.isFeasible){
                ofs << "* ";
            }
            ofs << bestSol.durationCost << endl;
                for (iDay = 0; iDay < HPGV::tDay; iDay++){
                    for (iVeh = 0; iVeh < HPGV::mVeh; iVeh++){
                        vod = iDay*HPGV::mVeh + iVeh;
                        Route::const_iterator uIter, endIter;
                        for (uIter = bestSol.m_route[vod].begin(), endIter = bestSol.m_route[vod].end(); uIter != endIter; ++uIter){
                            if (uIter == bestSol.m_route[vod].begin()){
                                ofs << (iDay+1) << "  " << (iVeh + 1) << "\t";
                                ofs << bestSol.m_data[vod]->cost << "\t" << bestSol.m_data[vod]->load << "\t";
                                ofs << "0(" << bestSol.m_data[vod]->timeLeaveDepot << ")  ";
                                ofs << (*uIter)->cus->id << "[" << bestSol.m_pattern[(*uIter)->cus->id - 1] << "]" << "(" << (*uIter)->timeStartService << ")  ";
                            }else{
                                ofs << (*uIter)->cus->id << "[" << bestSol.m_pattern[(*uIter)->cus->id - 1] << "]" << "(" << (*uIter)->timeStartService << ")  ";
                            }
                        }
                        if (bestSol.m_route[vod].size() > 0){
                            uIter = --bestSol.m_route[vod].end();
                            ofs << "0 (" << (*uIter)->timeDeparture + HPGV::gDistance[(*uIter)->cus->id][0] << ")\n";
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
void cacheDistances(){
    unsigned int i = 0;
    unsigned int j = 0;
    HPGV::gDistance.resize(HPGV::nCus+1);
    for (i = 0; i <= HPGV::nCus; i++){
        HPGV::gDistance[i].resize(HPGV::nCus+1);
        for (j = 0; j <= HPGV::nCus; j++){
            if (i == j){
                HPGV::gDistance[i][j] = 0;
            }else{
                if (i < j){
                    if (i == 0){
                        HPGV::gDistance[i][j] = distance(HPGV::gDepot, &HPGV::gArrC[j-1]);
                    }else{
                        HPGV::gDistance[i][j] = distance(&HPGV::gArrC[i-1], &HPGV::gArrC[j-1]);
                    }
                }else{
                    HPGV::gDistance[i][j] = HPGV::gDistance[j][i];
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
        cerr << "Usage:\n hgafinal [input] [settings] [params] [output]" << endl;
        exit(1);
    } else {
        // read input file
        getInputData(ifs, argv[1]);
        getParams(ifs, argv[2]);
        cacheDistances();
        HGAGenome genome(0);
        HPGradProjectGA ga(genome);
        HPGradSelector select;
        HPGradScaling scaling;

        ga.parameters(argv[2]); // read parameters from settings file
        HPGV::nKeep = (int)(ga.populationSize() - ga.nReplacement());
        HPGV::nPop = 2*((int)ga.nReplacement()/2);
        if (HPGV::nPop > 2){
            HPGV::nPop -= 2;
        }
        ga.minimize();          // minimize objective function
        ga.selector(select);
        ga.scaling(scaling);
        cout << "\nEvolving..." << endl;
        cpu_time();
        ga.evolve();
        double totalTime = cpu_time();

        if (HPGV::bestFeasibleCost == 0){
            bestSol = (HGAGenome&) (ga.statistics().bestIndividual());
        }
        cout << "\nHGA finished! Total time: " << totalTime << "(sec)\n";
        writeOutputData(ofs, argv[4], argv[1], totalTime);

        exit(1);
    }
    return 0;
}
