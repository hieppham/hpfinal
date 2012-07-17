/*
 * HGAGenome.cc
 *
 *  Created on: Jun 24, 2012
 *      Author: hieppham
 */

#include "HGAGenome.h"

extern vector<Customer> gArrC;
extern vector<vector<double> > gDistance;
/**
 * Some methods of class CustomerInDay
 */
CustomerInDay::CustomerInDay(int tc, unsigned int tv){
    cid = tc;
    vod = tv;
}
CustomerInDay::~CustomerInDay(){

}
/**
 * Initialize genome
 */
void HGAGenome::Initializer(GAGenome& g) {
    unsigned int ic = 1; // index of customer
    // TODO: complete code here.
    HGAGenome & hg = (HGAGenome &) g;
    vector<Customer*> refArr;
    refArr.resize(HPGV::nCus);
    hg.arrC = gArrC;
    hg.m_route.resize(HPGV::numRoute);
    hg.m_data.resize(HPGV::numRoute);
    for (ic = 0; ic < HPGV::nCus; ++ic) {
        int idx = GARandomInt(0, hg.arrC[ic].a - 1);
        hg.arrC[ic].randomAssignedPattern(idx);
        refArr[ic] = &(hg.arrC[ic]);
    }
    int randType = GARandomInt(1, 3);
    if (randType > 1){
        hg.clusterFirstInit(refArr);
    }else{
        hg.SolomonTONNInit(refArr);
    }


    hg.m_pattern.resize(HPGV::nCus);
    for (unsigned int i = 0; i < HPGV::nCus; i++){
        hg.m_pattern[i] = hg.arrC[i].pattern;
    }
    for (unsigned int vod = 0; vod < (HPGV::numRoute); vod++){
        for (Route::iterator rIter = hg.m_route[vod].begin(), endIter = hg.m_route[vod].end(); rIter != endIter; ++rIter){
            hg.m_tour.push_back(CidPtr(new CustomerInDay((*rIter)->cus->id, vod)));
        }
    }

    hg.updateTotalVio();
}

float HGAGenome::Evaluator(GAGenome& g) {
    HGAGenome & hgenome = (HGAGenome &) g;
    // TODO: Evaluator

    return 0;
}

int HGAGenome::Crossover(const GAGenome& a, const GAGenome& b, GAGenome* c, GAGenome* d) {
    int numOffsping = 0;


    HGAGenome& p1 = (HGAGenome &)a;
    HGAGenome& p2 = (HGAGenome &)b;
//    HGAGenome::printSolution(p1, "parent1.txt");
//    HGAGenome::printSolution(p2, "parent2.txt");

    if (c){
        HGAGenome& sis = (HGAGenome &) *c;
        sis.arrC = gArrC;
        HGAGenome::explorationCrossover(p1, p2, sis);
        numOffsping++;
    }
    // and second offspring is created similarly
    if (d){
        HGAGenome& bro = (HGAGenome &) *d;
        bro.arrC = gArrC;
        HGAGenome::exploitationCrossover(p1, p2, bro);
        numOffsping++;
    }

    return numOffsping;
}

int HGAGenome::explorationCrossover(const HGAGenome& p1, const HGAGenome& p2, HGAGenome& child){
    unsigned int iDay = 0;
    unsigned int iVeh = 0;
    unsigned int vod = 0;
    unsigned int tempSize = 0;
    double newStartTime = 0;
    vector<vector<int> > checkInherit(0, vector<int>(0));
    vector<int>::iterator findPos;

    VCus backupArr;
    VCus cloneArr;  // copy of current cluster
    VCus::iterator minVPos, skipCus, kIter;

    Route::iterator prevCus, nextCus, uIter, beforeIter;
    checkInherit.resize(HPGV::tDay);

    // STEP 1: Assign a pattern to each customer
    // inherit pattern assignments from parent P1
    tempSize = p1.m_tour.size();
    child.m_pattern.resize(HPGV::nCus);

    int cutPointOne = GARandomInt(0, tempSize - 2);
    int cutPointTwo = GARandomInt(cutPointOne + 1, tempSize - 1);

    for (int i = cutPointOne; i <= cutPointTwo; i++){
        int lcid = p1.m_tour[i]->cid - 1;
        iDay = (int) (p1.m_tour[i]->vod / HPGV::mVeh);
        child.m_pattern[lcid] |= (int) pow(2, (double)(HPGV::tDay - iDay - 1));
        child.arrC[lcid].pattern = child.m_pattern[lcid];
        child.arrC[lcid].token++;
    }
    // inherit pattern assignments from pattern P2
    for (iDay = 0; iDay < HPGV::tDay; iDay++){
        for (unsigned int i = 0; i < HPGV::nCus; i++){
            if (child.arrC[i].token < child.arrC[i].f){
                int flag = (int) pow(2, (double)(HPGV::tDay - iDay - 1));
                int flagBkp = flag;
                for (int j = 0; j < child.arrC[j].a; j++){
                    flag &= child.arrC[i].comb[j];
                    if (flag == 0){
                        continue;
                    }else{
                        flagBkp |= child.arrC[i].pattern;
                        flag = flagBkp;
                        flagBkp &= child.arrC[i].comb[j];
                        if (flagBkp != 0){
                            child.arrC[i].pattern = flag;
                            child.arrC[i].token = numberOfSetBits(flag);
                        }
                    }
                }
            }
        }
    }
    // complete pattern assignments
    for (unsigned int i = 0; i < HPGV::nCus; i++){
        if (!child.arrC[i].isServiced){
            // assign a random pattern that includes day(i)
            int ord = GARandomInt(0, child.arrC[i].a - 1);
            while(1){
                int flagBkp = child.arrC[i].pattern;
                flagBkp &= child.arrC[i].comb[ord];
                if (flagBkp !=0){
                    child.arrC[i].pattern = child.arrC[i].comb[ord];
                    child.m_pattern[i] = child.arrC[i].pattern;
                    child.arrC[i].isServiced = true;
                    break;
                }
                ord++;
                ord = ord % child.arrC[i].a;
            }
        }
    }

    // tempSize = HPGV::numRoute;
    // STEP 2 : Assign customers to routes
    // Copy routes from P1
    // initialize all routes
    child.m_route.resize(HPGV::numRoute);
    child.m_data.resize(HPGV::numRoute);
    for (unsigned int tm = 0; tm < HPGV::numRoute; tm++){
        child.m_data[tm] = RinfoPtr(new RouteInfo());
    }

    backupArr.resize(HPGV::nCus);
    for (unsigned int ic = 0; ic < HPGV::nCus; ic++){
        backupArr[ic] = &(child.arrC[ic]);
    }

    // The customers between the two cutting points determined in Step 1a are routed
    // as in the P1 parent and the corresponding sequences are copied into C
    for (int i = cutPointOne; i <= cutPointTwo; i++){
        int s2cid = p1.m_tour[i]->cid - 1;
        int s2vod = p1.m_tour[i]->vod;
        int tmpDay = s2vod/HPGV::mVeh;

        HGAGenome::pushbackRoute(child.m_route[s2vod], child.m_data[s2vod], backupArr[s2cid]);

        checkInherit[tmpDay].push_back(p1.m_tour[i]->cid);
    }

    // Assign remaining customers to routes
    for (iDay = 0; iDay < HPGV::tDay; iDay++) {
        cloneArr.clear();
        cloneArr = backupArr;
        // check patterns and remove all customers that should not be visited on current day
        // of course, we also remove all customers that has been satisfied
        for (skipCus = cloneArr.begin(); skipCus != cloneArr.end(); ++skipCus) {
            int flag = (int) pow(2, (double) (HPGV::tDay - iDay - 1));
            flag &= (*skipCus)->pattern;
            if (flag == 0) {
                cloneArr.erase(skipCus);
                skipCus--;
            }else{
                findPos = find(checkInherit[iDay].begin(), checkInherit[iDay].end(), (*skipCus)->id);
                if (findPos != checkInherit[iDay].end()){
                    cloneArr.erase(skipCus);
                    skipCus--;
                    checkInherit[iDay].erase(findPos);
                }
            }
        }
        // here we use cost insertion heuristic (Potvin - Rosseau, 1993)
        // sort(cloneArr.begin(), cloneArr.end(), compareEndingTime);
        HGAGenome::PRheuristic(child.m_route, child.m_data, cloneArr, iDay, true);
    }

    // complete tour genome
    for (vod = 0; vod < (HPGV::numRoute); vod++){
        for (Route::iterator rIter = child.m_route[vod].begin(), endIter = child.m_route[vod].end(); rIter != endIter; ++rIter){
            child.m_tour.push_back(CidPtr(new CustomerInDay((*rIter)->cus->id, vod)));
        }
    }

    child.updateTotalVio();
    HGAGenome::printSolution(child, "childExplor.txt");
    child._evaluated = gaFalse;

    return 1;
}

int HGAGenome::exploitationCrossover(const HGAGenome& p1, const HGAGenome& p2, HGAGenome& child){
    // TODO: exploitationCrossover
    child._evaluated = gaFalse;

    return 1;
}

int HGAGenome::Mutator(GAGenome& g, float pMut) {
    // TODO: Mutator
    int nMut = 0;
    HGAGenome & hgenome = (HGAGenome &) g;
    if(nMut) hgenome._evaluated = gaFalse;
    return nMut;
}
