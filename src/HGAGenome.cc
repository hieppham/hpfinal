/*
 * HGAGenome.cc
 *
 *  Created on: Jun 24, 2012
 *      Author: hieppham
 */

#include "HGAGenome.h"
extern HGAGenome bestSol;

/**
 * Some methods of class CustomerInDay
 */
CustomerInDay::CustomerInDay(int tc, unsigned int tv){
    cid = tc;
    vod = tv;
}
CustomerInDay::~CustomerInDay(){
    // cout << "destructing cid: " << cid << " - " << vod << "\n";
}
/**
 * Initialize genome
 */
void HGAGenome::Initializer(GAGenome& g) {
    unsigned int ic = 1; // index of customer
    HGAGenome & hg = (HGAGenome &) g;
    vector<Customer*> refArr;
    refArr.resize(HPGV::nCus);
    hg.arrC = HPGV::gArrC;
    hg.m_route.clear();
    hg.m_data.clear();
    hg.m_route.resize(HPGV::numRoute);
    hg.m_data.resize(HPGV::numRoute);
    for (ic = 0; ic < HPGV::nCus; ++ic) {
        int idx = GARandomInt(0, hg.arrC[ic].a - 1);
        hg.arrC[ic].randomAssignedPattern(idx);
        refArr[ic] = &(hg.arrC[ic]);
    }
    int randType = GARandomInt(1, 5);
    if (randType > 2){
        hg.clusterFirstInit(refArr);
    }else{
        if (randType == 1){
            hg.SolomonTONNInit(refArr);
        }else{
            hg.PRInit(refArr);
        }
    }

    hg.tourConstruct();
    HGAGenome::improveRoute(hg);

    hg.updateTotalVio();
}

float HGAGenome::Evaluator(GAGenome& g) {
    HGAGenome & hgenome = (HGAGenome &) g;
    float s = (float)(HGAGenome::calcObjectValue(hgenome));
    if (hgenome.isFeasible){
        cout << "***********************************************Acceptable with cost = " << hgenome.durationCost << endl;
        if (HPGV::bestFeasibleCost == 0){
            HPGV::bestFeasibleCost = hgenome.durationCost;
            bestSol = hgenome;
        }else{
            if (HPGV::bestFeasibleCost > hgenome.durationCost){
                HPGV::bestFeasibleCost = hgenome.durationCost;
                bestSol = hgenome;
            }
        }
    }
    return s;
}

int HGAGenome::Crossover(const GAGenome& a, const GAGenome& b, GAGenome* c, GAGenome* d) {
    int numOffsping = 0;


    HGAGenome& p1 = (HGAGenome &)a;
    HGAGenome& p2 = (HGAGenome &)b;
//    HGAGenome::printSolution(p1, "parent1.txt");
//    HGAGenome::printSolution(p2, "parent2.txt");

    if (c){
        HGAGenome& sis = (HGAGenome &) *c;
        sis.resetAll();
        HGAGenome::explorationCrossover(p1, p2, sis);
        numOffsping++;
    }
    // and second offspring is created similarly
    if (d){
        HGAGenome& bro = (HGAGenome &) *d;
        bro.resetAll();
        HGAGenome::exploitationCrossover(p1, p2, bro);
        numOffsping++;
    }

    return numOffsping;
}

int HGAGenome::Mutator(GAGenome& g, float pMut) {
    int nMut = 0;
    unsigned int oldPattern, newPattern, insertMask, removeMask;
    unsigned int iDay = 0;
    unsigned int iVeh = 0;
    unsigned int vod = 0;
    VCus tmpCus(0);

    HGAGenome & hg = (HGAGenome &) g;

    for (unsigned int i = 0; i < HPGV::nCus; i++){
        // changing the pattern assignment of a few customers, which are selected through a low probability
        if (GAFlipCoin(pMut)){
            // assign new pattern for selected customer
            // skip if this customer has only one pattern
            if (hg.arrC[i].a > 1){
                nMut++;
                Customer* mixer = &(hg.arrC[i]);
                oldPattern = hg.m_pattern[i];
                unsigned int ord = GARandomInt(0, hg.arrC[i].a - 1);
                while(1){
                    if (hg.arrC[i].comb[ord] != oldPattern){
                        newPattern = hg.arrC[i].comb[ord];
                        break;
                    }
                    ord++;
                    ord %= hg.arrC[i].a;
                }
                hg.m_pattern[i] = newPattern;
                insertMask = newPattern;
                removeMask = oldPattern;
                insertMask ^= (newPattern & oldPattern);
                removeMask ^= (newPattern & oldPattern);

                // update genome by insert customer into new routes and remove it from old routes
                for (iDay = 0; iDay < HPGV::tDay; iDay++){
                    int flagInsert = (int) pow(2, (double) (HPGV::tDay - iDay - 1));
                    int flagRemove = flagInsert;
                    flagInsert &= insertMask;
                    flagRemove &= removeMask;

                    if (flagInsert){
                        // check if customer has already existed or not
                        bool needServiced = true;
                        for (iVeh = 0; iVeh < HPGV::mVeh; iVeh++){
                            unsigned int currVod = iDay * HPGV::mVeh + iVeh;
                            if (HGAGenome::isInRoute(hg.m_route[currVod], mixer->id)){
                                needServiced = false;
                                break;
                            }
                        }
                        // here we insert customer into one route of this day
                        if (needServiced){
                            tmpCus.clear();
                            tmpCus.push_back(mixer);
                            HGAGenome::PRheuristic(hg.m_route, hg.m_data, tmpCus, iDay, false);
                        }
                    }else if (flagRemove){
                        for (iVeh = 0; iVeh < HPGV::mVeh; iVeh++){
                            vod = iDay * HPGV::mVeh + iVeh;
                            if (HGAGenome::isInRoute(hg.m_route[vod], mixer->id)){
                                HGAGenome::removeFromRoute(hg.m_route[vod], hg.m_data[vod], mixer->id);
                                break;
                            }
                        }
                    }
                }
            }
        }
    }

    hg.tourConstruct();
    hg.updateTotalVio();

    // HGAGenome::printSolution(hg, "Mutation.txt");
    // cout << "Mutation\n";

    if(nMut) hg._evaluated = gaFalse;
    return nMut;
}

/**
 * Education - improve quality of offsprings
 */
int HGAGenome::Education(GAGenome& g, const int CNG){
    HGAGenome & hg = (HGAGenome &) g;

    if ((CNG % 2) == 0){
        hg = HGAGenome::UTS(hg);
    }else{
        hg = HGAGenome::RVNS(hg);
    }
    HGAGenome::improveRoute(hg);

    return 0;
}
