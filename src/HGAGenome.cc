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
    if (randType == 1){
        hg.clusterFirstInit(refArr);
    }else{
        hg.SolomonTONNInit(refArr);
    }
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

    if(c){
        HGAGenome& sis = (HGAGenome &) *c;
        HGAGenome::explorationCrossover(p1, p2, sis);
        numOffsping++;
    }
    // and second offspring is created similarly
    if (d){
        HGAGenome& bro = (HGAGenome &) *d;
        HGAGenome::exploitationCrossover(p1, p2, bro);
        numOffsping++;
    }

    return numOffsping;
}

int HGAGenome::explorationCrossover(const HGAGenome& p1, const HGAGenome& p2, HGAGenome& child){
    // TODO: explorationCrossover
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
