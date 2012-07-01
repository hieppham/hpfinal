/*
 * HGAGenome.cc
 *
 *  Created on: Jun 24, 2012
 *      Author: hieppham
 */

#include "HGAGenome.h"

HGAGenome::HGAGenome(int initCost) : GAGenome(Initializer, Mutator) {
    evaluator(Evaluator);
    crossover(Crossover);
    durationCost = initCost;
}
HGAGenome::~HGAGenome() {
    // TODO Auto-generated destructor stub
}

void HGAGenome::copy(const GAGenome& g) {
    if (&g != this && sameClass(g)) {
        GAGenome::copy(g); // copy the base class part
        HGAGenome & hgenome = (HGAGenome &)g;
        m_route = hgenome.m_route;

        m_pattern = hgenome.m_pattern;
        m_tour = hgenome.m_tour;

        durationCost = hgenome.durationCost;
        isFeasible = hgenome.isFeasible;
        totalTimeVio = hgenome.totalTimeVio;
        totalCapacityVio = hgenome.totalCapacityVio;
        totalDurationVio = hgenome.totalDurationVio;
    }
}

GAGenome*
HGAGenome::clone(GAGenome::CloneMethod) const {
    return new HGAGenome(*this);
}

int HGAGenome::write(ostream & os) const {
    os << durationCost << "\n";
    return os.fail() ? 1 : 0;
}

/**
 * Initialize genome
 */
void HGAGenome::Initializer(GAGenome& g) {
    // TODO: complete code here.
    cout << "Init" << endl;
    exit(1);
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
