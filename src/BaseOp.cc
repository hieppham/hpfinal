#include "HGAGenome.h"

extern vector<Customer> gArrC;
extern vector<vector<double> > gDistance;

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
