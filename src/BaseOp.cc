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

/**
 * Insert a customer at last position of a route
 */
void HGAGenome::pushbackRoute(Route& mRoute, RinfoPtr& mRinfo, Customer* mCus){
    if (mRoute.empty()){
        mRoute.push_back(VertexPtr(new Vertex(mCus)));
        mRoute.back()->timeArrive = gDistance[0][mCus->id];
        mRinfo->load = mCus->q;
        mRinfo->cost = gDistance[0][mCus->id];
        mRoute.back()->timeStartService = max(mCus->e, gDistance[0][mCus->id]);
        mRoute.back()->timeWait = mCus->e - mRoute.back()->timeStartService;
        mRoute.back()->timeDeparture = mRoute.back()->timeStartService + mCus->d;
    }else{
        double newTimeArrive = mRoute.back()->timeDeparture + gDistance[mRoute.back()->cus->id][mCus->id];
        double newTimeStart = max(mCus->e, newTimeArrive);

        mRinfo->load += mCus->q;
        mRinfo->cost += gDistance[mRoute.back()->cus->id][mCus->id];
        mRoute.push_back(VertexPtr(new Vertex(mCus)));
        mRoute.back()->timeArrive = newTimeArrive;
        mRoute.back()->timeStartService = newTimeStart;
        mRoute.back()->timeWait = mCus->e - mRoute.back()->timeStartService;
        mRoute.back()->timeDeparture = newTimeStart + mCus->d;
    }
    if (mRoute.back()->timeStartService > mCus->l){
        mRinfo->timeVio += (mRoute.back()->timeStartService - mCus->l);
    }
}
/**
 * route construction using Solomon heuristic, the first variant (I1)
 */
void HGAGenome::SolomonI1(Route& mRoute, RinfoPtr& mRinfo, VCus& arrCus){
    double maxProfit = -1000000;    // negative infinite
    double newEarliestTime = 0;
    double newLastestTime = 0;
    double newProfit;
    unsigned int pTrace = 0;
    unsigned int currPos = 0;
    VCus::iterator minVPos;
    Route::iterator nextCus, prevCus, endOfRoute, uIter, beforeIter;

    for (VCus::iterator kIter = arrCus.begin(), endArr = arrCus.end(); kIter != endArr; ++kIter){
        // evaluate all possible move
        for (nextCus = mRoute.begin(), endOfRoute = mRoute.end(); nextCus != endOfRoute; ++nextCus){
             currPos = nextCus - mRoute.begin();
             bool legal = true;
             if (currPos == 0){
                 newEarliestTime = max((*kIter)->e, gDistance[0][(*kIter)->id]);
                 newLastestTime = min((*kIter)->l, gDepot->l);
                 (*nextCus)->newTimeArrive = newEarliestTime + (*kIter)->d
                         + gDistance[(*kIter)->id][(*nextCus)->cus->id];
                 (*nextCus)->newTimeStart = max((*nextCus)->newTimeArrive, (*nextCus)->cus->e);
                 // check this insertion is feasible or not
                 if ((*nextCus)->newTimeStart > (*nextCus)->cus->l){
                     legal = false;
                 }
                 (*nextCus)->pf = (*nextCus)->newTimeStart - (*nextCus)->timeStartService;
                 uIter = nextCus;
                 beforeIter = nextCus;
                 uIter++;
                 while((uIter != endOfRoute) && legal){
                     (*uIter)->pf = max(0, (*beforeIter)->pf);
                     if ((*uIter)->timeStartService + (*uIter)->pf > (*uIter)->cus->l){
                         legal = false;
                     }
                     uIter++;
                     beforeIter++;
                 }
                 if (legal){
                     // calculate c_2(i, u, j)
                     double c11 = gDistance[0][(*kIter)->id]
                                               + gDistance[(*kIter)->id][(*nextCus)->cus->id]
                                               - HPGV::I1Mu * gDistance[0][(*nextCus)->cus->id];
                     double c12 = (*nextCus)->newTimeStart - (*nextCus)->timeStartService;
                     newProfit = HPGV::I1Lambda * gDistance[0][(*kIter)->id]
                                           - HPGV::I1Alpha * c11
                                           - (1 - HPGV::I1Alpha) * c12;
                     if (newProfit > maxProfit){
                         maxProfit = newProfit;
                         pTrace = currPos;
                         minVPos = kIter;
                     }
                 }
             }else{
                 bool legal = true;
                 prevCus = nextCus;
                 prevCus--;
                 newEarliestTime = max((*kIter)->e,
                         (*prevCus)->timeStartService + (*prevCus)->cus->d
                         + gDistance[(*prevCus)->cus->id][(*kIter)->id]);
                 newLastestTime = min((*kIter)->l,
                         (*nextCus)->cus->l - (*kIter)->d - gDistance[(*nextCus)->cus->id][(*kIter)->id]);
                 (*nextCus)->newTimeArrive = newEarliestTime + (*kIter)->d
                         + gDistance[(*kIter)->id][(*nextCus)->cus->id];
                 (*nextCus)->newTimeStart = max((*nextCus)->newTimeArrive, (*nextCus)->cus->e);
                 // check this insertion is feasible or not
                 if ((*nextCus)->newTimeStart > (*nextCus)->cus->l){
                     legal = false;
                 }
                 (*nextCus)->pf = (*nextCus)->newTimeStart - (*nextCus)->timeStartService;
                 uIter = nextCus;
                 beforeIter = nextCus;
                 uIter++;
                 while ((uIter != endOfRoute) && (legal)){
                     (*uIter)->pf = max(0, (*beforeIter)->pf);
                     if ((*uIter)->timeStartService + (*uIter)->pf > (*uIter)->cus->l){
                         legal = false;
                     }
                     uIter++;
                     beforeIter++;
                 }
                 if (legal){
                     // calculate c_2(i, u, j)
                     double c11 = gDistance[(*prevCus)->cus->id][(*kIter)->id]
                                       + gDistance[(*kIter)->id][(*nextCus)->cus->id]
                                       - HPGV::I1Mu*gDistance[(*prevCus)->cus->id][(*nextCus)->cus->id];
                     double c12 = (*nextCus)->newTimeStart - (*nextCus)->timeStartService;
                     newProfit = HPGV::I1Lambda * gDistance[0][(*kIter)->id] - HPGV::I1Alpha * c11 - (1 - HPGV::I1Alpha) * c12;
                     if (newProfit > maxProfit){
                         maxProfit = newProfit;
                         pTrace = currPos;
                         minVPos = kIter;
                     }
                 }
             }
             // TODO: insert at the end of route

        }

    }
}
