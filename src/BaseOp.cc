#include "HGAGenome.h"

extern Customer* gDepot;
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
 * Insert a customer at a defined position in a route (pTrace = -1 means end of route)
 */
void HGAGenome::insertIntoRoute(Route& mRoute, RinfoPtr& mRinfo, Customer* mCus, unsigned int& pTrace){
    if (pTrace == -1){
        HGAGenome::pushbackRoute(mRoute, mRinfo, mCus);
    }else{
        Route::iterator pos = mRoute.begin();
        for (unsigned int i = 0; i < pTrace; i++){
            pos++;
        }
        mRoute.insert(pos, VertexPtr(new Vertex(mCus)));
        // TODO: update information for all customers
        Route::iterator endOfRoute = mRoute.end();
        for (Route::iterator uIter = pos; uIter != endOfRoute; ++uIter){
            if (pTrace == 0){
                (*uIter)->timeArrive = gDistance[0][(*uIter)->cus->id];
                (*uIter)->timeStartService = max((*uIter)->cus->e, gDistance[0][(*uIter)->cus->id]);
                (*uIter)->timeWait = (*uIter)->timeStartService - (*uIter)->timeArrive;
                (*uIter)->timeDeparture = (*uIter)->timeStartService + (*uIter)->cus->d;
                mRinfo->cost += gDistance[0][(*uIter)->cus->id];
                mRinfo->load += (*uIter)->cus->q;
            }else{
                Route::iterator beforeIter = uIter;
                beforeIter--;
                (*uIter)->timeArrive = (*beforeIter)->timeStartService + (*beforeIter)->cus->d
                                                        + gDistance[(*uIter)->cus->id][(*beforeIter)->cus->id];

                mRinfo->cost += gDistance[(*uIter)->cus->id][(*beforeIter)->cus->id];
                if ((*uIter)->timeArrive < (*uIter)->cus->e){
                    (*uIter)->timeStartService = (*uIter)->cus->e;
                    (*uIter)->timeWait = (*uIter)->cus->e - (*uIter)->timeArrive;
                }else{
                    (*uIter)->timeStartService = (*uIter)->timeArrive;
                    if ((*uIter)->timeArrive > (*uIter)->cus->l){
                        mRinfo->timeVio += ((*uIter)->timeArrive - (*uIter)->cus->l);
                    }
                }
                (*uIter)->timeDeparture = (*uIter)->timeStartService + (*uIter)->cus->d;
                mRinfo->load += (*uIter)->cus->q;
            }
        }
        // link the last customer with depot
        double newStartTime = mRoute.back()->timeDeparture + gDistance[0][mRoute.back()->cus->id];
        mRinfo->cost += gDistance[0][mRoute.back()->cus->id];

        if (newStartTime > gDepot->l){
            mRinfo->timeVio += (newStartTime - gDepot->l);
        }
    }
}
/**
 * route construction using Solomon heuristic, the first variant (I1)
 */
void HGAGenome::SolomonI1(Route& mRoute, RinfoPtr& mRinfo, VCus& arrCus){
    double maxProfit = NEGINF;    // negative infinite
    double newEarliestTime = 0;
    double newLastestTime = 0;
    double newProfit;
    unsigned int pTrace = 0;
    unsigned int currPos = 0;
    VCus::iterator minVPos;
    Route::iterator nextCus, prevCus, endOfRoute, uIter, beforeIter;

    while(1){
        maxProfit = NEGINF;
        for (VCus::iterator kIter = arrCus.begin(), endArr = arrCus.end(); kIter != endArr; ++kIter){
            // evaluate all possible move
            for (nextCus = mRoute.begin(), endOfRoute = mRoute.end(), currPos = 0; nextCus != endOfRoute; ++nextCus){
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
                 currPos++;

            }
            // insert at the end of route
            prevCus = endOfRoute;
            prevCus--;
            newEarliestTime = max((*kIter)->e, (*prevCus)->timeStartService
                                     + gDistance[(*prevCus)->cus->id][(*kIter)->id]);
            newLastestTime = min((*kIter)->l, gDepot->l);
            if (newEarliestTime < newLastestTime){
                double c11 = gDistance[(*prevCus)->cus->id][(*kIter)->id] + gDistance[(*kIter)->id][0] - HPGV::I1Mu * gDistance[(*prevCus)->cus->id][0];
                double c12 = newEarliestTime + (*kIter)->d + gDistance[(*kIter)->id][0] - (*prevCus)->timeStartService - (*prevCus)->cus->d - gDistance[(*prevCus)->cus->id][0];
                newProfit = HPGV::I1Lambda * gDistance[0][(*kIter)->id] - HPGV::I1Alpha * c11 - (1 - HPGV::I1Alpha) * c12;
                if (newProfit > maxProfit){
                    maxProfit = newProfit;
                    pTrace = -1;
                    minVPos = kIter;
                }
            }
        }
        // if found feasible customer and best move
        if (maxProfit != NEGINF){
            HGAGenome::insertIntoRoute(mRoute, mRinfo, *minVPos, pTrace);
            (*minVPos)->checkServiced();
            arrCus.erase(minVPos);
            if (arrCus.empty()){
                break;
            }
        }else{
            break;
        }
    }
    return;
}
