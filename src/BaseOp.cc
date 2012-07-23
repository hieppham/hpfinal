#include "HGAGenome.h"

extern Customer* gDepot;
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

        arrC = hgenome.arrC;
        m_route = hgenome.m_route;
        m_data = hgenome.m_data;

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
void HGAGenome::pushbackRoute(Route& mRoute, RinfoPtr& mRinfo, Customer*& mCus){
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

        mRoute.back()->timeWait = max(0, mCus->e - mRoute.back()->timeStartService);
        mRoute.back()->timeDeparture = newTimeStart + mCus->d;
    }
    if (mRoute.back()->timeStartService > mCus->l){
        mRinfo->timeVio += (mRoute.back()->timeStartService - mCus->l);
    }
}
/**
 * Insert a customer at a defined position in a route (pTrace = VERYLAST means end of route)
 */
void HGAGenome::insertIntoRoute(Route& mRoute, RinfoPtr& mRinfo, Customer*& mCus, unsigned int& pTrace){
    if (pTrace == VERYLAST){
        HGAGenome::pushbackRoute(mRoute, mRinfo, mCus);
    }else{
        Route::iterator pos = mRoute.begin();
        for (unsigned int i = 0; i < pTrace; i++){
            pos++;
        }
        mRoute.insert(pos, VertexPtr(new Vertex(mCus)));
        pos--;
        // TODO: update information for all customers
        HGAGenome::updateInfo(mRoute, mRinfo);
    }
}
/**
 * trace best position and insert customer into given route, using Potvin - Rosseau heuristic
 */
void HGAGenome::PRinsert(Route& mRoute, RinfoPtr& mRinfo, Customer*& mixer){
    double newEarliestTime = 0;
    double newLastestTime = 0;
    double newProfit = 0;
    double minProfit = POSINF;

    unsigned int pTrace = 0;
    unsigned int currPos = 0;
    Route::iterator nextCus, prevCus, uIter, beforeIter;

    if (mRoute.empty()){
        HGAGenome::pushbackRoute(mRoute, mRinfo, mixer);
    }else{
        // for (i-1, i) \in r
        if (mRoute.empty()){
            newEarliestTime = max(mixer->e, gDistance[0][mixer->id]);
            newLastestTime = mixer->l;

            newProfit = gDistance[0][mixer->id];
            if ((newEarliestTime < newLastestTime) && (newProfit < minProfit)){
                pTrace = 0;
                minProfit = newProfit;
            }
        }else{
            for (nextCus = mRoute.begin(), currPos = 0; nextCus != mRoute.end(); ++nextCus){
                // prevCus = nextCus - 1;
                // isFeasible(i, j)
                if (currPos == 0){
                    // i - 1 means the depot
                    newEarliestTime = max(mixer->e, mRinfo->timeLeaveDepot + gDistance[0][mixer->id]);
                    newLastestTime = min(mixer->l, gDepot->l);
                    newProfit = gDistance[0][mixer->id];
                    newProfit += gDistance[mixer->id][(*nextCus)->cus->id];
                    newProfit -= gDistance[(*nextCus)->cus->id][0];

                    if ((newEarliestTime < newLastestTime) && (newProfit < minProfit)){
                        pTrace = 0;
                        minProfit = newProfit;
                    }
                }else{
                    prevCus = nextCus;
                    prevCus--;
                    newEarliestTime = max(mixer->e, (*prevCus)->timeDeparture + gDistance[(*prevCus)->cus->id][mixer->id]);
                    newLastestTime = min(mixer->l, (*nextCus)->cus->l - mixer->d - gDistance[(*nextCus)->cus->id][mixer->id]);

                    newProfit = gDistance[(*prevCus)->cus->id][mixer->id] + gDistance[mixer->id][(*nextCus)->cus->id] - gDistance[(*nextCus)->cus->id][(*prevCus)->cus->id];

                    if ((newEarliestTime < newLastestTime) && (newProfit < minProfit)){
                        pTrace = currPos;
                        minProfit = newProfit;
                    }
                }
                currPos++;
            }
            // insert after the last item
            prevCus = mRoute.end();
            prevCus--;
            newEarliestTime = max(mixer->e, (*prevCus)->timeDeparture + gDistance[(*prevCus)->cus->id][mixer->id]);
            newLastestTime = min(mixer->l, gDepot->l);

            newProfit = gDistance[0][mixer->id] + gDistance[mixer->id][(*prevCus)->cus->id] - gDistance[(*prevCus)->cus->id][0];
            if ((newEarliestTime < newLastestTime) && (newProfit < minProfit)){
                pTrace = VERYLAST;
                minProfit = newProfit;
            }
        }
        // if not found the best move
        if (minProfit == POSINF){
            pTrace = VERYLAST;
        }
        // Insert (i*, j*) and Update(r*)
        HGAGenome::insertIntoRoute(mRoute, mRinfo, mixer, pTrace);
    }
    cout << "PRinsert\n";
}

void HGAGenome::updateInfo(Route& mRoute, RinfoPtr& mRinfo){
    mRinfo->resetAll();
    if (mRoute.empty()){
        return;
    }
    for (Route::iterator uIter = mRoute.begin(); uIter != mRoute.end(); ++uIter){
        if (uIter == mRoute.begin()){
            (*uIter)->timeArrive = gDistance[0][(*uIter)->cus->id];
            (*uIter)->timeStartService = max((*uIter)->cus->e, gDistance[0][(*uIter)->cus->id]);
            (*uIter)->timeWait = (*uIter)->timeStartService - (*uIter)->timeArrive;
            (*uIter)->timeDeparture = (*uIter)->timeStartService + (*uIter)->cus->d;
            mRinfo->cost += gDistance[0][(*uIter)->cus->id];
            mRinfo->load += (*uIter)->cus->q;
        }else{
            Route::iterator beforeIter = uIter;
            beforeIter--;
            (*uIter)->timeArrive = (*beforeIter)->timeDeparture + gDistance[(*uIter)->cus->id][(*beforeIter)->cus->id];

            mRinfo->cost += gDistance[(*uIter)->cus->id][(*beforeIter)->cus->id];
            (*uIter)->timeStartService = max((*uIter)->timeArrive, (*uIter)->cus->e);
            (*uIter)->timeWait = (*uIter)->timeStartService - (*uIter)->timeArrive;
            if ((*uIter)->timeStartService > (*uIter)->cus->l){
                mRinfo->timeVio += ((*uIter)->timeStartService - (*uIter)->cus->l);
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
    Route::iterator nextCus, prevCus, uIter, beforeIter;

    while(1){
        maxProfit = NEGINF;
        for (VCus::iterator kIter = arrCus.begin(), endArr = arrCus.end(); kIter != endArr; ++kIter){
            // evaluate all possible move
            for (nextCus = mRoute.begin(), currPos = 0; nextCus != mRoute.end(); ++nextCus){
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
                     while((uIter != mRoute.end()) && legal){
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
                     while ((uIter != mRoute.end()) && (legal)){
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
            prevCus = mRoute.end();
            prevCus--;
            newEarliestTime = max((*kIter)->e, (*prevCus)->timeDeparture + gDistance[(*prevCus)->cus->id][(*kIter)->id]);
            newLastestTime = min((*kIter)->l, gDepot->l);
            if (newEarliestTime < newLastestTime){
                double c11 = gDistance[(*prevCus)->cus->id][(*kIter)->id] + gDistance[(*kIter)->id][0] - HPGV::I1Mu * gDistance[(*prevCus)->cus->id][0];
                double c12 = newEarliestTime + (*kIter)->d + gDistance[(*kIter)->id][0] - (*prevCus)->timeStartService - (*prevCus)->cus->d - gDistance[(*prevCus)->cus->id][0];
                newProfit = HPGV::I1Lambda * gDistance[0][(*kIter)->id] - HPGV::I1Alpha * c11 - (1 - HPGV::I1Alpha) * c12;
                if (newProfit > maxProfit){
                    maxProfit = newProfit;
                    pTrace = VERYLAST;
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
/**
 * Potvin - Rosseau heuristic
 */
void HGAGenome::PRheuristic(vector<Route>& m_route, RouteData& m_data, VCus& arrCus, unsigned int& iDay, bool checkService = false){
    unsigned int iVeh = 0;
    unsigned int vod = 0;
    double newEarliestTime = 0;
    double newLastestTime = 0;

    unsigned int currPos = 0;
    unsigned int pTrace = 0;
    unsigned int rVod = 0;
    double newProfit = 0;
    double minProfit = POSINF;
    VCus::iterator minVPos, kIter;
    Route::iterator prevCus, nextCus, uIter, beforeIter;
    while(!arrCus.empty()){
        minProfit = POSINF;
        for (kIter = arrCus.begin(); kIter != arrCus.end(); ++kIter){
            for (iVeh = 0; iVeh < HPGV::mVeh; iVeh++){
                // for r \in R
                vod = iDay * HPGV::mVeh + iVeh;
                // for (i-1, i) \in r
                if (m_route[vod].empty()){
                    newEarliestTime = max((*kIter)->e, gDistance[0][(*kIter)->id]);
                    newLastestTime = (*kIter)->l;

                    newProfit = gDistance[0][(*kIter)->id];
                    if ((newEarliestTime < newLastestTime) && (newProfit < minProfit)){
                        pTrace = 0;
                        rVod = vod;
                        minVPos = kIter;
                        minProfit = newProfit;
                    }
                }else{
                    for (nextCus = m_route[vod].begin(), currPos = 0; nextCus != m_route[vod].end(); ++nextCus){
                        // prevCus = nextCus - 1;
                        // isFeasible(i, j)
                        if (currPos == 0){
                            // i - 1 means the depot
                            newEarliestTime = max((*kIter)->e, m_data[vod]->timeLeaveDepot + gDistance[0][(*kIter)->id]);
                            newLastestTime = min((*kIter)->l, gDepot->l);
                            newProfit = gDistance[0][(*kIter)->id];
                            newProfit += gDistance[(*kIter)->id][(*nextCus)->cus->id];
                            newProfit -= gDistance[(*nextCus)->cus->id][0];

                            if ((newEarliestTime < newLastestTime) && (newProfit < minProfit)){
                                rVod = vod;
                                pTrace = 0;
                                minVPos = kIter;
                                minProfit = newProfit;
                            }
                        }else{
                            prevCus = nextCus;
                            prevCus--;
                            newEarliestTime = max((*kIter)->e, (*prevCus)->timeDeparture + gDistance[(*prevCus)->cus->id][(*kIter)->id]);
                            newLastestTime = min((*kIter)->l, (*nextCus)->cus->l - (*kIter)->d - gDistance[(*nextCus)->cus->id][(*kIter)->id]);

                            newProfit = gDistance[(*prevCus)->cus->id][(*kIter)->id] + gDistance[(*kIter)->id][(*nextCus)->cus->id] - gDistance[(*nextCus)->cus->id][(*prevCus)->cus->id];

                            if ((newEarliestTime < newLastestTime) && (newProfit < minProfit)){
                                rVod = vod;
                                pTrace = currPos;
                                minVPos = kIter;
                                minProfit = newProfit;
                            }
                        }
                        currPos++;
                    }
                    // insert after the last item
                    prevCus = m_route[vod].end();
                    prevCus--;
                    newEarliestTime = max((*kIter)->e, (*prevCus)->timeDeparture + gDistance[(*prevCus)->cus->id][(*kIter)->id]);
                    newLastestTime = min((*kIter)->l, gDepot->l);

                    newProfit = gDistance[0][(*kIter)->id] + gDistance[(*kIter)->id][(*prevCus)->cus->id] - gDistance[(*prevCus)->cus->id][0];
                    if ((newEarliestTime < newLastestTime) && (newProfit < minProfit)){
                        rVod = vod;
                        pTrace = VERYLAST;
                        minVPos = kIter;
                        minProfit = newProfit;
                    }
                }
            }
        }// end for(kIter....)
        // if not found the best move
        if (minProfit == POSINF){
            minVPos = arrCus.begin();
            rVod = iDay * HPGV::mVeh + GARandomInt(0, HPGV::mVeh - 1);
            pTrace = VERYLAST;
        }
        // Insert (i*, j*) and Update(r*)
        HGAGenome::insertIntoRoute(m_route[rVod], m_data[rVod], *minVPos, pTrace);
        if (checkService){
            (*minVPos)->checkServiced();
        }
        // Now N = N\j*
        arrCus.erase(minVPos);
    }// end while
}
/**
 * Print out solution (for testing)
 */
void HGAGenome::printSolution(HGAGenome& hg, char* fileout){
    unsigned int iVeh = 0;
    unsigned int iDay = 0;
    unsigned int vod = 0;
    fstream ofs;
    ofs.open(fileout, ios::out|ios::app);
    if (ofs.is_open()) {
        while (ofs.good()) {
            if (!hg.isFeasible){
                ofs << "* ";
            }
            ofs << hg.durationCost << endl;
                for (iDay = 0; iDay < HPGV::tDay; iDay++){
                    for (iVeh = 0; iVeh < HPGV::mVeh; iVeh++){
                        vod = iDay*HPGV::mVeh + iVeh;
                        Route::iterator uIter, endIter;
                        for (uIter = hg.m_route[vod].begin(), endIter = hg.m_route[vod].end(); uIter != endIter; ++uIter){
                            if (uIter == hg.m_route[vod].begin()){
                                ofs << (iDay+1) << "  " << (iVeh + 1) << "\t";
                                ofs << hg.m_data[vod]->cost << "\t" << hg.m_data[vod]->load << "\t";
                                ofs << "0(" << hg.m_data[vod]->timeLeaveDepot << ")  ";
                                ofs << (*uIter)->cus->id << "[" << (*uIter)->cus->pattern << "]" << "(" << (*uIter)->timeStartService << ")  ";
                            }else{
                                ofs << (*uIter)->cus->id << "[" << (*uIter)->cus->pattern << "]" << "(" << (*uIter)->timeStartService << ")  ";
                            }
                            if (((*uIter)->cus->id < 1) || ((*uIter)->cus->id > HPGV::nCus)){
                                ofs << "####################################";
                            }
                        }
                        if (hg.m_route[vod].size() > 0){
                            uIter = --hg.m_route[vod].end();
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
}

void HGAGenome::removeFromRoute(Route& mRoute, RinfoPtr& mRinfo, int idToErase){
    if (mRoute.empty()){
        return;
    }
    for (Route::iterator uIter = mRoute.begin(), endIter = mRoute.end(); uIter != endIter; ++uIter){
        if ((*uIter)->cus->id == idToErase){
            mRoute.erase(uIter);
            break;
        }
    }
    // update information
    HGAGenome::updateInfo(mRoute, mRinfo);
}

void HGAGenome::tourConstruct(void){
    this->m_pattern.resize(HPGV::nCus);
    for (unsigned int i = 0; i < HPGV::nCus; i++){
        this->m_pattern[i] = this->arrC[i].pattern;
    }

    this->m_tour.clear();
    for (unsigned int vod = 0; vod < (HPGV::numRoute); vod++){
        for (Route::iterator rIter = this->m_route[vod].begin(), endIter = this->m_route[vod].end(); rIter != endIter; ++rIter){
            this->m_tour.push_back(CidPtr(new CustomerInDay((*rIter)->cus->id, vod)));
        }
    }
}

double HGAGenome::calcObjectValue(HGAGenome& hg){
    double score = 0;

    // hg.updateTotalVio();
    score = hg.durationCost;
    // calculate score after updating parameters
    if (HPGV::hPenalty != 0){
        double sumSq = (HPGV::avgQ * HPGV::avgQ) + (HPGV::avgD * HPGV::avgD) + (HPGV::avgW * HPGV::avgW);
        double aPen = HPGV::hPenalty * HPGV::avgQ / sumSq;
        double bPen = HPGV::hPenalty * HPGV::avgD / sumSq;
        double cPen = HPGV::hPenalty * HPGV::avgW / sumSq;
        score += aPen * hg.totalCapacityVio + bPen * hg.totalDurationVio + cPen * hg.totalTimeVio;
    }
    return score;
}

void HGAGenome::testRoute(Route& mRoute){
    bool flag = true;
    for (Route::iterator r = mRoute.begin(), e = mRoute.end(); r != e; ++r){
        if (((*r)->cus->id < 1) || ((*r)->cus->id > HPGV::nCus)){
            cout << "Error!!!\n";
            flag = false;
        }
    }
    if (flag){
        cout << "testRoute OK\n";
    }
}
