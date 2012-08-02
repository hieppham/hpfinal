#include "HGAGenome.h"

void HGAGenome::delayDeparture(Route& mRoute, RinfoPtr& mRinfo){
    Route::iterator prevCus, uIter;
    Route::reverse_iterator ruIter, rprevIter;

    mRinfo->resetAll();
    if (mRoute.empty()){
        return;
    }
    // Set D_0 = e_0
    mRinfo->timeLeaveDepot = HPGV::gDepot->e;

    // compute Ai, Bi, Wi, Di for each vertex in the route

    for (uIter = mRoute.begin(); uIter != mRoute.end(); ++uIter){
        if (uIter == mRoute.begin()){
            (*uIter)->timeArrive = mRinfo->timeLeaveDepot + HPGV::gDistance[0][(*uIter)->cus->id];
        }else{
            prevCus = uIter;
            prevCus--;
            (*uIter)->timeArrive = (*prevCus)->timeDeparture + HPGV::gDistance[(*prevCus)->cus->id][(*uIter)->cus->id];
        }
        (*uIter)->timeStartService = max((*uIter)->cus->e, (*uIter)->timeArrive);
        (*uIter)->timeWait = (*uIter)->timeStartService - (*uIter)->timeArrive;
        (*uIter)->timeDeparture = (*uIter)->timeStartService + (double)(*uIter)->cus->d;
    }

    // Compute FTS0 using backward recursion
    double sumWait = 0;
    for (ruIter = mRoute.rbegin(); ruIter != mRoute.rend(); ++ruIter){
        if (ruIter == mRoute.rbegin()){
            (*ruIter)->FTS = max(0, (*ruIter)->cus->l - (*ruIter)->timeStartService);
        }else{
            rprevIter = ruIter;
            rprevIter--;
            double  tmpF = max(0, (double)(*ruIter)->cus->l - (*ruIter)->timeStartService);
            (*ruIter)->FTS = min((*rprevIter)->FTS + (*rprevIter)->timeWait, tmpF);
        }
        sumWait += (*ruIter)->timeWait;
    }
    uIter = mRoute.begin();
    double tmpF01 = max(0, HPGV::gDepot->l - mRinfo->timeLeaveDepot);
    double tmpF02 = (*uIter)->FTS + (*uIter)->timeWait;
    mRinfo->FTS0 = min(tmpF01, tmpF02);

    // set D_0 = e_0 + min{F_0, \sum W_p}
    mRinfo->timeLeaveDepot = min(mRinfo->FTS0, sumWait);
    mRinfo->timeLeaveDepot += HPGV::gDepot->e;

    // update A_i, W_i, B_i, D_i for each vertex in the route
    // also evaluate violations in respect to capacity, duration and time windows constraints
    for (uIter = mRoute.begin(); uIter != mRoute.end(); ++uIter){
        if (uIter == mRoute.begin()){
            (*uIter)->timeArrive = mRinfo->timeLeaveDepot + HPGV::gDistance[0][(*uIter)->cus->id];
            mRinfo->cost += HPGV::gDistance[0][(*uIter)->cus->id];
        }else{
            prevCus = uIter;
            prevCus--;
            (*uIter)->timeArrive = (*prevCus)->timeDeparture + HPGV::gDistance[(*prevCus)->cus->id][(*uIter)->cus->id];
            mRinfo->cost += HPGV::gDistance[(*prevCus)->cus->id][(*uIter)->cus->id];
        }
        (*uIter)->timeStartService = max((double)(*uIter)->cus->e, (*uIter)->timeArrive);
        (*uIter)->timeWait = (*uIter)->timeStartService - (*uIter)->timeArrive;
        (*uIter)->timeDeparture = (*uIter)->timeStartService + (double)(*uIter)->cus->d;

        mRinfo->load += (double)(*uIter)->cus->q;
        if ((*uIter)->timeStartService > (double)(*uIter)->cus->l){
            mRinfo->timeVio += ((*uIter)->timeStartService - (double)(*uIter)->cus->l);
        }
    }
    uIter = --mRoute.end();
    double timeBack = HPGV::gDistance[0][(*uIter)->cus->id];
    mRinfo->cost += timeBack;
    timeBack += (*uIter)->timeDeparture;
    if (timeBack > HPGV::gDepot->l){
        mRinfo->timeVio += (timeBack - HPGV::gDepot->l);
    }
}
/**
 * Improve quality of solution using some local search procedures
 */
void HGAGenome::improveRoute(HGAGenome& hg){
    unsigned int iDay, iVehicle, vod;
    unsigned int i, seed1, seed2, rv1, rv2, vod1, vod2;

    // penal parameters
    double aPen, bPen, cPen;
    double sumSq = (HPGV::avgQ * HPGV::avgQ) + (HPGV::avgD * HPGV::avgD) + (HPGV::avgW * HPGV::avgW);
    if (sumSq == 0){
        aPen = 100;
        bPen = 100;
        cPen = 100;
    }else{
        aPen = HPGV::hPenalty * HPGV::avgQ / sumSq;
        bPen = HPGV::hPenalty * HPGV::avgD / sumSq;
        cPen = HPGV::hPenalty * HPGV::avgW / sumSq;
    }

    for (iDay = 0; iDay < HPGV::tDay; iDay++){
        // consider all possible pairs of routes in random order
        seed1 = GARandomInt(0, HPGV::mVeh - 2);
        seed2 = GARandomInt(seed1+1, HPGV::mVeh - 1);
        for (i = 0; i < HPGV::mVeh; i++){
            rv1 = (seed1 + i) % HPGV::mVeh;
            rv2 = (seed2 + i) % HPGV::mVeh;

            vod1 = iDay * HPGV::mVeh + rv1;
            vod2 = iDay * HPGV::mVeh + rv2;

            // apply inter-route operators : 2-opt*, lambda-interchange and CROSS-exchange
            if (hg.interRouteOpt(aPen, bPen, cPen, vod1, vod2)){
                break;
            }
        }

        // no further improvement can be found
        // locally improving each route of the current day in turn
        for (iVehicle = 0; iVehicle < HPGV::mVeh; iVehicle++){
            vod = iDay * HPGV::mVeh + iVehicle;
            for (int ilocal = 0; ilocal < 10; ilocal++){
                if (hg.intra2Opt(aPen, bPen, cPen, vod)){
                    break;
                }

                if (hg.intraOrOpt(aPen, bPen, cPen, vod)){
                    break;
                }
            }
        }
    }
}

void HGAGenome::apply2OptForAllRoutes(double aPen, double bPen, double cPen, HGAGenome& hg){
    unsigned int iDay, iVehicle, vod;

    for (iDay = 0; iDay < HPGV::tDay; iDay++){
        for (iVehicle = 0; iVehicle < HPGV::mVeh; iVehicle++){
            vod = iDay * HPGV::mVeh + iVehicle;
            hg.intra2Opt(aPen, bPen, cPen, vod);
        }
    }
    hg.tourConstruct();
}

void HGAGenome::apply2OptUntilFirstImprovement(double aPen, double bPen, double cPen, HGAGenome& hg){
    // update 310712: priotize stochastic strategy
    unsigned int tmpvod = GARandomInt(0, HPGV::numRoute - 1);

    for (unsigned i = 0; i < 20; i++){
        unsigned int vod = (tmpvod + i) % HPGV::numRoute;
        if (hg.intra2Opt(aPen, bPen, cPen, vod)){
            break;
        }
    }
    hg.tourConstruct();
}

void HGAGenome::apply2OptStarForAllRoutes(double aPen, double bPen, double cPen, HGAGenome& hg){
    // update 310712: priotize stochastic strategy
    unsigned int iDay, iVF, iVS, vod1, vod2;

    for (iDay = 0; iDay < HPGV::tDay; iDay++){
        for (iVF = 0; iVF < HPGV::mVeh - 1; iVF++){
            vod1 = iDay * HPGV::mVeh + iVF;
            for (iVS = iVF + 1; iVS < HPGV::mVeh; iVS++){
                vod2 = iDay * HPGV::mVeh + iVS;
                hg.inter2OptStar(aPen, bPen, cPen, vod1, vod2);
            }
        }
    }
    hg.tourConstruct();
}

void HGAGenome::apply2OptStarUntilFirstImprovement(double aPen, double bPen, double cPen, HGAGenome& hg){
    unsigned int iDay, iVF, iVS, vod1, vod2;
    bool improved = false;

    if (HPGV::mVeh < 2){
        return;
    }

    for (iDay = 0; iDay < HPGV::tDay; iDay++){
        iVF = GARandomInt(0, HPGV::mVeh - 2);
        iVS = GARandomInt(iVF + 1, HPGV::mVeh - 1);
        vod1 = iDay * HPGV::mVeh + iVF;
        vod2 = iDay * HPGV::mVeh + iVS;
        if (hg.inter2OptStar(aPen, bPen, cPen, vod1, vod2)){
            improved = true;
            break;
        }
    }
    if (improved){
        hg.tourConstruct();
    }
}

/**
 * intra Or-Opt operation
 */
bool HGAGenome::intraOrOpt(double aPen, double bPen, double cPen, unsigned int vod){
    bool improved = false;
    Route firstSeg, lastSeg, midSeg, tempSeg, newRoute, curRoute, bestRoute;

    unsigned int count = 0;
    unsigned int i;
    unsigned int rSize = this->m_route[vod].size();

    if (rSize <= 3){
        return false;
    }

    HGAGenome::delayDeparture(this->m_route[vod], this->m_data[vod]);
    double minCost = cost(aPen, bPen, cPen, HPGV::maxLoad, HPGV::maxDuration)(this->m_route[vod], this->m_data[vod]);

    bestRoute = this->m_route[vod];

    do {
        curRoute = this->m_route[vod];
        firstSeg.clear();
        midSeg.clear();
        lastSeg.clear();
        for (i = 0; i < count; i++){
            firstSeg.push_back(curRoute.front());
            curRoute.pop_front();
        }
        // middle segment contains 2 consecutive vertices
        midSeg.push_back(curRoute.front());
        curRoute.pop_front();
        midSeg.push_back(curRoute.front());
        curRoute.pop_front();

        lastSeg = curRoute;

        while (lastSeg.size() > 0){
            newRoute.clear();

            newRoute = firstSeg;

            newRoute.push_back(lastSeg.front());    // 0 -> ... -> i-1 -> i+2 -> j
            tempSeg = midSeg;
            // 0 -> ... -> i-1 -> i+2 -> j -> i -> i+1
            newRoute.splice(newRoute.end(), tempSeg, tempSeg.begin(), tempSeg.end());

            firstSeg.push_back(lastSeg.front());
            lastSeg.pop_front();

            tempSeg = lastSeg;
            // 0 -> ... -> i-1 -> i+2 -> j -> i -> i+1 -> j+1 -> ...
            newRoute.splice(newRoute.end(), tempSeg, tempSeg.begin(), tempSeg.end());

            HGAGenome::delayDeparture(newRoute, this->m_data[vod]);

            double newCost = cost(aPen, bPen, cPen, HPGV::maxLoad, HPGV::maxDuration)(newRoute, this->m_data[vod]);
            if (newCost < minCost){
                minCost = newCost;
                bestRoute = newRoute;
                improved = true;
            }
        }
        count++;
    } while (count < rSize-2);

    this->m_route[vod] = bestRoute;
    HGAGenome::delayDeparture(this->m_route[vod], this->m_data[vod]);

    return improved;
}

bool HGAGenome::stoIntraOrOpt(double aPen, double bPen, double cPen, unsigned int vod){
    bool improved = false;
    Route firstSeg, lastSeg, midSeg, tempSeg, newRoute, curRoute, bestRoute;

    unsigned int count = 0;
    unsigned int i;
    unsigned int rSize = this->m_route[vod].size();

    if (rSize <= 3){
        return false;
    }

    HGAGenome::delayDeparture(this->m_route[vod], this->m_data[vod]);
    double minCost = cost(aPen, bPen, cPen, HPGV::maxLoad, HPGV::maxDuration)(this->m_route[vod], this->m_data[vod]);

    bestRoute = this->m_route[vod];

    do {
        curRoute = this->m_route[vod];
        firstSeg.clear();
        midSeg.clear();
        lastSeg.clear();
        for (i = 0; i < count; i++){
            firstSeg.push_back(curRoute.front());
            curRoute.pop_front();
        }
        // middle segment contains 2 consecutive vertices
        midSeg.push_back(curRoute.front());
        curRoute.pop_front();
        midSeg.push_back(curRoute.front());
        curRoute.pop_front();

        lastSeg = curRoute;

        while ((lastSeg.size() > 0) && (!improved)){
            newRoute.clear();

            newRoute = firstSeg;

            newRoute.push_back(lastSeg.front());    // 0 -> ... -> i-1 -> i+2 -> j
            tempSeg = midSeg;
            // 0 -> ... -> i-1 -> i+2 -> j -> i -> i+1
            newRoute.splice(newRoute.end(), tempSeg, tempSeg.begin(), tempSeg.end());

            firstSeg.push_back(lastSeg.front());
            lastSeg.pop_front();

            tempSeg = lastSeg;
            // 0 -> ... -> i-1 -> i+2 -> j -> i -> i+1 -> j+1 -> ...
            newRoute.splice(newRoute.end(), tempSeg, tempSeg.begin(), tempSeg.end());

            HGAGenome::delayDeparture(newRoute, this->m_data[vod]);

            double newCost = cost(aPen, bPen, cPen, HPGV::maxLoad, HPGV::maxDuration)(newRoute, this->m_data[vod]);
            if (newCost < minCost){
                minCost = newCost;
                bestRoute = newRoute;
                improved = true;
            }
        }
        count++;
    } while ((count < rSize-2) && (!improved));

    this->m_route[vod] = bestRoute;
    HGAGenome::delayDeparture(this->m_route[vod], this->m_data[vod]);

    return improved;
}

bool HGAGenome::intra2Opt(double aPen, double bPen, double cPen, unsigned int vod){
    bool improved = false;

    Route firstSeg, lastSeg, midSeg, tempSeg, newRoute, curRoute, bestRoute;
    // Route::iterator uIter;
    // RinfoPtr rInfo = (RinfoPtr)(new RouteInfo());

    unsigned int count = 0;
    unsigned int i, lastSize;
    unsigned int rSize = this->m_route[vod].size();

    if (rSize <= 2){
        // this route has 0, 2 or 3 edge(s)
        return false;
    }

    HGAGenome::delayDeparture(this->m_route[vod], this->m_data[vod]);
    double minCost = cost(aPen, bPen, cPen, HPGV::maxLoad, HPGV::maxDuration)(this->m_route[vod], this->m_data[vod]);

    bestRoute = this->m_route[vod];

    do {
        if (count == 0){
            lastSize = 1;
        }else{
            lastSize = 0;
        }
        curRoute = this->m_route[vod];
        firstSeg.clear();
        midSeg.clear();
        lastSeg.clear();
        for (i = 0; i < count; i++){
            firstSeg.push_back(curRoute.front());
            curRoute.pop_front();
        }
        midSeg.push_back(curRoute.front());
        curRoute.pop_front();
        lastSeg = curRoute;

        while (lastSeg.size() > lastSize){
            newRoute.clear();

            newRoute = firstSeg;

            newRoute.push_back(lastSeg.front());    // 0 -> i
            tempSeg = midSeg;
            // 0 -> i -> i-1 -> ... ->2
            newRoute.splice(newRoute.end(), tempSeg, tempSeg.begin(), tempSeg.end());

            midSeg.push_front(lastSeg.front());
            lastSeg.pop_front();

            tempSeg = lastSeg;
            // 0 -> i -> i-1 -> ... 2 -> 1 -> i+1 -> ...
            newRoute.splice(newRoute.end(), tempSeg, tempSeg.begin(), tempSeg.end());
            HGAGenome::delayDeparture(newRoute, this->m_data[vod]);
            double newCost = cost(aPen, bPen, cPen, HPGV::maxLoad, HPGV::maxDuration)(newRoute, this->m_data[vod]);
            if (newCost < minCost){
                minCost = newCost;
                bestRoute = newRoute;
                improved = true;
            }
        }
        count++;
    } while (count < rSize);

    this->m_route[vod] = bestRoute;
    HGAGenome::delayDeparture(this->m_route[vod], this->m_data[vod]);
    return improved;
}
/**
 * another version of 2-Opt operator, that prefers stochastic strategy
 */
bool HGAGenome::stoIntra2Opt(double aPen, double bPen, double cPen, unsigned int vod){
    bool improved = false;

    Route firstSeg, lastSeg, midSeg, tempSeg, newRoute, curRoute, bestRoute;
    // Route::iterator uIter;
    // RinfoPtr rInfo = (RinfoPtr)(new RouteInfo());

    unsigned int count = 0;
    unsigned int i, lastSize;
    unsigned int rSize = this->m_route[vod].size();

    if (rSize <= 2){
        // this route has 0, 2 or 3 edge(s)
        return false;
    }

    HGAGenome::delayDeparture(this->m_route[vod], this->m_data[vod]);
    double minCost = cost(aPen, bPen, cPen, HPGV::maxLoad, HPGV::maxDuration)(this->m_route[vod], this->m_data[vod]);

    bestRoute = this->m_route[vod];

    do {
        if (count == 0){
            lastSize = 1;
        }else{
            lastSize = 0;
        }
        curRoute = this->m_route[vod];
        firstSeg.clear();
        midSeg.clear();
        lastSeg.clear();
        for (i = 0; i < count; i++){
            firstSeg.push_back(curRoute.front());
            curRoute.pop_front();
        }
        midSeg.push_back(curRoute.front());
        curRoute.pop_front();
        lastSeg = curRoute;

        while ((lastSeg.size() > lastSize) && (!improved)){
            newRoute.clear();

            newRoute = firstSeg;

            newRoute.push_back(lastSeg.front());    // 0 -> i
            tempSeg = midSeg;
            // 0 -> i -> i-1 -> ... ->2
            newRoute.splice(newRoute.end(), tempSeg, tempSeg.begin(), tempSeg.end());

            midSeg.push_front(lastSeg.front());
            lastSeg.pop_front();

            tempSeg = lastSeg;
            // 0 -> i -> i-1 -> ... 2 -> 1 -> i+1 -> ...
            newRoute.splice(newRoute.end(), tempSeg, tempSeg.begin(), tempSeg.end());
            HGAGenome::delayDeparture(newRoute, this->m_data[vod]);
            double newCost = cost(aPen, bPen, cPen, HPGV::maxLoad, HPGV::maxDuration)(newRoute, this->m_data[vod]);
            if (newCost < minCost){
                minCost = newCost;
                bestRoute = newRoute;
                improved = true;
            }
        }
        count++;
    } while ((count < rSize) && (!improved));

    this->m_route[vod] = bestRoute;
    HGAGenome::delayDeparture(this->m_route[vod], this->m_data[vod]);
    return improved;
}

bool HGAGenome::interRouteOpt(double aPen, double bPen, double cPen, unsigned int vod1, unsigned int vod2){
    bool improved = false;

    int rOrder = GARandomInt(1, 2);
    if (rOrder == 1){
        improved = this->inter2OptStar(aPen, bPen, cPen, vod1, vod2);
    }else if (rOrder == 2){
        improved = this->interCrossExchange(aPen, bPen, cPen, vod1, vod2);
    }
    return improved;
}

bool HGAGenome::inter2OptStar(double aPen, double bPen, double cPen, unsigned int vod1, unsigned int vod2){
    bool improved = false;

    Route firstHead, secondHead, firstTail, secondTail, tempSeg;
    Route curFirst, curSecond, newFirst, newSecond, bestFirst, bestSecond;

    unsigned int i, j, k;
    unsigned int firstSize = this->m_route[vod1].size();
    unsigned int secondSize = this->m_route[vod2].size();

    double minCost;

    bestFirst = this->m_route[vod1];
    bestSecond = this->m_route[vod2];

    HGAGenome::delayDeparture(this->m_route[vod1], this->m_data[vod1]);
    HGAGenome::delayDeparture(this->m_route[vod2], this->m_data[vod2]);
    minCost = cost(aPen, bPen, cPen, HPGV::maxLoad, HPGV::maxDuration)(this->m_route[vod1], this->m_data[vod1]);
    minCost += cost(aPen, bPen, cPen, HPGV::maxLoad, HPGV::maxDuration)(this->m_route[vod2], this->m_data[vod2]);

    // update 310712: prioritize stochastic strategy
    int tryCounter = 0;
    while ((tryCounter < 10) && (!improved)){
        // splice first route
        curFirst = this->m_route[vod1];
        i = (unsigned int) GARandomInt(0, firstSize - 2);

        firstHead.clear();
        for (k = 0; k < i; k++){
            firstHead.push_back(curFirst.front());
            curFirst.pop_front();
        }
        firstTail = curFirst;

        // splice second route
        curSecond = this->m_route[vod2];
        j = (unsigned int) GARandomInt(0, secondSize - 2);
        secondHead.clear();
        for (k = 0; k < j; k++){
            secondHead.push_back(curSecond.front());
            curSecond.pop_front();
        }
        secondTail = curSecond;

        // create 2 new routes
        newFirst.clear();
        newSecond.clear();

        newFirst = firstHead;
        tempSeg = secondTail;
        newFirst.splice(newFirst.end(), tempSeg, tempSeg.begin(), tempSeg.end());

        newSecond = secondHead;
        tempSeg = firstTail;
        newSecond.splice(newSecond.end(), tempSeg, tempSeg.begin(), tempSeg.end());

        HGAGenome::delayDeparture(newFirst, this->m_data[vod1]);
        HGAGenome::delayDeparture(newSecond, this->m_data[vod2]);
        double newCost = cost(aPen, bPen, cPen, HPGV::maxLoad, HPGV::maxDuration)(newFirst, this->m_data[vod1]);
        newCost += cost(aPen, bPen, cPen, HPGV::maxLoad, HPGV::maxDuration)(newSecond, this->m_data[vod2]);

        if (newCost < minCost){
            minCost = newCost;
            bestFirst = newFirst;
            bestSecond = newSecond;
            improved = true;
        }
        tryCounter++;
    }

    this->m_route[vod1] = bestFirst;
    this->m_route[vod2] = bestSecond;
    HGAGenome::delayDeparture(this->m_route[vod1], this->m_data[vod1]);
    HGAGenome::delayDeparture(this->m_route[vod2], this->m_data[vod2]);

    return improved;
}

bool HGAGenome::interCrossExchange(double aPen, double bPen, double cPen, unsigned int vod1, unsigned int vod2){
    bool improved = false;

    Route firstHead, secondHead, firstMid, secondMid, firstTail, secondTail, tempSeg;
    Route curFirst, curSecond, newFirst, newSecond, bestFirst, bestSecond;

    unsigned int i, j, k, mF, mS;
    unsigned int firstSize = this->m_route[vod1].size();
    unsigned int secondSize = this->m_route[vod2].size();

    if ((firstSize  <= 2) || (secondSize <= 2)){
        return false;
    }

    double minCost;

    bestFirst = this->m_route[vod1];
    bestSecond = this->m_route[vod2];

    HGAGenome::delayDeparture(this->m_route[vod1], this->m_data[vod1]);
    HGAGenome::delayDeparture(this->m_route[vod2], this->m_data[vod2]);
    minCost = cost(aPen, bPen, cPen, HPGV::maxLoad, HPGV::maxDuration)(this->m_route[vod1], this->m_data[vod1]);
    minCost += cost(aPen, bPen, cPen, HPGV::maxLoad, HPGV::maxDuration)(this->m_route[vod2], this->m_data[vod2]);

    // update 310712: prioritize stochastic strategy
    int tryCounter = 0;
    while ((tryCounter < 10) && (!improved)){
        // splice first route
        curFirst = this->m_route[vod1];

        i = (unsigned int) GARandomInt(0, firstSize - 2);
        mF = (unsigned int) GARandomInt(i + 1, firstSize - 1);

        firstHead.clear();
        for (k = 0; k < i; k++){
            firstHead.push_back(curFirst.front());
            curFirst.pop_front();
        }
        firstTail = curFirst;
        firstMid.clear();
        for (k = i; k <= mF; k++){
            firstMid.push_back(firstTail.front());
            firstTail.pop_front();
        }

        // splice second route
        curSecond = this->m_route[vod2];

        j = (unsigned int) GARandomInt(0, secondSize - 2);
        mS = (unsigned int) GARandomInt(j + 1, secondSize - 1);

        secondHead.clear();
        for (k = 0; k < j; k++){
            secondHead.push_back(curSecond.front());
            curSecond.pop_front();
        }
        secondTail = curSecond;
        secondMid.clear();
        for (k = j; k <= mS; k++){
            secondMid.push_back(secondTail.front());
            secondTail.pop_front();
        }

        // All segments are ready, here we create and test 2 new routes
        newFirst.clear();
        newSecond.clear();

        newFirst = firstHead;
        tempSeg = secondMid;
        newFirst.splice(newFirst.end(), tempSeg, tempSeg.begin(), tempSeg.end());
        tempSeg = firstTail;
        newFirst.splice(newFirst.end(), tempSeg, tempSeg.begin(), tempSeg.end());

        newSecond = secondHead;
        tempSeg = firstMid;
        newSecond.splice(newSecond.end(), tempSeg, tempSeg.begin(), tempSeg.end());
        tempSeg = secondTail;
        newSecond.splice(newSecond.end(), tempSeg, tempSeg.begin(), tempSeg.end());

        HGAGenome::delayDeparture(newFirst, this->m_data[vod1]);
        HGAGenome::delayDeparture(newSecond, this->m_data[vod2]);
        double newCost = cost(aPen, bPen, cPen, HPGV::maxLoad, HPGV::maxDuration)(newFirst, this->m_data[vod1]);
        newCost += cost(aPen, bPen, cPen, HPGV::maxLoad, HPGV::maxDuration)(newSecond, this->m_data[vod2]);

        if (newCost < minCost){
            minCost = newCost;
            bestFirst = newFirst;
            bestSecond = newSecond;
            improved = true;
        }
        tryCounter++;
    }

    this->m_route[vod1] = bestFirst;
    this->m_route[vod2] = bestSecond;
    HGAGenome::delayDeparture(this->m_route[vod1], this->m_data[vod1]);
    HGAGenome::delayDeparture(this->m_route[vod2], this->m_data[vod2]);

    return improved;
}
