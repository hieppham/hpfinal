#include "HGAGenome.h"

extern Customer* gDepot;
extern vector<vector<double> > gDistance;

void HGAGenome::delayDeparture(Route& mRoute, RinfoPtr& mRinfo){
    Route::iterator prevCus, uIter;
    Route::reverse_iterator ruIter, rprevIter;

    mRinfo->resetAll();
    if (mRoute.empty()){
        return;
    }
    // Set D_0 = e_0
    mRinfo->timeLeaveDepot = gDepot->e;

    // compute Ai, Bi, Wi, Di for each vertex in the route

    for (uIter = mRoute.begin(); uIter != mRoute.end(); ++uIter){
        if (uIter == mRoute.begin()){
            (*uIter)->timeArrive = mRinfo->timeLeaveDepot + gDistance[0][(*uIter)->cus->id];
        }else{
            prevCus = uIter;
            prevCus--;
            (*uIter)->timeArrive = (*prevCus)->timeDeparture + gDistance[(*prevCus)->cus->id][(*uIter)->cus->id];
        }
        (*uIter)->timeStartService = max((*uIter)->cus->e, (*uIter)->timeArrive);
        (*uIter)->timeWait = (*uIter)->timeStartService - (*uIter)->timeArrive;
        (*uIter)->timeDeparture = (*uIter)->timeStartService + (*uIter)->cus->d;
    }

    // Compute FTS0 using backward recursion
    double sumWait = 0;
    for (ruIter = mRoute.rbegin(); ruIter != mRoute.rend(); ++ruIter){
        if (ruIter == mRoute.rbegin()){
            (*ruIter)->FTS = max(0, (*ruIter)->cus->l - (*ruIter)->timeStartService);
        }else{
            rprevIter = ruIter;
            rprevIter--;
            double  tmpF = max(0, (*ruIter)->cus->l - (*ruIter)->timeStartService);
            (*ruIter)->FTS = min((*rprevIter)->FTS + (*rprevIter)->timeWait, tmpF);
        }
        sumWait += (*ruIter)->timeWait;
    }
    uIter = mRoute.begin();
    double tmpF01 = max(0, gDepot->l - mRinfo->timeLeaveDepot);
    double tmpF02 = (*uIter)->FTS + (*uIter)->timeWait;
    mRinfo->FTS0 = min(tmpF01, tmpF02);

    // set D_0 = e_0 + min{F_0, \sum W_p}
    mRinfo->timeLeaveDepot = min(mRinfo->FTS0, sumWait);
    mRinfo->timeLeaveDepot += gDepot->e;

    // update A_i, W_i, B_i, D_i for each vertex in the route
    // also evaluate violations in respect to capacity, duration and time windows constraints
    for (uIter = mRoute.begin(); uIter != mRoute.end(); ++uIter){
        if (uIter == mRoute.begin()){
            (*uIter)->timeArrive = mRinfo->timeLeaveDepot + gDistance[0][(*uIter)->cus->id];
            mRinfo->cost += gDistance[0][(*uIter)->cus->id];
        }else{
            prevCus = uIter;
            prevCus--;
            (*uIter)->timeArrive = (*prevCus)->timeDeparture + gDistance[(*prevCus)->cus->id][(*uIter)->cus->id];
            mRinfo->cost += gDistance[(*prevCus)->cus->id][(*uIter)->cus->id];
        }
        (*uIter)->timeStartService = max((*uIter)->cus->e, (*uIter)->timeArrive);
        (*uIter)->timeWait = (*uIter)->timeStartService - (*uIter)->timeArrive;
        (*uIter)->timeDeparture = (*uIter)->timeStartService + (*uIter)->cus->d;

        mRinfo->load += (*uIter)->cus->q;
        if ((*uIter)->timeStartService > (*uIter)->cus->l){
            mRinfo->timeVio += ((*uIter)->timeStartService - (*uIter)->cus->l);
        }
    }
    uIter = --mRoute.end();
    double timeBack = gDistance[0][(*uIter)->cus->id];
    mRinfo->cost += timeBack;
    timeBack += (*uIter)->timeDeparture;
    if (timeBack > gDepot->l){
        mRinfo->timeVio += (timeBack - gDepot->l);
    }
}
/**
 * Improve quality of solution using some local search procedures
 */
void HGAGenome::improveRoute(HGAGenome& hg){
    // TODO: improveRoute
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

void HGAGenome::apply2OptStarForAllRoutes(double aPen, double bPen, double cPen, HGAGenome& hg){
    unsigned int iDay, iVF, iVS, vod1, vod2;

    for (iDay = 0; iDay < HPGV::tDay; iDay++){
        for (iVF = 0; iVF < HPGV::mVeh - 1; iVF++){
            vod1 = iDay * HPGV::mVeh + iVF;
            for (iVS = iVF+1; iVS < HPGV::mVeh; iVS++){
                vod2 = iDay * HPGV::mVeh + iVS;
                hg.inter2OptStar(aPen, bPen, cPen, vod1, vod2);
            }
        }
    }
    hg.tourConstruct();
}

bool HGAGenome::intraOrOpt(double aPel, double bPel, double cPel, unsigned int vod){
    bool improved = false;
    Route firstSeg, lastSeg, midSeg, tempSeg, newRoute, curRoute, bestRoute;

    unsigned int count = 0;
    unsigned int i;
    unsigned int rSize = this->m_route[vod].size();

    if (rSize <= 3){
        return false;
    }

    HGAGenome::delayDeparture(this->m_route[vod], this->m_data[vod]);
    double minCost = cost(aPel, bPel, cPel, HPGV::maxLoad, HPGV::maxDuration)(this->m_route[vod], this->m_data[vod]);

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

        while(lastSeg.size() > 0){
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

            double newCost = cost(aPel, bPel, cPel, HPGV::maxLoad, HPGV::maxDuration)(newRoute, this->m_data[vod]);
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

bool HGAGenome::intra2Opt(double aPel, double bPel, double cPel, unsigned int vod){
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
    double minCost = cost(aPel, bPel, cPel, HPGV::maxLoad, HPGV::maxDuration)(this->m_route[vod], this->m_data[vod]);

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

        while(lastSeg.size() > lastSize){
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
            double newCost = cost(aPel, bPel, cPel, HPGV::maxLoad, HPGV::maxDuration)(newRoute, this->m_data[vod]);
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

bool HGAGenome::interRouteOpt(double aPel, double bPel, double cPel, unsigned int vod1, unsigned int vod2){
    bool improved = false;

    int rOrder = GARandomInt(1, 2);
    if (rOrder == 1){
        improved = this->inter2OptStar(aPel, bPel, cPel, vod1, vod2);
    }else if (rOrder == 2){
        improved = this->interCrossExchange(aPel, bPel, cPel, vod1, vod2);
    }
    // TODO: complete your code here asap. (Add more operators).
    return improved;
}

bool HGAGenome::inter2OptStar(double aPel, double bPel, double cPel, unsigned int vod1, unsigned int vod2){
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
    minCost = cost(aPel, bPel, cPel, HPGV::maxLoad, HPGV::maxDuration)(this->m_route[vod1], this->m_data[vod1]);
    minCost += cost(aPel, bPel, cPel, HPGV::maxLoad, HPGV::maxDuration)(this->m_route[vod2], this->m_data[vod2]);

    for (i = 0; i <= firstSize; i++){
        curFirst = this->m_route[vod1];

        firstHead.clear();
        for (k = 0; k < i; k++){
            firstHead.push_back(curFirst.front());
            curFirst.pop_front();
        }
        firstTail = curFirst;

        for (j = 0; j <= secondSize; j++){
            curSecond = this->m_route[vod2];

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
            double newCost = cost(aPel, bPel, cPel, HPGV::maxLoad, HPGV::maxDuration)(newFirst, this->m_data[vod1]);
            newCost += cost(aPel, bPel, cPel, HPGV::maxLoad, HPGV::maxDuration)(newSecond, this->m_data[vod2]);

            if (newCost < minCost){
                minCost = newCost;
                bestFirst = newFirst;
                bestSecond = newSecond;
                improved = true;
            }
        }
    }

    this->m_route[vod1] = bestFirst;
    this->m_route[vod2] = bestSecond;
    HGAGenome::delayDeparture(this->m_route[vod1], this->m_data[vod1]);
    HGAGenome::delayDeparture(this->m_route[vod2], this->m_data[vod2]);

    return improved;
}

bool HGAGenome::interCrossExchange(double aPel, double bPel, double cPel, unsigned int vod1, unsigned int vod2){
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
    delayDeparture(this->m_route[vod2], this->m_data[vod2]);
    minCost = cost(aPel, bPel, cPel, HPGV::maxLoad, HPGV::maxDuration)(this->m_route[vod1], this->m_data[vod1]);
    minCost += cost(aPel, bPel, cPel, HPGV::maxLoad, HPGV::maxDuration)(this->m_route[vod2], this->m_data[vod2]);

    for (i = 0; i < firstSize - 1; i++){
        curFirst = this->m_route[vod1];

        firstHead.clear();
        for (k = 0; k < i; k++){
            firstHead.push_back(curFirst.front());
            curFirst.pop_front();
        }
        mF = 1;
        firstTail = curFirst;
        firstMid.clear();
        while ((mF < 4) && (!firstTail.empty())){
            firstMid.push_back(firstTail.front());
            firstTail.pop_front();

            for (j = 0; j < secondSize - 1; j++){
                curSecond = this->m_route[vod2];

                secondHead.clear();
                for (k = 0; k < j; k++){
                    secondHead.push_back(curSecond.front());
                    curSecond.pop_front();
                }
                mS = 1;
                secondTail = curSecond;
                secondMid.clear();

                while ((mS < 4) && (!secondTail.empty())){
                    secondMid.push_back(secondTail.front());
                    secondTail.pop_front();
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
                    double newCost = cost(aPel, bPel, cPel, HPGV::maxLoad, HPGV::maxDuration)(newFirst, this->m_data[vod1]);
                    newCost += cost(aPel, bPel, cPel, HPGV::maxLoad, HPGV::maxDuration)(newSecond, this->m_data[vod2]);

                    if (newCost < minCost){
                        minCost = newCost;
                        bestFirst = newFirst;
                        bestSecond = newSecond;
                        improved = true;
                    }

                    mS++;
                }
            }
            mF++;
        }

    }

    this->m_route[vod1] = bestFirst;
    this->m_route[vod2] = bestSecond;
    HGAGenome::delayDeparture(this->m_route[vod1], this->m_data[vod1]);
    HGAGenome::delayDeparture(this->m_route[vod2], this->m_data[vod2]);

    return improved;
}
