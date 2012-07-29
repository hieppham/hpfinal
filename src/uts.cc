#include "HGAGenome.h"


inline double utsFit(HGAGenome& hg, double& aQ, double& bD, double& cW){
    return (hg.durationCost + aQ*hg.totalCapacityVio + bD*hg.totalDurationVio + cW*hg.totalTimeVio);
}
double utsPenalty(HGAGenome& hs, HGAGenome& hg, vector<vector<int> >& pFreq){
    double result = 0;
    if (hg.durationCost  < hs.durationCost){
        return result;
    }
    cusinday::iterator iterTour;
    for (iterTour = hg.m_tour.begin(); iterTour != hg.m_tour.end(); ++iterTour){
        result += pFreq[(*iterTour)->vod][(*iterTour)->cid];
    }
    result *= HPGV::utsLambda * hg.durationCost;
    result *= sqrt(double(HPGV::nCus * HPGV::mVeh));
    return result;
}
/**
 * update tabu list
 */
void updateTabuMap(TabuMap& tabuMap, int& tabuLength){
    TabuMap::iterator pos;
    for (pos = tabuMap.begin(); pos != tabuMap.end(); ){
        (*pos).second++;
        if ((*pos).second > tabuLength){
            tabuMap.erase(pos++);
        }else{
            ++pos;
        }
    }
}

/**
 * check tabu status
 */
bool isTabu(unsigned int cid, unsigned int vod, TabuMap& tabuMap){
    unsigned int key = (unsigned int)POSINF * cid + vod;
    TabuMap::iterator pos = tabuMap.find(key);
    return (pos != tabuMap.end());
}

void HGAGenome::tourUpdate(vector<vector<int> >& pFreq){
    unsigned int iDay = 0;
    unsigned int iVehicle = 0;
    unsigned int vod;

    this->m_pattern.resize(HPGV::nCus);
    for (unsigned int i = 0; i < HPGV::nCus; i++){
        this->m_pattern[i] = this->arrC[i].pattern;
    }

    this->m_tour.clear();
    for (iDay = 0; iDay < HPGV::tDay; iDay++) {
        for (iVehicle = 0; iVehicle < HPGV::mVeh; iVehicle++) {
            vod = iDay * HPGV::mVeh + iVehicle;
            for (Route::iterator iterRoute = this->m_route[vod].begin(); iterRoute != this->m_route[vod].end(); ++iterRoute) {
                this->m_tour.push_back(CidPtr(new CustomerInDay((*iterRoute)->cus->id, vod)));
                pFreq[vod][(unsigned int)(*iterRoute)->cus->id]++;
            }
        }
    }
}

HGAGenome HGAGenome::UTS(HGAGenome& hg){
    double aQ = 1;
    double bD = 1;
    double cW = 1;

    TabuMap g_tabu;
    // number of times attribute (i, k) has been added to the solution during the search process
    vector<vector<int> > pFreq;

    HGAGenome hs(hg);
    HGAGenome hs2(hg);
    HGAGenome hs1(hg);
    bool flag = false;

    // Tabu Length \theta = 1.5lg(n)
    int tabuLength = (int) (1.5 * log((double) HPGV::nCus));
    // initialize array pFreq
    pFreq.resize(HPGV::numRoute);
    for (unsigned int ip = 0; ip < pFreq.size(); ip++){
        pFreq[ip].resize(HPGV::nCus);
        for (unsigned int ic = 0; ic < HPGV::nCus; ic++){
            pFreq[ip][ic] = 0;
        }
    }

    // main loop
    for (int it = 0; it < HPGV::utsIter; it++){
        // Generate N(s)
        HGAGenome bestNeighbor(hs);

        for (int uni = 0; uni < 50; uni++){
            if (GAFlipCoin(0.5)){
                if (HGAGenome::UTSNeighborByRouting(bestNeighbor, g_tabu, pFreq, it, tabuLength, aQ, bD, cW)){
                    // HGAGenome::printSolution(bestNeighbor, "utsRouting.txt");
                    break;
                }
            }else{
                if (HGAGenome::UTSNeighborByPattern(bestNeighbor, g_tabu, pFreq, it, tabuLength, aQ, bD, cW)){
                    // HGAGenome::printSolution(bestNeighbor, "utsPattern.txt");
                    break;
                }
            }
        }

        // if _s_ is feasible and c(_s_) < c(s1)
        if ((bestNeighbor).isFeasible && (bestNeighbor).durationCost < hs1.durationCost){
            hs1 = bestNeighbor;
            flag = 1;
        }else if (utsFit(bestNeighbor, aQ, bD, cW) < utsFit(hs2, aQ, bD, cW)){
            hs2 = bestNeighbor;
        }
        // Compute q(_s_), d(_s_), w(_s_) and update penalty parameters
        if ((bestNeighbor).totalCapacityVio > 0){
            aQ /= 1.5;
        }else{
            aQ *= 1.5;
        }

        if ((bestNeighbor).totalDurationVio > 0){
            bD /= 1.5;
        }else{
            bD *= 1.5;
        }

        if ((bestNeighbor).totalTimeVio > 0){
            cW /= 1.5;
        }else{
            cW *= 1.5;
        }

        hs = bestNeighbor;
    }

    if (!hg.isFeasible){
        hs1.durationCost = POSINF;
    }
    if (flag){
        return hs1;
    }else{
        return hs2;
    }
}

/**
 * find a neighbor against pattern
 */
bool HGAGenome::UTSNeighborByPattern(HGAGenome& hg, TabuMap& g_tabu, vector<vector<int> >& pFreq,
        int& currUtsIter, int& tabuLength, double& aQ, double& bD, double& cW){
    bool foundNewCid = false;
    unsigned int oldPattern, newPattern, insertMask, removeMask;

    unsigned int iDay = 0;
    unsigned int iVeh = 0;
    unsigned int currVod = 0;
    unsigned int newVod = 0;
    unsigned int newVeh = 0;

    VCus tmpCus(0);

    int maxTries = 50;    // avoid infinite loop
    int counter = 0;

    while ((!foundNewCid) && (counter < maxTries)){
        // select customer that has more than 1 pattern
        int rc = GARandomInt(0, HPGV::nCus - 1);
        if (hg.arrC[rc].a > 1){
            Customer* mixer = &(hg.arrC[rc]);

            oldPattern = hg.m_pattern[rc];
            int ord = GARandomInt(0, hg.arrC[rc].a - 1);
            while (1){
                if (hg.arrC[rc].comb[ord] != oldPattern){
                    newPattern = hg.arrC[rc].comb[ord];
                    break;
                }
                ord++;
                ord %= hg.arrC[rc].a;
            }
            hg.arrC[rc].pattern = newPattern;
            insertMask = newPattern;
            removeMask = oldPattern;
            insertMask ^= (newPattern & oldPattern);
            removeMask ^= (newPattern & oldPattern);

            for (iDay = 0; iDay < HPGV::tDay; iDay++){
                int flagInsert = (int) pow(2, (double) (HPGV::tDay - iDay - 1));
                flagInsert &= insertMask;
                if (flagInsert){
                    int randVeh = GARandomInt(0, HPGV::mVeh - 1);
                    for (unsigned int iv = 0; iv < HPGV::mVeh; iv++){
                        newVeh = (randVeh + iv) % HPGV::mVeh;
                        newVod = iDay * HPGV::mVeh + newVeh;
                        if (!isTabu(mixer->id, newVod, g_tabu)){
                            foundNewCid = true;
                            break;
                        }
                    }
                }
                if (foundNewCid){
                    break;
                }
            }
            // update genome by insert customer into new routes and remove it from old routes
            if (foundNewCid){
                for (iDay = 0; iDay < HPGV::tDay; iDay++){
                    int flagInsert = (int) pow(2, (double)(HPGV::tDay - iDay -1));
                    int flagRemove = flagInsert;
                    flagInsert &= insertMask;
                    flagRemove &= removeMask;

                    if (flagInsert){
                        // check if this customer has been serviced on current day
                        bool needServiced = true;

                        for (iVeh = 0; iVeh < HPGV::mVeh; iVeh++){
                            // for r \in R
                            currVod = iDay * HPGV::mVeh + iVeh;
                            if (HGAGenome::isInRoute(hg.m_route[currVod], mixer->id)){
                                needServiced = false;
                                break;
                            }
                        }
                        // insert into one route of this day
                        if (needServiced){
                            int randVeh = GARandomInt(0, HPGV::mVeh - 1);
                            for (unsigned int iv = 0; iv < HPGV::mVeh; iv++){
                                newVeh = (randVeh + iv) % HPGV::mVeh;
                                newVod = iDay * HPGV::mVeh + newVeh;
                                if (!isTabu(mixer->id, newVod, g_tabu)){
                                    HGAGenome::PRinsert(hg.m_route[newVod], hg.m_data[newVod], mixer);
                                    needServiced = false;
                                    break;
                                }
                            }
                        }
                        if (needServiced){
                            newVeh = GARandomInt(0, HPGV::mVeh - 1);
                            newVod = iDay * HPGV::mVeh + newVeh;
                            HGAGenome::PRinsert(hg.m_route[newVod], hg.m_data[newVod], mixer);
                        }
                        // HGAGenome::printSolution(hg, "utsPattern.txt");
                    }else if (flagRemove){
                        // remove customer from current day
                        for (iVeh = 0; iVeh < HPGV::mVeh; iVeh++){
                            // for r \in R
                            currVod = iDay * HPGV::mVeh + iVeh;
                            if (hg.m_route[currVod].empty()){
                                continue;
                            }else{
                                if (HGAGenome::isInRoute(hg.m_route[currVod], mixer->id)){
                                    HGAGenome::removeFromRoute(hg.m_route[currVod], hg.m_data[currVod], mixer->id);

                                    // HGAGenome::testRoute(hg.m_route[currVod]);
                                    // HGAGenome::printSolution(hg, "utsPattern.txt");

                                    unsigned int key = (unsigned int)POSINF * (mixer->id) + currVod;
                                    TabuMap::iterator pos = g_tabu.find(key);
                                    if (pos != g_tabu.end()){
                                        pos->second++;
                                    }else{
                                        g_tabu.insert(make_pair(key, 0));
                                    }

                                    break;
                                }
                            }
                        }
                    }
                }
                hg.tourUpdate(pFreq);
                hg.updateTotalVio();
                // update tabu list
                updateTabuMap(g_tabu, tabuLength);
            }
        }
        counter++;
    }
    return foundNewCid;
}

/**
 * find a neighbor against routing
 */
bool HGAGenome::UTSNeighborByRouting(HGAGenome& hg, TabuMap& g_tabu, vector<vector<int> >& pFreq,
        int& currUtsIter, int& tabuLength, double& aQ, double& bD, double& cW){
    bool foundNewCid = false;
    unsigned int cid, vod, currDay, currVeh, newVeh;

    unsigned int newVod, newRouteSize;
    Route::iterator prevCus, bestPos, uIter, saveIter;
    Route::reverse_iterator ruIter, rprevIter;

    int maxTries = 50;    // avoid infinite loop
    int counter = 0;

    while ((!foundNewCid) && (counter < maxTries)){
        int rc = GARandomInt(0, hg.m_tour.size() - 1);
        cid = (unsigned int)hg.m_tour[rc]->cid;
        vod = (unsigned int) hg.m_tour[rc]->vod;
        currDay = vod / HPGV::mVeh;
        currVeh = vod % HPGV::mVeh;
        int randVeh = GARandomInt(0, HPGV::mVeh - 1);

        for (unsigned int iv = 1; iv < HPGV::mVeh; iv++){
            newVeh = (randVeh + iv) % HPGV::mVeh;
            if (newVeh == currVeh){
                continue;
            }
            newVod = currDay * HPGV::mVeh + newVeh;
            if (!isTabu(cid, newVod, g_tabu)){
                foundNewCid = true;
                break;
            }
        }
        // can not found any route for moving this customer
        counter++;
    }
    if (foundNewCid){
        // now we apply routing modification for this customer

        // insert into new route (another vehicle, same day)
        newRouteSize = hg.m_route[newVod].size();
        // initialize changes in violations of load, duration and time windows by inserting
        // new customer at the end of route
        Route tempRoute;
        Route saveRoute = hg.m_route[newVod];

        Customer* choice = &(hg.arrC[cid - 1]);
        double minObj = 0;
        tempRoute.clear();

        VertexPtr vChoice = VertexPtr(new Vertex(choice));

        hg.m_route[newVod].push_back(vChoice);
        tempRoute = hg.m_route[newVod];

        Route bkpVodRoute = hg.m_route[newVod];
        HGAGenome::updateInfo(tempRoute, hg.m_data[newVod]);

        minObj = hg.m_data[newVod]->cost;

        if (hg.m_data[newVod]->load > HPGV::maxLoad){
            minObj += aQ * (hg.m_data[newVod]->load - HPGV::maxLoad);
        }
        if (hg.m_data[newVod]->cost > HPGV::maxDuration){
            minObj += bD * (hg.m_data[newVod]->cost - HPGV::maxDuration);
        }
        minObj += cW * (hg.m_data[newVod]->timeVio);

        for (unsigned int i = 0; i < newRouteSize; i++){
            hg.m_route[newVod] = saveRoute;
            tempRoute.clear();
            for (unsigned int j = 0; j < i; j++){
                tempRoute.push_back(hg.m_route[newVod].front());
                hg.m_route[newVod].pop_front();
            }
            tempRoute.push_back(vChoice);
            tempRoute.splice(tempRoute.end(), hg.m_route[newVod], hg.m_route[newVod].begin(), hg.m_route[newVod].end());

            HGAGenome::updateInfo(tempRoute, hg.m_data[newVod]);

            double newObjVal = hg.m_data[newVod]->cost;
            if (hg.m_data[newVod]->load > HPGV::maxLoad){
                newObjVal += aQ * (hg.m_data[newVod]->load - HPGV::maxLoad);
            }
            if (hg.m_data[newVod]->cost > HPGV::maxDuration){
                newObjVal += bD * (hg.m_data[newVod]->cost - HPGV::maxDuration);
            }
            newObjVal += cW * (hg.m_data[newVod]->timeVio);

            if (newObjVal < minObj){
                bkpVodRoute = tempRoute;
                minObj = newObjVal;
            }
        }

        // here we found best neighbor
        hg.m_route[newVod] = bkpVodRoute;

        // remove customer from current route

        for (Route::iterator rIter = hg.m_route[vod].begin(), endR = hg.m_route[vod].end(); rIter != endR; ++rIter){
            if ((*rIter)->cus->id == cid){
                hg.m_route[vod].erase(rIter);
                break;
            }
        }
        HGAGenome::updateInfo(hg.m_route[vod], hg.m_data[vod]);
        HGAGenome::updateInfo(hg.m_route[newVod], hg.m_data[newVod]);

        hg.tourUpdate(pFreq);
        hg.updateTotalVio();

        // update tabu list
        unsigned int key = (unsigned int)POSINF * cid + vod;
        g_tabu.insert(make_pair(key, 0));
        updateTabuMap(g_tabu, tabuLength);
    }

    return foundNewCid;
}
