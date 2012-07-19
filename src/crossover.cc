#include "HGAGenome.h"

int HGAGenome::explorationCrossover(const HGAGenome& p1, const HGAGenome& p2, HGAGenome& child){
    unsigned int iDay = 0;
    unsigned int tempSize = 0;
    vector<vector<int> > checkInherit(0, vector<int>(0));
    vector<int>::iterator findPos;

    VCus backupArr(0);
    VCus cloneArr(0);  // copy of current cluster
    VCus::iterator minVPos, skipCus, kIter;

    Route::iterator prevCus, nextCus, uIter, beforeIter;
    checkInherit.resize(HPGV::tDay);

    // STEP 1: Assign a pattern to each customer
    // inherit pattern assignments from parent P1
    tempSize = p1.m_tour.size();
    child.m_pattern.resize(HPGV::nCus);

    int cutPointOne = GARandomInt(0, tempSize - 2);
    int cutPointTwo = GARandomInt(cutPointOne + 1, tempSize - 1);

    for (int i = cutPointOne; i <= cutPointTwo; i++){
        int lcid = p1.m_tour[i]->cid - 1;
        iDay = (int) (p1.m_tour[i]->vod / HPGV::mVeh);
        child.m_pattern[lcid] |= (int) pow(2, (double)(HPGV::tDay - iDay - 1));
        child.arrC[lcid].pattern = child.m_pattern[lcid];
        child.arrC[lcid].token++;
    }
    // inherit pattern assignments from pattern P2
    for (iDay = 0; iDay < HPGV::tDay; iDay++){
        for (unsigned int i = 0; i < HPGV::nCus; i++){
            if (child.arrC[i].token < child.arrC[i].f){
                int flag = (int) pow(2, (double)(HPGV::tDay - iDay - 1));
                int flagBkp = flag;
                for (int j = 0; j < child.arrC[j].a; j++){
                    flag &= child.arrC[i].comb[j];
                    if (flag == 0){
                        continue;
                    }else{
                        flagBkp |= child.arrC[i].pattern;
                        flag = flagBkp;
                        flagBkp &= child.arrC[i].comb[j];
                        if (flagBkp != 0){
                            child.arrC[i].pattern = flag;
                            child.arrC[i].token = numberOfSetBits(flag);
                            child.arrC[i].isServiced = true;
                        }
                    }
                }
            }else{
                child.arrC[i].isServiced = true;
            }
        }
    }
    // complete pattern assignments
    for (unsigned int i = 0; i < HPGV::nCus; i++){
        if (!child.arrC[i].isServiced){
            // assign a random pattern that includes day(i)
            int ord = GARandomInt(0, child.arrC[i].a - 1);
            while (1){
                int flagBkp = child.arrC[i].pattern;
                flagBkp &= child.arrC[i].comb[ord];
                if (flagBkp !=0){
                    child.arrC[i].pattern = child.arrC[i].comb[ord];
                    child.m_pattern[i] = child.arrC[i].pattern;
                    child.arrC[i].isServiced = true;
                    break;
                }
                ord++;
                ord = ord % child.arrC[i].a;
            }
        }
    }

    // tempSize = HPGV::numRoute;
    // STEP 2 : Assign customers to routes
    // Copy routes from P1
    // initialize all routes
    child.m_route.resize(HPGV::numRoute);
    child.m_data.resize(HPGV::numRoute);
    for (unsigned int tm = 0; tm < HPGV::numRoute; tm++){
        child.m_data[tm] = RinfoPtr(new RouteInfo());
    }

    backupArr.resize(HPGV::nCus);
    for (unsigned int ic = 0; ic < HPGV::nCus; ic++){
        backupArr[ic] = &(child.arrC[ic]);
    }

    // The customers between the two cutting points determined in Step 1a are routed
    // as in the P1 parent and the corresponding sequences are copied into C
    for (int i = cutPointOne; i <= cutPointTwo; i++){
        int s2cid = p1.m_tour[i]->cid - 1;
        int s2vod = p1.m_tour[i]->vod;
        int tmpDay = s2vod/HPGV::mVeh;

        HGAGenome::pushbackRoute(child.m_route[s2vod], child.m_data[s2vod], backupArr[s2cid]);

        checkInherit[tmpDay].push_back(p1.m_tour[i]->cid);
    }

    // Assign remaining customers to routes
    for (iDay = 0; iDay < HPGV::tDay; iDay++) {
        cloneArr.clear();
        cloneArr = backupArr;
        // check patterns and remove all customers that should not be visited on current day
        // of course, we also remove all customers that has been satisfied
        for (skipCus = cloneArr.begin(); skipCus != cloneArr.end(); ++skipCus) {
            int flag = (int) pow(2, (double) (HPGV::tDay - iDay - 1));
            flag &= (*skipCus)->pattern;
            if (flag == 0) {
                cloneArr.erase(skipCus);
                skipCus--;
            }else{
                findPos = find(checkInherit[iDay].begin(), checkInherit[iDay].end(), (*skipCus)->id);
                if (findPos != checkInherit[iDay].end()){
                    cloneArr.erase(skipCus);
                    skipCus--;
                    checkInherit[iDay].erase(findPos);
                }
            }
        }
        // here we use cost insertion heuristic (Potvin - Rosseau, 1993)
        // sort(cloneArr.begin(), cloneArr.end(), compareEndingTime);
        HGAGenome::PRheuristic(child.m_route, child.m_data, cloneArr, iDay, true);
    }

    // complete tour genome
    child.tourConstruct();
    child.updateTotalVio();
    // HGAGenome::printSolution(child, "childExplor.txt");
    child._evaluated = gaFalse;

    return 1;
}

int HGAGenome::exploitationCrossover(const HGAGenome& p1, const HGAGenome& p2, HGAGenome& child){
    unsigned int iDay = 0;
    unsigned int iVeh = 0;
    unsigned int vod = 0;

    cusinday tmpTour(0);

    VCus backupArr(0);

    // initialize all routes
    child.m_pattern.resize(HPGV::nCus);
    child.m_route.resize(HPGV::numRoute);
    child.m_data.resize(HPGV::numRoute);
    for (unsigned int tm = 0; tm < HPGV::numRoute; tm++){
        child.m_data[tm] = RinfoPtr(new RouteInfo());
    }

    backupArr.resize(HPGV::nCus);
    for (unsigned int ic = 0; ic < HPGV::nCus; ic++){
        backupArr[ic] = &(child.arrC[ic]);
    }

    // treats each day t in the planning horizon
    for (iDay = 0; iDay < HPGV::tDay; iDay++){
        if (GAFlipCoin(0.5)){
            tmpTour = p1.m_tour;
        }else{
            tmpTour = p2.m_tour;
        }
        int bitFlag = (int) pow(2, (double) (HPGV::tDay - iDay - 1));
        for (iVeh = 0; iVeh < HPGV::mVeh; iVeh++){
            vod = iDay * HPGV::mVeh + iVeh;
            for (cusinday::iterator cIter = tmpTour.begin(), endIter = tmpTour.end(); cIter != endIter; ++cIter){
                if ((*cIter)->vod == (int)vod){
                    int bcid = (*cIter)->cid - 1;

                    HGAGenome::pushbackRoute(child.m_route[vod], child.m_data[vod], backupArr[bcid]);

                    child.m_pattern[bcid] |= bitFlag;
                    child.arrC[bcid].pattern = child.m_pattern[bcid];
                    child.arrC[bcid].token++;
                }else if ((*cIter)->vod > (int)vod){
                    break;
                }
            }
        }
    }

    child.modifyEachCustomer();

    child.tourConstruct();
    child.updateTotalVio();
    // HGAGenome::printSolution(child, "childExploit.txt");
    child._evaluated = gaFalse;

    return 1;
}

void HGAGenome::modifyEachCustomer(void){
    unsigned int iDay = 0;
    unsigned int iVeh = 0;
    unsigned int vod = 0;

    int oldPattern, newPattern, insertMask, removeMask;
    vector<int>::iterator findPos;
    Customer* mixer;
    VCus tmpCus(0);

    Route::iterator prevCustomer, nextCustomer, uIter, beforeIter, endIter;

    for (unsigned int i = 0; i < HPGV::nCus; i++){
        mixer = &(this->arrC[i]);
        // assign a random pattern that includes day(i)
        int ord = GARandomInt(0, this->arrC[i].a - 1);
        oldPattern = this->arrC[i].pattern;
        if (oldPattern != 0){
            while(1){
                int flagBkp = oldPattern;
                flagBkp &= this->arrC[i].comb[ord];
                if (flagBkp !=0){
                    this->arrC[i].pattern = this->arrC[i].comb[ord];
                    this->m_pattern[i] = this->arrC[i].pattern;
                    this->arrC[i].isServiced = true;
                    break;
                }
                ord++;
                ord = ord % this->arrC[i].a;
            }
        }else{
            this->arrC[i].pattern = this->arrC[i].comb[ord];
        }

        newPattern = this->arrC[i].pattern;
        insertMask = newPattern;
        removeMask = oldPattern;
        insertMask ^= (newPattern & oldPattern);
        removeMask ^= (newPattern & oldPattern);

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
                    if (this->m_route[currVod].empty()){
                        continue;
                    }else{
                        for (uIter = this->m_route[currVod].begin(), endIter = this->m_route[currVod].end(); uIter != endIter; ++uIter){
                            if ((*uIter)->cus->id == mixer->id){
                                needServiced = false;
                                break;
                            }
                        }
                    }
                }
                // here we insert customer into one route of this day
                if (needServiced){
                    tmpCus.clear();
                    tmpCus.push_back(mixer);
                    HGAGenome::PRheuristic(this->m_route, this->m_data, tmpCus, iDay, false);
                }
            }else if (flagRemove){
                for (iVeh = 0; iVeh < HPGV::mVeh; iVeh++){
                    vod = iDay * HPGV::mVeh + iVeh;
                    if (this->m_route[vod].empty()){
                        continue;
                    }else{
                        for (uIter = this->m_route[vod].begin(), endIter = this->m_route[vod].end(); uIter != endIter; ++uIter){
                            if ((*uIter)->cus->id == mixer->id){
                                HGAGenome::removeFromRoute(this->m_route[vod], this->m_data[vod], mixer->id);
                                break;
                            }
                        }
                    }
                }
            }
        }
    }
}
