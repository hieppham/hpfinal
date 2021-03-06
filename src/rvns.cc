#include "HGAGenome.h"

HGAGenome HGAGenome::RVNS(HGAGenome& hgenome){
    HGAGenome bestNeighbor(hgenome);    // s*
    HGAGenome orig(hgenome);
    HGAGenome hgs(hgenome);     // s'

    double T = 10;      // initial temperature
    double dT = 100;    // parameter for linear cooling scheme
    vector<unsigned int> ar;
    unsigned int curIter = 1;
    unsigned int curNeighbor, k;
    double curObjVal, recordVal;
    double probAccept = 0;

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

    curObjVal = HGAGenome::calcObjectValue(orig);
    recordVal = curObjVal;

    ar.clear();
    ar.resize(18);
    for (unsigned int i = 1; i <= 18; i++){
        ar[i - 1] = i;
    }

    do {
        random_shuffle(ar.begin(), ar.end());

        curNeighbor = 1;
        while ((curNeighbor <= 18) && (curIter <= HPGV::rvnsIter)){
            k = ar[curNeighbor - 1];
            hgs = HGAGenome::Shaking(hgenome, k, HPGV::rvnsPRev);
            // HGAGenome::apply2OptForAllRoutes(aPen, bPen, cPen, hgs);
            // update 310712: priotize stochastic strategy
            HGAGenome::apply2OptUntilFirstImprovement(aPen, bPen, cPen, hgs);
            double newObjVal = HGAGenome::calcObjectValue(hgs);

            // if (s' better than s)
            if (newObjVal < curObjVal){
                // update 310712: priotize stochastic strategy
                // HGAGenome::apply2OptStarForAllRoutes(aPen, bPen, cPen, hgs);
                HGAGenome::apply2OptStarUntilFirstImprovement(aPen, bPen, cPen, hgs);
                // s <- s'
                orig = hgs;
                curObjVal = newObjVal;

                // if (s better than s*)
                if (curObjVal < recordVal){
                    // update record
                    bestNeighbor = orig;
                    recordVal = curObjVal;
                }
                curNeighbor = 1;
                random_shuffle(ar.begin(), ar.end());
            }else{
                // probability for accepting inferior solution
                probAccept = exp((curObjVal - newObjVal) / T);
                if (GAFlipCoin(float(probAccept))){
                    // s <- s'
                    orig = hgs;
                    curObjVal = newObjVal;

                    curNeighbor = 1;
                    random_shuffle(ar.begin(), ar.end());
                }else{
                    curNeighbor++;
                }
            }
            curIter++;
            // Update temperature T if needed
            T = T * (1 - dT/HPGV::rvnsTmax);
        }

    } while (curIter <= HPGV::rvnsIter);

    return bestNeighbor;
}

HGAGenome HGAGenome::Shaking(HGAGenome& hgenome, unsigned int k, double pRev){
    HGAGenome hgs(hgenome);

    if ((k >= 1) && (k <= 6)){
        hgs = HGAGenome::ShakingPattern(hgenome, k, pRev);
    }else if ((k >= 7) && (k <= 12)){
        hgs = HGAGenome::ShakingMoveSegment(hgenome, k - 6, pRev);
    }else if ((k >= 13) && (k <= 18)){
        hgs = HGAGenome::ShakingExchangeSegments(hgenome, k - 12, pRev);
    }

    return hgs;
}

HGAGenome HGAGenome::ShakingPattern(HGAGenome& hgenome, unsigned int k, double pRev){
    HGAGenome hg(hgenome);
    unsigned int oldPattern, newPattern, insertMask, removeMask;

    unsigned int iDay = 0;
    unsigned int iVeh = 0;
    unsigned int vod = 0;
    VCus tmpCus(0);

    // TODO: remove after debugging
    // cout << HPGV::genCounter << " - shaking Pattern..." << flush;

    // move to neighbor by changing pattern up to k time(s)
    unsigned int counter = 0;

    do {
        // select customer that has more than 1 pattern
        unsigned int trc = GARandomInt(0, HPGV::numOfMP - 1);
        unsigned int rc = HPGV::multipattern[trc] - 1;
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
        hg.m_pattern[rc] = newPattern;
        hg.arrC[rc].pattern = newPattern;

        insertMask = newPattern;
        removeMask = oldPattern;
        insertMask ^= (newPattern & oldPattern);
        removeMask ^= (newPattern & oldPattern);

        // update genome by insert customer into new routes and remove it from old routes
        for (iDay = 0; iDay < HPGV::tDay; iDay++){
            int flagInsert = (int) pow(2, (double) (HPGV::tDay - iDay -1));
            int flagRemove = flagInsert;
            flagInsert &= insertMask;
            flagRemove &= removeMask;

            if (flagInsert){
                // check if this customer has been serviced on current day
                bool needServiced = true;
                for (iVeh = 0; iVeh < HPGV::mVeh; iVeh++){
                    // for r \in R
                    unsigned currVod = iDay * HPGV::mVeh + iVeh;
                    if (HGAGenome::isInRoute(hg.m_route[currVod], mixer->id)){
                        needServiced = false;
                        break;
                    }
                }
                // insert into one route of this day
                if (needServiced){
                    tmpCus.clear();
                    tmpCus.push_back(mixer);
                    HGAGenome::PRheuristic(hg.m_route, hg.m_data, tmpCus, iDay, false);
                }
            } else if (flagRemove){
                // remove customer from current day
                for (iVeh = 0; iVeh < HPGV::mVeh; iVeh++){
                    vod = iDay * HPGV::mVeh + iVeh;
                    if (HGAGenome::isInRoute(hg.m_route[vod], mixer->id)){
                        HGAGenome::removeFromRoute(hg.m_route[vod], hg.m_data[vod], mixer->id);
                        break;
                    }
                }
            }
        }
        counter++;
    } while (counter < k);

    hg.tourConstruct();
    // HGAGenome::printSolution(hg, "ShakingPattern.txt");
    // TODO: remove after debugging
    // cout << "end" << endl;

    return hg;
}

HGAGenome HGAGenome::ShakingMoveSegment(HGAGenome& hgenome, unsigned int k, double pRev){
    HGAGenome hg(hgenome);
    Route firstHead, movSeg, firstTail, secondHead;
    // move to neighbor by moving a random segment of maximal length k of a route
    // to another route on the same day
    unsigned int iDay, iVF, iVS, vod1, vod2, fSize, sSize;

    int maxTries = 10;
    int tryCounter = 0;

    bool flagStop = true;

    // TODO: remove after debugging
    // cout << HPGV::genCounter << " - rvns shaking MoveSegment..." << flush;

    while (tryCounter < maxTries){
        iDay = GARandomInt(0, HPGV::tDay - 1);
        iVF = GARandomInt(0, HPGV::mVeh - 2);

        iVS = GARandomInt(iVF + 1, HPGV::mVeh - 1);

        vod1 = iDay * HPGV::mVeh + iVF;
        vod2 = iDay * HPGV::mVeh + iVS;
        fSize = hg.m_route[vod1].size();
        sSize = hg.m_route[vod2].size();
        if ((fSize >= 1) && (sSize >= 1)){
            flagStop = false;
            break;
        }
        tryCounter++;
    }

    if (flagStop){
        return hg;
    }
    // select a random segment of maximal length k from first route

    unsigned int pivot1 = GARandomInt(1, fSize);
    unsigned int pivot2 = GARandomInt(1, sSize);

    firstHead.clear();
    for (unsigned int i = 0; i < pivot1 - 1; i++){
        firstHead.push_back(hg.m_route[vod1].front());
        hg.m_route[vod1].pop_front();
    }
    unsigned int counter = 0;
    movSeg.clear();
    if (k != 6){
        do{
            movSeg.push_back(hg.m_route[vod1].front());
            hg.m_route[vod1].pop_front();
            counter++;
        } while ((counter < k) && (!hg.m_route[vod1].empty()));
    }else{
        do{
            movSeg.push_back(hg.m_route[vod1].front());
            hg.m_route[vod1].pop_front();
            counter++;
        } while (!hg.m_route[vod1].empty());
    }

    // re-connect first tour
    firstTail.clear();
    firstTail = hg.m_route[vod1];

    hg.m_route[vod1].clear();
    hg.m_route[vod1] = firstHead;

    hg.m_route[vod1].splice(hg.m_route[vod1].end(), firstTail, firstTail.begin(), firstTail.end());

    // reverse the segment with a small probability pRev
    if (GAFlipCoin((float)pRev)){
        movSeg.reverse();
    }

    // insert new segment into second route
    secondHead.clear();
    for (unsigned int j = 0; j < pivot2 - 1; j++){
        secondHead.push_back(hg.m_route[vod2].front());
        hg.m_route[vod2].pop_front();
    }
    secondHead.splice(secondHead.end(), movSeg, movSeg.begin(), movSeg.end());
    secondHead.splice(secondHead.end(), hg.m_route[vod2], hg.m_route[vod2].begin(), hg.m_route[vod2].end());

    hg.m_route[vod2].clear();
    hg.m_route[vod2] = secondHead;

    HGAGenome::delayDeparture(hg.m_route[vod1], hg.m_data[vod1]);
    HGAGenome::delayDeparture(hg.m_route[vod2], hg.m_data[vod2]);

    hg.tourConstruct();
    hg.updateTotalVio();
    // HGAGenome::printSolution(hg, "ShakingMoveSegment.txt");
    // TODO: remove after debugging
    // cout << "end" << endl;

    return hg;
}

HGAGenome HGAGenome::ShakingExchangeSegments(HGAGenome& hgenome, unsigned int k, double pRev){
    HGAGenome hg(hgenome);
    // move to neighbor by exchange 2 random segments of maximum length k of 2 routes
    // on the same day
    Route firstHead, firstSeg, secondSeg, secondHead;
    unsigned int iDay, iVF, iVS, vod1, vod2, fSize, sSize;

    int maxTries = 10;
    int tryCounter = 0;
    bool flagStop = true;

    // TODO: remove after debugging
    // cout << HPGV::genCounter << " - rvns shaking ExchangeSegments..." << flush;

    while (tryCounter < maxTries){
        iDay = GARandomInt(0, HPGV::tDay - 1);
        iVF = GARandomInt(0, HPGV::mVeh - 2);
        iVS = GARandomInt(iVF + 1, HPGV::mVeh - 1);

        vod1 = iDay*HPGV::mVeh + iVF;
        vod2 = iDay*HPGV::mVeh + iVS;

        fSize = hg.m_route[vod1].size();
        sSize = hg.m_route[vod2].size();

        if ((fSize >= 1) && (sSize >= 1)){
            flagStop = false;
            break;
        }

        tryCounter++;
    }

    if (flagStop){
        return hg;
    }

    unsigned int pivot1 = GARandomInt(1, fSize);
    unsigned int pivot2 = GARandomInt(1, sSize);

    // select a random segment of maximal length k from first route
    firstHead.clear();
    for (unsigned int i = 0; i < pivot1 - 1; i++){
        firstHead.push_back(hg.m_route[vod1].front());
        hg.m_route[vod1].pop_front();
    }
    unsigned int counter = 0;
    firstSeg.clear();
    if (k != 6){
        do{
            firstSeg.push_back(hg.m_route[vod1].front());
            hg.m_route[vod1].pop_front();
            counter++;
        } while ((counter < k) && (!hg.m_route[vod1].empty()));
    }else{
        do{
            firstSeg.push_back(hg.m_route[vod1].front());
            hg.m_route[vod1].pop_front();
            counter++;
        } while (!hg.m_route[vod1].empty());
    }

    // select a random segment of maximal length k from second route
    for (unsigned int m = 0; m < pivot2 - 1; m++){
        secondHead.push_back(hg.m_route[vod2].front());
        hg.m_route[vod2].pop_front();
    }
    counter = 0;
    secondSeg.clear();
    if (k != 6){
        do{
            secondSeg.push_back(hg.m_route[vod2].front());
            hg.m_route[vod2].pop_front();
            counter++;
        } while ((counter < k) && (!hg.m_route[vod2].empty()));
    }else{
        do{
            secondSeg.push_back(hg.m_route[vod2].front());
            hg.m_route[vod2].pop_front();
            counter++;
        } while (!hg.m_route[vod2].empty());
    }

    // reverse 2 segments with a small probability pRev
    if (GAFlipCoin((float)pRev)){
        firstSeg.reverse();
        secondSeg.reverse();
    }

    // insert new segment into second route
    secondHead.splice(secondHead.end(), firstSeg, firstSeg.begin(), firstSeg.end());
    secondHead.splice(secondHead.end(), hg.m_route[vod2], hg.m_route[vod2].begin(), hg.m_route[vod2].end());

    hg.m_route[vod2].clear();
    hg.m_route[vod2] = secondHead;

    // insert new segment into first route
    firstHead.splice(firstHead.end(), secondSeg, secondSeg.begin(), secondSeg.end());
    firstHead.splice(firstHead.end(), hg.m_route[vod1], hg.m_route[vod1].begin(), hg.m_route[vod1].end());

    hg.m_route[vod1].clear();
    hg.m_route[vod1] = firstHead;

    HGAGenome::delayDeparture(hg.m_route[vod1], hg.m_data[vod1]);
    HGAGenome::delayDeparture(hg.m_route[vod2], hg.m_data[vod2]);

    hg.tourConstruct();
    hg.updateTotalVio();
    // HGAGenome::printSolution(hg, "ShakingExchangeSegments.txt");
    // TODO: remove after debugging
    // cout << "end" << endl;

    return hg;
}
