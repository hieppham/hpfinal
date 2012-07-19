#include "HGAGenome.h"

HGAGenome HGAGenome::RVNS(HGAGenome& hgenome){
    // TODO: RVNS
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
            HGAGenome::apply2OptForAllRoutes(aPen, bPen, cPen, hgs);
            double newObjVal = HGAGenome::calcObjectValue(hgs);

            // if (s' better than s)
            if (newObjVal < curObjVal){
                HGAGenome::apply2OptStarForAllRoutes(aPen, bPen, cPen, hgs);
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
    // TODO: ShakingPattern
    return hg;
}

HGAGenome HGAGenome::ShakingMoveSegment(HGAGenome& hgenome, unsigned int k, double pRev){
    HGAGenome hg(hgenome);
    // TODO: ShakingMoveSegment
    return hg;
}

HGAGenome HGAGenome::ShakingExchangeSegments(HGAGenome& hgenome, unsigned int k, double pRev){
    HGAGenome hg(hgenome);
    // TODO: ShakingExchangeSegments
    return hg;
}
