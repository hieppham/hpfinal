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
bool isTabu(int cid, unsigned int vod, TabuMap& tabuMap){
    CustomerInDay* key = new CustomerInDay(cid, vod);
    TabuMap::iterator pos = tabuMap.find(key);
    delete key;
    return (pos != tabuMap.end());
}

HGAGenome HGAGenome::UTS(HGAGenome& hg){
    // TODO: UTS
    double aQ = 1;
    double bD = 1;
    double cW = 1;

    TabuMap g_tabu;
    // number of times attribute (i, k) has been added to the solution during the search process
    vector<vector<int> > pFreq;

    HGAGenome hs = hg;
    HGAGenome hs2 = hg;
    HGAGenome hs1 = hg;
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
        HGAGenome* bestNeighbor = new HGAGenome(hs);

        while(1){
            if (GAFlipCoin(0.5)){
                if (HGAGenome::UTSNeighborByRouting(*bestNeighbor, g_tabu, pFreq, it, tabuLength, aQ, bD, cW)){
                    break;
                }
            }else{
                if (HGAGenome::UTSNeighborByPattern(*bestNeighbor, g_tabu, pFreq, it, tabuLength, aQ, bD, cW)){
                    break;
                }
            }
        }

        // if _s_ is feasible and c(_s_) < c(s1)
        if ((*bestNeighbor).isFeasible && (*bestNeighbor).durationCost < hs1.durationCost){
            hs1 = *bestNeighbor;
            flag = 1;
        }else if (utsFit(*bestNeighbor, aQ, bD, cW) < utsFit(hs2, aQ, bD, cW)){
            hs2 = *bestNeighbor;
        }
        // Compute q(_s_), d(_s_), w(_s_) and update penalty parameters
        if ((*bestNeighbor).totalCapacityVio > 0){
            aQ /= 1.5;
        }else{
            aQ *= 1.5;
        }

        if ((*bestNeighbor).totalDurationVio > 0){
            bD /= 1.5;
        }else{
            bD *= 1.5;
        }

        if ((*bestNeighbor).totalTimeVio > 0){
            cW /= 1.5;
        }else{
            cW *= 1.5;
        }

        hs = *bestNeighbor;
        delete bestNeighbor;
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
    // TODO: UTSNeighborByPattern
    return foundNewCid;
}

/**
 * find a neighbor against routing
 */
bool HGAGenome::UTSNeighborByRouting(HGAGenome& hg, TabuMap& g_tabu, vector<vector<int> >& pFreq,
        int& currUtsIter, int& tabuLength, double& aQ, double& bD, double& cW){
    bool foundNewCid = false;
    // TODO: UTSNeighborByRouting
    return foundNewCid;
}
