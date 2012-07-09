#include "HGAGenome.h"

extern vector<Customer> gArrC;
extern vector<vector<double> > gDistance;
extern Customer* gDepot;
/**
 * compare the angles that 2 customers make with the depot
 * ====================================
 */
bool compareAngular(const Customer* c1, const Customer* c2) {
    double areaTrapezoid = (c1->x - gDepot->x) * (c2->y - gDepot->y)
            - (c1->y - gDepot->y) * (c2->x - gDepot->x);
    return (areaTrapezoid > 0);
}
bool compareEndingTime(Customer* c1, Customer* c2) {
    return ((c1->l - c2->l) > 0);
}
void HGAGenome::initRouting(unsigned int& iDay, Route& mRoute, VCus& mCluster, RinfoPtr& mRinfo){
    // TODO: initRouting (treat one route with a cluster of customers)
    // check patterns and remove all customers that should not be visited on current day
    // of course, we also remove all customers that has been satisfied
    for (VCus::iterator skipCus = mCluster.begin(), endCluster = mCluster.end(); skipCus != endCluster; ++skipCus){
        if ((*skipCus)->isServiced) {
            mCluster.erase(skipCus);
            skipCus--;
            continue;
        }
        int flag = (int) pow(2, (double) (HPGV::tDay - iDay - 1));
        flag &= (*skipCus)->pattern;
        if (flag == 0) {
            mCluster.erase(skipCus);
            skipCus--;
        }
    }
    if (mCluster.empty()){
        return;
    }
    sort(mCluster.begin(), mCluster.end(), compareEndingTime);
    HGAGenome::pushbackRoute(mRoute, mRinfo, mCluster.back());
    mCluster.pop_back();

    // sequentially add remaining customers
    HGAGenome::SolomonI1(mRoute, mRinfo, mCluster);
}
/**
 * Initialize solution using "cluster first - route second"
 */
void HGAGenome::clusterFirstInit(VCus& refArr){
    // TODO: clusterFirstInit
    unsigned int clusterSize = 1;
    unsigned int iDay = 0;
    unsigned int iVeh = 0;
    unsigned int vod = 0;
    vector<VCus> cluster(1, VCus(0));
    vector<VCus> tempCluster(1, VCus(0));
    // clustering
    // sort customers in increasing order of the angle they make with the depot
    sort(refArr.begin(), refArr.end(), compareAngular);
    if ((HPGV::nCus % HPGV::mVeh) == 0) {
        clusterSize = (int) (HPGV::nCus / HPGV::mVeh);
    } else {
        clusterSize = (int) (HPGV::nCus / HPGV::mVeh) + 1;
    }

    int pivot = GARandomInt(0, HPGV::nCus - 1);
    cluster.resize(HPGV::mVeh);
    for (unsigned int i = 0; i < HPGV::mVeh; i++){
        cluster[i].resize(clusterSize);
        for(unsigned int j = 0; j < clusterSize; j++){
            cluster[i][j] = refArr[pivot];
            pivot = (1 + pivot) % HPGV::nCus;
        }
    }
    // start routing here
    for (iDay = 0; iDay < HPGV::tDay; iDay++){
        // solve VRPTW for each day
        tempCluster.clear();
        tempCluster = cluster;
        for (iVeh = 0; iVeh < HPGV::mVeh; iVeh++){
            vod = iDay * HPGV::mVeh + iVeh;
            this->m_data[vod] = (RinfoPtr)(new RouteInfo());
            this->m_route[vod].clear();
            HGAGenome::initRouting(iDay, this->m_route[vod], tempCluster[iVeh], this->m_data[vod]);
        }
    }
}
void HGAGenome::SolomonTONNInit(VCus& refArr){
    // TODO: SolomonTONNInit
    cout << "solomon" << endl;
    exit(1);
}
