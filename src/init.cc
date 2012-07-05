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
void HGAGenome::initRouting(Route& mRoute, vector<Customer*>& mCluster, RinfoPtr& mRinfo){
    // TODO: initRouting (treat one route with a cluster of customers)
}
void HGAGenome::clusterFirstInit(vector<Customer*>& refArr){
    // TODO: clusterFirstInit
    unsigned int clusterSize = 1;
    unsigned int iDay = 0;
    unsigned int iVeh = 0;
    unsigned int vod = 0;
    vector<vector<Customer*> > cluster(1, vector<Customer*>(0));
    vector<vector<Customer*> > tempCluster(1, vector<Customer*>(0));
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
    for (int i = 0; i < HPGV::mVeh; i++){
        cluster[i].resize(clusterSize);
        for(int j = 0; j < clusterSize; j++){
            cluster[i][j] = refArr[pivot];
            pivot = (1 + pivot) % HPGV::nCus;
        }
    }
    // start routing here
    // solve VRPTW for each day
    tempCluster.clear();
    tempCluster = cluster;
    for (iVeh = 0; iVeh < HPGV::mVeh; iVeh++){
        vod = iDay * HPGV::mVeh + iVeh;
        this->m_data[vod] = (RinfoPtr)(new RouteInfo());
        this->m_route[vod].clear();
        HGAGenome::initRouting(this->m_route[vod], tempCluster[iDay], this->m_data[vod]);
    }
}
void HGAGenome::SolomonTONNInit(vector<Customer*>& refArr){
    // TODO: SolomonTONNInit
    cout << "solomon" << endl;
    exit(1);
}
