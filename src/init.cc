#include "HGAGenome.h"
#include "hpfinal.h"

/**
 * compare the angles that 2 customers make with the depot
 * ====================================
 */
bool compareAngular(const Customer* c1, const Customer* c2) {
    double areaTrapezoid = (c1->x - HPGV::gDepot->x) * (c2->y - HPGV::gDepot->y)
            - (c1->y - HPGV::gDepot->y) * (c2->x - HPGV::gDepot->x);
    return (areaTrapezoid > 0);
}
struct angularLess {
    bool operator ()(Customer* const& c1, Customer* const& c2) const {

        double areaTrapezoid = (c1->x - HPGV::gDepot->x) * (c2->y - HPGV::gDepot->y)
                    - (c1->y - HPGV::gDepot->y) * (c2->x - HPGV::gDepot->x);
        if (areaTrapezoid > 0) return true;

        return false;
    }
};
/**
 * compareEndingTime: deprecated since July 25th, 2012
 */
bool compareEndingTime(Customer* c1, Customer* c2) {
    return ((c1->l - c2->l) > 0);
}
struct endingTimeLess {
    bool operator ()(Customer* const& c1, Customer* const& c2) const {

        if (c1->l > c2->l) return true;

        return false;
    }
};

/**
 * compare the metrics from 2 customers to the depot
 * ====================================
 * deprecated since July 25th, 2012
 */
bool compareMetric(Customer* c1, Customer* c2) {
    double firstMetric = metric(HPGV::gDepot, c1, 0, HPGV::gDistance);
    double secondMetric = metric(HPGV::gDepot, c2, 0, HPGV::gDistance);
    return (firstMetric - secondMetric > 0);
}
struct metricLess {
    bool operator ()(Customer* const& c1, Customer* const& c2) const {
        double firstMetric = metric(HPGV::gDepot, c1, 0, HPGV::gDistance);
        double secondMetric = metric(HPGV::gDepot, c2, 0, HPGV::gDistance);
        if (firstMetric > secondMetric) return true;

        return false;
    }
};

bool compareMetricDynamic(Customer* c1, Customer* c2) {
    double firstMetric = metric(HPGV::gDynamic, c1, HPGV::gDynamicStart, HPGV::gDistance);
    double secondMetric = metric(HPGV::gDynamic, c2, HPGV::gDynamicStart, HPGV::gDistance);
    return (firstMetric - secondMetric > 0);
}
struct metricDynamicLess {
    bool operator ()(Customer* const& c1, Customer* const& c2) const {
        double firstMetric = metric(HPGV::gDynamic, c1, HPGV::gDynamicStart, HPGV::gDistance);
        double secondMetric = metric(HPGV::gDynamic, c2, HPGV::gDynamicStart, HPGV::gDistance);
        if (firstMetric > secondMetric) return true;

        return false;
    }
};
/**
 * function initCluster() : begin inserting customers into routes in "cluster first - route second" schema
 */
void HGAGenome::initCluster(unsigned int& iDay, Route& mRoute, VCus& mCluster, RinfoPtr& mRinfo){
    // check patterns and remove all customers that should not be visited on current day
    // of course, we also remove all customers that has been satisfied
    for (VCus::iterator skipCus = mCluster.begin(); skipCus != mCluster.end(); ++skipCus){
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
    sort(mCluster.begin(), mCluster.end(), endingTimeLess());
    HGAGenome::pushbackRoute(mRoute, mRinfo, mCluster.back());
    mCluster.pop_back();

    // sequentially add remaining customers
    HGAGenome::SolomonI1(mRoute, mRinfo, mCluster);
}
/**
 * function initSolomon()
 */
void HGAGenome::initSolomon(unsigned int& iDay, Route& mRoute, VCus& mArr, RinfoPtr& mRinfo){
    if (mArr.empty()){
        // all customers have been serviced that day
        return;
    }
    // start new route, so we find the "closest" customer to the depot
    sort(mArr.begin(), mArr.end(), endingTimeLess());

    HGAGenome::pushbackRoute(mRoute, mRinfo, mArr.back());
    mArr.back()->checkServiced();
    mArr.pop_back();

    while(1){
        HPGV::gDynamic = mRoute.back()->cus;
        HPGV::gDynamicStart = mRoute.back()->timeStartService;
        if (mArr.empty()){
            break;
        }
        sort(mArr.begin(), mArr.end(), metricDynamicLess());
        double newDuration = mRinfo->cost + HPGV::gDistance[HPGV::gDynamic->id][mArr.back()->id];
        double newStartTime = mRoute.back()->timeDeparture + HPGV::gDistance[HPGV::gDynamic->id][mArr.back()->id];
        int newLoad = mRinfo->load + mArr.back()->q;
        if ((newDuration > HPGV::maxDuration) || (newStartTime + mArr.back()->d > mArr.back()->l) || (newLoad > HPGV::maxLoad)){
            break;
        }else{
            HGAGenome::pushbackRoute(mRoute, mRinfo, mArr.back());
            mArr.back()->checkServiced();
            mArr.pop_back();
            continue;
        }
    }
    return;
}
/**
 * Initialize solution using "cluster first - route second"
 */
void HGAGenome::clusterFirstInit(VCus& refArr){
    unsigned int clusterSize = 1;
    unsigned int iDay = 0;
    unsigned int iVeh = 0;
    unsigned int vod = 0;
    vector<VCus> cluster(1, VCus(0));
    vector<VCus> tempCluster(1, VCus(0));
    // clustering
    // sort customers in increasing order of the angle they make with the depot
    sort(refArr.begin(), refArr.end(), angularLess());
    if ((HPGV::nCus % HPGV::mVeh) == 0) {
        clusterSize = (int) (HPGV::nCus / HPGV::mVeh);
    }else{
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
            RinfoPtr nDat(new RouteInfo());
            this->m_data[vod] = nDat;
            this->m_route[vod].clear();
            HGAGenome::initCluster(iDay, this->m_route[vod], tempCluster[iVeh], this->m_data[vod]);
        }
        // All routes have been initialized, so we should check if any customer left
        VCus remainClus(0);
        for (iVeh = 0; iVeh < HPGV::mVeh; iVeh++){
            if (!tempCluster[iVeh].empty()){
                remainClus.insert(remainClus.end(), tempCluster[iVeh].begin(), tempCluster[iVeh].end());
            }
        }
        if (!remainClus.empty()){
            HGAGenome::PRheuristic(this->m_route, this->m_data, remainClus, iDay, true);
        }
        // finish routing and calculate duration cost
        for (iVeh = 0; iVeh < HPGV::mVeh; iVeh++){
            vod = iDay * HPGV::mVeh + iVeh;
            HGAGenome::delayDeparture(this->m_route[vod], this->m_data[vod]);
        }
    }
}
void HGAGenome::SolomonTONNInit(VCus& refArr){
    unsigned int iDay = 0;
    unsigned int iVeh = 0;
    unsigned int vod = 0;
    VCus cloneArr(0);

    for (iDay = 0; iDay < HPGV::tDay; iDay++){
        // solve VRPTW for each day
        cloneArr.clear();
        cloneArr = refArr;
        VCus::iterator skipCus;
        // check patterns and remove all customers that should not be visited on current day
        // of course, we also remove all customers that has been satisfied
        for (skipCus = cloneArr.begin(); skipCus != cloneArr.end();
                ++skipCus) {
            if ((*skipCus)->isServiced) {
                cloneArr.erase(skipCus);
                skipCus--;
                continue;
            }
            int flag = (int) pow(2, (double) (HPGV::tDay - iDay - 1));
            flag &= (*skipCus)->pattern;
            if (flag == 0) {
                cloneArr.erase(skipCus);
                skipCus--;
            }
        }
        for (iVeh = 0; iVeh < HPGV::mVeh; iVeh++){
            vod = iDay * HPGV::mVeh + iVeh;
            RinfoPtr nDat(new RouteInfo());
            this->m_data[vod] = nDat;
            this->m_route[vod].clear();
            HGAGenome::initSolomon(iDay, this->m_route[vod], cloneArr, this->m_data[vod]);
        }
        // All routes have been initialized, so we should check if any customer left
        if (!cloneArr.empty()){
            HGAGenome::PRheuristic(this->m_route, this->m_data, cloneArr, iDay, true);
        }
        // finish routing and calculate duration cost
        for (iVeh = 0; iVeh < HPGV::mVeh; iVeh++){
            vod = iDay * HPGV::mVeh + iVeh;
            HGAGenome::delayDeparture(this->m_route[vod], this->m_data[vod]);
        }
    }
}

void HGAGenome::PRInit(VCus& refArr){
    unsigned int iDay = 0;
    unsigned int iVeh = 0;
    unsigned int vod = 0;
    VCus cloneArr(0);

    for (iDay = 0; iDay < HPGV::tDay; iDay++){
        // solve VRPTW for each day
        cloneArr.clear();
        cloneArr = refArr;
        VCus::iterator skipCus;
        // check patterns and remove all customers that should not be visited on current day
        // of course, we also remove all customers that has been satisfied
        for (skipCus = cloneArr.begin(); skipCus != cloneArr.end();
                ++skipCus) {
            if ((*skipCus)->isServiced) {
                cloneArr.erase(skipCus);
                skipCus--;
                continue;
            }
            int flag = (int) pow(2, (double) (HPGV::tDay - iDay - 1));
            flag &= (*skipCus)->pattern;
            if (flag == 0) {
                cloneArr.erase(skipCus);
                skipCus--;
            }
        }
        for (iVeh = 0; iVeh < HPGV::mVeh; iVeh++){
            vod = iDay * HPGV::mVeh + iVeh;
            RinfoPtr nDat(new RouteInfo());
            this->m_data[vod] = nDat;
            this->m_route[vod].clear();
            // HGAGenome::initSolomon(iDay, this->m_route[vod], cloneArr, this->m_data[vod]);
        }
        // All routes have been initialized, so we should check if any customer left
        HGAGenome::PRheuristic(this->m_route, this->m_data, cloneArr, iDay, true);
        // finish routing and calculate duration cost
        for (iVeh = 0; iVeh < HPGV::mVeh; iVeh++){
            vod = iDay * HPGV::mVeh + iVeh;
            HGAGenome::delayDeparture(this->m_route[vod], this->m_data[vod]);
        }
    }
}
/**
 * update violations after initializing
 */
void HGAGenome::updateTotalVio(void){
    unsigned int iDay = 0;
    unsigned int iVeh = 0;
    unsigned int vod = 0;
    // reset all value
    this->durationCost = 0;
    this->totalCapacityVio = 0;
    this->totalDurationVio = 0;
    this->totalTimeVio = 0;

    for (iDay = 0; iDay < HPGV::tDay; iDay++){
        for (iVeh = 0; iVeh < HPGV::mVeh; iVeh++){
            vod = iDay * HPGV::mVeh + iVeh;
            // total travel cost
            this->durationCost += this->m_data[vod]->cost;
            // total violation of duration
            if (this->m_data[vod]->cost > HPGV::maxDuration){
                this->totalDurationVio += (this->m_data[vod]->cost - HPGV::maxDuration);
            }
            // total violation of capacity
            if (this->m_data[vod]->load > HPGV::maxLoad){
                this->totalCapacityVio += (this->m_data[vod]->load - HPGV::maxLoad);
            }
            // total violation of time windows
            this->totalTimeVio += this->m_data[vod]->timeVio;
        }
    }
    this->isFeasible = (this->totalCapacityVio == 0) && (this->totalDurationVio == 0) && (this->totalTimeVio == 0);
}
