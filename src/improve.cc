#include "HGAGenome.h"

extern Customer* gDepot;
extern vector<Customer> gArrC;
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
    // TODO: 2-Opt for all routes
}

void HGAGenome::apply2OptStarForAllRoutes(double aPen, double bPen, double cPen, HGAGenome& hg){
    // TODO: 2-Opt* for all routes
}
