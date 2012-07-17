/*
 * RouteInfo.h
 *
 *  Created on: Jun 23, 2012
 *      Author: hieppham
 */

#include <list>
#ifndef ROUTEINFO_H_
#define ROUTEINFO_H_

class RouteInfo {
public:
    double cost;
    double load;
    double timeVio;
    double timeLeaveDepot;
    double FTS0;
public:
    RouteInfo();
    virtual ~RouteInfo();
    void resetAll(void);
};

#endif /* ROUTEINFO_H_ */
