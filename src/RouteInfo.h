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
    RouteInfo(const RouteInfo & orig){
        copy(orig);
    }
    virtual ~RouteInfo();
    virtual void copy(const RouteInfo & orig);
    RouteInfo operator=(const RouteInfo & arg) {
        copy(arg);
        return *this;
    }
    void resetAll(void);
};

#endif /* ROUTEINFO_H_ */
