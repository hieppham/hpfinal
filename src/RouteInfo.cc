/*
 * RouteInfo.cpp
 *
 *  Created on: Jun 23, 2012
 *      Author: hieppham
 */

#include "RouteInfo.h"

RouteInfo::RouteInfo() {
    cost = 0;
    load = 0;
    timeVio = 0;
    timeLeaveDepot = 0;

}

RouteInfo::~RouteInfo() {
}

void RouteInfo::resetAll(void){
    this->cost = 0;
    this->load = 0;
    this->timeLeaveDepot = 0;
    this->timeVio = 0;
    this->FTS0 = 0;
}
