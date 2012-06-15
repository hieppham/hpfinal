/*
 * Customer.h
 *
 *  Created on: Mar 26, 2012
 *      Author: hieppham
 */
#ifndef CUSTOMER_H_
#define CUSTOMER_H_

#include <iostream>
#include <vector>
using namespace std;

class Customer {
public:
    int id; // customer number
    double x;  // x coordinate
    double y;  // y coordinate
    int d;  // service duration
    int q;  // demand
    int f;  // frequency of visit
    int a;  // number of possible visit combinations
    vector<int> comb;
    int e;  // beginning of time window (earliest time for start of service)
    int l;  // end of time window (latest time for start of service)
    int pattern; // current pattern assigned
    int token; // number of visits until now
    bool isServiced; // this flag is used to mark customer as serviced or not

public:
    Customer();
    Customer(int, double, double, int, int, int, int);
    void setTime(int te, int tl);
    void checkServiced();
    void randomAssignedPattern(int);
    ~Customer();
};

class Vertex{
public:
    Customer* cus;
    double timeArrive;
    double timeWait;
    double timeStartService;
    double timeDeparture;
    // temporary values for calculating insertion cost
    double newTimeArrive;
    double newTimeStart;

    // push forward
    double pf;
    double FTS;
public:
    Vertex(Customer*);
    ~Vertex();
};
#endif /* CUSTOMER_H_ */
