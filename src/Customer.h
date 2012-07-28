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
    unsigned int id; // customer number
    double x;  // x coordinate
    double y;  // y coordinate
    unsigned int d;  // service duration
    unsigned int q;  // demand
    unsigned int f;  // frequency of visit
    unsigned int a;  // number of possible visit combinations
    vector<unsigned int> comb;
    unsigned int e;  // beginning of time window (earliest time for start of service)
    unsigned int l;  // end of time window (latest time for start of service)
    unsigned int pattern; // current pattern assigned
    unsigned int token; // number of visits until now
    bool isServiced; // this flag is used to mark customer as serviced or not

public:
    Customer();
    Customer(const Customer & orig);
    Customer(unsigned int, double, double, unsigned int, unsigned int, unsigned int, unsigned int);
    void setTime(unsigned int te, unsigned int tl);
    void checkServiced();
    void randomAssignedPattern(unsigned int);
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
    Vertex(Customer*&);
    ~Vertex();
};
#endif /* CUSTOMER_H_ */
