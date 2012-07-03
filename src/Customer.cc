/*
 * Customer.cc
 *
 *  Created on: Mar 26, 2012
 *      Author: hieppham
 */

#include "Customer.h"

Customer::Customer(){
    this->token = 0;
    this->isServiced = false;
}
Customer::Customer(int tid, double tx, double ty, int td, int tq, int tf, int ta) {
    // TODO Auto-generated constructor stub
    this->id = tid;
    this->x = tx;
    this->y = ty;
    this->d = td;
    this->q = tq;
    this->f = tf;
    this->a = ta;

    this->token = 0;
    this->pattern = 0;
    this->isServiced = false;

}
void Customer::checkServiced(){
    this->token++;
    if (this->token >= this->f){
        this->isServiced = true;
    }
}
void Customer::setTime(int te, int tl){
    this->e = te;
    this->l = tl;
}
void Customer::randomAssignedPattern(int idx){
    this->pattern = this->comb[idx];
}
Customer::~Customer() {
    // TODO Auto-generated destructor stub
}

Vertex::Vertex(Customer* c){
    this->cus = c;
    this->timeArrive = 0;
    this->timeStartService = 0;
    this->timeWait = 0;
    this->newTimeArrive = 0;
    this->newTimeStart = 0;

    this->pf = 0;
    this->FTS = 0;
}
Vertex::~Vertex(){
}
