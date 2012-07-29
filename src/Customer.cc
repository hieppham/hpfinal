/*
 * Customer.cc
 *
 *  Created on: Mar 26, 2012
 *      Author: hieppham
 */

#include "Customer.h"
#include "hpfinal.h"

Customer::Customer(){
    this->id = 0;
    this->x = 0;
    this->y = 0;
    this->d = 0;
    this->q = 0;
    this->f = 0;
    this->a = 0;

    this->e = 0;
    this->l = 0;

    this->pattern = 0;
    this->token = 0;
    this->isServiced = false;
}

Customer::Customer(const Customer & orig){
    this->id = orig.id;
    this->x = orig.x;
    this->y = orig.y;
    this->d = orig.d;
    this->q = orig.q;
    this->f = orig.f;
    this->a = orig.a;
    this->comb = orig.comb;

    this->e = orig.e;
    this->l = orig.l;

    this->pattern = orig.pattern;
    this->token = numberOfSetBits(this->pattern);
    this->isServiced = orig.isServiced;
}
Customer::Customer(unsigned int tid, double tx, double ty, unsigned int td, unsigned int tq, unsigned int tf, unsigned int ta) {
    this->id = tid;
    this->x = tx;
    this->y = ty;
    this->d = td;
    this->q = tq;
    this->f = tf;
    this->a = ta;

    this->e = 0;
    this->l = 0;

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
void Customer::setTime(unsigned int te, unsigned int tl){
    this->e = te;
    this->l = tl;
}
void Customer::randomAssignedPattern(unsigned int idx){
    this->pattern = this->comb[idx];
}
Customer::~Customer() {
}

Vertex::Vertex(Customer*& c){
    this->cus = new Customer(*c);
    this->timeArrive = 0;
    this->timeStartService = 0;
    this->timeWait = 0;
    this->newTimeArrive = 0;
    this->newTimeStart = 0;
    this->timeDeparture = 0;

    this->pf = 0;
    this->FTS = 0;
}
Vertex::~Vertex(){
    delete cus;
}
