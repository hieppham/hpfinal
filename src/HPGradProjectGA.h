/*
 * HPGradProjectGA.h
 *
 *  Created on: Jun 29, 2012
 *      Author: hieppham
 */

#ifndef HPGRADPROJECTGA_H_
#define HPGRADPROJECTGA_H_

#include "hpfinal.h"
/**
 * Define my own GA here
 */
class HPGradProjectGA : public GASteadyStateGA {
public:
    GADefineIdentity("HPGradProjectGA", 246);
    HPGradProjectGA(const GAGenome& g) : GASteadyStateGA(g) {}
    virtual ~HPGradProjectGA() {}
    virtual void step();
    HPGradProjectGA & operator++() { step(); return *this; }
};

/**
 * Customized selector
 */
class HPGradSelector : public GASelectionScheme{
public:
    GADefineIdentity("HPGradSelector", 273);
    HPGradSelector(int w=GASelectionScheme::SCALED) : GASelectionScheme(w){}
    HPGradSelector(const HPGradSelector& orig) { copy(orig); }
    HPGradSelector& operator=(const HPGradSelector& orig){
        if(&orig != this) copy(orig); return *this;
    }
    virtual ~HPGradSelector() { }
    virtual GASelectionScheme* clone() const { return new HPGradSelector; }
    virtual GAGenome& select() const;
};

/**
 * Customized scaling scheme
 */
class HPGradScaling : public GAScalingScheme{
public:
    GADefineIdentity("HPGradScaling", 286);
    HPGradScaling(){}
    HPGradScaling(const HPGradScaling & arg){copy(arg);}
    HPGradScaling & operator=(const GAScalingScheme & arg){
        copy(arg);
        return *this;
    }
    virtual ~HPGradScaling(){}
    virtual GAScalingScheme * clone() const{
        return new HPGradScaling(*this);
    }
    virtual void evaluate(const GAPopulation & p);
    virtual void copy(const GAScalingScheme & arg){
        if(&arg != this && sameClass(arg)){
            GAScalingScheme::copy(arg);
        }
    }
};
#endif /* HPGRADPROJECTGA_H_ */
