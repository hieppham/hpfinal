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
#endif /* HPGRADPROJECTGA_H_ */
