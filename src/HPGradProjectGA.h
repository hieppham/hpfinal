/*
 * HPGradProjectGA.h
 *
 *  Created on: Jun 29, 2012
 *      Author: hieppham
 */

#ifndef HPGRADPROJECTGA_H_
#define HPGRADPROJECTGA_H_

#include "hpfinal.h"

class HPGradProjectGA : public GASteadyStateGA {
public:
    GADefineIdentity("HPGradProjectGA", 246);
    HPGradProjectGA(const GAGenome& g) : GASteadyStateGA(g) {}
    virtual ~HPGradProjectGA() {}
    virtual void step();
    HPGradProjectGA & operator++() { step(); return *this; }
};

#endif /* HPGRADPROJECTGA_H_ */
