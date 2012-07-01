/*
 * HPGradProjectGA.cc
 *
 *  Created on: Jun 29, 2012
 *      Author: hieppham
 */

#include "HPGradProjectGA.h"

void
HPGradProjectGA::step()
{
    int i, mut, c1, c2;
    GAGenome *mom, *dad;          // tmp holders for selected genomes
    int CNG = stats.generation();

    // Generate the individuals in the temporary population from individuals in
    // the main population.

    for(i=0; i < (int) HPGV::nPop; i+=2){   // takes care of odd population
        mom = &(pop->select());
        dad = &(pop->select());
        stats.numsel += 2;      // keep track of number of selections

        c1 = c2 = 0;
        if(GAFlipCoin(pCrossover())){
            stats.numcro += (*scross)(*mom, *dad, &tmpPop->individual(i),
                    &tmpPop->individual(i+1));
            c1 = c2 = 1;
        }
        else{
            tmpPop->individual( i ).copy(*mom);
            tmpPop->individual(i+1).copy(*dad);
        }
        stats.nummut += (mut = tmpPop->individual( i ).mutate(pMutation()));
        if(mut > 0) c1 = 1;
        stats.nummut += (mut = tmpPop->individual(i+1).mutate(pMutation()));
        if(mut > 0) c2 = 1;
        // After crossover and mutation, we educate all offsprings
//        HGAGenome::Education(tmpPop->individual( i ), CNG);
//        HGAGenome::Education(tmpPop->individual( i+1 ), CNG);

        stats.numeval += c1 + c2;
    }
    // TODO: fix here!
//    for (i = nPop; i < tmpPop->size(); i++){
//        mom = &(pop->select());
//        stats.numsel += 1;
//        tmpPop->individual(i).copy(*mom);
//        HGAGenome::UTS((HGAGenome &)tmpPop->individual(i));
//        HGAGenome::improveRoute((HGAGenome&)tmpPop->individual(i));
//        stats.numeval += 1;
//    }

    for(i=0; i<tmpPop->size(); i++)
            pop->remove(GAPopulation::WORST, GAPopulation::SCALED);
    // Replace the worst genomes in the main population with all of the individuals
    // we just created.  Notice that we invoke the population's add member with a
    // genome pointer rather than reference.  This way we don't force a clone of
    // the genome - we just let the population take over.  Then we take it back by
    // doing a remove then a replace in the tmp population.

    for(i=0; i<tmpPop->size(); i++)
        pop->add(&tmpPop->individual(i));
    pop->evaluate(gaTrue);      // get info about current pop for next time

    /*
    // update penalty parameters
    HGAGenome & best = (HGAGenome &) pop->best();
    if (best.isFeasible){
        if (HPGV::bestFeasiblecost == 0){
            HPGV::bestFeasibleCost = best.durationCost;
            bestSol = best;
        }else{
            if (HPGV::bestFeasibleCost > best.durationCost){
                HPGV::bestFeasibleCost = best.durationCost;
                bestSol = best;
            }
        }
    }
    if (HPGV::bestFeasibleCost != 0){
        HPGV::hPenalty = HPGV::bestFeasibleCost;
    }else{
        HGAGenome & worst = (HGAGenome &) pop->worst();
        HPGV::hPenalty = worst.getTotalCost();
    }
    // redefine average values
    aveQ = 0;
    aveD = 0;
    aveW = 0;
    for (int i = 0; i < pop->size(); i++){
        HGAGenome& tmpHg = (HGAGenome &) (pop->individual(i));
        aveQ += tmpHg.totalCapacityVio;
        aveD += tmpHg.totalDurationVio;
        aveW += tmpHg.totalTimeVio;
    }
    aveQ /= (pop->size());
    aveD /= (pop->size());
    aveW /= (pop->size());
    */

    pop->scale(gaTrue);         // remind the population to do its scaling

    stats.numrep += tmpPop->size();

    stats.update(*pop);       // update the statistics by one generation
}