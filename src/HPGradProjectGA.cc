/*
 * HPGradProjectGA.cc
 *
 *  Created on: Jun 29, 2012
 *      Author: hieppham
 */

#include "HPGradProjectGA.h"

extern HGAGenome bestSol;

void
HPGradProjectGA::step()
{
    int i, mut, c1, c2;
    GAGenome *mom, *dad;          // tmp holders for selected genomes
    int CNG = stats.generation();
    HPGV::genCounter = CNG;
    cout << CNG << endl;
    // Generate the individuals in the temporary population from individuals in
    // the main population.

    for(i = 0; i < (int) HPGV::nPop; i += 2){   // takes care of odd population
        mom = &(pop->select());
        dad = &(pop->select());
        stats.numsel += 2;      // keep track of number of selections

        c1 = c2 = 0;

        // TODO: remove after debugging
        // cout << HPGV::genCounter << " - costs of parents: " << ((HGAGenome &) *mom).durationCost << " - " << ((HGAGenome &) *dad).durationCost << endl;
        /*
        if(GAFlipCoin(pCrossover())){
            stats.numcro += (*scross)(*mom, *dad, &tmpPop->individual(i),
                    &tmpPop->individual(i+1));
            c1 = c2 = 1;
        }
        else{
            tmpPop->individual( i ).copy(*mom);
            tmpPop->individual(i+1).copy(*dad);
        }
        */
        stats.numcro += (*scross)(*mom, *dad, &tmpPop->individual(i),
                            &tmpPop->individual(i+1));
        c1 = c2 = 1;
        stats.nummut += (mut = tmpPop->individual( i ).mutate(pMutation()));
        if(mut > 0) c1 = 1;
        stats.nummut += (mut = tmpPop->individual(i+1).mutate(pMutation()));
        if(mut > 0) c2 = 1;
        // After crossover and mutation, we educate all offsprings
        HGAGenome::Education(tmpPop->individual( i ), CNG);
        HGAGenome::Education(tmpPop->individual( i+1 ), CNG);

        // TODO: remove after debugging
        // cout << " ## costs of children: " << ((HGAGenome &) tmpPop->individual(i)).durationCost << " - " << ((HGAGenome &) tmpPop->individual(i + 1)).durationCost << endl;

        stats.numeval += c1 + c2;
    }

    for (i = HPGV::nPop; i < tmpPop->size(); i++){
        // mom = &(pop->select());
        // stats.numsel += 1;
        // tmpPop->individual(i).copy(*mom);
        // HGAGenome::UTS((HGAGenome&)tmpPop->individual(i));
        // HGAGenome::improveRoute((HGAGenome&)tmpPop->individual(i));

        HGAGenome::Initializer(tmpPop->individual(i));
        HGAGenome::Education(tmpPop->individual( i ), CNG);

        stats.numeval += 1;
    }

    for(i=0; i<tmpPop->size(); i++)
            pop->remove(GAPopulation::WORST, GAPopulation::SCALED);
    // Replace the worst genomes in the main population with all of the individuals
    // we just created.  Notice that we invoke the population's add member with a
    // genome pointer rather than reference.  This way we don't force a clone of
    // the genome - we just let the population take over.  Then we take it back by
    // doing a remove then a replace in the tmp population.

    for(i=0; i<tmpPop->size(); i++)
        pop->add(tmpPop->individual(i));
    pop->evaluate(gaTrue);      // get info about current pop for next time

    // update penalty parameters
    /*
    HGAGenome & best = (HGAGenome &) pop->best();
    if (best.isFeasible){
        cout << "***********************************************Acceptable with cost = " << best.durationCost << endl;
        if (HPGV::bestFeasibleCost == 0){
            HPGV::bestFeasibleCost = best.durationCost;
            bestSol = best;
        }else{
            if (HPGV::bestFeasibleCost > best.durationCost){
                HPGV::bestFeasibleCost = best.durationCost;
                bestSol = best;
            }
        }
    }
    */
    if (HPGV::bestFeasibleCost != 0){
        HPGV::hPenalty = HPGV::bestFeasibleCost;
    }else{
        HGAGenome & worst = (HGAGenome &) pop->worst();
        HPGV::hPenalty = worst.durationCost;
    }
    // redefine average values
    HPGV::avgQ = 0;
    HPGV::avgD = 0;
    HPGV::avgW = 0;
    for (int i = 0; i < pop->size(); i++){
        HGAGenome& tmpHg = (HGAGenome &) (pop->individual(i));
        cout << HPGV::genCounter << " - " << tmpHg.durationCost << endl;

        HPGV::avgQ += tmpHg.totalCapacityVio;
        HPGV::avgD += tmpHg.totalDurationVio;
        HPGV::avgW += tmpHg.totalTimeVio;
    }
    HPGV::avgQ /= (pop->size());
    HPGV::avgD /= (pop->size());
    HPGV::avgW /= (pop->size());

    pop->scale(gaTrue);         // remind the population to do its scaling

    stats.numrep += tmpPop->size();

//    for(i=0; i<tmpPop->size(); i++)
//        pop->remove(GAPopulation::WORST, GAPopulation::SCALED);

    stats.update(*pop);       // update the statistics by one generation
}

/**
 * Main method of selector
 */
GAGenome& HPGradSelector::select() const {
    return pop->best(GARandomInt(0, pop->size()/2), GAPopulation::SCALED);
}

/**
 * Main scaling method
 */
void HPGradScaling::evaluate(const GAPopulation & p){
    float f;
    for(int i=0; i<p.size(); i++){
        HGAGenome & hg = (HGAGenome &)p.individual(i);
        if (hg.isFeasible){
            f = p.individual(i).score();
        }else{
            f = p.individual(i).score() + HPGV::bestFeasibleCost;
        }
        p.individual(i).fitness(f);
    }
}
