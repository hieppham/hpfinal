/*
 * HGAGenome.h
 *
 *  Created on: Jun 24, 2012
 *      Author: hieppham
 */

#ifndef HGAGENOME_H_
#define HGAGENOME_H_

#include <ga/ga.h>
#include <ga/std_stream.h>
#include <boost/shared_ptr.hpp>
#include "hpfinal.h"

class CustomerInDay{
public:
    short cid;
    short vod;
public:
    CustomerInDay(int, unsigned int);
    ~CustomerInDay();
};

struct ltcid
{
  bool operator()(const CustomerInDay* c1, const CustomerInDay* c2) const
  {
      double tmp = c1->vod - c2->vod;
      if (tmp != 0){
          return (tmp < 0);
      }else{
          return ((c1->cid - c2->cid) < 0);
      }
  }
};

typedef boost::shared_ptr<CustomerInDay> CidPtr;
typedef std::vector<CidPtr> cusinday;

typedef boost::shared_ptr<RouteInfo> RinfoPtr;
typedef std::vector<RinfoPtr> RouteData;

typedef boost::shared_ptr<Vertex> VertexPtr;
typedef std::list<VertexPtr> Route;

class HGAGenome: public GAGenome {
public:
    GADefineIdentity("HGAGenome", 251);

    static void Initializer(GAGenome&);
    static int Mutator(GAGenome&, float);
    static double Comparator(const GAGenome&, const GAGenome&);
    static float Evaluator(GAGenome&);

    static int Crossover(const GAGenome&, const GAGenome&, GAGenome*, GAGenome*);
    static int explorationCrossover(const HGAGenome&, const HGAGenome&, HGAGenome&);
    static int exploitationCrossover(const HGAGenome&, const HGAGenome&, HGAGenome&);
    static int Education(const GAGenome&, const int);
public:
    void clusterFirstInit(vector<Customer*>&);
    void SolomonTONNInit(vector<Customer*>&);
public:
    HGAGenome(int);
    HGAGenome(const HGAGenome & orig) {
        copy(orig);
    }
    HGAGenome operator=(const GAGenome & arg) {
        copy(arg);
        return *this;
    }
    virtual ~HGAGenome();
    virtual GAGenome *clone(GAGenome::CloneMethod) const;
    virtual void copy(const GAGenome & c);
    // virtual int equal(const GAGenome& g) const;
    // virtual int read(istream & is);
    virtual int write(ostream & os) const;
public:
    bool isFeasible;
    double durationCost;
    double totalTimeVio;
    double totalCapacityVio;
    double totalDurationVio;

    RouteData m_data;
    vector<Customer> arrC;
    vector<Route> m_route;
    vector<int> m_pattern;
    cusinday m_tour;
};

#endif /* HGAGENOME_H_ */
