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
      short tmp = c1->vod - c2->vod;
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

typedef vector<Customer*> VCus;

typedef map<unsigned int, int> TabuMap;

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
    static int Education(GAGenome&, const int);

    void resetAll();
public:
    void clusterFirstInit(VCus&);
    void SolomonTONNInit(VCus&);
    void PRInit(VCus&);
    static void initCluster(unsigned int&, Route&, VCus&, RinfoPtr&);
    static void initSolomon(unsigned int&, Route&, VCus&, RinfoPtr&);

    static void pushbackRoute(Route&, RinfoPtr&, Customer*&);
    static void insertIntoRoute(Route&, RinfoPtr&, Customer*&, unsigned int&);
    static void PRinsert(Route&, RinfoPtr&, Customer*&);
    static void removeFromRoute(Route&, RinfoPtr&, unsigned int);
    static void updateInfo(Route&, RinfoPtr&);
    static void SolomonI1(Route&, RinfoPtr&, VCus&);
    static void PRheuristic(vector<Route>&, RouteData&, VCus&, unsigned int&, bool);

    static void testRoute(Route&);
    static bool isInRoute(Route&, unsigned int);
    void updateTotalVio(void);
    void tourConstruct(void);

    void modifyEachCustomer(void);
    static void delayDeparture(Route&, RinfoPtr&);
    static double calcObjectValue(HGAGenome&);
public:
    static HGAGenome UTS(HGAGenome&);
    static bool UTSNeighborByPattern(HGAGenome&, TabuMap&, vector<vector<unsigned int> >&, int&, int&, double&, double&, double&);
    static bool UTSNeighborByRouting(HGAGenome&, TabuMap&, vector<vector<unsigned int> >&, int&, int&, double&, double&, double&);
    void tourUpdate(vector<vector<unsigned int> >&);

    static HGAGenome RVNS(HGAGenome&);
    static HGAGenome Shaking(HGAGenome&, unsigned int, double);
    static HGAGenome ShakingPattern(HGAGenome&, unsigned int, double);
    static HGAGenome ShakingMoveSegment(HGAGenome&, unsigned int, double);
    static HGAGenome ShakingExchangeSegments(HGAGenome&, unsigned int, double);

    static void improveRoute(HGAGenome&);

    static void apply2OptUntilFirstImprovement(double, double, double, HGAGenome&);
    static void apply2OptForAllRoutes(double, double, double, HGAGenome&);
    static void apply2OptStarForAllRoutes(double, double, double, HGAGenome&);
    static void apply2OptStarUntilFirstImprovement(double, double, double, HGAGenome&);

    bool inter2OptStar(double, double, double, unsigned int, unsigned int);
    bool interCrossExchange(double, double, double, unsigned int, unsigned int);
    bool interRouteOpt(double, double, double, unsigned int, unsigned int);

    bool stoIntraOrOpt(double, double, double, unsigned int);
    bool stoIntra2Opt(double, double, double, unsigned int);

    bool intraOrOpt(double, double, double, unsigned int);
    bool intra2Opt(double, double, double, unsigned int);
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
    static void printSolution(HGAGenome&, char*);
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

class cost{
public:
    double aQ;
    double bD;
    double cW;

    double maxLoad;
    double maxDuration;
public:
    cost(double a, double b, double c, double mL, double mD): aQ(a), bD(b), cW(c), maxLoad(mL), maxDuration(mD){}
    double operator()(Route& r, RinfoPtr& rInfo){
        double costVal = 0;
        Route::iterator uIter, prevIter;
        // consider duplication of calculating duration here.

        costVal = rInfo->cost;
        if (rInfo->load > maxLoad){
            costVal += (rInfo->load - maxLoad) * aQ;
        }
        if (rInfo->cost > maxDuration){
            costVal += (rInfo->cost - maxDuration) * bD;
        }
        if (rInfo->timeVio > 0){
            costVal += rInfo->timeVio * cW;
        }
        return costVal;
    }
};
#endif /* HGAGENOME_H_ */
