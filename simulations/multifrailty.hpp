#include <vector>
#include <map>
#include <queue>
#include <utility>
#include <list>
#include <random>
#include <iterator>
#include <algorithm>
#include <unistd.h>
#include <chrono>
#include "jsonutils.hpp"

using uint = unsigned int;
using std::vector;
#pragma once

#ifndef MULTIFRAILTY_H
#define MULTIFRAILTY_H

struct edge {
  uint target;
  uint status;
  edge(uint t):target(t){}
};

struct node {
  uint idx;
  double restoration_time;
  node(uint idx, double t):idx(idx),restoration_time(t){}
    
};

enum NodeStatus {Operational, PrimaryFail, SecondaryFail, PermanentFail};
enum NetworkTopology {Lattice,ER, SF};
enum DynamicsType {Percolation, Watts};
enum RestorationMode {RestoreRate,RestoreTime};
enum DamageProcess {Constant, Multiplicative, Stochastic, StochasticCorrelated};
class PriorityQueueNodeCompare
{
public:
    bool operator() (node a, node b)
    {
        return a.restoration_time > b.restoration_time;
    }
};

struct beta_distribution{
        std::gamma_distribution<> X;
        std::gamma_distribution<> Y;
        
        beta_distribution(double alpha, double beta){
            X = std::gamma_distribution<>(alpha,1);
            Y = std::gamma_distribution<>(beta,1);
        }
        double operator() (std::mt19937& gen){
            double x = X(gen);
            double y = Y(gen);
            return x / (x+y);
        }
};

class Multifrailty
{
public:
    Multifrailty();
    Multifrailty(uint N);
    void make_nets();
    
    void do_disabling_round(uint net_idx);
    void do_enabling_round(uint net_idx);
    void do_enabling_round_restore_rate_mode(uint net_idx);
    void do_enabling_round_restore_time_mode(uint net_idx);
    void do_secondary_disables(uint net_idx);
    void do_secondary_disables_percolation(uint net_idx);
    void do_secondary_disables_watts(uint net_idx);
    void disable_node(uint node_idx, uint net_idx);
    void do_enable_disable_round_rate_mode(uint net_idx);
    void restore_all();
    
    void set_node_status(uint net_idx, uint node_idx, NodeStatus new_status,double restore_time=0);
    void assess_functionality(uint nx);
    double assess_neighborhood_functionality(uint nx, uint node_idx,bool include_self=true);
    double determine_restore_time(uint ndet_idx,uint node_idx);
    void assess_total_functionality();
    //vector<int> get_degree_distribution(uint nx){ return get_degree_distribution(adjlist[0]);}
    void evolve(double until_t,bool keep_history=false);
    void evolve_continuous(double until_t);
    void evolve_one_step();    
    void do_local_attack(uint net_idx, uint root_node, uint depth);
    uint do_random_attack(uint net_idx, double damage_fraction, double tmax, double stop_above, double stop_below);

    
    //GETTERS//
    const double get_t(){return t;}
    const NodeStatus get_node_status(uint net_idx, uint node_idx);
    const vector<double> get_total_functionality(){ return is_smoothed_functionality_available? smoothed_functionality : functionality; }
    const vector<vector<double>> get_functionality_history(){ return func_history;}
    const vector<vector<double>> get_repair_history(){ return repair_history;}
    const vector<vector<double>> get_primary_fails_history(){ return primaryfails_history;}
    const vector<vector<double>> get_newprimary_fails_history(){ return newprimaryfails_history;}
    const vector<double> get_primary_fails(){ return primaryfails; }
    double get_disable_rate(uint net_idx);
    const vector<int> get_status_as_int(uint nx);
    const double get_impacted_restore_rate(uint net_idx,uint node_idx);
    const vector<double> get_all_neighborhood_functionality(uint nx);
    const std::vector<std::vector<uint>> get_adjlist(uint nx){return adjlist[nx];};
    
    //SETTERS//
    void set_dt(double dt_){dt=dt_;}
    void set_disable_rate(double rate0,double rate1){disable_rate[0]=rate0;disable_rate[1]=rate1;}
    void set_disable_rate(double same_rate){set_disable_rate(same_rate,same_rate);}
    void set_disable_behavior(std::string argv_input);
    void set_restore_rate(double rate0,double rate1){restore_rate[0]=rate0;restore_rate[1]=rate1;}
    void set_restore_rate(double same_rate){set_restore_rate(same_rate,same_rate);}
    void set_impacting_factor(double factor0,double factor1){impacting_factor[0]=factor0;impacting_factor[1]=factor1;}
    void set_impacting_factor(double same_factor){set_impacting_factor(same_factor,same_factor);}
    void set_restoration_mode(RestorationMode new_mode){restoration_mode = new_mode;}
    void set_sfgamma(double sfgamma0,double sfgamma1){sfgamma[0]=sfgamma0; sfgamma[1]=sfgamma1; }
    void set_kavg(double kavg0,double kavg1){kavg[0]=kavg0; kavg[1]=kavg1; }
    void set_short_restore_time(double t0,double t1){short_restore_time[0]=t0; short_restore_time[1]=t1; }
    void set_long_restore_time(double t0,double t1){long_restore_time[0]=t0; long_restore_time[1]=t1; }
    void set_topology(NetworkTopology topo0,NetworkTopology topo1){topology[0]=topo0; topology[1] = topo1;}
    void set_dynamics_type(DynamicsType type0,DynamicsType type1){dynamics_type[0]=type0; dynamics_type[1] = type1; }
    void set_damage_process_type(DamageProcess type0,DamageProcess type1){damage_process[0]=type0; damage_process[1] = type1; }
    void set_include_self_in_neighborhood_functionality(int yesno){include_self_in_neighborhood_functionality = static_cast<bool>(yesno);}
    void set_include_neighbors_in_neighborhood_functionality(int yesno){include_neighbors_in_neighborhood_functionality = static_cast<bool>(yesno);}
    void set_impervious_rate(double r){impervious_rate = r;}
    void set_do_enable_disable_together(bool x){do_enable_disable_together=x;}
    void set_use_repair_capacity(bool x){use_repair_capacity = x;}
    void set_repair_capacity(double this_repair_capacity){repair_capacity[0]=this_repair_capacity; repair_capacity[1]=this_repair_capacity;}
    void set_status_from_int_vector(uint net_idx, std::vector<int> status);
    void set_stochastic_parameters(std::vector<double> stochastic_parameters){ this->stochastic_parameters = stochastic_parameters;}
    
    

private:
    uint N;
    double t,dt;
    std::vector<std::vector<NodeStatus>> status;
    std::vector<std::priority_queue<node,std::vector<node>,PriorityQueueNodeCompare>> restore_queue;
    std::mt19937 gen;
    std::vector<uint> nets{0,1}; 
    std::vector<double> disable_rate;
    std::vector<double> stochastic_disable_rate;
    std::vector<double> impacting_factor;
    std::vector<double> restore_rate;
    std::vector<double> stochastic_parameters;
    std::vector<std::vector<std::vector<uint>>> adjlist;
    std::vector<int> restore_time;
    std::vector<double> functionality;
    std::vector<std::vector<double>> func_history;
    std::vector<std::vector<double>> repair_history;
    std::vector<std::vector<double>> primaryfails_history;
    std::vector<std::vector<double>> newprimaryfails_history;
    std::vector<double> primaryfails;
    std::vector<double> repaired_this_round;
    std::vector<double> disabled_this_round;
    std::vector<double> repaired_last_round;
    std::vector<double> repair_capacity;
    std::vector<double> sfgamma;
    std::vector<double> kavg;
    std::vector<double> short_restore_time;
    std::vector<double> long_restore_time;
    std::vector<NetworkTopology> topology;
    std::vector<DynamicsType> dynamics_type;
    std::vector<DamageProcess> damage_process;
    RestorationMode restoration_mode;
    std::vector<double> smoothed_functionality;
    bool is_smoothed_functionality_available;
    double impervious_rate;
    bool print_degree_distribution,include_neighbors_in_neighborhood_functionality,include_self_in_neighborhood_functionality,do_enable_disable_together,use_repair_capacity;
};




template <typename Stream> void write_mf_history_json( Multifrailty  mf, Stream & stream){
        stream <<"{";

    //  if(!isfirst)
    //    stream<<",\n";
    stream << "\"functionality_history\" : " << std::endl;
    jsonArrayofArrays(mf.get_functionality_history(),stream);
    stream << ",\n\"repair_history\" :"<< std::endl;
    jsonArrayofArrays(mf.get_repair_history(),stream);
    stream << ",\n\"primary_failures_history\" :"<< std::endl;
    jsonArrayofArrays(mf.get_primary_fails_history(),stream);
    stream << ",\n\"newprimary_failures_history\" :"<< std::endl;
    jsonArrayofArrays(mf.get_newprimary_fails_history(),stream);        //stream << restore_rate << "\t" << func[0] << "\t" << func[1] << std::endl;
    //isfirst=false;
    
    //for(auto f1f2 : mf.get_functionality_history()){
    //    stream << "\t" << f1f2[0] << "\t" <<f1f2[1];
    //}
    //        stream <<std::endl;

    //}
    stream<<"}\n";
}

template <typename Stream> void write_mf_history_csv( Multifrailty  mf, Stream & stream){
    stream << "functionality_history_0,";
    stream << "functionality_history_1,";
    stream << "repair_history_0,";
    stream << "repair_history_1,";
    stream << "primary_failures_history_0,";
    stream << "primary_failures_history_1,";
    stream << "newprimary_failures_history_0,";
    stream << "newprimary_failures_history_1";
    auto fv = mf.get_functionality_history();
    auto rv = mf.get_repair_history();
    auto pv = mf.get_primary_fails_history();
    auto npv = mf.get_newprimary_fails_history(); 
    
    auto M = fv.size();
    for (auto i=M*0; i<M; i++){
        stream << "\n";
        stream << fv[i][0] << ",";
        stream << fv[i][1] << ",";
        stream << rv[i][0] << ",";
        stream << rv[i][1] << ",";
        stream << pv[i][0] << ",";
        stream << pv[i][1] << ",";
        stream << npv[i][0] << ",";
        stream << npv[i][1] ;
    }
}

        




#endif // MULTIFRAILTY_H
