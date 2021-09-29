#include <iostream>
#include "multifrailty.hpp"

int main(int argc, char **argv) {

    if(argc < 5){
        std::cout << "usage:\n"<<argv[0]<<" N restore_rate damage_fraction tmax impacting_factor capacity\n";
        return 1;
    }
    uint N = atoi(argv[1]);
    double restore_rate = atof(argv[2]);
    double damage_fraction_input = atof(argv[3]);
    double tmax = atof(argv[4]);
    double impacting_factor_max = atof(argv[5]);
    double capacity = atof(argv[6]);
    Multifrailty mf(N);
    mf.set_topology(NetworkTopology::ER,NetworkTopology::ER);
    mf.set_kavg(3,3);
    mf.make_nets();
    double dt=1;
    mf.set_dt(dt);
    mf.do_secondary_disables(0);
    mf.do_secondary_disables(1);
    mf.assess_total_functionality();
    auto func = mf.get_total_functionality();
    //std::cout << func[0] << "\t" << func[1] << std::endl;
    auto initial_func = func;
    
    int t=0;
    double op1,op2;
    std::vector < std::vector<double> > history;
    mf.set_disable_rate(0);
    mf.set_restore_rate(restore_rate);
    mf.set_impacting_factor(impacting_factor_max);
    mf.set_dynamics_type(DynamicsType::Percolation,DynamicsType::Percolation);
    double restore_start = 1;
    int restore_nsteps = 6;
    double restore_step_size = restore_start / restore_nsteps;
    double impacting_factor_step=0.499999999;
    //std::cout <<"{";
    bool isfirst=true;
    
    double damage_rate = 0.05*restore_rate;
    //for(double restore_rate = 0; restore_rate<=0.8; restore_rate+=restore_step_size ){
    mf.set_do_enable_disable_together(1);
    if(capacity <= 0){
        mf.set_use_repair_capacity(0);
    } else {
        mf.set_use_repair_capacity(1);
        mf.set_repair_capacity(capacity*N);
    }
    
    
        mf.restore_all();
        
        
        mf.set_disable_rate(damage_rate);
        mf.set_restore_rate(restore_rate);
        mf.evolve(tmax/2);
        auto init_func = mf.get_total_functionality();
    
        auto state_0 = mf.get_status_as_int(0);
        auto state_1 = mf.get_status_as_int(1);
        
        if (init_func[0] < 0.01 * N){
            std::cerr << "Spontaneous collapse, no critical transition\n";
            return 0;
        }
        
        auto net0func = mf.do_random_attack(0,1,tmax,init_func[0]*0.9,0.01*N);
        
        if(net0func > 0.01 * N){
            std::cerr << "Collapse impossible, no critical transition\n";
            return 0;
        }
        auto longest_history = mf.get_functionality_history();
        double damage_fraction_max = 1;
        double damage_fraction_min = 0;
        double damage_fraction = (damage_fraction_max + damage_fraction_min) / 2;
        double precision = 1e-5;
        while (damage_fraction_max - damage_fraction_min > precision){
            mf.set_status_from_int_vector(0,state_0);
            mf.set_status_from_int_vector(1,state_1);
            
            damage_fraction = (damage_fraction_max + damage_fraction_min) / 2;
            std::cerr << damage_fraction << "\t";
            net0func = mf.do_random_attack(0,damage_fraction,tmax,init_func[0]*0.9,0.01*N);
            std::cerr << net0func << "\t" << mf.get_functionality_history().size() << std::endl;
            if(net0func > 0.01 * N){
                damage_fraction_min = damage_fraction;
                longest_history = mf.get_functionality_history();
                
            } else {
                damage_fraction_max = damage_fraction;
            }
        

        }        
      //  if(!isfirst)
        //    std::cout<<",\n";
        //std::cout << "\"" <<damage_fraction << "\": ";
        //jsonArrayofArrays(mf.get_functionality_history(),std::cout);
        //std::cout << restore_rate << "\t" << func[0] << "\t" << func[1] << std::endl;
        //isfirst=false;
        std::cout << damage_fraction_min;
        std::cout <<  "\t" << init_func[0] << "\t" <<init_func[1];
        for(auto f1f2 : longest_history){
            std::cout << "\t" << f1f2[0] << "\t" <<f1f2[1];
        }
                std::cout <<std::endl;

        
        //std::cout<<"}\n";
        

    return 0;
}

