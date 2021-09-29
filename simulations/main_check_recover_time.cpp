#include <iostream>
#include "multifrailty.hpp"

int main(int argc, char **argv) {

    if(argc < 5){
        std::cout << "usage:\n"<<argv[0]<<" N restore_rate damage_fraction tmax impacting_factor CSM\n";
        return 1;
    }
    uint N = atoi(argv[1]);
    double restore_rate = atof(argv[2]);
    double damage_fraction_max = atof(argv[3]);
    double tmax = atof(argv[4]);
    double impacting_factor_max = atof(argv[5]);
    double capacity = 0;
    
    Multifrailty mf(N);
    mf.set_topology(NetworkTopology::ER,NetworkTopology::ER);
    mf.set_kavg(5,5);
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
    double damage_fraction = 0;
    double damage_fraction_min = 0.63;
    double damage_fraction_step = (damage_fraction_max - damage_fraction_min)/ restore_nsteps;
    double damage_rate = 0.25*restore_rate;
    //for(double restore_rate = 0; restore_rate<=0.8; restore_rate+=restore_step_size ){
    mf.set_do_enable_disable_together(1);
    if( (std::string(argv[6]) == "C") )
        mf.set_damage_process_type(DamageProcess::Constant,DamageProcess::Constant);
    if( (std::string(argv[6]) == "M") )
        mf.set_damage_process_type(DamageProcess::Multiplicative,DamageProcess::Multiplicative);
    if( (std::string(argv[6]) == "S") )
        mf.set_damage_process_type(DamageProcess::Stochastic,DamageProcess::Stochastic);
    if( (std::string(argv[6]) == "SC") )
        mf.set_damage_process_type(DamageProcess::StochasticCorrelated,DamageProcess::StochasticCorrelated);
        
    
    if(capacity <= 0){
        mf.set_use_repair_capacity(0);
    } else {
        mf.set_use_repair_capacity(1);
        mf.set_repair_capacity(capacity*N);
    }
    
    
      for(damage_fraction = damage_fraction_min; damage_fraction<damage_fraction_max; damage_fraction+=damage_fraction_step){
        mf.restore_all();
        
        
        mf.set_disable_rate(damage_rate);
        mf.set_restore_rate(restore_rate);
        mf.evolve(tmax/2);
        func = mf.get_total_functionality();
        
        
        mf.set_disable_rate(damage_fraction);
        mf.set_restore_rate(0);
        mf.set_dt(1);
        mf.evolve_one_step();
        
        mf.set_dt(dt);
        mf.set_disable_rate(damage_rate);
        mf.set_restore_rate(restore_rate);
        mf.evolve(tmax,true);
        
        

        func = mf.get_total_functionality();
      //  if(!isfirst)
        //    std::cout<<",\n";
        //std::cout << "\"" <<damage_fraction << "\": ";
        //jsonArrayofArrays(mf.get_functionality_history(),std::cout);
        //std::cout << restore_rate << "\t" << func[0] << "\t" << func[1] << std::endl;
        //isfirst=false;
        std::cout << damage_fraction;
        for(auto f1f2 : mf.get_functionality_history()){
            std::cout << "\t" << f1f2[0] << "\t" <<f1f2[1];
        }
                std::cout <<std::endl;

        }
        //std::cout<<"}\n";
        

    return 0;
}
