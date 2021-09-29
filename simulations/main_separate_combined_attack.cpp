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
    
    mf.set_disable_behavior(std::string(argv[6]));
    
    
    double restore_start = 1;
    int restore_nsteps = 6;
    double restore_step_size = restore_start / restore_nsteps;
    double impacting_factor_step=0.499999999;
    bool isfirst=true;
    double damage_fraction = 0;
    double damage_fraction_min = 0.63;
    double damage_fraction_step = (damage_fraction_max - damage_fraction_min)/ restore_nsteps;
    double damage_rate = 0.25*restore_rate;
    //for(double restore_rate = 0; restore_rate<=0.8; restore_rate+=restore_step_size ){
    mf.set_do_enable_disable_together(1);
    if(capacity <= 0){
        mf.set_use_repair_capacity(0);
    } else {
        mf.set_use_repair_capacity(1);
        mf.set_repair_capacity(capacity*N);
    }
    
    
      //for(damage_fraction = damage_fraction_min; damage_fraction<damage_fraction_max; damage_fraction+=damage_fraction_step){
        damage_fraction = damage_fraction_max;
        mf.restore_all();
        
        
        mf.set_disable_rate(damage_rate);
        mf.set_restore_rate(restore_rate);
        mf.evolve(10*tmax);
        func = mf.get_total_functionality();
        
        
        //Attack Network 1
        mf.set_disable_rate(damage_fraction,damage_rate);
        mf.set_restore_rate(0);
        mf.set_dt(1);
        mf.evolve_one_step();

        //Recover
        mf.set_dt(dt);
        mf.set_disable_rate(damage_rate);
        mf.set_restore_rate(restore_rate);
        mf.evolve(tmax,true);
        
        //Attack Network 2
        mf.set_disable_rate(damage_rate,damage_fraction);
        mf.set_restore_rate(0);
        mf.set_dt(1);
        mf.evolve_one_step();

        //Recover
        mf.set_dt(dt);
        mf.set_disable_rate(damage_rate);
        mf.set_restore_rate(restore_rate);
        mf.evolve(tmax,true);
        
        //Attack Both
        mf.set_disable_rate(damage_fraction);
        mf.set_restore_rate(0);
        mf.set_dt(1);
        mf.evolve_one_step();

        //Recover
        mf.set_dt(dt);
        mf.set_disable_rate(damage_rate);
        mf.set_restore_rate(restore_rate);
        mf.evolve(tmax,true);
        
                

        func = mf.get_total_functionality();
        std::cout <<"{";

        //  if(!isfirst)
        //    std::cout<<",\n";
        std::cout << "\"damage_fraction\" : " << damage_fraction << "," << std::endl;
        std::cout << "\"functionality_history\" : " << std::endl;
        jsonArrayofArrays(mf.get_functionality_history(),std::cout);
        std::cout << ",\n\"repair_history\" :"<< std::endl;
        jsonArrayofArrays(mf.get_repair_history(),std::cout);
        std::cout << ",\n\"primary_failures_history\" :"<< std::endl;
        jsonArrayofArrays(mf.get_primary_fails_history(),std::cout);
        std::cout << ",\n\"newprimary_failures_history\" :"<< std::endl;
        jsonArrayofArrays(mf.get_newprimary_fails_history(),std::cout);        //std::cout << restore_rate << "\t" << func[0] << "\t" << func[1] << std::endl;
        //isfirst=false;
        
        //for(auto f1f2 : mf.get_functionality_history()){
        //    std::cout << "\t" << f1f2[0] << "\t" <<f1f2[1];
        //}
        //        std::cout <<std::endl;

        //}
        std::cout<<"}\n";
        

    return 0;
}
