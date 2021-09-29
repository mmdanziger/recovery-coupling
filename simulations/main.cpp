#include <iostream>
#include "multifrailty.hpp"
int main(int argc, char **argv) {

    if(argc < 5){
        std::cout << "usage:\n"<<argv[0]<<" N k disable_rate tmax impacting_factor\n";
        return 1;
    }
    uint N = atoi(argv[1]);
    double k = atof(argv[2]);
    double disable_rate = atof(argv[3]);
    double tmax = atof(argv[4]);
    double impacting_factor = atof(argv[5]);
    
    Multifrailty mf(N);
    mf.set_topology(NetworkTopology::ER,NetworkTopology::ER);
    mf.set_kavg(k,k);
    mf.make_nets();
    mf.set_include_neighbors_in_neighborhood_functionality(1);
    mf.set_include_self_in_neighborhood_functionality(0);
    mf.set_use_repair_capacity(1);
    mf.set_repair_capacity(.1*N);
    mf.set_dt(1);
    mf.do_secondary_disables(0);
    mf.do_secondary_disables(1);
    mf.assess_total_functionality();
    auto func = mf.get_total_functionality();
    auto pfails = mf.get_primary_fails();

    //std::cout << func[0] << "\t" << func[1] << std::endl;
    
    int t=0;
    double op1,op2;
    std::vector < std::vector<double> > history,pfailhistory;
    double restore_rate=0.02;
    mf.set_disable_rate(disable_rate);
    mf.set_restore_rate(restore_rate);
    mf.set_impacting_factor(impacting_factor);
    mf.set_dynamics_type(DynamicsType::Watts,DynamicsType::Watts);
    double restore_start = 0.8;
    double restore_min = 0.01;
    int restore_nsteps = 1000;
    double restore_step_size = (restore_start - restore_min) / restore_nsteps;
    mf.restore_all();
    mf.set_do_enable_disable_together(true);
    restore_rate=restore_start;
    restore_step_size*=-1;
    while(true){
        std::cout<< mf.get_disable_rate(0) << std::endl;
    }
    while(true){
    
        mf.set_disable_rate(disable_rate);
        mf.set_restore_rate(restore_rate);
        
        mf.evolve(tmax,true);
        mf.assess_total_functionality();
        func = mf.get_total_functionality();
        pfails = mf.get_primary_fails();

        history.push_back(func);
        pfailhistory.push_back(pfails);
        //1/(1 - restore_rate*(1 - 1/disable_rate))
        std::cout << restore_rate << "\t" << disable_rate << "\t" << func[0] << "\t" << pfails[0] << std::endl;
        restore_rate += restore_step_size;
        if(restore_rate < restore_min){
            restore_step_size*=-1;
            restore_rate += restore_step_size;
        }
        if(restore_rate > restore_start){
            break;
        }
        
    }
    
        

    return 0;
}
