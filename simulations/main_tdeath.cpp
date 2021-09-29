#include <iostream>
#include "multifrailty.hpp"
int main(int argc, char **argv) {

    uint N = 5000;
    Multifrailty mf(N);
    
    int t=0,tmax=2000;
    double op1,op2;
    double t_short=1,t_long=1;
    double kavg=atof(argv[1]);
    double dt=.01;
    double disable_rate = 0.5;
    double pc_theory = 1.0 / kavg;
    double pc_theory_id = 2.4554 / kavg;
    double tdeath_theory = log( pc_theory ) / -disable_rate;
    double tdeath_theory_id = log( pc_theory_id ) / (-disable_rate*disable_rate);
    
    t_short = tdeath_theory*0.05;
    t_long=tdeath_theory_id*atof(argv[2]);;
    tmax = tdeath_theory_id * 20;
    mf.set_topology(NetworkTopology::ER,NetworkTopology::ER);
    mf.set_dynamics_type(DynamicsType::Percolation,DynamicsType::Percolation);
    mf.set_kavg(kavg,kavg);
    mf.make_nets();

    mf.set_include_neighbors_in_neighborhood_functionality(0);
    mf.set_dt(dt);
    mf.do_secondary_disables(0);
    mf.do_secondary_disables(1);
    mf.assess_total_functionality();
    auto func = mf.get_total_functionality();
    //std::cout << func[0] << "\t" << func[1] << std::endl;
    auto pf = mf.get_primary_fails();
    mf.set_disable_rate(disable_rate,disable_rate);
    mf.set_short_restore_time(t_short,t_short);
    mf.set_long_restore_time(t_long,t_long);
    while(mf.get_t() < tmax){
        
        mf.evolve_one_step();
        func = mf.get_total_functionality();
        //pf = mf.get_primary_fails();
        
        //std::cout <<  pf[0] << "\t" << func[0] << std::endl;
        if(func[0]<0.01*N)
            break;
        if(mf.get_t()>=tmax){
            std::cerr << "Made it to tmax...\n";
        }
    }
    std::cout << mf.get_t() << std::endl;
//     mf.assess_total_functionality();
//     func = mf.get_total_functionality();
//     pf = mf.get_primary_fails();
    //std::cout << mf.get_t() << "\t" << tdeath_theory << std::endl;
    
    //<<":\t" << func[0] << "\t" << func[1] << std::endl;
    //std::cout << mf.get_status_as_int()
    return 0;
}
