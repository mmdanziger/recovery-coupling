#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "multifrailty.hpp"
//#include "jsonutils.hpp"


int main(int argc, char **argv) {

    if(argc < 6){
        std::cout << "usage:\n"<<argv[0]<<" N restore_rate damage_rate tmax impacting_factor damage_process_code shape scale identifier_string (optional)\n";
        return 1;
    }
    srand(time(0));
    uint N = atoi(argv[1]);
    double restore_rate = atof(argv[2]);
    double damage_rate = atof(argv[3]);
    double tmax = atof(argv[4]);
    double impacting_factor = atof(argv[5]);
    auto damage_process_code = std::string(argv[6]);
    double shape = atof(argv[7]);
    double scale = atof(argv[8]);
    bool oneway_dependency = false;
    std::string identifier = argc > 8 ? std::string(argv[9]) : std::to_string(rand()%999999);
    std::ostringstream ofname;
    ofname << "multifrailty_detailed_output_N" << N << "_" << "gr_"<<restore_rate <<"gd_"<<damage_rate <<"tmax_"<<tmax<<"alpha_"<<impacting_factor<< "damage_"<< damage_process_code <<"shape_"<<shape <<"scale_" << scale <<"id_"<<identifier <<".csv";
    Multifrailty mf(N);
    mf.set_topology(NetworkTopology::ER,NetworkTopology::ER);
    mf.set_kavg(5,5);
    mf.make_nets();
    mf.set_dt(1);
    mf.set_stochastic_parameters(std::vector<double>{shape,scale});
    mf.do_secondary_disables(0);
    mf.do_secondary_disables(1);
    mf.assess_total_functionality();
    mf.set_disable_behavior(damage_process_code);
    mf.set_restore_rate(restore_rate);
    if(oneway_dependency){
        mf.set_impacting_factor(impacting_factor, 0);
    }
    else { 
        mf.set_impacting_factor(impacting_factor, impacting_factor);
    }
    mf.set_dynamics_type(DynamicsType::Percolation,DynamicsType::Percolation);
    mf.set_do_enable_disable_together(true);
    mf.restore_all();
   
    mf.set_disable_rate(damage_rate);
    mf.set_restore_rate(restore_rate);
    mf.evolve(tmax,false);

    auto of = std::ofstream(ofname.str());
    write_mf_history_csv(mf, of);


    return 0;
}
