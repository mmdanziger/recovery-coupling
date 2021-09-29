#include "multifrailty.hpp"
#include "adjlistgen/adjlistgen.hpp"

Multifrailty::Multifrailty()
{
    gen.seed(time(0));
}

Multifrailty::Multifrailty(uint N) : N(N)
{

    status.resize(2);
    status[0].resize(N, NodeStatus::Operational);
    status[1].resize(N, NodeStatus::Operational);
    restore_queue.resize(2);
    disable_rate.resize(2);
    stochastic_disable_rate.resize(2);
    functionality.resize(2);
    restore_rate.resize(2);
    impacting_factor.resize(2);
    smoothed_functionality.resize(2);
    repaired_this_round.resize(2);
    disabled_this_round.resize(2);
    repaired_last_round.resize(2);
    repair_capacity.resize(2);
    is_smoothed_functionality_available = false;
    restoration_mode = RestorationMode::RestoreRate;
    sfgamma.resize(2, 3);
    kavg.resize(2, 4);
    topology.resize(2);
    primaryfails.resize(2);
    short_restore_time.resize(2);
    long_restore_time.resize(2);
    dynamics_type.resize(2);
    damage_process.resize(2, DamageProcess::Constant);
    auto microseconds_since_epoch = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    gen.seed(microseconds_since_epoch * getpid());
    dt = 0.1;
    print_degree_distribution = false;
    include_self_in_neighborhood_functionality = true;
    include_neighbors_in_neighborhood_functionality = true;
    impervious_rate = 0;
    do_enable_disable_together = false;
    //    make_nets();
}

void Multifrailty::make_nets()
{
    adjlist.clear();
    adjlist.resize(2);
    adjlist_t newadjlist;
    for (auto nx : nets)
    {
        switch (topology[nx])
        {
        case NetworkTopology::ER:
            newadjlist = generate_random_adjlist(N, kavg[nx], gen);
            break;
        case NetworkTopology::SF:
            newadjlist = generate_scale_free_adjlist(N, sfgamma[nx], gen);
            break;
        case NetworkTopology::Lattice:
            newadjlist = generate_lattice_adjlist(N);
            break;
        default:
            std::cerr << "Topology not supported.\n";
            break;
        }
        //         for(int i=0; i<N; i++){
        //             std::cout << "\n" << i << ":";
        //             for(auto j : newadjlist[i]){
        //                 std::cout <<" " << j;
        //             }
        //         }
        adjlist[nx] = newadjlist;
    }
    if (print_degree_distribution)
    {
        auto dd = get_degree_distribution(adjlist[0]);
        for (int i = 0; i < dd.size(); i++)
        {

            std::cerr << i << " " << dd[i] << std::endl;
        }
    }
}

void Multifrailty::set_disable_behavior(std::string argv_input)
{
    if (argv_input == "C")
        set_damage_process_type(DamageProcess::Constant, DamageProcess::Constant);
    if (argv_input == "M")
        set_damage_process_type(DamageProcess::Multiplicative, DamageProcess::Multiplicative);
    if (argv_input == "S")
        set_damage_process_type(DamageProcess::Stochastic, DamageProcess::Stochastic);
    if (argv_input == "SC")
        set_damage_process_type(DamageProcess::StochasticCorrelated, DamageProcess::StochasticCorrelated);
    //std::cout << "Using damage mode " << damage_process[0] << std::endl;
}

void Multifrailty::restore_all()
{
    std::fill(status[0].begin(), status[0].end(), NodeStatus::Operational);
    std::fill(status[1].begin(), status[1].end(), NodeStatus::Operational);
    do_secondary_disables(0);
    do_secondary_disables(1);
}

void Multifrailty::evolve_continuous(double until_t)
{
    std::uniform_real_distribution<double> real(0, 1);
    std::uniform_int_distribution<int> randint(0, N - 1);
    repaired_last_round = repaired_this_round;
    repaired_this_round[0] = 0;
    repaired_this_round[1] = 0;
    disabled_this_round[0] = 0;
    disabled_this_round[1] = 0;

    t = 0;
    while (t < until_t)
    {
        int node_idx = randint(gen);
        int net_idx = real(gen) > 0.5;
        if (get_node_status(net_idx, node_idx) == NodeStatus::PrimaryFail)
        {
            if ((real(gen) < get_impacted_restore_rate(net_idx, node_idx) * dt))
            {
                status[net_idx][node_idx] = NodeStatus::Operational;
                repaired_this_round[net_idx]++;
            }
        }
        else
        { //not primaryfail, eligible for failure
            if (real(gen) < get_disable_rate(net_idx) * dt)
            {
                disable_node(net_idx, node_idx);
            }
        }
        auto effective_dt = dt / (2 * N);
        if (static_cast<int>(t + effective_dt) - static_cast<int>(t))
        {
            repaired_this_round[0] = 0;
            repaired_this_round[1] = 0;

            if (impacting_factor[0] > 0 || impacting_factor[1] > 0)
            {
                do_secondary_disables(0);
                do_secondary_disables(1);
            }
            assess_functionality(0);
            assess_functionality(1);
        }
        t += effective_dt;
    }
}

void Multifrailty::evolve(double until_t, bool keep_history)
{

    t = 0;
    is_smoothed_functionality_available = false;
    if (!keep_history)
    {
        func_history.clear();
    }
    smoothed_functionality[0] = 0;
    smoothed_functionality[1] = 0;
    double smoothing_norm = 0;
    auto shape = stochastic_parameters[0];
    auto scale = stochastic_parameters[1];
    std::weibull_distribution<> weibull(shape, scale);
    while (t < until_t)
    {
        if (do_enable_disable_together)
        {
            for (auto netidx : nets)
            {
                do_enable_disable_round_rate_mode(netidx);
                do_secondary_disables(netidx);
                assess_functionality(netidx);
            }
        }
        else
        {
            for (auto netidx : nets)
            {
                do_enabling_round(netidx);
            }
            for (auto netidx : nets)
            {
                do_disabling_round(netidx);
                assess_functionality(netidx);
            }
        }
        if (t > 0.9 * until_t)
        {
            smoothed_functionality[0] += functionality[0];
            smoothed_functionality[1] += functionality[1];
            smoothing_norm += 1;
        }
        if ((abs(dt - 1) < 1e-6) || (static_cast<int>(t + dt) - static_cast<int>(t)))
        {

            repair_history.push_back(repaired_this_round);
            func_history.push_back(functionality);
            primaryfails_history.push_back(primaryfails);
            newprimaryfails_history.push_back(disabled_this_round);
            //std::cerr << disable_rate[0] << "\t" << disabled_this_round[0] << std::endl;
            double lambda_memory = 0;
            if (damage_process[0] == DamageProcess::Stochastic)
            {
                stochastic_disable_rate[0] = lambda_memory * stochastic_disable_rate[0] + (1 - lambda_memory) * weibull(gen) / N;
                stochastic_disable_rate[1] = lambda_memory * stochastic_disable_rate[1] + (1 - lambda_memory) * weibull(gen) / N;
            }
            else if (damage_process[0] == DamageProcess::StochasticCorrelated)
            {
                stochastic_disable_rate[0] = lambda_memory * stochastic_disable_rate[0] + (1 - lambda_memory) * weibull(gen) / N;
                stochastic_disable_rate[1] = stochastic_disable_rate[0];
            }
            else if (damage_process[0] == DamageProcess::Multiplicative)
            {
                for (auto netidx : nets)
                {
                    double mu = (1 - static_cast<double>(functionality[netidx]) / (N)) * disable_rate[netidx];
                    double sigma = std::sqrt(disable_rate[netidx] * (1 - disable_rate[netidx]));
                    std::normal_distribution<> normal(mu, sigma);
                    stochastic_disable_rate[netidx] = normal(gen);
                }
            }
            repaired_last_round = repaired_this_round;
            repaired_this_round[0] = 0;
            repaired_this_round[1] = 0;
            disabled_this_round[0] = 0;
            disabled_this_round[1] = 0;
        }
        //   std::cout << primaryfails[0]<<"\t"<< functionality[0]   <<"\n";
        //if(functionality[0] < 0.01*N && functionality[1] < 0.01*N)
        //  break;
        t += dt;
    }

    smoothed_functionality[0] /= smoothing_norm;
    smoothed_functionality[1] /= smoothing_norm;
    if (smoothing_norm > 5)
    {
        is_smoothed_functionality_available = true;
    }
}

void Multifrailty::evolve_one_step()
{
    is_smoothed_functionality_available = false;
    for (auto netidx : nets)
    {
        do_enabling_round(netidx);
    }
    for (auto netidx : nets)
    {
        do_disabling_round(netidx);
        assess_functionality(netidx);
    }
    //  std::cout << primaryfails[0]<<"\t"<< functionality[0]   <<"\n";

    t += dt;
}

void Multifrailty::do_disabling_round(uint net_idx)
{
    std::uniform_real_distribution<double> real(0, 1);
    for (auto node_idx = 0; node_idx < N; node_idx++)
    {
        if (real(gen) < get_disable_rate(net_idx) * dt && get_node_status(net_idx, node_idx) != NodeStatus::PrimaryFail)
        {
            disable_node(net_idx, node_idx);
        }
    }
    do_secondary_disables(net_idx);
}

void Multifrailty::do_secondary_disables(uint nx)
{
    switch (dynamics_type[nx])
    {
    case DynamicsType::Percolation:
        do_secondary_disables_percolation(nx);
        break;
    case DynamicsType::Watts:
        do_secondary_disables_watts(nx);
        break;
    default:
        throw(1);
        break;
    }
}

void Multifrailty::do_secondary_disables_percolation(uint nx)
{
    vector<int> components(N, -1);
    std::queue<uint> node_queue;
    uint compsize = 0, maxcompsize = 0, maxcompidx = -1;
    for (int i = 0; i < N; i++)
    {
        if (status[nx][i] == NodeStatus::SecondaryFail)
            status[nx][i] = NodeStatus::Operational;
    }
    for (uint i = 0; i < N; i++)
    {
        compsize = 0;
        if (status[nx][i] == NodeStatus::Operational && components[i] == -1)
        {
            node_queue.push(i);
            components[i] = i;
            compsize++; //count root on algorithm start
            while (!node_queue.empty())
            {
                int j = node_queue.front();
                node_queue.pop();
                for (auto t : adjlist[nx][j])
                {
                    if (components[t] == -1 && status[nx][t] == NodeStatus::Operational)
                    {
                        node_queue.push(t);
                        components[t] = i;
                        compsize++;
                    }
                }
            }
        }
        if (compsize > maxcompsize)
        {
            maxcompsize = compsize;
            maxcompidx = i;
        }
    }
    int secondaryfailcount = 0;
    std::uniform_real_distribution<double> real(0, 1);
    for (uint i = 0; i < N; i++)
    {
        if (components[i] != maxcompidx && status[nx][i] == NodeStatus::Operational)
        {
            status[nx][i] = NodeStatus::SecondaryFail;
            secondaryfailcount++;
        }
    }
    if (secondaryfailcount > 0)
    {
        //std::cout << t << " secondaryfailcount: "<<secondaryfailcount<<std::endl;
    }
}

void Multifrailty::do_secondary_disables_watts(uint net_idx)
{
    bool changes = false;
    for (int i = 0; i < N; i++)
    {
        if (status[net_idx][i] == NodeStatus::SecondaryFail)
            status[net_idx][i] = NodeStatus::Operational;
    }
    do
    {
        changes = false;
        for (int i = 0; i < N; i++)
        {
            if (status[net_idx][i] == Operational && assess_neighborhood_functionality(net_idx, i, false) < 0.5)
            {
                status[net_idx][i] = NodeStatus::SecondaryFail;
                changes = true;
            }
        }

    } while (changes);
}

void Multifrailty::set_node_status(uint net_idx, uint node_idx, NodeStatus new_status, double restore_time)
{
}

double Multifrailty::determine_restore_time(uint net_idx, uint node_idx)
{
    int nx = (net_idx + 1) % 2;
    return assess_neighborhood_functionality(nx, node_idx) > 0.5 ? short_restore_time[nx] : long_restore_time[nx]; //never! bwahahahahah!!!
}

const NodeStatus Multifrailty::get_node_status(uint net_idx, uint node_idx)
{
    return status[net_idx][node_idx];
}

const vector<int> Multifrailty::get_status_as_int(uint nx)
{
    vector<int> output(N);
    for (uint i = 0; i < N; i++)
    {
        switch (status[nx][i])
        {
        case NodeStatus::Operational:
            output[i] = 0;
            break;
        case NodeStatus::PrimaryFail:
            output[i] = -1;
            break;
        case NodeStatus::SecondaryFail:
            output[i] = -2;
            break;
        case NodeStatus::PermanentFail:
            output[i] = -3;
            break;
        default:
            break;
        }
    }
    return output;
}

void Multifrailty::set_status_from_int_vector(uint net_idx, std::vector<int> new_status)
{
    for (uint i = 0; i < N; i++)
    {
        switch (new_status[i])
        {
        case 0:
            status[net_idx][i] = NodeStatus::Operational;
            break;
        case -1:
            status[net_idx][i] = NodeStatus::PrimaryFail;
            break;
        case -2:
            status[net_idx][i] = NodeStatus::SecondaryFail;
            break;
        case -3:
            status[net_idx][i] = NodeStatus::PermanentFail;
            break;
        default:
            break;
        }
    }
}

void Multifrailty::disable_node(uint net_idx, uint node_idx)
{
    double restore_time = determine_restore_time(net_idx, node_idx);
    status[net_idx][node_idx] = NodeStatus::PrimaryFail;
    if (restoration_mode == RestorationMode::RestoreTime)
        restore_queue[net_idx].push(node(node_idx, t + restore_time));
    disabled_this_round[net_idx]++;
}

void Multifrailty::do_enabling_round(uint net_idx)
{
    switch (restoration_mode)
    {
    case RestoreTime:
        do_enabling_round_restore_time_mode(net_idx);
        break;
    case RestoreRate:
        do_enabling_round_restore_rate_mode(net_idx);
        break;
    default:
        std::cerr << "Restoration mode not supported." << std::endl;
    }
}

void Multifrailty::do_enabling_round_restore_time_mode(uint net_idx)
{
    int i = 0;
    while (!restore_queue[net_idx].empty())
    {
        auto node = restore_queue[net_idx].top();
        if (node.restoration_time < t)
        {
            status[net_idx][node.idx] = NodeStatus::Operational;
            restore_queue[net_idx].pop();
            i++;
        }
        else
        {
            break;
        }
    }
}

const vector<double> Multifrailty::get_all_neighborhood_functionality(uint nx)
{
    std::vector<double> out(N);
    for (int i = 0; i < N; i++)
    {
        out[i] = assess_neighborhood_functionality(nx, i, true);
    }
    return out;
}

double Multifrailty::get_disable_rate(uint net_idx)
{

    switch (damage_process[net_idx])
    {
    case DamageProcess::Constant:
        return disable_rate[net_idx];
        break;
    case DamageProcess::Stochastic:
    case DamageProcess::StochasticCorrelated:
    case DamageProcess::Multiplicative:
        return stochastic_disable_rate[net_idx];
        break;
    }

    return disable_rate[net_idx];

    double alpha = 0.85;
    double beta = 14;
    beta_distribution B(alpha, beta);

    return B(gen);
}

const double Multifrailty::get_impacted_restore_rate(uint net_idx, uint node_idx)
{

    uint other_net_idx = (net_idx + 1) % 2;
    double slow_factor = 1.0;
    if (impacting_factor[net_idx] > 0)
    {
        //slow_factor = pow(assess_neighborhood_functionality(other_net_idx,node_idx,true),impacting_factor[net_idx]);
        slow_factor = 1 - impacting_factor[net_idx] * (1 - assess_neighborhood_functionality(other_net_idx, node_idx, true));
    }
    if (use_repair_capacity)
    {
        double w = t - static_cast<int>(t);
        double repaired_this_round_avg = (1 - w) * repaired_last_round[net_idx] + (w)*repaired_this_round[net_idx];
        slow_factor *= (1 - (repaired_this_round_avg / repair_capacity[net_idx]));
        //std::cerr<<(1 - (repaired_this_round_avg / repair_capacity[net_idx]))<<std::endl;
    }
    return restore_rate[net_idx] * slow_factor;
    if (status[other_net_idx][node_idx] == NodeStatus::Operational)
        return restore_rate[net_idx];
    return impacting_factor[net_idx] * restore_rate[net_idx];
}

void Multifrailty::do_enabling_round_restore_rate_mode(uint net_idx)
{
    std::uniform_real_distribution<double> real(0, 1);
    for (auto node_idx = 0; node_idx < N; node_idx++)
    {
        if ((get_node_status(net_idx, node_idx) != NodeStatus::Operational) && (real(gen) < get_impacted_restore_rate(net_idx, node_idx) * dt))
        {
            status[net_idx][node_idx] = NodeStatus::Operational;
        }
        if (impervious_rate > 0)
        {
            if ((get_node_status(net_idx, node_idx) != NodeStatus::Operational) && (real(gen) < impervious_rate * dt) && (real(gen) < restore_rate[net_idx] * dt))
            {
                status[net_idx][node_idx] = NodeStatus::Operational;
            }
        }
    }
}

void Multifrailty::do_enable_disable_round_rate_mode(uint net_idx)
{
    std::uniform_real_distribution<double> real(0, 1);
    std::uniform_int_distribution<int> randint(0, N - 1);

    for (auto node_idx = 0; node_idx < N; node_idx++)
    {
        //if( get_node_status(net_idx,node_idx)  == NodeStatus::PrimaryFail){
        if (get_node_status(net_idx, node_idx) != NodeStatus::Operational)
        {
            if ((real(gen) < get_impacted_restore_rate(net_idx, node_idx) * dt))
            {
                status[net_idx][node_idx] = NodeStatus::Operational;
                repaired_this_round[net_idx]++;
            }
        }
        else
        { //not primaryfail, eligible for failure
            if (real(gen) < get_disable_rate(net_idx) * dt)
            {
                disable_node(net_idx, node_idx);
            }
        }
        t += (dt / (2 * N));
    }
}

void Multifrailty::assess_functionality(uint nx)
{
    functionality[nx] = 0;
    primaryfails[nx] = 0;
    for (uint i = 0; i < N; i++)
    {
        switch (status[nx][i])
        {
        case NodeStatus::Operational:
            functionality[nx]++;
            break;
        case NodeStatus::PrimaryFail:
            primaryfails[nx]++;
            break;
        default:
            break; //don't bother counting other states
        }
    }
}

void Multifrailty::assess_total_functionality()
{
    for (auto nx : nets)
        assess_functionality(nx);
}

double Multifrailty::assess_neighborhood_functionality(uint nx, uint node_idx, bool include_self)
{

    double num = 0, denom = 0;
    if (include_self_in_neighborhood_functionality && include_self)
    {
        denom++;
        if (status[nx][node_idx] == NodeStatus::Operational)
            num++;
    }
    if (include_neighbors_in_neighborhood_functionality)
    {
        for (auto j : adjlist[nx][node_idx])
        {
            denom++;
            if (status[nx][j] == NodeStatus::Operational)
                num++;
        }
    }

    return denom > 0 ? num / denom : 0;
}

uint Multifrailty::do_random_attack(uint net_idx, double damage_fraction, double tmax, double stop_above, double stop_below)
{
    func_history.clear();
    auto orig_damage_rate = disable_rate[0];
    auto orig_restore_rate = restore_rate[0];
    auto orig_dt = dt;
    set_disable_rate(damage_fraction);
    set_restore_rate(0);
    set_dt(1);
    evolve_one_step();

    auto func = get_total_functionality();

    set_dt(orig_dt);
    set_disable_rate(orig_damage_rate);
    set_restore_rate(orig_restore_rate);
    for (int i = 0; i < 100; i++)
    {
        evolve(tmax, true);
        func = get_total_functionality();
        if (func[0] < stop_below || func[0] > stop_above)
        {
            return func[0];
        }
    }
    return func[0];
}

void Multifrailty::do_local_attack(uint net_idx, uint root_node, uint attack_depth)
{
    //std::cout << "starting local attack at root "<<root_node<<" of depth " <<attack_depth << std::endl;
    disable_node(net_idx, root_node);
    std::queue<int> new_roots;
    int depth = 0;
    new_roots.push(root_node);
    new_roots.push(-1);
    std::map<uint, bool> visited;
    visited.insert(std::make_pair(root_node, true));
    while (!new_roots.empty())
    {
        int this_root = new_roots.front();
        new_roots.pop();
        if (this_root == -1)
        {
            depth++;
            if (depth > attack_depth)
                break;
            this_root = new_roots.front();
            if (this_root == -1) //crawled whole net
                break;
            new_roots.pop();
            new_roots.push(-1); //depth_marker
        }
        for (auto j : adjlist[net_idx][this_root])
        {
            if (visited.find(j) != visited.end()) //if visited
                continue;
            disable_node(net_idx, j);
            new_roots.push(j);
            visited.insert(std::make_pair(j, true));
            //              std::cout << j << "\t";
        }
    }
    //std::cout << std::endl;
}
