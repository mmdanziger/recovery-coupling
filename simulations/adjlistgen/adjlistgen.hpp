#ifndef ADJLIST_H
#define ADJLIST_H

#include <vector>
#include <algorithm>
#include <random>
#include <utility>
#include <iostream>
#include <exception>
#include <string>
#include <sstream>
#include <fstream>


using std::vector;
using std::pair;
using std::cout;
using std::endl;
using std::exception;
using uint_t = unsigned int;
using adjlist_t = std::vector<std::vector<uint_t>>;
using uvector = std::vector<uint_t>;
typedef std::pair<int,std::pair<int,int> > DistancePair;

// class AdjListGen{
// private:
//   std::mt19937 gen;
//   
// public:
//   void generate_double_random(long N, double kavg1, double kavg2, uvector & flatAdjList, uvector & offsetList);
//   AdjListGen(){gen.seed(time(0));}
//   
//   
//   
// };

class WeightedListGen{
private:
	std::mt19937 gen;
	float get_weight(int s, int t);
	vector< vector < float> > weight_list;
	adjlist_t adjacency_list;
	int N,num_links;
	std::string ifname;
public:
	WeightedListGen():num_links(0){gen.seed(time(0));}
	WeightedListGen(std::string input_file_name):num_links(0),ifname(input_file_name){gen.seed(time(0));}
	void create_weighted_adjacency_list();
	void load_weighted_adjacency_list();
	void copy_weighted_adjacency_list(uvector& flatAdjList, uvector& offsetList, vector<float>& flatWeightList);

};

void WeightedListGen::load_weighted_adjacency_list(){

	std::ifstream ifs(ifname);

	for (std::string line; std::getline(ifs, line); )
	{

		std::istringstream iss(line);
		int a, b;
		float c;
		if (!(iss >> a >> b >> c) || iss.get() != EOF) { /* error! die? */ }
		//iss >> a >> b >> c;
		if(a >= adjacency_list.size() || b >= adjacency_list.size() ){
			adjacency_list.resize(std::max(a,b)+1);
			weight_list.resize(std::max(a,b)+1);
		}
		adjacency_list[a].push_back(b);
		weight_list[a].push_back(c);
		num_links++;
	}
	N = adjacency_list.size() / 2 ;
//	std::cout << "Loaded " << adjacency_list.size() << " nodes with " << num_links << " links.\n";
}

void WeightedListGen::copy_weighted_adjacency_list(uvector& flatAdjList, uvector& offsetList, vector<float>& flatWeightList){


	flatAdjList.resize( num_links );
	flatWeightList.resize( num_links);
	offsetList.resize(adjacency_list.size() + 1);

	int link_count = 0;
	for(int i = 0; i< adjacency_list.size(); ++i){
		offsetList[ i ] = link_count;
		for(int j=0; j<adjacency_list[i].size(); j++){
			flatAdjList[ link_count ] = adjacency_list[i][j];
			flatWeightList[ link_count ] = weight_list[i][j];
			link_count++;
		}
	}
	offsetList[adjacency_list.size()] = link_count;


}

void WeightedListGen::create_weighted_adjacency_list(){
	adjacency_list.resize(N);
	weight_list.resize(N);
	for(int i=0; i<N; ++i){


	}

}


void generate_double_random(long int N, double kavg1, double kavg2, uvector& flatAdjList, uvector& offsetList)
{
std::mt19937 gen;
  
gen.seed(time(0));

for(int net_idx=0; net_idx<2; ++net_idx){
	adjlist_t adjacency_list(N);
	uint_t s, t, num_links = 0, link_count = 0;
	vector<double> kbar = {kavg1,kavg2};
	std::uniform_int_distribution<int> randint(0,N-1);
	while (num_links < N * kbar[net_idx] / 2) {
		do {
			s = randint(gen);
			t = randint(gen);
		} while (s == t);
		if (std::find(adjacency_list[s].begin(), adjacency_list[s].end(), t)
				== adjacency_list[s].end()) {
			adjacency_list[s].push_back(t);
			adjacency_list[t].push_back(s);
			num_links++;
		}
	}

uvector ordering(N);
for (uint_t i = 0; i < N; i++) {
	ordering[i] = i;
}
uvector inverse_ordering(N);
for (uint_t i = 0; i < N; i++) {
	inverse_ordering[ordering[i]] = i;
}
//std::cout << "Added " <<num_links << " links\n";
uint lastNumLinks=0, offsetCorrection=0;
if (net_idx == 0){
	flatAdjList.resize( num_links * 2);
	offsetList.resize( N + 1);
} else {
	lastNumLinks = flatAdjList.size();
	offsetCorrection = N;
	flatAdjList.resize(lastNumLinks + num_links * 2);
	offsetList.resize(2*N + 1);
}
uint_t block_number;
for (uint_t i = 0; i < N; i++) {
	block_number = inverse_ordering[i];
	offsetList[offsetCorrection + i] = lastNumLinks + link_count;
//	std::cout << "Offset["<<i<<"] = "<<link_count<<"\n";
	for (uint_t j = 0; j < adjacency_list[block_number].size(); j++) {
		flatAdjList[link_count + lastNumLinks] = ordering[adjacency_list[block_number][j]] + offsetCorrection;
		link_count++;

	}
}
offsetList[N + offsetCorrection] = link_count + lastNumLinks;

}

}

adjlist_t generate_random_adjlist(long int N, double kavg, std::mt19937 &gen)
{



	adjlist_t adjacency_list(N);
	uint_t s, t, num_links = 0, link_count = 0;
	std::uniform_int_distribution<int> randint(0,N-1);
	while (num_links < N * kavg / 2) {
		do {
			s = randint(gen);
			t = randint(gen);
		} while (s == t);
		if (std::find(adjacency_list[s].begin(), adjacency_list[s].end(), t)
				== adjacency_list[s].end()) {
			adjacency_list[s].push_back(t);
			adjacency_list[t].push_back(s);
			num_links++;
		}
	}
	return adjacency_list;

}

adjlist_t generate_lattice_adjlist(long int N,bool periodic_bc=true){

    adjlist_t adjacency_list(N);
    uint L = int(sqrt(N));
    if(L*L != N){
        throw(1);
    }
    if(periodic_bc){
        for(int i=0; i<N; i++){
            if((i+1)%L == 0){
                adjacency_list[i].push_back(i+1-L);
                adjacency_list[i+1-L].push_back(i);
            }
            else{
                adjacency_list[i].push_back(i+1);
                adjacency_list[i+1].push_back(i);
            }
            
            adjacency_list[i].push_back((i+L)%N);
            adjacency_list[(i+L)%N].push_back(i);
            
        }
        
    } else {
        for(int i=0; i<N; i++){
            if((i+1)%L != 0){
                adjacency_list[i].push_back(i+1);
                adjacency_list[i+1].push_back(i); 
            }
            if (i/L < L-1){
                adjacency_list[i].push_back(i+L);
                adjacency_list[i+L].push_back(i);
            }
        }
        
    }
    
    return adjacency_list;
    
    
}




vector< pair<int, pair<int,int> > > calculate_distances(uint L)
{
    
    int uniquePairs = static_cast<int>(L*(L+1)/2) -1; //-1 to disallow (0,0) -> force at least (0,1)
     vector< pair<int, pair<int,int> > > pencils;
    pencils.resize( uniquePairs );
    //cout << "Resized pencils" <<endl;
    int m=0;
    for (int i = 1; i<L;++i){//start at ONE: no dist=0 links
        for (int j=0; j<=i;++j){
            try{
                if (m >= uniquePairs){
                    cout <<"I'm sorry Dave. I can't let "<<m<<" do that to ("<<i<<","<<j<<")..."<<endl;
                    return pencils;
                }
           pencils[m++] = std::make_pair(i*i+j*j, std::make_pair(i,j));
            }catch(exception &e){
                cout << "Pencil construction failed" <<endl;
                cout <<e.what() <<endl;
            }
        }
    }
    std::sort(pencils.begin(), pencils.end(), [](const DistancePair& dpA, const DistancePair& dpB){ return dpA.first < dpB.first;});
    
    return pencils;

}

int randsign(std::mt19937& gen){
 return gen() > (gen.max() / 2) ? 1 : -1; 
}

int draw_target(int sNode, int L, std::mt19937& gen, double zeta, vector< pair<int, pair<int,int> > >& pencils, int periodic_bc = 1 )
{
  std::exponential_distribution<double> exp(1 / zeta);
  double rexp,rexp2;
  do{
      rexp = exp(gen);
      rexp2 = rexp*rexp;}while(rexp > L/2 );
  auto it = std::upper_bound(pencils.cbegin(), pencils.cend(), rexp2, [](const int A, const DistancePair B){return A<B.first;});
  if(it!=pencils.cbegin() && ((*it).first + (*(it-1)).first > 2*rexp2 )  ){
    it--;
    if(it!=pencils.cbegin() && (*(it-1)).first == (*(it)).first){
      if (randsign(gen) > 0)
	it--;
    }
  } else if ((it+1) != pencils.cend() && (*(it+1)).first == (*(it)).first){
      if (randsign(gen) > 0)
	it++;
  }
  int dx=randsign(gen)*(*it).second.first;
  int dy=randsign(gen)*(*it).second.second;
  if (randsign(gen)>0){//switch i and j with 50% chance because only one pair appears in the pencils list
    std::swap(dx,dy);
    }
  int tNode_x = sNode%L+dx;
  int tNode_y = sNode/L+dy;
  if(periodic_bc){
      tNode_x = (tNode_x +L)%L;
      tNode_y = (tNode_y +L)%L;
  } else {
      if (tNode_y >= L || tNode_y <0)
          tNode_y-=(2*dy);
      if (tNode_x >= L || tNode_x <0)
          tNode_x-=(2*dx);
      if (tNode_y >= L || tNode_y <0 || tNode_x >= L || tNode_x <0) //if you still can't find something, return failure
          return -1;
  } 
  return tNode_y*L+tNode_x;
}


void generate_double_exponential(long int N, double kavg1, double kavg2, double zeta1, double zeta2, uvector& flatAdjList, uvector& offsetList)
{
std::mt19937 gen;
  
gen.seed(time(0));
uint_t L = sqrt(N);
auto pencils = calculate_distances(L);

for(int net_idx=0; net_idx<2; ++net_idx){
	adjlist_t adjacency_list(N);
	uint_t s, t, num_links = 0, link_count = 0;
	vector<double> kbar = {kavg1,kavg2};
	vector<double> zeta = {zeta1,zeta2};
	
	std::uniform_int_distribution<int> randint(0,N-1);
	while (num_links < N * kbar[net_idx] / 2) {
		do {
			s = randint(gen);
			t = draw_target(s,L,gen,zeta[net_idx],pencils);
		} while (s == t || t>=N);
		
		if (std::find(adjacency_list[s].begin(), adjacency_list[s].end(), t)
				== adjacency_list[s].end()) {
			adjacency_list[s].push_back(t);
			adjacency_list[t].push_back(s);
			num_links++;
		}
	}

uvector ordering(N);
for (uint_t i = 0; i < N; i++) {
	ordering[i] = i;
}
uvector inverse_ordering(N);
for (uint_t i = 0; i < N; i++) {
	inverse_ordering[ordering[i]] = i;
}
//std::cout << "Added " <<num_links << " links\n";
uint lastNumLinks=0, offsetCorrection=0;
if (net_idx == 0){
	flatAdjList.resize( num_links * 2);
	offsetList.resize( N + 1);
} else {
	lastNumLinks = flatAdjList.size();
	offsetCorrection = N;
	flatAdjList.resize(lastNumLinks + num_links * 2);
	offsetList.resize(2*N + 1);
}
uint_t block_number;
for (uint_t i = 0; i < N; i++) {
	block_number = inverse_ordering[i];
	offsetList[offsetCorrection + i] = lastNumLinks + link_count;
//	std::cout << "Offset["<<i<<"] = "<<link_count<<"\n";
	for (uint_t j = 0; j < adjacency_list[block_number].size(); j++) {
		flatAdjList[link_count + lastNumLinks] = ordering[adjacency_list[block_number][j]] + offsetCorrection;
		link_count++;

	}
}
offsetList[N + offsetCorrection] = link_count + lastNumLinks;

}

}

struct inverse_sample_sf{
	double kmaxgamma,kmingamma,gamma,deltakmaxkmingamma,oneOverOneMinusGamma;
	std::mt19937 gen;
	std::uniform_real_distribution<double> unireal;
	inverse_sample_sf(double kmin, double kmax, double gamma, std::mt19937 &gen) : gamma(gamma),gen(gen),unireal(0.0,1.0){
		kmaxgamma = pow(kmax,1-gamma);
		kmingamma = pow(kmin,1-gamma);
		deltakmaxkmingamma = kmaxgamma-kmingamma;
		oneOverOneMinusGamma = 1 / (1 - gamma);

	}
	int operator()()
	{
		return static_cast<int>(round(pow(unireal(gen) * deltakmaxkmingamma + kmingamma , oneOverOneMinusGamma)));
	}
};



adjlist_t generate_scale_free_adjlist(long int N, double sfgamma, std::mt19937 &gen){

    double kmin=1,kmax=sqrt(N);

	vector<int> degree_credits(N);
	inverse_sample_sf net_sampler(kmin,kmax,sfgamma,gen);
	std::generate(degree_credits.begin(), degree_credits.begin() + N, net_sampler);
	
	
	long Nstubs = std::accumulate(degree_credits.begin(), degree_credits.begin()+N, 0);
	if (Nstubs%2!=0){
		degree_credits[gen()%N]++;
		Nstubs++;
	}
	//std::cout << "generating net with " << Nstubs << " stubs.\n";
	std::vector<int> stubs( Nstubs);
	int up_to=0,copies=0;
	for(int i=0; i<N; i++){
		copies =  degree_credits[i];
		if(copies + up_to > Nstubs){
			std::cerr << "Overrunning stubs at node "<< i << " of " << N << "...expect problems\n";
			std::cerr << "Trying to access "<< copies + up_to << " >= " << Nstubs << "...expect problems\n";
		}
		std::fill_n(stubs.begin() +up_to, copies, i);
		up_to+=copies;
	}
	std::shuffle(stubs.begin(), stubs.end(), gen);
	adjlist_t adjacency_list(N);
	uint_t s, t, num_links = 0, link_count = 0;
	std::vector<int> fails;
	for(uint_t i=0; i<stubs.size(); i+=2){
	    int tries=0;
	    while(true){
		tries+=1;
		if(tries > 20){
			fails.push_back(i);
		}
		std::uniform_int_distribution<uint_t> randint(i+2,stubs.size()); //only used if s==t or s and t already linked
		s = stubs[i];
		t = stubs[i+1];
		if (s==t){
			std::swap(stubs[i+1],stubs[randint(gen)]);
			continue;
		}
		if (std::find(adjacency_list[s].begin(), adjacency_list[s].end(), t) == adjacency_list[s].end()) {
		adjacency_list[s].push_back(t);
		adjacency_list[t].push_back(s);
		num_links++;
		break;
		} else {
			std::swap(stubs[i+1],stubs[randint(gen)]);
			continue;
		}

	    }
  
	}
	if(!fails.empty())
	  std::cerr << "Failed to create " << fails.size() << " out of " << stubs.size() / 2 << " links.\n";
    return adjacency_list;
}

vector<int> get_degree_distribution(adjlist_t& adjlist, int normalized=0){
        vector<int> degree_distribution;
        int k;
        double total_degree=0;
        for( auto node_it = adjlist.begin(); node_it != adjlist.end(); node_it++){
              k = node_it->size();
              if (k >= degree_distribution.size()){
                  degree_distribution.resize(k+1);
              }
              degree_distribution[k]++;

              total_degree+=k;
        }
        if(normalized){
            std::for_each(degree_distribution.begin(), degree_distribution.end(), [&](int &n){ n/=total_degree; });
        }
        return degree_distribution;
}


void generate_double_scale_free(long int N, double sfgamma1, double sfgamma2,  uvector& flatAdjList, uvector& offsetList){
	double kmin=1,kmax=sqrt(N);

	vector<int> degree_credits(2*N);
	std::mt19937 gen;
	gen.seed(time(0));
	inverse_sample_sf net1_sampler(kmin,kmax,sfgamma1,gen);
	std::generate(degree_credits.begin(), degree_credits.begin() + N, net1_sampler);
	inverse_sample_sf net2_sampler(kmin,kmax,sfgamma2,gen);
	std::generate(degree_credits.begin()+N, degree_credits.begin() + 2*N, net2_sampler);
	
	
for(int net_idx=0; net_idx<2; ++net_idx){

	
	
	long Nstubs = std::accumulate(degree_credits.begin() + net_idx*N, degree_credits.begin()+N*(1 + net_idx), 0);
	if (Nstubs%2!=0){
		degree_credits[N*net_idx + gen()%N]++;
		Nstubs++;
	}
	std::vector<int> stubs( Nstubs);
	int up_to=0,copies=0;
	for(int i=0; i<N; i++){
		copies =  degree_credits[N*net_idx + i];
		if(copies + up_to > Nstubs){
			std::cerr << "Overrunning stubs at node "<< i << " of " << N << "...expect problems\n";
			std::cerr << "Trying to access "<< copies + up_to << " >= " << Nstubs << "...expect problems\n";
		}
		std::fill_n(stubs.begin() +up_to, copies, i);
		up_to+=copies;
	}
	std::shuffle(stubs.begin(), stubs.end(), gen);
	adjlist_t adjacency_list(N);
	uint_t s, t, num_links = 0, link_count = 0;
	std::vector<int> fails;
	for(uint_t i=0; i<stubs.size(); i+=2){
	    int tries=0;
	    while(true){
		tries+=1;
		if(tries > 20){
			fails.push_back(i);
		}
		std::uniform_int_distribution<uint_t> randint(i+2,stubs.size()); //only used if s==t or s and t already linked
		s = stubs[i];
		t = stubs[i+1];
		if (s==t){
			std::swap(stubs[i+1],stubs[randint(gen)]);
			continue;
		}
		if (std::find(adjacency_list[s].begin(), adjacency_list[s].end(), t) == adjacency_list[s].end()) {
		adjacency_list[s].push_back(t);
		adjacency_list[t].push_back(s);
		num_links++;
		break;
		} else {
			std::swap(stubs[i+1],stubs[randint(gen)]);
			continue;
		}

	    }
  
	}
	if(!fails.empty())
	  std::cerr << "Failed to create " << fails.size() << " out of " << stubs.size() / 2 << " links.\n";

	uvector ordering(N);
	for (uint_t i = 0; i < N; i++) {
		ordering[i] = i;
	}
	uvector inverse_ordering(N);
	for (uint_t i = 0; i < N; i++) {
		inverse_ordering[ordering[i]] = i;
	}
	//std::cout << "Added " <<num_links << " links\n";
	uint lastNumLinks=0, offsetCorrection=0;
	if (net_idx == 0){
		flatAdjList.resize( num_links * 2);
		offsetList.resize( N + 1);
	} else {
		lastNumLinks = flatAdjList.size();
		offsetCorrection = N;
		flatAdjList.resize(lastNumLinks + num_links * 2);
		offsetList.resize(2*N + 1);
	}
	uint_t block_number;
	for (uint_t i = 0; i < N; i++) {
		block_number = inverse_ordering[i];
		offsetList[offsetCorrection + i] = lastNumLinks + link_count;
	//	std::cout << "Offset["<<i<<"] = "<<link_count<<"\n";
		for (uint_t j = 0; j < adjacency_list[block_number].size(); j++) {
			flatAdjList[link_count + lastNumLinks] = ordering[adjacency_list[block_number][j]] + offsetCorrection;
			link_count++;

		}
	}
offsetList[N + offsetCorrection] = link_count + lastNumLinks;
}
}

template <typename Stream> void write_edge_list(const adjlist_t & adjlist, Stream & stream){
    stream << "#" << adjlist.size() << std::endl;
    for(uint s=0; s<adjlist.size(); s++){
        for(auto t: adjlist[s]){
                if(s<t)//don't repeat
                    stream << s << "\t" << t << std::endl;
        }
    }
    
}

#endif
