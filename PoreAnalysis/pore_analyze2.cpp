#include <iostream>
#include <fstream>
#include <vector>
#include <unistd.h>
#include <cstdlib>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <iomanip> 

#define DX 0.2

struct VertexProperty{
     unsigned int idx;
};
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, VertexProperty> Graph;
//struct NodeIdxMap {
//	NodeIdxMap(const std::vector<VertexProperty>& nodes): nodes_(nodes) {}
//	NodeIdxMap& operator=(const NodeIdxMap&)= delete ;
//	NodeIdxMap(const NodeIdxMap&) = delete;
//
//	typedef boost::readble_property_map_tag category;
//	typedef VertexProperty key_type;
//	typedef int value_type;
//
//	int operator[](const VertexProperty& node) const {
//		return node.idx;
//	}
//	private:
//	     const std::vector<VertexProperty>& nodes_;
//};

void idx2xyz(unsigned long int idx,std::vector<int> dims, std::vector<int> &xyz){
    // idx  = x + y*dims[0] + z*dims[0]*dims[1]
    xyz[0] = idx%dims[0];
    xyz[1] = (idx/dims[0])%dims[1];
    xyz[2] = idx/(dims[0]*dims[1]);
}

unsigned long int xyz2idx(std::vector<int> xyz, std::vector<int> dims, bool reflect){
    unsigned long int D=1, idx=0;
    int i;
  //  std::cout<<std::endl;
    for (i=0; i<3; i++){
        if (reflect==i){
	    idx+= (unsigned long int)xyz[i]*D;
//	    std::cout << "a. i="<< i << " idx = " << idx << " D = " << D << " xyz[i] = " << xyz[i]*D <<std::endl;
	} else{
	    idx+= (unsigned long int)((xyz[i]+dims[i])%dims[i])*D;
//	    std::cout << "b. i="<< i << " idx = " << idx << " D = " << D<< " xyz[i] = " << xyz[i]%dims[i] << std::endl;
	}
	D*=(unsigned long int)dims[i];
    }
    //std::cout<<std::endl;
    return idx;
}

void readGridFromFile(const std::string &filename, std::vector<bool> &grid, std::vector<int> &dims, double cutoff, int reflect, int fac) {

    std::ifstream file(filename);
    int i; 
    std::vector<int> xyz={0,0,0}, old_dim={1,1,1};
    unsigned long int Nmax,cnt=0,diff,idx;
    std::string word;
 
    if (!file.is_open()) {
        std::cerr << "Error opening file" << std::endl;
        return;
    }

    // Skip header lines
    std::string line;
    for (i = 0; i < 3; i++) { // Adjust this if the number of header lines changes
        std::getline(file, line);
    }
    for (i=0;i<3;i++){
         file >> word;
         dims.push_back(std::stoi(word)+1);
    }
    std::getline(file,line);
    std::getline(file,line);
    std::cout << "Dims is : "  << dims[0] << " " << dims[1] << " " << dims[2] << std::endl;
    old_dim[0]=dims[0];
    old_dim[1]=dims[1];
    old_dim[2]=dims[2];

    Nmax = dims[0]*dims[1]*dims[2];
    if (reflect >=0 ){
         Nmax/=2;
	 dims[reflect]/=2;
	 std::cout << "Reflected in direction " << reflect << " Dims is now : " << dims[0] << " " << dims[1] << " " << dims[2] << std::endl; 
	 std::cout << "Reflected in direction " << reflect << "old  Dims is now : " << old_dim[0] << " " << old_dim[1] << " " << old_dim[2] << std::endl; 
    }
    if (fac>1){
         dims[0]/=fac;
         dims[1]/=fac;
         dims[2]/=fac;
	 Nmax=dims[0]*dims[1]*dims[2];
	 std::cout << "Precision scaled by  " << fac << " Dims is now : " << dims[0] << " " << dims[1] << " " << dims[2] << std::endl; 

    }
    for (cnt=0;cnt<Nmax; cnt++){
	 grid.push_back(false);
    }
    cnt=0;
   // Read grid data
   while (file >> word){
        idx2xyz(cnt, old_dim,xyz) ;
        if (xyz[0]/fac >=dims[0] || xyz[1]/fac >=dims[1] || xyz[2]/fac >=dims[2]){
             //std::cout << "cnt " << cnt << " xyz " << xyz[0] << " " << xyz[1] << " "<< xyz[2] << std::endl;
	     cnt++;
             continue;
	}	

        if (std::stod(word)>cutoff){
	     idx= xyz[0]/fac + xyz[1]/fac*dims[0] + (xyz[2]/fac)*dims[0]*dims[1];
             grid[idx]=true;
        }
	cnt++;
   }

    file.close();
    // Example .grd format
    // Poreblazer nitrogen accessible grid: field = r
    // (1p,e12.5)
    //      48.047      48.047     169.360      90.000      90.000      90.000
    //         239         239         844
    // 1  0          239  0          239  0          844
    // 0.00000E+00
    // 0.00000E+00
    // 0.00000E+00
    // 0.00000E+00
    // 0.13747E+00
    // 0.29951E+00
}


void grid2graph(Graph &poreG, std::vector<bool> grid, std::vector<int> dims, int reflect, std::map<int,int> &poreG_map){

    unsigned long int cnt=0,cnt2=0, idx1, idx2;
    int i,j;
    
    std::vector<std::vector<int>> direc;
    direc.push_back({1, 0, 0});
    direc.push_back({0, 1, 0});
    direc.push_back({0, 0, 1});
    direc.push_back({-1, 0, 0});
    direc.push_back({0, -1, 0});
    direc.push_back({0, 0, -1});
    
    //std::cout << "direc size " << direc.size() << std::endl;
    std::vector<int> xyz={0,0,0}, foo_xyz={0,0,0};

    int max_z=0;
    for (const bool& ele: grid){
         if (ele){
 	     idx2xyz(cnt, dims,xyz) ;
             Graph::vertex_descriptor v1 = boost::add_vertex(poreG);
             poreG[v1].idx=cnt;
	     poreG_map[cnt]=cnt2;
	     if (xyz[2]>max_z){
                  max_z=xyz[2];
	     }
             cnt2++;
         }
         cnt++;
    }
    std::cout << "max_z " << max_z << ", "<< cnt2 << "  nodes, ";
    cnt=0;
    // Add periodic edges
    Graph::vertex_iterator v, vend;
    Graph::vertex_descriptor node1, node2; 
    max_z=0;
    bool cont=false;
    for (boost::tie(v,vend) = vertices(poreG); v != vend ; ++v){
//	 std::cout << "HERE" << std::endl;
//	 std::cout << poreG[*v].idx << std::endl;
 	 idx1 = poreG[*v].idx;
 	 idx2xyz(idx1, dims,xyz) ;
	 if (xyz[2]>max_z){
              max_z=xyz[2];
	 }
 
//	 std::cout << " idx1 = " << idx1 << ", xyz = " << xyz[0] << " " << xyz[1] << " " << xyz[2] << std::endl; 
 	 for (i=0;i<direc.size();i++){//different directions
	       cont=false;
               for (j=0; j< 3 ; j++){
      	          foo_xyz[j]=xyz[j]+direc[i][j];
		  if (reflect!=j){ //Apply PBC as long as the direction is not reflected
	             foo_xyz[j]=(foo_xyz[j]+dims[j])%dims[j];
		  }
		  if (foo_xyz[j]<0){
	             cont=true;
		  }
 	      }
	      if (cont){
		      continue;
	      }

//	      std::cout << " foo_xyz = " <<  foo_xyz[0] << " " << foo_xyz[1] << " " << foo_xyz[2]; 
//	      
//

 	      idx2 = xyz2idx(foo_xyz, dims, reflect); //Applies PBC in non-reflect directions
//              std::cout << " idx2 = " << idx2 << std::endl;
//            

               if (poreG_map.find(idx2) != poreG_map.end() ){ // found
                   boost::add_edge( *v, poreG_map[idx2], poreG);
		   cnt+=1;
 	      }
 
          }
    }
    std::cout << cnt << " edges" << std::endl;
    std::cout << "In graph  max_z = " << max_z << std::endl; 
}

void print_graph(Graph G){
    
    Graph::vertex_iterator v, vend;
    for (boost::tie(v,vend) = vertices(G); v != vend ; ++v){

         Graph::adjacency_iterator ai, a_end;
	 std::cout << "Node " << *v  << " (" << G[*v].idx << ") : " ;
         for (boost::tie(ai,a_end) = boost::adjacent_vertices(*v, G); ai != a_end; ++ai){
	     std::cout << *ai << " " ; 
	 }
    }

}

void find_connected_pores(Graph poreG, std::vector<int> dims, std::vector<int> &node2comp, std::vector<std::vector<Graph::vertex_descriptor>> &org_comp, int* num_comp){

//        std::vector<int> components(boost::num_vertices(poreG));
        std::size_t size=boost::num_vertices(poreG);
	node2comp.resize(size);

	*num_comp = boost::connected_components(poreG, &node2comp[0]);

	std::cout << *num_comp << " found " << std::endl;
        int i;

	for (i=0;i< *num_comp ; i++){//Initialize
	    org_comp.push_back({});
	}

        Graph::vertex_iterator v, vend;
        for (boost::tie(v,vend) = boost::vertices(poreG); v != vend ; ++v){
	    org_comp[(int)node2comp[*v]].push_back(*v);
	}

}


void find_surface_graph(Graph poreG, std::map<int,int> poreG_map, Graph &surfG, std::map<int,int> &poreS_map, std::vector<int> dims, std::vector<std::vector<Graph::vertex_descriptor>> Vnodes , std::vector<std::vector<Graph::vertex_descriptor>> &Snodes){

    Graph::vertex_iterator v, vend, v1;
    Graph::vertex_descriptor node1, node2; 
    unsigned long int cnt=0,cnt2=0,foo=0,cnt3=0;
    int max_z=0;
    //std::cout << " In surf: dims = " << dims[0] << " " << dims[1] << " " << dims[2] << std::endl;
    std::vector<int> xyz={0,0,0}; 
    //for (boost::tie(v,vend) = vertices(poreG); v != vend ; ++v){
    
    for (auto const comp: Vnodes){
       Snodes.push_back({});
       for (auto const node1: comp){   
         if (boost::degree(node1, poreG)!=6){
             Graph::vertex_descriptor v1 = boost::add_vertex(surfG);
             surfG[v1].idx=poreG[node1].idx;
 	     idx2xyz(poreG[node1].idx, dims,xyz) ;
	     Snodes[cnt3].push_back(v1);
//             if (xyz[2]>max_z){
//		  max_z=xyz[2];
//		  if (max_z>=dims[2]){
//		      std::cout<< "Error xyz " << xyz[0] << " " << xyz[1] << " " << xyz[2] << " > " << dims[2] << " ; idx = " << poreG[*v].idx << std::endl;
//		      exit(-1);
//		  }
//	     }
	     
	     poreS_map[poreG[node1].idx]=cnt2;
             cnt2++;
         } else if (boost::degree(node1,poreG)>6){
             std::cout << "Error degree = " << boost::degree(node1,poreG) << std::endl; 
	 }

         cnt++;
       }
       cnt3++;
    }
//    std::cout << "max_z = " << max_z << ", " << cnt2 << " nodes, " ;
    // Add edges
    for (boost::tie(v,vend) = vertices(surfG); v != vend ; ++v){
         Graph::adjacency_iterator ai, a_end;
	 
         for (boost::tie(ai,a_end) = boost::adjacent_vertices(poreG_map[surfG[*v].idx], poreG); ai != a_end; ++ai){
	      if (boost::degree(*ai,poreG)!=6){
		  boost::add_edge(*v,poreS_map[poreG[*ai].idx],surfG);
		  foo++;
	      }
	 }    
    }
    //std::cout << foo <<  " edges" << std::endl;

}

void create_sid(Graph::vertex_iterator v,std::map<Graph::vertex_descriptor,std::vector<int>>* pbc, std::map<int,std::vector<Graph::vertex_descriptor>>* sid_to_nodes, std::map<Graph::vertex_descriptor,int>* node_2_sid, int sid){
	node_2_sid->insert({*v,sid});
	pbc->insert({*v,{0,0,0}});
	sid_to_nodes->insert({sid,{*v}});
}

void add_to_sid(Graph::vertex_descriptor v,std::map<Graph::vertex_descriptor,std::vector<int>>* pbc, std::map<int,std::vector<Graph::vertex_descriptor>>* sid_to_nodes, std::map<Graph::vertex_descriptor,int>* node_2_sid, int sid, std::vector<int> pbc_diff){
	node_2_sid->insert({v,sid});
	(*sid_to_nodes)[sid].push_back(v);
        for (int i=0;i<3;i++){
   	    (*pbc)[v][i]+=pbc_diff[i];
	}
}

void transfer_sid(int sid1, int sid2, std::map<Graph::vertex_descriptor,std::vector<int>>* pbc, std::map<int,std::vector<Graph::vertex_descriptor>>* sid_to_nodes, std::map<Graph::vertex_descriptor,int>* node_2_sid, std::vector<int> pbc_diff){ 
	// Transfer from sid1 -> sid2;
	for( auto n: (*sid_to_nodes)[sid1]){
	    (*sid_to_nodes)[sid2].push_back(n);
	    (*node_2_sid)[n]=sid2;
            for (int i=0;i<3;i++){
   	        (*pbc)[n][i]+=pbc_diff[i];
	    }

	}
        (*sid_to_nodes).erase(sid1);
}


void find_pbc(Graph G, std::vector<int> dims, std::map<Graph::vertex_descriptor,std::vector<int>>* pbc, int reflect){
     	
    Graph::vertex_iterator v, vend, v1;
    Graph::vertex_descriptor node1, node2; 
    int j, sid_cnt=0;
    //std::map<Graph::vertex_descriptor,std::vector<int>> pbc; //PBC vector for subcomponent termed sub_pbc
    std::map<int,std::vector<Graph::vertex_descriptor>> sid_to_nodes; //All nodes for a given sub_pbc
    std::map<Graph::vertex_descriptor,int> node_2_sid; //sub_pbc id of a fiven node
    std::vector<int> pbc_diff={0,0,0};


    for (boost::tie(v,vend) = vertices(G); v != vend ; ++v){
	pbc->insert({*v,{0,0,0}});
    }
    return ;
    std::vector<int> xyz1={0,0,0}, xyz2={0,0,0};
    for (boost::tie(v,vend) = boost::vertices(G); v != vend ; ++v){
	 if ( node_2_sid.find(*v) == node_2_sid.end()){ //Initialize new sub_pbc
              create_sid(v, pbc, &sid_to_nodes, &node_2_sid, sid_cnt);
	      sid_cnt++;
	 }	 

 	 idx2xyz(G[*v1].idx, dims,xyz1) ; //pbc known
	 Graph::adjacency_iterator ai, a_end;
	
         for (boost::tie(ai,a_end) = boost::adjacent_vertices(*v, G); ai != a_end; ++ai){
 	      idx2xyz(G[*ai].idx, dims,xyz2) ;
	      for (j=0;j<3;j++){
		  if (xyz2[j]-xyz1[j] > dims[j]/2 && reflect!=j){
                       pbc_diff[j]= -1;
		  }else if (xyz2[j]-xyz1[j] < -dims[j]/2 && reflect!=j){
                       pbc_diff[j]= 1;
		  } 
	      }

	      if ( node_2_sid.find(*ai) != node_2_sid.end()){ //sid already assigned
		  if (node_2_sid[*ai]<node_2_sid[*v]){// Add all nodes in sid of v1 to ai
		      pbc_diff[0]*=-1;
		      pbc_diff[1]*=-1;
		      pbc_diff[2]*=-1;
                      transfer_sid(node_2_sid[*v], node_2_sid[*ai], pbc, &sid_to_nodes, &node_2_sid, pbc_diff);
		       
    		  } else if (node_2_sid[*ai]>node_2_sid[*v]){ // Add all nodes of ai to v1 
                      transfer_sid(node_2_sid[*ai], node_2_sid[*v], pbc, &sid_to_nodes, &node_2_sid, pbc_diff);
		  }	
	      }   else{ // Add node ai to sid of v1
	
                  add_to_sid(*ai, pbc, &sid_to_nodes, &node_2_sid, node_2_sid[*v], pbc_diff);
	      } 
	 }
    }

}

void find_spanning_pores(Graph G, std::vector<int> dims, std::vector<int> node2comp, std::vector<std::vector<Graph::vertex_descriptor>> org_comp, int num_comp, std::vector<std::vector<bool>> &spanning_components, std::vector<std::vector<bool>> &water_accessible, std::vector<std::vector<bool>> &water_inaccessible){


        int i,j,k, cnt,streak,streak_max;
        std::vector<int> xyz={0,0,0}, foo_xyz={0,0,0};
	int max_d = std::max(dims[0],std::max(dims[1],dims[2]));
	std::vector<std::vector<std::vector<bool>>> touch(num_comp, std::vector<std::vector<bool>>(3,std::vector<bool>(max_d,false)));
	for (i=0; i<num_comp;i++){
	    for (auto node: org_comp[i]){
		idx2xyz(G[node].idx, dims, xyz);
    	        for (j=0; j<3; j++){
	            touch[i][j][xyz[j]]=true;
    	   	}
	    }
	    spanning_components.push_back({false, false, false});
	    water_accessible.push_back({false, false, false});
	    water_inaccessible.push_back({false, false, false});
	}


        for (i=0;i< num_comp;i++){
	    for (j=0; j<3; j++){
		streak=0;
		streak_max=0;
		for (k=0;k<dims[j];k++){
		    if (touch[i][j][k]){
			 streak++;
			 streak_max=std::max(streak_max,streak);
		    } else {
			 streak=0;
		    }
                    	    

		}
  //              std::cout<< "Direction " << j << " Pore " << i  << " max_streak " << streak_max ;
	        if (streak_max>0.95*(float)dims[j]){
		    spanning_components[i][j]=true;
		}
                else if (touch[i][j][0] or touch[i][j][dims[j]]){
                    water_accessible[i][j]=true;
                }
                else {
                    water_inaccessible[i][j]=true;
                }
	        //std::cout << std::endl; 	
	    }

	    if (spanning_components[i][0] ||  spanning_components[i][1] || spanning_components[i][2]){ //CHANGE this condition for spanning a specific direction
	        std::cout<< "Pore " << i << " spanning in : x " << spanning_components[i][0] << " y " << spanning_components[i][1] << " z " << spanning_components[i][2] << std::endl; 
	    }

	}
        return ;

}

void gen_spanning_pore_nodes(Graph G, std::vector<int> dims, std::vector<std::vector<Graph::vertex_descriptor>> &nodes,std::vector<std::vector<Graph::vertex_descriptor>> org_comp, int num_comp, std::vector<std::vector<bool>> spanning_components, int direction ){

        int i;	
	for (i=0; i<num_comp; i++){
	    if (direction>0){
	        if (spanning_components[i][direction]){
		    nodes.push_back(org_comp[i]);
	        }
	    }
	    else {
	        if (spanning_components[i][0] || spanning_components[i][1] || spanning_components[i][2] ){
		    nodes.push_back(org_comp[i]);
	        }
	    }
	}

}

void gen_inacc_pore_nodes(Graph G, std::vector<int> dims, std::vector<std::vector<Graph::vertex_descriptor>> &nodes,std::vector<std::vector<Graph::vertex_descriptor>> org_comp, int num_comp, std::vector<std::vector<bool>> spanning_components, int direction ){

        int i;	
	for (i=0; i<num_comp; i++){
	    if (direction>0){
	        if (spanning_components[i][direction]){
		    nodes.push_back(org_comp[i]);
	        }
	    }
	    else {
	        if (spanning_components[i][0] && spanning_components[i][1] && spanning_components[i][2] ){
		    nodes.push_back(org_comp[i]);
	        }
	    }
	}

}



int write_pores(Graph G,std::vector<int> dims,  std::vector<std::vector<Graph::vertex_descriptor>> nodes, std::string outname, std::string header,std::map<Graph::vertex_descriptor,std::vector<int>> pbc, int reflect, int fac ){

     std::ofstream ofile(outname);
     if (!ofile.is_open()){
	std::cerr << outname << " already open. Unable to write" << std::endl;
	return 1;
     }
     ofile << header << " . Created by pore_analysis.cpp #Subhamoy Mahajan" << std::endl;
     unsigned long int N=0;
     std::vector<int> xyz={0,0,0};
     for( auto ci : nodes){
	  N+=ci.size();
     } 
     ofile << N << std::endl;
     int resid=1,aid=1,bid=0,i;
     std::stringstream foo;
     std::string resname, aname;
     int pbc_foo;
     for(auto ci: nodes){
	 for (Graph::vertex_descriptor ai : ci){
 	      idx2xyz(G[ai].idx, dims,xyz) ;
	      foo.str("");
	      if (resid<65535){
 	          foo << "P" << std::left << std::setw(4) << std::hex << resid ;
	      } else if (resid<65535*2){
 	          foo << "Q" << std::left << std::setw(4) << std::hex << resid%65535 ;
	      }	 else if (resid<65535*3){
 	          foo << "R" << std::left << std::setw(4) << std::hex << resid%65535 ;
	      }	  else if (resid<65535*4){
 	          foo << "T" << std::left << std::setw(4) << std::hex << resid%65535 ;
	      }	   
	      resname=foo.str();
	      foo.str("");
	      foo << "S" << std::left << std::setw(4) << bid ;
	      aname=foo.str();

	      ofile << std::setw(5) << resid  <<  resname << std::setw(5) << aname << std::setw(5) << aid ;
	      for (i=0;i<3;i++){
		   pbc_foo=0;
		   if (pbc.size()>0){
			pbc_foo=pbc[ai][i];
		   }
		      
   	           ofile << std::setw(8) << std::setprecision(3) << std::fixed << 0.1*DX*(float)((xyz[i]+dims[i]*pbc_foo)*fac)+1E-6;
		  
	      }
	      ofile << std::endl;
	      if (pbc.size()==0){
		  pbc_foo=0;
	      }
	      else if (pbc[ai][reflect]!=0){
		  std::cout << "PBC not behaving " << std::endl;
		  exit(-1);
              }
	      if (aid==99999){
		  bid++;
	      }
	      aid=(aid+1)%100000;
	      
	 }
	 resid=(resid+1)%100000;
     }
     ofile << 0.1*DX*(float)(dims[0]*fac) << " ";
     ofile << 0.1*DX*(float)(dims[1]*fac) << " ";
     ofile << 0.1*DX*(float)(dims[2]*fac) << " " << std::endl;
     ofile.close();
     return 0;
}

void get_vol_from_surf(Graph poreG, std::vector<int> dims, std::vector<int>pnode2comp, std::vector<std::vector<Graph::vertex_descriptor>> porg_comp, int num_pcomp, std::vector<std::vector<Graph::vertex_descriptor>> nodes, std::vector<std::vector<Graph::vertex_descriptor>> &Vnodes  ){

     int cid;
     for (auto const comp: nodes){
	  cid=pnode2comp[comp[0]];
	  std::cout << " Surface node = " << comp[0] << "  Component idx  = " << cid << " no. of nodes in volume " << porg_comp[cid].size() << std::endl;
          Vnodes.push_back(porg_comp[cid]);
     } 
}


int main(int argc, char *argv[]) {
    std::string filename = "nitrogen_network.grd";
    std::vector<bool> grid;
    std::vector<int> dims;
    double cutoff=1.4;
    std::string outname = "output.dat";
    int reflect = -1; //None
    // Parse command-line options using getopt
    int option;
    while ((option = getopt(argc, argv, "f:c:o:r:")) != -1) {
        switch (option) {
            case 'f':
                filename = optarg;
                break;
            case 'c':
                cutoff = std::stod(optarg);
                break;
            case 'o':
                outname = optarg;
                break;
            case 'r':
                reflect=std::stoi(optarg);
		break;
            case '?':
		if (optopt == 'f' || optopt == 'c'){
			std::cerr << "Option -" << static_cast<char>(optopt) << "requires an argument." << std::endl;
		} else {
			std::cerr << "Unknown option: -" <<  static_cast<char>(optopt) << std::endl;
		}
            default:
                std::cerr << "Usage: " << argv[0] << " -f <filename> -c <cutoff> -o <outname>" << std::endl;
                return 1;
        }
    }
    std::cout << "Reading GRD file : ";       
    int fac = 2; 
    readGridFromFile(filename, grid, dims, cutoff, reflect, fac);
    std::cout << "Done" << std::endl;
    
    std::cout << "Generating Graph : ";        
    Graph poreG, surfG;
    std::map<int,int> poreG_map, surfG_map;
    grid2graph(poreG, grid, dims, reflect, poreG_map);
    std::cout << "Done" << std::endl;
    int num_comp, num_comp2;
    std::vector<int> node2comp;
    std::vector<std::vector<Graph::vertex_descriptor>> components;
    std::cout <<"Finding connected volume components (Unique pores): "; 
    find_connected_pores(poreG, dims, node2comp, components, &num_comp);
    std::cout <<"Done" << std::endl;

    std::vector<std::vector<bool>> spanning_components, water_accessible, water_inaccessible;  

    std::cout <<"Finding pores spanning/water accessible/water inaccessible the boundaries: "; 

    find_spanning_pores(poreG, dims, node2comp, components, num_comp, spanning_components, water_accessible, water_inaccessible);
    std::vector<std::vector<Graph::vertex_descriptor>> Vnodes, Snodes;
    gen_spanning_pore_nodes(poreG, dims, Vnodes, components, num_comp, spanning_components, reflect);
    std::map<Graph::vertex_descriptor,std::vector<int>> vol_pbc, surf_pbc;
    write_pores(poreG, dims, Vnodes, "vol_mem_spanning.gro", " ", vol_pbc, reflect, fac );
    find_surface_graph(poreG, poreG_map, surfG, surfG_map,  dims, Vnodes, Snodes);
    write_pores(surfG, dims, Snodes, "surf_mem_spanning.gro", " ", surf_pbc, reflect, fac );


    std::vector<std::vector<Graph::vertex_descriptor>> Vnodes2, Snodes2;
    if (reflect>0){
        gen_spanning_pore_nodes(poreG, dims, Vnodes2, components, num_comp, water_accessible, reflect);
    }
    //std::cout <<"Done: "<< Vnodes2.size() << std::endl;
    std::map<Graph::vertex_descriptor,std::vector<int>> vol_pbc2, surf_pbc2;
    write_pores(poreG, dims, Vnodes2, "vol_wat_acc.gro", " ", vol_pbc2, reflect, fac );
    find_surface_graph(poreG, poreG_map, surfG, surfG_map,  dims, Vnodes2, Snodes2);
    write_pores(surfG, dims, Snodes2, "surf_wat_acc.gro", " ", surf_pbc2, reflect, fac );


    std::vector<std::vector<Graph::vertex_descriptor>> Vnodes3, Snodes3;
    gen_inacc_pore_nodes(poreG, dims, Vnodes3, components, num_comp, water_inaccessible, reflect);
    //std::cout <<"Done: "<< Vnodes3.size() << std::endl;
    std::map<Graph::vertex_descriptor,std::vector<int>> vol_pbc3, surf_pbc3;
    write_pores(poreG, dims, Vnodes3, "vol_wat_inacc.gro", " ", vol_pbc3, reflect, fac );
    find_surface_graph(poreG, poreG_map, surfG, surfG_map,  dims, Vnodes3, Snodes3);
    write_pores(surfG, dims, Snodes3, "surf_wat_inacc.gro", " ", surf_pbc3, reflect, fac );


    return 0;
}
