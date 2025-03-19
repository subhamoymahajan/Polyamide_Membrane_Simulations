#include <iostream>
#include <fstream>
#include <vector>
#include <unistd.h>
#include <cstdlib>
#include <iomanip> 

#define DX 0.2


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

void calc_poro_z(std::vector<bool> grid, std::vector<int> dims, std::string outname) {
     int i,j,k,poro;
     unsigned long int idx, A;
     A=dims[0]*dims[1];
     std::ofstream file(outname);
     file << "#z (nm), porosity" <<std::endl;
     for (k=0; k<dims[2];k++){
	     poro=0;
	     for (i=0; i<dims[0]; i++){
		     for(j=0;j<dims[0]; j++){
			     idx=i+j*dims[0]+k*dims[1]*dims[0];
			     if (grid[idx]){
				     poro++;
			     }

		     }
	     }
   	     file <<  std::setprecision(3) << std::fixed << 0.1*DX*(float)k << "," ;
   	     file <<  std::setprecision(6) << std::fixed << poro/(float)A << std::endl ;
     }
     file.close();
}



int main(int argc, char *argv[]) {
    std::string filename = "nitrogen_network.grd";
    std::vector<bool> grid;
    std::vector<int> dims;
    double cutoff=0.00001;
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
    int fac = 1; 
    readGridFromFile(filename, grid, dims, cutoff, reflect, fac);
    std::cout << "Done" << std::endl;

    calc_poro_z(grid, dims, "poro_z.dat");
    
    return 0;
}
