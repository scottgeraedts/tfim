
#include "version.h"
#include "tfim.h"
#include <Eigen/Core>

using namespace std;

int main(){

#ifdef USE_COMPLEX
	MatrixTFIM< complex<double> > A(1);
#else
	MatrixTFIM<double> A(1);
#endif
//	ifstream datfile;
//	datfile.open("temp");
//	int size=4096;
//	vector<double> energies(size);
//	for(int i=0;i<size;i++){ datfile>>energies[i]; }
//	vector<double> s,raw_spacings;
//	int meshes[]={10,20,50,100,200,500,1000,2000,5000,10000};
//	
//	for(int m=0;m<10;m++){
//		s=unfoldE(energies,meshes[m]);
//		cout<<"s r: "<<meshes[m]<<" "<<compute_r(s)<<endl;
//		raw_spacings=spacings(s);
//		for(int i=0;i<size-1;i++) cout<<meshes[m]<<" "<<raw_spacings[i]<<" #s spacings"<<endl;
//	}
	
}


