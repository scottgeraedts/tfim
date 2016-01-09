
#include "version.h"
#include "tfim.h"
#include <Eigen/Core>

using namespace std;

int main(){

	//read parameters from a file
	ifstream infile;
	infile.open("params");
	int Lx=value_from_file(infile,3);
//	cout<<this->nrows()<<endl;
	double Jx=value_from_file(infile,1.);
	double hx=value_from_file(infile,1.);
	double hy=value_from_file(infile,1.);
	double hz=value_from_file(infile,1.);
	#ifndef USE_COMPLEX
	if(hy!=0){
		cout<<"you can't use hy!=0 if you don't enable complex numbers!"<<endl;
		exit(0);
	}
	#endif	
	int seed=value_from_file(infile,1);
	infile.close();

	//setup random coefficients
	MTRand ran;
	ran.seed(seed);
	vector<double> alphax=vector<double>(Lx,0);
	vector<double> alphay=vector<double>(Lx,0);
	vector<double> alphaz=vector<double>(Lx,0);
//	ifstream datfile;
//	datfile.open("nicrans");
	for(int i=0; i<Lx; i++){
		alphax[i]= 2.*hx*(ran.rand()  -0.5);
		alphay[i]= 2.*hy*(ran.rand()  -0.5);
		alphaz[i]= 2.*hz*(ran.rand()  -0.5);
//getting random alphas from a file file sent by nicolas
//		datfile>>alphax[i]>>alphay[i]>>alphaz[i];
//		cout<<alphax[i]<<" "<<alphay[i]<<" "<<alphaz[i]<<endl;
	}
//	datfile.close();

#ifdef USE_COMPLEX
	MatrixTFIM< complex<double> > A(Lx,1,Jx,Jx,Jx,alphax,alphay,alphaz,-1);
#else
	MatrixTFIM<double> A(Lx,1,Jx,Jx,Jx,alphax,alphay,alphaz,-1);
#endif

	bool sparseSolve=true;
	int N_output_states,start,end;
	A.makeDense();
	if(!sparseSolve){	
		N_output_states=A.nrows();
		A.EigenDenseEigs();
		start=A.nrows()/3; end=2*A.nrows()/3;
	}
	else{
		N_output_states=200;
		N_output_states=A.eigenvalues(N_output_states,0.1);
		start=0; end=N_output_states;
	}
	A.entanglement_spacings(start,end);
}


