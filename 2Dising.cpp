
#include "version.h"
#include "tfim.h"
#include <Eigen/Core>

using namespace std;

int main(){

	//read parameters from a file
	ifstream infile;
	infile.open("params2d");
	int Lx=value_from_file(infile,3);
	int Ly=value_from_file(infile,3);
	double Jz=value_from_file(infile,1.);
	double hx=value_from_file(infile,1.);
	double hz_const=value_from_file(infile,1.);
	double hz_var=value_from_file(infile,1.);
	#ifndef USE_COMPLEX
//	if(hy!=0){
//		cout<<"you can't use hy!=0 if you don't enable complex numbers!"<<endl;
//		exit(0);
//	}
	#endif	
	int seed=value_from_file(infile,1);
	infile.close();

	//setup random coefficients
	MTRand ran;
	ran.seed(seed);
	vector<double> alphax=vector<double>(Lx*Ly,0);
	vector<double> alphay=vector<double>(Lx*Ly,0);
	vector<double> alphaz=vector<double>(Lx*Ly,0);
//	ifstream datfile;
//	datfile.open("nicrans");
	for(int i=0; i<Lx*Ly; i++){
		alphax[i]= hx; //2.*hx*(ran.rand()  -0.5);
		alphay[i]= 0; //2.*hy*(ran.rand()  -0.5);
		alphaz[i]= hz_const+hz_var*2*(ran.rand()-0.5);
//getting random alphas from a file file sent by nicolas
//		datfile>>alphax[i]>>alphay[i]>>alphaz[i];
//		cout<<alphax[i]<<" "<<alphay[i]<<" "<<alphaz[i]<<endl;
	}
//	datfile.close();

#ifdef USE_COMPLEX
	MatrixTFIM< complex<double> > A(Lx,Ly,0,0,Jz,alphax,alphay,alphaz,-1);
#else
	MatrixTFIM<double> A(Lx,Ly,0,0,Jz,alphax,alphay,alphaz,-1);
#endif

	bool sparseSolve=false;
	if(Lx*Ly>0) sparseSolve=true;
	int N_output_states,start,end;
	if(!sparseSolve){	
		N_output_states=A.nrows();
		A.makeDense();
		A.EigenDenseEigs();
		start=A.nrows()/3; end=2*A.nrows()/3;
	}
	else{
		N_output_states=200;
		clock_t t=clock();
		double target=A.find_middle();
		t=clock()-t;
		cout<<"middle time"<<(float)t/CLOCKS_PER_SEC<<endl;
		t=clock();
		N_output_states=A.eigenvalues(N_output_states,target);
		t=clock()-t;
		cout<<"diag time"<<(float)t/CLOCKS_PER_SEC<<endl;
		start=0; end=N_output_states;
	}
	for(int i=0;i<N_output_states;i++) cout<<A.eigvals[i]<<endl;
	A.energy_spacings();
}


