
#include "version.h"
#include "tfim.h"
#include <string>
#include <bitset>
#include <Eigen/Core>

using namespace std;
//return integers that make nice entanglement cuts
int cuts(int Lx,int Ly){

	bitset<30> out;
//4x3:
// 0 0 0 0
// 1 1 1 0
// 1 1 1 0 ==0,1,2,4,5,6
	if(Lx*Ly==12) //careful! this won't work for 6x2!
		out=bitset<30>(string("000001110111"));
//5x3:
// 0 0 0 0 0
// 1 1 1 1 0
// 1 1 1 1 0 ==0,1,2,4,5,6
	if(Lx*Ly==15)
		out=bitset<30>(string("000000111101111"));

	if(Lx*Ly==16)
		out=bitset<30>(string("0000000011111111"));
	return out.to_ulong();	
}
int main(){

	clock_t CPUtime1=clock();
	time_t walltime1=time(NULL);
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
	if(Lx*Ly>12) sparseSolve=true;
	int N_output_states,start,end;
	if(!sparseSolve){	
		N_output_states=A.nrows();
		A.makeDense();
		A.EigenDenseEigs();
		start=A.nrows()/3; end=2*A.nrows()/3;
	}
	else{
		N_output_states=1000;
	//	clock_t CPUtime2=clock();
	//	time_t walltime2=time(NULL);
	//	double target=A.find_middle();
	//	CPUtime2=clock()-CPUtime2;
		double target=0.5;
		start=0; end=N_output_states;
//		walltime2=time(NULL)-walltime2;
//		cout<<"time to find middlish energies: "<<(float)CPUtime2/CLOCKS_PER_SEC<<" CPU time and "<<walltime2<<" walltime"<<endl;
		N_output_states=A.eigenvalues(N_output_states,target);
		start=0; end=N_output_states;
	}
//	for(int i=0;i<N_output_states;i++) cout<<A.eigvals[i]<<endl;
	A.energy_spacings();
	A.entanglement_spacings(start,end,cuts(Lx,Ly));
	CPUtime1=clock()-CPUtime1;
	walltime1=time(NULL)-walltime1;
	cout<<"total time: "<<(float)CPUtime1/CLOCKS_PER_SEC<<" CPU time and "<<walltime1<<" walltime"<<endl;
}


