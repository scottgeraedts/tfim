
#include "version.h"
#include "pardiso_wrapper.h"
#include "tfim.h"
#include <Eigen/Core>
#include <iomanip>

using namespace std;
ArpackError::ErrorCode ArpackError::code = NO_ERRORS;

int main(){

	clock_t CPUtime=clock();
	//time_t walltime=time(NULL);
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
	int NROD=value_from_file(infile,1);
	infile.close();

	//setup random coefficients

for(int j=seed;j<=seed+NROD;j++){
	MTRand ran;
	ran.seed(j);
	vector<double> alphax=vector<double>(Lx,0);
	vector<double> alphay=vector<double>(Lx,0);
	vector<double> alphaz=vector<double>(Lx,0);
	ofstream datfile;
	datfile.open("nicrans");
	for(int i=0; i<Lx; i++){
		alphax[i]= 2.*hx*(ran.rand()  -0.5);
		alphay[i]= 2.*hy*(ran.rand()  -0.5);
		alphaz[i]= 2.*hz*(ran.rand()  -0.5);
//getting random alphas from a file file sent by nicolas
		datfile<<setprecision(14)<<alphax[i]<<" "<<alphay[i]<<" "<<alphaz[i]<<endl;
//		cout<<alphax[i]<<" "<<alphay[i]<<" "<<alphaz[i]<<endl;
	}
	datfile.close();

#ifdef USE_COMPLEX
	MatrixTFIM< complex<double> > A(Lx,1,Jx,Jx,Jx,alphax,alphay,alphaz,-1);
#else
	MatrixTFIM<double> A(Lx,1,Jx,Jx,Jx,alphax,alphay,alphaz,-1);//careful!
	//MatrixTFIM<complex<double> >A(Lx,1,Jx,Jx,Jx,alphax,alphay,alphaz,-1);//careful!
#endif

	bool sparseSolve=false;
	if(Lx>=16) sparseSolve=true;
	int N_output_states,start,end;
	if(!sparseSolve){	
		A.makeDense();
		N_output_states=A.nrows();
		A.EigenDenseEigs();
		start=A.nrows()/3; end=2*A.nrows()/3;
		//start=0; end=A.nrows();
	}
	else{
		N_output_states=1000;
//		N_output_states=A.eigenvalues(N_output_states,0.1);
		CPUtime=clock();	
		A.makeSparse();
		cout<<"time to make matrix: "<<(float)(clock()-CPUtime)/(CLOCKS_PER_SEC)<<endl;
		CPUtime=clock();
		//SuperLU_Wrapper<complex<double> >mat2(A.nrows());
		Pardiso_Wrapper<double> mat2(A.nrows());
		mat2.CSR_from_Sparse(A.sparse);
		cout<<"time to setup LU stuff: "<<(float)(clock()-CPUtime)/(CLOCKS_PER_SEC)<<endl;
		CPUtime=clock();
		mat2.eigenvalues(N_output_states,0.0);
		cout<<"time to get eigenvalues: "<<(float)(clock()-CPUtime)/(CLOCKS_PER_SEC)<<endl;
		//A.eigvecs=vector< vector <complex< double> > >(N_output_states);
		A.eigvecs=vector< vector <double> >(N_output_states);
		for(int i=0;i<N_output_states;i++) A.eigvecs[i]=mat2.eigvecs[i];
		start=0; end=N_output_states;

	}
	//A.energy_spacings();
	//start=0; end=A.nrows();
	//for(int i=start;i<end;i++)
	//	A.entanglement_spacings(i,i+1,A.rangeToBitstring(0,Lx/2),i);
	
	stringstream folder;
	folder.str("");
	if(NROD>1) folder<<j;
	else folder<<".";
	//full cut EE
	for(int cut=2;cut<=Lx/2;cut++)
		A.entanglement_spacings(start,end,A.rangeToBitstring(0,cut),cut,folder.str());

	//one charge sector
//	for(int cut=2;cut<=Lx/2;cut++)
//		A.entanglement_spacings(start,end,A.rangeToBitstring(0,cut),cut,folder.str(),cut/2,cut/2+1);

//	A.plot_states(start,end);
//	unsigned int everyother=0;
//	for(int i=0; i<Lx; i+=2)
//		everyother=everyother | 1<<i;
//	A.entanglement_spacings(start,end,everyother);
	//cout<<"total time is"<<(float)(clock()-CPUtime)/CLOCKS_PER_SEC<<" CPU time and "<<time(NULL)-walltime<<" walltime"<<endl;
}
}


