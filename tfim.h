
#ifndef TFIM_H
#define TFIM_H

//#include "blas1c.h"
//#include "lapackc.h"
#include "MersenneTwister.h"
#include "utils.h"
#include <vector>
#include "matprod.h"
#include "wavefunction.h"
#include "setup.h"
using namespace std;

template<class ART>
class MatrixTFIM: 
	public MatrixWithProduct<ART>,
	public Wavefunction<ART> {

 private:

	vector<double> alphax,alphay,alphaz;
	double Jx,Jy,Jz, hx,hy,hz;
	int nStates;
	int next(int);

 public:

	vector<int> states;
	int bitflip(int,int);
	int bittest(int,int);
	void MultMv(ART* v, ART* w);
//	Eigen::Matrix<ART,-1,1> MultEigen(Eigen::Matrix<ART,-1,1>);

	void make_disorder(int seed);
	void make_states(int charge);
  MatrixTFIM(int);
  	int L(){return nStates;}

}; // MatrixTFIM.

template<class ART>
int MatrixTFIM<ART>::bittest(int state,int bit){	
	if (state & 1<<bit) return 1;
	else return 0;
}

template<class ART>
int MatrixTFIM<ART>::next(int i){
	if(i==nStates-1) return 0;
	else return i+1;
}
template<class ART>
void MatrixTFIM<ART>::make_states(int charge=-1){
	for(int i=0;i<1<<nStates;i++)
		if(charge==-1 || count_bits(i)==charge) states.push_back(i);
	this->setrows(states.size());
}
//  Matrix-vector multiplication w <- M*v.
//you might think (from Glen Evenbly's website) that you need some factors of 1/2, but apparently you do not
template<class ART>
void MatrixTFIM<ART>::MultMv(ART* v, ART* w){
	int countJ=0;
	int sign;
	double countH;
	for(int in=0; in<this->ncols(); in++) w[in]=0.;
	for(int in=0; in<this->ncols(); in++){
		//off-diagonal elements
		countH=0.;
		for(int i=0; i<nStates; i++){
			if(bittest(states[in],i)) sign=1;
			else sign=-1;

			//hx, hy
			//w[lookup_flipped(states[in],i,states)]+=(alphax[i]*hx + sign*alphay[i]*hy*complex<double>(0,1.))*v[in];
			w[states[in]^1<<i]+=0.5*(alphax[i]*hx + sign*alphay[i]*hy*complex<double>(0,1.))*v[in];
			//w[states[in]^1<<i]+=0.5*(alphax[i]*hx )*v[in];
			//hz
			countH+=0.5*sign*alphaz[i];
		}
		
		//diagonal elements
		countJ=0;
		for(int i=0; i<nStates;i++){
			if(bittest(states[in],i) == bittest(states[in],next(i)) ) sign=1;
			else sign=-1;
			//Jx, Jy
//			if(sign==-1) w[lookup_flipped(states[in],i,next(i),states)]+=(Jx+Jy)*v[in];
			w[ states[in] ^ ( (1<<i) + (1<<next(i)) ) ]+=(Jx-sign*Jy)*v[in]*0.25;
			//Jz
			countJ+=sign;			
			
		}
		w[in]+=(countJ*Jz*0.25+countH*hz)*v[in];	
	}
} //  MultMv.

//template<class ART>
//Eigen::Matrix<ART,-1,1> MatrixTFIM<ART>::MultEigen(Eigen::Matrix<ART,-1,1> v){
//	int countJ=0;
//	Eigen::Matrix<ART,-1,1> w=Eigen::Matrix<ART,-1,1>::Zero(this->ncols());
//	for(int in=0; in<this->ncols(); in++){
//		//off-diagonal elements
//		for(int i=0; i<nStates; i++){
//			w(bitflip(in,i))+=alpha[i]*v(in);
//		}
//		//diagonal elements
//		countJ=0;
//		for(int i=0; i<nStates-1;i++){
//			if(bittest(in,i) == bittest(in,i+1) ) countJ--;
//			else countJ++;
//		}
////		if (bittest(in, nStates-1) == bittest(in,0) ) countJ--;
////		else countJ++;
//		w(in)+=countJ*(0.5*Jz)*v(in);	
//	}
//	return w;
//} //  MultMv.


////initializes the onsite disorder
template<class ART>
void MatrixTFIM<ART>::make_disorder(int seed){
	MTRand ran;
	ran.seed(seed);
	alphax=vector<double>(nStates,0);
	alphay=vector<double>(nStates,0);
	alphaz=vector<double>(nStates,0);
	ifstream datfile;
	datfile.open("nicrans");
	for(int i=0; i<nStates; i++){
//		alphax[i]= 2.*(ran.rand()  -0.5);
//		alphay[i]= 2.*(ran.rand()  -0.5);
//		alphaz[i]= 2.*(ran.rand()  -0.5);
//getting random alphas from a file file sent by nicolas
		datfile>>alphax[i]>>alphay[i]>>alphaz[i];
		cout<<alphax[i]<<" "<<alphay[i]<<" "<<alphaz[i]<<endl;
	}
	datfile.close();
}

//template<class ART>
//inline MatrixTFIM<ART>::MatrixTFIM(int x): MatrixWithProduct<ART>(){

//	nStates=12;
//	make_states(-1);
//	Jz=1; Jx=1; Jy=1;
//	hx=1; hy=1; hz=1;
//	make_disorder(nStates);
//	//sparse
//	this->eigenvalues(10,1.);
//	for(int i=0;i<10;i++) cout<<this->eigvals[this->lowlevpos[i]]<<endl;
//	cout<<endl;
//	//dense
////	this->makeDense();
////	this->EigenDenseEigs();
////	for(int i=0;i<this->nrows();i++) cout<<this->eigvals[i]<<endl;
//	
//}	
template<class ART>
inline MatrixTFIM<ART>::MatrixTFIM(int x): MatrixWithProduct<ART>()
// Constructor

{
	//initialize
	ifstream infile;
	infile.open("params");
	nStates=value_from_file(infile,3);
	int charge=-1;
	make_states(charge);
//	cout<<this->nrows()<<endl;
	Jz=value_from_file(infile,1.);
	Jx=Jz; Jy=Jz;
	hx=value_from_file(infile,1.);
	hy=value_from_file(infile,1.);
	hz=value_from_file(infile,1.);
	int seed=value_from_file(infile,1);
	infile.close();
	this->init_wavefunction(nStates);

	vector<double> s;
	make_disorder(seed);

	this->makeDense();
	this->EigenDenseEigs();
	
//	cout<<this->EigenDense<<endl;

//entanglement testing
//	int rhosize;
//	Eigen::Matrix<ART,-1,-1> rho;
//	Eigen::SelfAdjointEigenSolver<Eigen::Matrix<ART,-1,-1> > rs;	
//	double vn=0;
//	for(int nspins=1;nspins<=4;nspins++){
//		rhosize=this->ee_setup(nspins,nStates,states);
//		rho=Eigen::Matrix<ART,-1,-1>::Zero(rhosize,rhosize);
//		this->ee_compute_rho(this->eigvecs[111],rho,states);
//		rs.compute(rho);
//		cout<<"nspins: "<<nspins<<"-----"<<endl;
//		cout<<rs.eigenvalues()<<endl;
//		vn=0;
//		for(int i=0;i<rhosize;i++) 		
//			vn+=-rs.eigenvalues()(i)*log(rs.eigenvalues()(i));
//		cout<<"entropy: "<<vn<<endl;
//	}

//	for(int i=0;i<this->nrows();i++) cout<<this->eigvals[i]<<endl;
	vector<double> EE_levels,s_spacings;
	vector<double> EE_levels_all,s_spacings_all;
	vector< vector<double> > EE_levels_storage;
	int rhosize=this->ee_setup(0,nStates/2,states);
	Eigen::Matrix<ART,-1,-1> rho=Eigen::Matrix<ART,-1,-1>::Zero(rhosize,rhosize);
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix<ART,-1,-1> > rs;	
	ofstream Lout,Sout;
	Lout.open("EE_levels");
	Sout.open("spacings");
	int start=this->nrows()/3,end=2*this->nrows()/3;
	for(int i=start;i<end;i++){
		rho=Eigen::Matrix<ART,-1,-1>::Zero(rhosize,rhosize);
		this->ee_compute_rho(this->eigvecs[i],rho,states);
		rs.compute(rho);
		cout<<"raw eigenvalues "<<i<<endl;
		cout<<rs.eigenvalues()<<endl;
		EE_levels.clear();
		for(int j=0;j<rhosize;j++) 
			if(rs.eigenvalues()(j)>0.) EE_levels.push_back(-log(rs.eigenvalues()(j)));
			else cout<<"error in eigenvalue! "<<rs.eigenvalues()(j)<<endl;
		sort(EE_levels.begin(),EE_levels.end());
		EE_levels_all.insert(EE_levels_all.end(),EE_levels.begin(),EE_levels.end());
		EE_levels_storage.push_back(EE_levels);
	}
	sort(EE_levels_all.begin(),EE_levels_all.end());
	vector<double> energy_grid=make_grid(EE_levels_all,50);
	vector<double> integrated_DOS=make_DOS(EE_levels_all,energy_grid);
	for(int i=start;i<end;i++){
		s=make_S(EE_levels_storage[i-start],energy_grid,integrated_DOS);
		s_spacings=spacings(s);
		s_spacings_all.insert(s_spacings_all.end(),s_spacings.begin(),s_spacings.end());
	}
	for(int j=0;j<(signed)EE_levels_all.size();j++) Lout<<EE_levels_all[j]<<" "<<endl;
	for(int j=0;j<(signed)s_spacings_all.size();j++) Sout<<s_spacings_all[j]<<endl;
	Lout.close();
	Sout.close();
}
#endif 


