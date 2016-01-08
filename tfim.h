
#ifndef TFIM_H
#define TFIM_H

//#include "blas1c.h"
//#include "lapackc.h"
#include "MersenneTwister.h"
#include "utils.h"
#include <vector>
#include "matprod.h"
#include "wavefunction.h"
using namespace std;

template<class ART>
class MatrixTFIM: 
	public MatrixWithProduct<ART>,
	public Wavefunction<ART> {

 private:

	double Jx,Jy,Jz;
	vector<double> alphax,alphay,alphaz;
	double hx,hy,hz;
	int Lx,Ly,nStates;
	int next(int i, int dir);

 public:

	vector<int> states;
	int bitflip(int,int);
	int bittest(int,int);
	void MultMv(ART* v, ART* w);
//	Eigen::Matrix<ART,-1,1> MultEigen(Eigen::Matrix<ART,-1,1>);

	void make_disorder(int seed);
	void make_states(int charge);
  	int nSpins(){return nStates;}
  	
	MatrixTFIM(int _Lx, int _Ly, double _Jx, double _Jy, double _Jz, vector<double> _alphax, vector<double> _alphay, vector<double> _alphaz, int charge);
	void entanglement_spacings(int start, int end);	
	
}; // MatrixTFIM.

template<class ART>
int MatrixTFIM<ART>::bittest(int state,int bit){	
	if (state & 1<<bit) return 1;
	else return 0;
}

template<class ART>
int MatrixTFIM<ART>::next(int i, int dir){
	if(dir==0){
		if(i%Lx==Lx-1) return i-(Lx-1);
		else return i+1;
	}else if(dir==1){
		if((i/Lx)%Ly==Ly-1) return i-Lx*(Ly-1);
		else return i+Lx;
	}else{
		cout<<"error in next!"<<endl;
		return -1;
	}
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
			//w[lookup_flipped(states[in],i,states)]+=(alphax[i] + sign*alphay[i]*complex<double>(0,1.))*v[in];
			#ifdef USE_COMPLEX
			w[states[in]^1<<i]+=0.5*(alphax[i] + sign*alphay[i]*complex<double>(0,1.))*v[in];
			#else
			w[states[in]^1<<i]+=0.5*(alphax[i] )*v[in];
			#endif
			//hz
			countH+=0.5*sign*alphaz[i];
		}
		
		//diagonal elements
		countJ=0;
		for(int i=0; i<nStates;i++){
			if(bittest(states[in],i) == bittest(states[in],next(i,0)) ) sign=1;
			else sign=-1;
			//Jx, Jy
			w[ states[in] ^ ( (1<<i) + (1<<next(i,0)) ) ]+=(Jx-sign*Jy)*v[in]*0.25;
			//Jz
			countJ+=sign;
			if(Ly>1){
				if(bittest(states[in],i) == bittest(states[in],next(i,1)) ) sign=1;
				else sign=-1;
				//Jx, Jy
				w[ states[in] ^ ( (1<<i) + (1<<next(i,1)) ) ]+=(Jx-sign*Jy)*v[in]*0.25;
				//Jz
				countJ+=sign;			
			}			
			
		}
		w[in]+=(countJ*Jz*0.25+countH)*v[in];	
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
		alphax[i]= 2.*(ran.rand()  -0.5);
		alphay[i]= 2.*(ran.rand()  -0.5);
		alphaz[i]= 2.*(ran.rand()  -0.5);
//getting random alphas from a file file sent by nicolas
//		datfile>>alphax[i]>>alphay[i]>>alphaz[i];
//		cout<<alphax[i]<<" "<<alphay[i]<<" "<<alphaz[i]<<endl;
	}
	datfile.close();
}

template<class ART>
void MatrixTFIM<ART>::entanglement_spacings(int start, int end){

	vector<double> s;
	vector<double> EE_levels,s_spacings;
	vector<double> EE_levels_all,s_spacings_all;
	vector< vector<double> > EE_levels_storage;
	int rhosize=this->ee_setup(0,nStates/2,states);
	Eigen::Matrix<ART,-1,-1> rho=Eigen::Matrix<ART,-1,-1>::Zero(rhosize,rhosize);
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix<ART,-1,-1> > rs;	
	ofstream Lout,Sout,rout;
	Lout.open("EE_levels");
	Sout.open("spacings");
	rout.open("r");

	for(int i=start;i<end;i++){
		rho=Eigen::Matrix<ART,-1,-1>::Zero(rhosize,rhosize);
		this->ee_compute_rho(this->eigvecs[i],rho,states);
		rs.compute(rho);
//		cout<<"raw eigenvalues "<<i<<endl;
//		cout<<rs.eigenvalues()<<endl;
		EE_levels.clear();
		for(int j=0;j<rhosize;j++) 
			if(rs.eigenvalues()(j)>0.) EE_levels.push_back(-log(rs.eigenvalues()(j)));
//			else cout<<"error in eigenvalue! "<<rs.eigenvalues()(j)<<endl;
		sort(EE_levels.begin(),EE_levels.end());
		EE_levels_all.insert(EE_levels_all.end(),EE_levels.begin(),EE_levels.end());
		EE_levels_storage.push_back(EE_levels);
	}
	sort(EE_levels_all.begin(),EE_levels_all.end());
	vector<double> energy_grid=make_grid(EE_levels_all,200);
	vector<double> integrated_DOS=make_DOS(EE_levels_all,energy_grid);
	for(int i=start;i<end;i++){
		s=make_S(EE_levels_storage[i-start],energy_grid,integrated_DOS);
//		for(int j=0;j<(signed)s.size();j++) cout<<s[j]<<endl;
		rout<<compute_r(s)<<endl;
		s_spacings=spacings(s);
		s_spacings_all.insert(s_spacings_all.end(),s_spacings.begin(),s_spacings.end());
	}
	for(int j=0;j<(signed)EE_levels_all.size();j++) Lout<<EE_levels_all[j]<<" "<<endl;
	for(int j=0;j<(signed)s_spacings_all.size();j++) Sout<<s_spacings_all[j]<<endl;
	rout.close();
	Lout.close();
	Sout.close();	

}
template<class ART>
inline MatrixTFIM<ART>::MatrixTFIM(int _Lx, int _Ly, double _Jx, double _Jy, double _Jz, vector<double> _alphax, vector<double> _alphay, vector<double> _alphaz, int charge): MatrixWithProduct<ART>()
// Constructor

{
	Lx=_Lx; Ly=_Ly; nStates=Lx*Ly;
	Jx=_Jx; Jy=_Jy; Jz=_Jz;
	alphax=_alphax; alphay=_alphay; alphaz=_alphaz;

	make_states(charge);

	this->init_wavefunction(nStates);

}
#endif 


