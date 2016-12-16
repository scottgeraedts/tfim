
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

 public:

	int Lx,Ly,Nspins;
	int next(int i, int dir);
	vector<int> states;
	int bitflip(int,int);
	void MultMv(ART* v, ART* w);
	void makeSparse();
//	Eigen::Matrix<ART,-1,1> MultEigen(Eigen::Matrix<ART,-1,1>);

	void make_disorder(int seed);
	void make_states(int charge);
  	int nSpins(){return Nspins;}
  	
	MatrixTFIM(int _Lx, int _Ly, double _Jx, double _Jy, double _Jz, vector<double> _alphax, vector<double> _alphay, vector<double> _alphaz, int charge);
	MatrixTFIM(){ };
	void plot_states(int start, int end);
	void entanglement_spacings(int start, int end, int to_trunc,double label, string folder, int kmin, int kmax);	
	void energy_spacings(double label);
	
}; // MatrixTFIM.

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
	for(int i=0;i<1<<Nspins;i++)
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
		for(int i=0; i<Nspins; i++){
			if(bittest(states[in],i)) sign=1;
			else sign=-1;

			//hx, hy
			if(alphax[i]!=0 || alphay[i]!=0){
				#ifdef USE_COMPLEX
				w[lookup_flipped(in,states, 1, i)]+=0.5*(alphax[i] + sign*alphay[i]*complex<double>(0,1.))*v[in];
				//w[states[in]^1<<i]+=0.5*(alphax[i] + sign*alphay[i]*complex<double>(0,1.))*v[in];
				#else
				w[lookup_flipped(in,states,1, i)]+=0.5*(alphax[i] )*v[in];
				//w[states[in]^1<<i]+=0.5*(alphax[i] )*v[in];
				#endif
			}
			//hz
			countH+=0.5*sign*alphaz[i];
		}
		
		//diagonal elements
		countJ=0;
		for(int i=0; i<Nspins;i++){
			if(bittest(states[in],i) == bittest(states[in],next(i,0)) ) sign=1;
			else sign=-1;
			//Jx, Jy
			if(sign==-1 && (Jx!=0 || Jy!=0)) w[ lookup_flipped(in , states ,2, i,next(i,0))  ]+=(Jx+Jy)*v[in]*0.25;
			if(sign==1 && (Jx!=0 || Jy!=0) && Jx!=Jy ) w[ lookup_flipped(in , states ,2, i,next(i,0))  ]+=(Jx-Jy)*v[in]*0.25;
			//if(Jx!=0 || Jy!=0) w[ states[in] ^ ( (1<<i) + (1<<next(i,0)) ) ]+=(Jx-sign*Jy)*v[in]*0.25;
			//Jz
			countJ+=sign;
			if(Ly>1){
				if(bittest(states[in],i) == bittest(states[in],next(i,1)) ) sign=1;
				else sign=-1;
				//Jx, Jy
				if(Jx!=0 || Jy!=0) w[ states[in] ^ ( (1<<i) + (1<<next(i,1)) ) ]+=(Jx-sign*Jy)*v[in]*0.25;
				//Jz
				countJ+=sign;			
			}			
			
		}
		w[in]+=(countJ*Jz*0.25+countH)*v[in];	
	}
} //  MultMv.

//similar to above but builds a sparse hamiltonian, hopefully way faster than above
template<class ART>
void MatrixTFIM<ART>::makeSparse(){
	int countJ=0;
	int sign;
	double countH;
	vector<Eigen::Triplet<ART> > triplets;
	for(int in=0; in<this->ncols(); in++){
		//off-diagonal elements
		countH=0.;
		for(int i=0; i<Nspins; i++){
			if(bittest(states[in],i)) sign=1;
			else sign=-1;

			//hx, hy
			if(alphax[i]!=0 || alphay[i]!=0){
				#ifdef USE_COMPLEX
				triplets.push_back(Eigen::Triplet<ART>(in,lookup_flipped(in,states, 1, i),0.5*(alphax[i] + sign*alphay[i]*complex<double>(0,1.))));
				//w[states[in]^1<<i]+=0.5*(alphax[i] + sign*alphay[i]*complex<double>(0,1.))*v[in];
				#else
				triplets.push_back(Eigen::Triplet<ART>(in,lookup_flipped(in,states,1, i),0.5*(alphax[i] )));
				//w[states[in]^1<<i]+=0.5*(alphax[i] )*v[in];
				#endif
			}
			//hz
			countH+=0.5*sign*alphaz[i];
		}
		
		//diagonal elements
		countJ=0;
		for(int i=0; i<Nspins;i++){
			if(bittest(states[in],i) == bittest(states[in],next(i,0)) ) sign=1;
			else sign=-1;
			//Jx, Jy
			if(sign==-1 && (Jx!=0 || Jy!=0)) triplets.push_back(Eigen::Triplet<ART>(in,lookup_flipped(in , states ,2, i,next(i,0)),(Jx+Jy)*0.25));
			if(sign==1 && (Jx!=0 || Jy!=0) && Jx!=Jy ) triplets.push_back(Eigen::Triplet<ART>(in,lookup_flipped(in , states ,2, i,next(i,0)) , (Jx-Jy)*0.25) );
			//if(Jx!=0 || Jy!=0) w[ states[in] ^ ( (1<<i) + (1<<next(i,0)) ) ]+=(Jx-sign*Jy)*v[in]*0.25;
			//Jz
			countJ+=sign;
			if(Ly>1){
				if(bittest(states[in],i) == bittest(states[in],next(i,1)) ) sign=1;
				else sign=-1;
				//Jx, Jy
		//		if(Jx!=0 || Jy!=0) w[ states[in] ^ ( (1<<i) + (1<<next(i,1)) ) ]+=(Jx-sign*Jy)*v[in]*0.25;
				//Jz
				countJ+=sign;			
			}			
			
		}
		triplets.push_back(Eigen::Triplet<ART>(in,in,(countJ*Jz*0.25+countH)));	
	}
	this->sparse=Eigen::SparseMatrix<ART>(this->nrows(),this->nrows());
	this->sparse.setFromTriplets(triplets.begin(),triplets.end());
	this->sparse.makeCompressed();
} //  makeSparse.

//template<class ART>
//Eigen::Matrix<ART,-1,1> MatrixTFIM<ART>::MultEigen(Eigen::Matrix<ART,-1,1> v){
//	int countJ=0;
//	Eigen::Matrix<ART,-1,1> w=Eigen::Matrix<ART,-1,1>::Zero(this->ncols());
//	for(int in=0; in<this->ncols(); in++){
//		//off-diagonal elements
//		for(int i=0; i<Nspins; i++){
//			w(bitflip(in,i))+=alpha[i]*v(in);
//		}
//		//diagonal elements
//		countJ=0;
//		for(int i=0; i<Nspins-1;i++){
//			if(bittest(in,i) == bittest(in,i+1) ) countJ--;
//			else countJ++;
//		}
////		if (bittest(in, Nspins-1) == bittest(in,0) ) countJ--;
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
	alphax=vector<double>(Nspins,0);
	alphay=vector<double>(Nspins,0);
	alphaz=vector<double>(Nspins,0);
	ifstream datfile;
	datfile.open("nicrans");
	for(int i=0; i<Nspins; i++){
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
void MatrixTFIM<ART>::energy_spacings(double label=-100){
	sort(this->eigvals.begin(),this->eigvals.end());
	ofstream energies;
	energies.open("energies");
	for(int i=0;i<(signed)this->eigvals.size();i++) energies<<this->eigvals[i]<<endl;
	energies.close();
	vector<double> s=unfoldE(this->eigvals,100);
	ofstream sout,rout;
	vector<double> energy_spacings=spacings(s);
	stringstream filename;
	filename<<"energy_spacings";
	if(label!=-100) filename<<label;
	sout.open(filename.str().c_str());
	for(int i=0;i<energy_spacings.size();i++) sout<<energy_spacings[i]<<endl;
//	for(int i=0;i<s.size();i++) sout<<s[i]<<endl;
	sout.close();
	filename.str("");
	filename<<"energy_r";
	if(label!=-100) filename<<label;
	rout.open(filename.str().c_str());
	rout<<compute_r(s)<<endl;
	rout.close();
}
template<class ART>
void MatrixTFIM<ART>::entanglement_spacings(int start, int end, int to_trunc,double label=-100,string folder=".",int kmin_t=-1, int kmax_t=-0){

	vector<double> s;
	vector<double> EE_levels,s_spacings;
	vector<double> EE_levels_all,s_spacings_all;
	vector< vector<double> > EE_levels_storage;
//	int rhosize=this->ee_setup(0,Nspins/2,states);
//	Eigen::Matrix<ART,-1,-1> rho=Eigen::Matrix<ART,-1,-1>::Zero(rhosize,rhosize);
//	Eigen::SelfAdjointEigenSolver<Eigen::Matrix<ART,-1,-1> > rs;	
	ofstream Lout,Sout,rout,entropy_out;
	stringstream filename;
	Lout.open("EE_levels");
	filename<<folder<<"/spacings";
	if(label!=-100) filename<<label;
	Sout.open(filename.str().c_str());
	filename.str("");
	filename<<folder<<"/r";
	if(label!=-100) filename<<label;
	rout.open(filename.str().c_str());
	filename.str("");
	filename<<folder<<"/entropy";
	if(label!=-100) filename<<label;
	entropy_out.open(filename.str().c_str());

	double entropy_sum;
	vector<double> from_svd;
	int kmin,kmax;
	if(this->nrows()==1<<Nspins){
		kmin=-1; kmax=0;
	}else{
		kmin=kmin_t; kmax=kmax_t;
//		kmin=4; kmax=5;
	}
	for(int i=start;i<end;i++){
//		rho=Eigen::Matrix<ART,-1,-1>::Zero(rhosize,rhosize);
//		this->ee_compute_rho(this->eigvecs[i],rho,states);
//		rs.compute(rho);
		entropy_sum=0;
		for(int k=kmin;k<kmax;k++){
			from_svd=this->entanglement_spectrum_SVD(this->eigvecs[i],states,to_trunc,k);
//		cout<<"raw eigenvalues "<<i<<endl;
//		for(int j=0;j<rhosize;j++) cout<<rs.eigenvalues()(j)<<" "<<from_svd[j]<<endl;
			EE_levels.clear();
			for(int j=0;j<(signed)from_svd.size();j++) 
				if(from_svd[j]>0.) EE_levels.push_back(-log(from_svd[j]));
//			else cout<<"error in eigenvalue! "<<rs.eigenvalues()(j)<<endl;
			sort(EE_levels.begin(),EE_levels.end());
			for(unsigned int j=0;j<EE_levels.size();j++) entropy_sum+=exp(-EE_levels[j])*EE_levels[j];
			EE_levels_all.insert(EE_levels_all.end(),EE_levels.begin(),EE_levels.end());
			EE_levels_storage.push_back(EE_levels);
		}
		entropy_out<<entropy_sum<<endl;
	}
	entropy_out.close();
	sort(EE_levels_all.begin(),EE_levels_all.end());
	vector<double> energy_grid=make_grid(EE_levels_all,200);
	vector<double> integrated_DOS=make_DOS(EE_levels_all,energy_grid);

	//print DOS
	filename.str("");
	if(kmin_t==-1) filename<<folder<<"/dos";
	else filename<<folder<<"/dos_single";
	if(label!=-100) filename<<label;
	ofstream dosout;
	dosout.open(filename.str().c_str());
	for(int i=0;i<(signed)integrated_DOS.size();i++){
		dosout<<energy_grid[i]<<" ";
		if(i==integrated_DOS.size()-1) dosout<<0.<<" ";
		else dosout<<(integrated_DOS[i+1]-integrated_DOS[i])/(energy_grid[i+1]-energy_grid[i])<<" ";
		dosout<<integrated_DOS[i]<<endl;
	}
	dosout.close();
/*
	for(auto it=EE_levels_storage.begin();it!=EE_levels_storage.end();++it){
		cout<<it->size()<<endl;
		s=make_S(*it,energy_grid,integrated_DOS);
//		for(int j=0;j<(signed)s.size();j++) cout<<s[j]<<endl;
		rout<<compute_r(s)<<endl;
		s_spacings=spacings(s);
		s_spacings_all.insert(s_spacings_all.end(),s_spacings.begin(),s_spacings.end());
	}
	for(int j=0;j<(signed)EE_levels_storage.size();j++) 
		for(int k=0;k<(signed)EE_levels_storage[j].size();k++) Lout<<EE_levels_storage[j][k]<<" "<<endl;
	for(int j=0;j<(signed)s_spacings_all.size();j++) Sout<<s_spacings_all[j]<<endl;
*/	rout.close();
	Lout.close();
	Sout.close();	

}

template<class ART>
void MatrixTFIM<ART>::plot_states(int start, int end){
	ofstream sout("eigenstates");
	ART max;
	int maxindex;
	vector<int> largest;
	for(int i=start; i<end; i++){
		//find the largest element in the eigenvectro
		max=0;
		for(int j=0;j<this->nrows();j++){
			if(abs(this->eigvecs[i][j])>max){
				max=abs(this->eigvecs[i][j]); 
				maxindex=j;
			}
		}
		largest=bitset_to_pos(this->states[maxindex],Nspins);
		for(int j=0;j<this->nrows();j++)
			sout<<distance_calc(largest,bitset_to_pos(this->states[j],Nspins))<<" "<<this->eigvecs[i][j]<<" "<<(bitset<12>)this->states[j]<<" "<<(bitset<12>)this->states[maxindex]<<endl;
	}
	sout.close();

}
template<class ART>
inline MatrixTFIM<ART>::MatrixTFIM(int _Lx, int _Ly, double _Jx, double _Jy, double _Jz, vector<double> _alphax, vector<double> _alphay, vector<double> _alphaz, int charge): MatrixWithProduct<ART>()
// Constructor

{
	Lx=_Lx; Ly=_Ly; Nspins=Lx*Ly;
	Jx=_Jx; Jy=_Jy; Jz=_Jz;
	alphax=_alphax; alphay=_alphay; alphaz=_alphaz;

	make_states(charge);

	this->init_wavefunction(Nspins);

}
#endif 


