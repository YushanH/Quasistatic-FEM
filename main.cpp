#include <iostream>
#include <string>
#include <Eigen/Dense>
#include <ctime>
#include <time.h>
#include <stdio.h>
#include "definitions.h"
#include "grid.h"
#include "mesh.h"
#include "io.h"
#include "krylov.h"
#include "mathTools.h"

#include "quasistaticFEM.h"

using namespace TGSL;

void squareDomain(const Grid& grid, IV& tri_mesh, TVP2& x){
	x.resize(grid.size());	//N*N
	tri_mesh.resize(6*(grid.gridN[0]-1)*(grid.gridN[0]-1));
	for(size_t ii=0;ii<grid.size();ii++)
		x[ii]=grid.node(ii);

	int count=0;
	for(size_t i=0;i<grid.gridN[0]-1;i++){
		for(size_t j=0;j<grid.gridN[1]-1;j++){
			size_t ij=grid.flatIndex({int(i),int(j)}); 	//i*N+j
			size_t ij1=grid.flatIndex({int(i),int(j+1)});	//northeast neighbors 
			size_t i1j=grid.flatIndex({int(i+1),int(j)});
			size_t i1j1=grid.flatIndex({int(i+1),int(j+1)});
			Vector3I tri0={int(ij),int(i1j1),int(ij1)}; //triangle mesh, counterclock wise
			Vector3I tri1={int(ij),int(i1j),int(i1j1)};
			tri_mesh[6*count]=int(ij);tri_mesh[6*count+1]=int(i1j1);tri_mesh[6*count+2]=int(ij1);	//load triangle indices into mesh (3*n*n)
			tri_mesh[6*count+3]=int(ij);tri_mesh[6*count+4]=int(i1j);tri_mesh[6*count+5]=int(i1j1);
			count++;
		}
	}
}

int main(void){

	int N=50;
	Index gridN={N,N};
	T dx=T(1)/T(N-1);

	Grid grid(gridN,dx,{T(0),T(0)});	//Initialize grid: N*N points with dx starting at (0,0)
	/* 		3*3 grid:
	2 5 8
	1 4 7
	0 3 6
	*/

	grid.setInterp(Grid::linear);	//Use linear interpolation
	IV tri_mesh;	//Integer vector for mesh indices
	TVP2 X,x;		//(double, double) vector
	squareDomain(grid,tri_mesh,X);
	x.resize(X.size());
	TV dPdF(tri_mesh.size()*d*d*d*d);

	TVI2 boundary_mesh;
	computeBoundaryMesh(tri_mesh,boundary_mesh);

	std::vector<std::vector<int>> incident_elements;
	computeIncidentElements(3,tri_mesh,incident_elements);		// e.g. incident_element[i]: vector of indices of points that is adjacent to point i


	T R_0=T(1);
	Vector g={T(0),T(-10)};
	T mu=T(10);
	T lambda=T(100);

	//initialize the left and right of the domain as Dirichlet
	IV dirichlet_nodes(2*N);
	for(sz n=0;n<N;n++){
		//dirichlet_nodes[2*n]=grid.gridN[1]*n;
		//dirichlet_nodes[2*n+1]=grid.gridN[1]*n + (grid.gridN[1]-1);
		dirichlet_nodes[2*n]=n; 
		dirichlet_nodes[2*n+1]=grid.gridN[1]*(grid.gridN[0]-1) + n;
	}
	
	//define elastic potential and stress
	//linear elasticity model
	/*auto psi = [&mu,&lambda](const Eigen::Matrix<T,d,d>& F){
		Eigen::Matrix<T,d,d> I,epsilon;I.setIdentity();
		epsilon=T(.5)*(F+F.transpose())-I;

		T total=T(0);
		for(sz alpha=0;alpha<d;alpha++){
			for(sz beta=0;beta<d;beta++){
				total+=mu*epsilon(alpha,beta)*epsilon(alpha,beta);
			}
		}
	
		total += T(.5)*lambda*(epsilon.trace())*(epsilon.trace());
		return total;
	};
	auto P = [&mu,&lambda](const Eigen::Matrix<T,d,d>& F,Eigen::Matrix<T,d,d>& P_out){
		Eigen::Matrix<T,d,d> I,epsilon;I.setIdentity();
		epsilon=T(.5)*(F+F.transpose())-I;
		P_out=T(2)*mu*epsilon+lambda*epsilon.trace()*I;
	};
	auto dP = [&mu,&lambda](const Eigen::Matrix<T,d,d>& F,const Eigen::Matrix<T,d,d>& dF,Eigen::Matrix<T,d,d>& dP_out){
		Eigen::Matrix<T,d,d> I,epsilon;I.setIdentity();
		epsilon=T(.5)*(dF+dF.transpose());
		dP_out=T(2)*mu*epsilon+lambda*epsilon.trace()*I;
	};*/


	//St-Venant Kirchoff model
	/*auto psi = [&mu,&lambda](const Eigen::Matrix<T,d,d>& F){
		Eigen::Matrix<T,d,d> I,E;I.setIdentity();
		E=T(.5)*((F.transpose())*F-I);

		T total=T(0);
		for(sz alpha=0;alpha<d;alpha++){
			for(sz beta=0;beta<d;beta++){
				total+=mu*E(alpha,beta)*E(alpha,beta);
			}
		}
	
		total += T(.5)*lambda*(E.trace())*(E.trace());
		return total;
	};
	auto P = [&mu,&lambda](const Eigen::Matrix<T,d,d>& F,Eigen::Matrix<T,d,d>& P_out){
		Eigen::Matrix<T,d,d> I,E;I.setIdentity();
		E=T(.5)*((F.transpose())*F-I);
		P_out=F*(T(2)*mu*E+lambda*E.trace()*I);
	};
	auto dP = [&mu,&lambda](const Eigen::Matrix<T,d,d>& F,const Eigen::Matrix<T,d,d>& dF,Eigen::Matrix<T,d,d>& dP_out){
		Eigen::Matrix<T,d,d> I,E,dE;I.setIdentity();
		E=T(.5)*((F.transpose())*F-I);
		dE=T(.5)*((dF.transpose())*F+(F.transpose())*dF);
		dP_out=dF*(T(2)*mu*E+lambda*E.trace()*I) + F*(T(2)*mu*dE+lambda*dE.trace()*I);
	};*/

	auto psi = [&mu,&lambda](const Eigen::Matrix<T,d,d>& F,const sz ie){
		Eigen::Matrix<T,d,d> R,S,E;
		T c,s,s00,s10,s11;
		polarDecomposition(F(0,0),F(1,0),F(0,1),F(1,1),c,s,s00,s10,s11);
		/*
		C = [c -s] 
			[s  c]		
		*/
		R(0,0)=c;R(0,1)=-s;
		R(1,0)=s;R(1,1)=c;
		S(0,0)=s00;S(0,1)=s10;
		S(1,0)=s10;S(1,1)=s11;
		E=F-R;
		T J = F.determinant();

		T total=T(0);
		for(sz alpha=0;alpha<d;alpha++){
			for(sz beta=0;beta<d;beta++){
				total+=mu*E(alpha,beta)*E(alpha,beta);
			}
		}
	
		total += T(.5)*lambda*(J-T(1))*(J-T(1));
		return total;
	};
	auto P = [&mu,&lambda](const Eigen::Matrix<T,d,d>& F,Eigen::Matrix<T,d,d>& P_out,const sz ie){
		Eigen::Matrix<T,d,d> R,S,E;
		T c,s,s00,s10,s11;
		polarDecomposition(F(0,0),F(1,0),F(0,1),F(1,1),c,s,s00,s10,s11);
		R(0,0)=c;R(0,1)=-s;
		R(1,0)=s;R(1,1)=c;
		S(0,0)=s00;S(0,1)=s10;
		S(1,0)=s10;S(1,1)=s11;
		T J = F.determinant();
		Eigen::Matrix<T,d,d> JFinvT;
		JFinvT(0,0)=F(1,1);JFinvT(0,1)=-F(1,0);
		JFinvT(1,0)=-F(0,1);JFinvT(1,1)=F(0,0);

		P_out=T(2)*mu*(F-R)+lambda*(J-T(1))*JFinvT;
	};
	auto dP = [&mu,&lambda](const Eigen::Matrix<T,d,d>& F,const Eigen::Matrix<T,d,d>& dF,Eigen::Matrix<T,d,d>& dP_out,const sz ie){
		Eigen::Matrix<T,d,d> R,S,e,RTdF,dR;
		e(0,0)=T(0);e(0,1)=-T(1);
		e(1,0)=T(1);e(1,1)=T(0);

		T c,s,s00,s10,s11;
		polarDecomposition(F(0,0),F(1,0),F(0,1),F(1,1),c,s,s00,s10,s11);		
		R(0,0)=c;R(0,1)=-s;
		R(1,0)=s;R(1,1)=c;
		S(0,0)=s00;S(0,1)=s10;
		S(1,0)=s10;S(1,1)=s11;

		RTdF=R.transpose()*dF;
		T contract=(RTdF(1,0)-RTdF(0,1))/S.trace();
		dR=contract*R*e;
		
		T J = F.determinant();
		Eigen::Matrix<T,d,d> JFinvT,dJFinvT;
		JFinvT(0,0)=F(1,1);JFinvT(0,1)=-F(1,0);
		JFinvT(1,0)=-F(0,1);JFinvT(1,1)=F(0,0);
		dJFinvT(0,0)=dF(1,1);dJFinvT(0,1)=-dF(1,0);
		dJFinvT(1,0)=-dF(0,1);dJFinvT(1,1)=dF(0,0);

		T dJ=T(0);
		for(sz alpha=0;alpha<d;alpha++){
			for(sz beta=0;beta<d;beta++){
				dJ+=JFinvT(alpha,beta)*dF(alpha,beta);
			}
		}

		dP_out=T(2)*mu*(dF-dR) + lambda*(dJ*JFinvT+(J-T(1))*dJFinvT);
	};
	auto dPDefinite = [&mu,&lambda,&dPdF](const Eigen::Matrix<T,d,d>& F,const Eigen::Matrix<T,d,d>& dF,Eigen::Matrix<T,d,d>& dP_out,const sz ie){
		

		Eigen::Matrix<int,2,2> flat;
		flat(0,0)=0;flat(0,1)=3;
		flat(1,0)=2;flat(1,1)=1;

		Eigen::Matrix<T,4,1> dF_flat,dP_flat;
		for(int alpha=0;alpha<d;alpha++){
			for(int beta=0;beta<d;beta++){
				dF_flat(flat(alpha,beta))=dF(alpha,beta);
			}
		}

		Eigen::Matrix<T,4,4> A;

		for(sz i=0;i<4;i++)
			for(sz j=0;j<4;j++)
				A(i,j)=dPdF[16*ie + 4*i + j];

		dP_flat=A*dF_flat;

		for(int alpha=0;alpha<d;alpha++){
			for(int beta=0;beta<d;beta++){
				dP_out(alpha,beta)=dP_flat(flat(alpha,beta));
			}
		}
	};

	quasistaticFEM qs(X,tri_mesh,x, incident_elements);
	x=X;
	
	//gravity
	auto addExternalForce =[&R_0,&g,&tri_mesh,&qs](TVP& f){
		for(sz e=0;e<tri_mesh.size()/(d+1);e++){
			T measure=qs.measure[e];
			for(sz ie=0;ie<(d+1);ie++){
				for(sz c=0;c<d;c++){
					f[tri_mesh[3*e+ie]][c]+=R_0*g[c]*measure/T(d+1);
				}
			}
		}
	};
	
	// phi_D on dirchlet boundary
	for(sz i=0;i<x.size();i++){
		x[i][0]*=T(1.15);
	}	
	//qs.refinementTestdP(P,dPDefinite,10,dx);

	auto computedPdF = [&qs,&dPdF,&mu,&lambda](const sz ie){
		// &dPdF flattened Ne*4*4
		Eigen::Matrix<T,d,d> F=qs.F(ie);

		Eigen::Matrix<T,d,d> U,sigma,V;
		T cu,su,cv,sv,s0,s1;
		SVD(F(0,0),F(1,0),F(0,1),F(1,1),cu,su,cv,sv,s0,s1);
		U(0,0)=cu;U(0,1)=-su;
		U(1,0)=su;U(1,1)=cu;
		V(0,0)=cv;V(0,1)=-sv;
		V(1,0)=sv;V(1,1)=cv;
		sigma(0,0)=s0;sigma(0,1)=T(0);
		sigma(1,0)=T(0);sigma(1,1)=s1;
		T tr_S=s0+s1;

		Eigen::Matrix<int,2,2> flat;
		flat(0,0)=0;flat(0,1)=3;
		flat(1,0)=2;flat(1,1)=1;

		Eigen::Matrix<T,4,4> block_dPdF,B,A;
		Eigen::Matrix<T,4,1> dF_flat,dP_flat;
		for(int alpha=0;alpha<d;alpha++){
			for(int beta=0;beta<d;beta++){
				for(int gamma=0;gamma<d;gamma++){
					for(int delta=0;delta<d;delta++){
						B(flat(alpha,beta),flat(gamma,delta))=U(gamma,alpha)*V(delta,beta);
					}
				}
			}
		}

		T a = T(2)*mu+lambda*s1*s1, b = lambda*(T(2)*s1*s0-T(1)), c = T(2)*mu+lambda*s0*s0, f = T(2)*mu*(T(1)-T(1)/(s0+s1)), e = T(2)*mu/(s0+s1) - lambda*(s0*s1-T(1));

		block_dPdF.setZero();
		Eigen::Matrix<T,2,2> b1,b2,eig_vect,eig_val;
		b1(0,0)=a;b1(0,1)=b;
		b1(1,0)=b;b1(1,1)=c;
		b2(0,0)=f;b2(0,1)=e;
		b2(1,0)=e;b2(1,1)=f;

		Jacobi(b1,eig_val,eig_vect);	// Eigenvalue decomposition by (sym + schur decomp)

		if(eig_val(0,0)<T(0))
			eig_val(0,0)=T(0);
		if(eig_val(1,1)<T(0))
			eig_val(1,1)=T(0);	//correct to be pos.definite
		b1=eig_vect*eig_val*(eig_vect.transpose());

		Jacobi(b2,eig_val,eig_vect);	

		if(eig_val(0,0)<T(0))
			eig_val(0,0)=T(0);
		if(eig_val(1,1)<T(0))
			eig_val(1,1)=T(0);
		b2=eig_vect*eig_val*(eig_vect.transpose());

		block_dPdF(0,0)=b1(0,0);block_dPdF(0,1)=b1(0,1);
		block_dPdF(1,0)=b1(1,0);block_dPdF(1,1)=b1(1,1);
		block_dPdF(2,2)=b2(0,0);block_dPdF(2,3)=b2(0,1);
		block_dPdF(3,2)=b2(1,0);block_dPdF(3,3)=b2(1,1);

		A=(B.transpose())*(block_dPdF*B);

		for(sz i=0;i<4;i++)
			for(sz j=0;j<4;j++)
				dPdF[16*ie + 4*i + j]=A(i,j);	
	};

	auto preIt = [&tri_mesh,&computedPdF](const sz it){
		for(sz e=0;e<sz(tri_mesh.size())/sz(d+1);e++)
			computedPdF(e);
	};

	auto postIt = [&tri_mesh,&x](const sz it){
		//"output/sim_frame_" + std::to_string(it+1) + ".obj"
		std::string filename = "output_sim_frame_" + std::to_string(it+1) + ".obj";
		writeOBJ(x,tri_mesh,filename);
		//std::cout << filename << std::endl;
	};
	time_t time_start = time(NULL);
	qs.solve(preIt,postIt,P,dPDefinite,addExternalForce,dirichlet_nodes,10);
	time_t time_end = time(NULL);
	std::cout << "N = " << N << ", Time used = " << difftime(time_end, time_start) << std::endl;
	return 1;
}