#pragma once

#include <Eigen/Dense>

#include "definitions.h"

namespace TGSL{
struct quasistaticFEM{
	const TVP& X;
	const IV& mesh;
	TV Dm_inverse;	// flatten, dim = d*d*Ne
	TV measure;	// Area of D_m
	TVP& x;
	const std::vector<std::vector<int>> incident_elements;

	quasistaticFEM(const TVP& X_in,const IV& mesh_in,TVP& x_in, const std::vector<std::vector<int>> incident_elements_in):X(X_in),mesh(mesh_in),x(x_in), incident_elements(incident_elements_in){
		measure.resize(mesh.size()/(d+1)); // dim = Ne
		Dm_inverse.resize(d*d*mesh.size()/(d+1));

		for(sz e=0;e<mesh.size()/(d+1);e++){
			Eigen::Matrix<T,d,d> Dm=Ds(e,X);
			
			T d_factorial=T(1);
			for(sz c=1;c<d;c++)
				d_factorial*=T(c+1);

			measure[e]=(T(1)/T(d_factorial))*Dm.determinant(); // D_m Area = 1/d! * det(D_m)
			
			Eigen::Matrix<T,d,d> Dm_inv=Dm.inverse();
			for(sz r=0;r<d;r++){
				for(sz c=0;c<d;c++){
					Dm_inverse[(d*d)*e + d*r + c]=Dm_inv(r,c);
				}
			}
		}
	}

	Eigen::Matrix<T,d,d> Dme_inv(const sz e)const{
		Eigen::Matrix<T,d,d> Dm_inv;
		for(sz r=0;r<d;r++){
			for(sz c=0;c<d;c++){
				Dm_inv(r,c)=Dm_inverse[(d*d)*e + d*r + c];
			}
		}
		return Dm_inv;
	}

	Eigen::Matrix<T,d,d> Ds(const sz e,const TVP& u) const {
		Eigen::Matrix<T,d,d> result;
		for(sz i=0;i<d;i++){
			for(sz c=0;c<d;c++){
				result(c,i)=u[mesh[3*e+i+1]][c]-u[mesh[3*e]][c];
			}
		}
		return result;
	}

	Eigen::Matrix<T,d,d> F(const sz e) const {
		return Ds(e,x)*Dme_inv(e);
	}

	Eigen::Matrix<T,d,d> dF(const sz e,const TVP& dx) const {
		// F = D_s*D_m^{-1}, dF = d(D_s)*D_m^{-1}, D_m constant
		return Ds(e,dx)*Dme_inv(e);
	}

	template<typename Func>
	T PE(Func psi)const {
		T total=T(0);
		for(sz e=0;e<mesh.size()/(d+1);e++){
			total += measure[e]*psi(F(e),e);
		}
		return total;
	}

/*
	template<typename Func>
	void addInternalForce(Func P,TVP& q)const{
		//(PARALLEL) add internal force

		TVP f_hat(mesh.size());
		tbb::parallel_for(tbb::blocked_range<sz>(0, mesh.size()/(d+1)), [&](const tbb::blocked_range<sz>& range){
			for(sz e = range.begin(); e < range.end(); e++){
				Eigen::Matrix<T,d,d> Fe=F(e);	
				Eigen::Matrix<T,d,d> Pe;
				P(Fe,Pe,e);
				Eigen::Matrix<T,d,d> g=-measure[e]*Pe*(Dme_inv(e).transpose());
				for(sz ie=1;ie<(d+1);ie++){
					for(sz c=0;c<d;c++){
						f_hat[3*e+ie][c]+=g(c,ie-1); //f_1_e, f_2_e
					}
				}
				for(sz c=0;c<d;c++){
					f_hat[3*e][c] = -(g(c,0)+g(c,1)); // f_0_e = -f_1_e - f_2_e
				}
			}
		});
		tbb:parallel_for(tbb::blocked_range<sz>(0, x.size()), [&](const tbb::blocked_range<sz>& range){
			for(sz i = range.begin(); i < range.end(); i++){
				for (sz e = 0; e < incident_elements[i].size(); e++){
					for (sz c=0; c<d; c++){
						q[i][c] += f_hat[incident_elements[i][e]][c];
					}
				}
			}
		});
	}

	template<typename Func>
	void addNegativeInternalForceDifferential(Func dP, const TVP& dx, TVP& ndf)const{
		//(PARALLEL) add internal force
		
		TVP df_hat(mesh.size());
		tbb::parallel_for(tbb::blocked_range<sz>(0, mesh.size()/(d+1)), [&](const tbb::blocked_range<sz>& range){
			for(sz e = range.begin(); e < range.end(); e++){
				Eigen::Matrix<T,d,d> Fe=F(e),dFe=dF(e,dx);	
				Eigen::Matrix<T,d,d> dPe;
				dP(Fe,dFe,dPe,e);
				Eigen::Matrix<T,d,d> dg=measure[e]*dPe*(Dme_inv(e).transpose());				
				for(sz ie=1;ie<(d+1);ie++){
					for(sz c=0;c<d;c++){
						df_hat[3*e+ie][c]+=dg(c,ie-1); //f_1_e, f_2_e
					}
				}
				for(sz c=0;c<d;c++){
					df_hat[3*e][c] = -(dg(c,0)+dg(c,1)); // f_0_e = -f_1_e - f_2_e
				}
			}
		});
		tbb:parallel_for(tbb::blocked_range<sz>(0, x.size()), [&](const tbb::blocked_range<sz>& range){
			for(sz i = range.begin(); i < range.end(); i++){
				for (sz e = 0; e < incident_elements[i].size(); e++){
					for (sz c=0; c<d; c++){
						ndf[i][c] += df_hat[incident_elements[i][e]][c];
					}
				}
			}
		});
	}
*/
	
	template<typename Func>
	void addInternalForce(Func P,TVP& q)const{
		//add internal force
		for(sz e=0;e<mesh.size()/(d+1);e++){
			Eigen::Matrix<T,d,d> Fe=F(e);	
			Eigen::Matrix<T,d,d> Pe;
			P(Fe,Pe,e);
			Eigen::Matrix<T,d,d> g=-measure[e]*Pe*(Dme_inv(e).transpose());
			for(sz ie=1;ie<(d+1);ie++){
				for(sz c=0;c<d;c++){
					q[mesh[3*e+ie]][c]+=g(c,ie-1); //f_1_e, f_2_e
				}
			}
			for(sz c=0;c<d;c++){
				q[mesh[3*e]][c]-=(g(c,0)+g(c,1)); // f_0_e = -f_1_e - f_2_e
			}
		}
	}

	template<typename Func>
	void addNegativeInternalForceDifferential(Func dP, const TVP& dx, TVP& ndf)const{
		//add internal force
		for(sz e=0;e<mesh.size()/(d+1);e++){
			Eigen::Matrix<T,d,d> Fe=F(e),dFe=dF(e,dx);	
			Eigen::Matrix<T,d,d> dPe;
			dP(Fe,dFe,dPe,e);
			Eigen::Matrix<T,d,d> dg=measure[e]*dPe*(Dme_inv(e).transpose());
			for(sz ie=1;ie<(d+1);ie++){
				for(sz c=0;c<d;c++){
					ndf[mesh[3*e+ie]][c]+=dg(c,ie-1);
				}
			}
			for(sz c=0;c<d;c++){
				ndf[mesh[3*e]][c]-=(dg(c,0)+dg(c,1));
			}
		}
	}

	

	template<typename Func1, typename Func2>
	void refinementTestdP(Func1 P,Func2 dP,const sz num_ref,const T rand_scale)const {
		TVP ndf(x.size()),x_save(x.size()),delta_x(x.size());
		x_save=x;
		
		srand(1);
		for(sz i=0;i<x.size();i++){
			for(sz c=0;c<d;c++){
				delta_x[i][c]=T(.5)*rand_scale*(T(2)*T(rand())/T(RAND_MAX) - 1);
				ndf[i][c]=T(0);
			}
		}

		addNegativeInternalForceDifferential(dP, delta_x, ndf);

		TV lg_eps(num_ref),lg_err(num_ref);
		for(sz r=0;r<num_ref;r++){
			T eps=T(.5)/T(r+1);
			lg_eps[r]=log(eps);

			TVP fp(x.size());
			for(sz i=0;i<x.size();i++){
				for(sz c=0;c<d;c++){
					x[i][c]=x_save[i][c]+T(.5)*eps*delta_x[i][c];
					fp[i][c]=T(0);
				}
			}
			addInternalForce(P,fp);

			TVP fm(x.size());
			for(sz i=0;i<x.size();i++){
				for(sz c=0;c<d;c++){
					x[i][c]=x_save[i][c]-T(.5)*eps*delta_x[i][c];
					fm[i][c]=T(0);
				}
			}
			addInternalForce(P,fm);

			T n_2=T(0);
			for(sz i=0;i<x.size();i++){
				for(sz c=0;c<d;c++){
					T ec=ndf[i][c]+(fp[i][c]-fm[i][c])/eps;
					n_2+=ec*ec;
				}
			}
			n_2=sqrt(n_2);

			lg_err[r]=log(n_2);

		}

		x=x_save;

		writeCSV(lg_eps,lg_err,"loglog.csv");
		
	}

	template<typename Func1, typename Func2>
	void refinementTestP(Func1 psi,Func2 P,const sz num_ref,const T rand_scale)const {
		TVP force(x.size()),x_save(x.size()),delta_x(x.size());
		x_save=x;
		
		srand(1);
		for(sz i=0;i<x.size();i++){
			for(sz c=0;c<d;c++){
				delta_x[i][c]=T(.5)*rand_scale*(T(2)*T(rand())/T(RAND_MAX) - 1);
				force[i][c]=T(0);
			}
		}

		addInternalForce(P,force);
		
		T dpsi=T(0);
		for(sz i=0;i<x.size();i++){
			for(sz c=0;c<d;c++){
				dpsi+=force[i][c]*delta_x[i][c];		
			}
		}

		TV lg_eps(num_ref),lg_err(num_ref);
		for(sz r=0;r<num_ref;r++){
			T eps=T(.5)/T(r+1);
			lg_eps[r]=log(eps);

			for(sz i=0;i<x.size();i++){
				for(sz c=0;c<d;c++){
					x[i][c]=x_save[i][c]+T(.5)*eps*delta_x[i][c];
				}
			}
			T psi_plus=PE(psi);

			for(sz i=0;i<x.size();i++){
				for(sz c=0;c<d;c++){
					x[i][c]=x_save[i][c]-T(.5)*eps*delta_x[i][c];
				}
			}
			T psi_minus=PE(psi);

			lg_err[r]=log(abs((psi_plus-psi_minus)/eps + dpsi));

		}

		x=x_save;

		writeCSV(lg_eps,lg_err,"loglog.csv");
	}

	void zeroDirichlet(const IV& dirichlet_nodes,TVP& u)const{
		for(sz i=0;i<dirichlet_nodes.size();i++){
			for(sz c=0;c<d;c++){
				u[dirichlet_nodes[i]][c]=T(0);
			}
		}
	}

	template<typename Func0, typename Func1,typename Func2, typename Func3, typename Func4>
	void solve(Func0 preItLambda, Func1 postItLambda, Func2 P,Func3 dP,Func4 addExternalForce,const IV& dirichlet_nodes,const sz max_it,const T tol=T(1e-10))const{
		TGSLAssert(x.size()==X.size(),"quasistaticFEM: deformed and undeformed vertices have mismatched sizes.");
		TVP q(x.size()),delta_x(x.size());

		postItLambda(0);	

		for(sz it=0;it<max_it;it++){
			preItLambda(it);	//compute dPdF

			//set q to zero, dim(q) = N*N*d
			for(sz i=0;i<q.size();i++){
				for(sz c=0;c<d;c++){
					q[i][c]=T(0);
					delta_x[i][c]=T(0);
				}
			}
			//add forces for residual
			addExternalForce(q);
			addInternalForce(P,q);
			zeroDirichlet(dirichlet_nodes,q);

			auto multiplyA = [&dirichlet_nodes,&dP,this](TVP& y, const TVP& x){
				for(sz i=0;i<y.size();i++){
					for(sz c=0;c<d;c++){
						y[i][c]=T(0);
					}
				}
				addNegativeInternalForceDifferential(dP, x, y);
				zeroDirichlet(dirichlet_nodes,y);
			};

			auto dotProduct = [&dirichlet_nodes](const TVP& y, const TVP& x){
				T result=T(0);
				for(sz i=0;i<x.size();i++){
					for(sz c=0;c<d;c++){
						result+=y[i][c]*x[i][c];
					}
				}
				return result;
			};
			
			auto AXPY = [&dirichlet_nodes](TVP& x, const T s, const TVP& y){
				for(sz i=0;i<x.size();i++){
					for(sz c=0;c<d;c++){
						x[i][c]=x[i][c] + s*y[i][c];
					}
				}
			};
			
			auto scale = [&dirichlet_nodes](TVP& x, const T s){
				for(sz i=0;i<x.size();i++){
					for(sz c=0;c<d;c++){
						x[i][c]=s*x[i][c];
					}
				}
			};
			
			auto set = [&dirichlet_nodes](TVP& x, const TVP& y){
				x=y;
			};

			auto residual = [&q,&multiplyA,&AXPY,&dotProduct,&dirichlet_nodes,this](TVP& x){
				TVP r(x.size());
				multiplyA(r,x);
				AXPY(r,T(-1),q);
				zeroDirichlet(dirichlet_nodes,r);
				T norm_2=dotProduct(r,r);
				return sqrt(norm_2);
			};

			T res=residual(delta_x);
			std::cout<<"Newton residual = " << res << std::endl;
			if(res<tol){
				std::cout << "Newton converged in " << it << " iterations." << std::endl;
				return;
			}
			lanczosCG<T>(multiplyA, dotProduct, AXPY, scale, set, delta_x,q,x.size());	//A(delta_x) = q
			std::cout<<"Linearized system residual after CG = " << residual(delta_x) << std::endl;

			for(sz i=0;i<delta_x.size();i++){
				for(sz c=0;c<d;c++){
					x[i][c]+=delta_x[i][c];
				}
			}

			postItLambda(it+1);
		}
	}
};
}