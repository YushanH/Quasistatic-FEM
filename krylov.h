#include <math.h>
#include <iostream>

#include "definitions.h"

namespace TGSL{
	template<typename T,typename Func1,typename Func2,typename Func3,typename Func4,typename Func5,typename Vector>
	void lanczosCG(Func1 multiplyA, Func2 dotProduct, Func3 AXPY, Func4 scale, Func5 set, Vector& x, const Vector& b,const int max_it){
		Vector v(x.size()),q(x.size()),q_old(x.size()),c(x.size());
		T beta=sqrt(T(dotProduct(b,b)));
		set(q,b);
		scale(q,T(1)/beta);
		multiplyA(v,q);
		T alpha=dotProduct(v,q);	
		T d=alpha;
		T p=beta/d;
		set(c,q);
		set(x,c);
		scale(x,p);
		AXPY(v,-alpha,q);

		for(int it=1;it<max_it;it++){
			beta=sqrt(T(dotProduct(v,v)));
			T mu=beta/d;
			set(q_old,q);
			set(q,v);
			scale(q,T(1)/beta);
			multiplyA(v,q);
			alpha=dotProduct(v,q);
			T d_new=alpha-mu*mu*d;
			p=-p*d*mu/d_new;
			scale(c,-mu);
			AXPY(c,T(1),q);
			AXPY(x,p,c);
			AXPY(v,-alpha,q);
			AXPY(v,-beta,q_old);
			d=d_new;
		}
	}

	template<typename T,typename Func1,typename Func2,typename Func3,typename Func4,typename Func5, typename Func6, typename Func7,typename Func8, typename Func9, typename Vector>
	void lanczosPCG(Func1 multiplyA, Func2 multiplyPrecond, Func3 dotProduct, Func4 set, Func5 setScaled, Func6 scaleAndAdd, Func7 addScaled, Func8 addScaled2, Func9 residual, Vector& x, const Vector& b,const int max_it,const T tol=1e-3){
		//multiplyA(a,b): a <- A * b
		//multiplyPrecond(a,b): a <- M * b	
		//set(a,b): a <- b
		//setScaled(a,b,s): a <- s * b
		//scaleAndAdd(a,s,b): a <- s * a + b
		//addScaled(a,b,s): a <- a + s * b
		//addScaled2(a,b1,s1,b2,s2); a <- a + s1 * b1 + s2 * b2
		//residual(r,x,b): r <- |b - A * x|
		
		setScaled(x,x,T(0));
		T r;residual(r,x,b);
		std::cout << "First PCG residual = " << r << std::endl;
		if(r<tol){
			std::cout <<"Converged without PCG iteration"<< std::endl;
			return;
		}
		Vector w_bar(x.size()),w_tilde(x.size()),q_bar(x.size()),q_bar_old(x.size()),p_bar(x.size());
		setScaled(q_bar_old,q_bar_old,T(0));
		set(w_tilde,b);
		multiplyPrecond(w_bar,w_tilde);
		T mag=sqrt(T(dotProduct(w_bar,w_tilde)));
		setScaled(q_bar,w_bar,T(1)/mag);//q_bar_0
		multiplyA(w_tilde,q_bar);
		T alpha=T(dotProduct(q_bar,w_tilde));
		T beta=T(0);
		T d=alpha;
		T t=mag;
		set(p_bar,q_bar);
		setScaled(x,p_bar,t/d);

		for(int it=1;it<max_it;it++){
			residual(r,x,b);
			if(r<tol){
				std::cout << "PCG residual = " << r << std::endl;
				std::cout << "PCG converged in " << it << " iterations." << std::endl;
				return;
			}
			
			multiplyPrecond(w_bar,w_tilde);
			T beta_new_squared = T(dotProduct(w_bar,w_tilde)) - alpha*alpha - beta*beta;
			if(beta_new_squared<=T(0)){
				std::cout << "PCG residual = " << r << std::endl;
				std::cout << "PCG converged to beta < 0 in " << it << " iterations." << std::endl;
				return;
			}
			
			addScaled2(w_bar,q_bar,-alpha,q_bar_old,-beta);
			set(q_bar_old,q_bar);
			beta=sqrt(beta_new_squared);
			setScaled(q_bar,w_bar,T(1)/beta);//q_bar_it
			multiplyA(w_tilde,q_bar);
			alpha=T(dotProduct(q_bar,w_tilde));
			T mu=beta/d;
			t=-mu*t;
			d=alpha-d*mu*mu;
			scaleAndAdd(p_bar,-mu,q_bar);
			addScaled(x,p_bar,t/d);
			
		}
		std::cout << "PCG used max iterations = " << max_it << std::endl;
	}

	struct LanczosPCG{
		TV w_bar;
		TV w_tilde;
		TV q_bar;
		TV q_bar_old;
		TV p_bar;

		void resize(const size_t N){
			w_bar.resize(N);
			w_tilde.resize(N);
			q_bar.resize(N);
			q_bar_old.resize(N);
			p_bar.resize(N);
		}

		template<typename Func1,typename Func2,typename Func3,typename Func4,typename Func5, typename Func6, typename Func7,typename Func8, typename Func9>
		void pcgSolve(Func1 multiplyA, Func2 multiplyPrecond, Func3 dotProduct, Func4 set, Func5 setScaled, Func6 scaleAndAdd, Func7 addScaled, Func8 addScaled2, Func9 residual, TV& x, const TV& b,const int max_it,const T tol=1e-3){
			//multiplyA(a,b): a <- A * b
			//multiplyPrecond(a,b): a <- M * b	
			//set(a,b): a <- b
			//setScaled(a,b,s): a <- s * b
			//scaleAndAdd(a,s,b): a <- s * a + b
			//addScaled(a,b,s): a <- a + s * b
			//addScaled2(a,b1,s1,b2,s2); a <- a + s1 * b1 + s2 * b2
			//residual(r,x,b): r <- |b - A * x|
			
			setScaled(x,x,T(0));
			T r;residual(r,x,b);
			std::cout << "First PCG residual = " << r << std::endl;
			if(r<tol){
				exitCG(setScaled,std::string("Converged without PCG iteration, exiting CG"));
				return;
			}
	
			setScaled(q_bar_old,q_bar_old,T(0));
			set(w_tilde,b);
			multiplyPrecond(w_bar,w_tilde);
			T mag=sqrt(T(dotProduct(w_bar,w_tilde)));
			setScaled(q_bar,w_bar,T(1)/mag);//q_bar_0
			multiplyA(w_tilde,q_bar);
			T alpha=T(dotProduct(q_bar,w_tilde));
			T beta=T(0);
			T d=alpha;
			T t=mag;
			set(p_bar,q_bar);
			setScaled(x,p_bar,t/d);

			for(int it=1;it<max_it;it++){
				residual(r,x,b);
				if(r<tol){
					std::cout << "PCG residual = " << r << std::endl;
					std::cout << "PCG converged in " << it << " iterations." << std::endl;
					exitCG(setScaled,std::string("Exiting CG"));
					return;
				}
				
				multiplyPrecond(w_bar,w_tilde);
				T beta_new_squared = T(dotProduct(w_bar,w_tilde)) - alpha*alpha - beta*beta;
				if(beta_new_squared<=T(0)){
					std::cout << "PCG residual = " << r << std::endl;
					std::cout << "PCG converged to beta < 0 in " << it << " iterations." << std::endl;
					exitCG(setScaled,std::string("Exiting CG"));
					return;
				}
				
				addScaled2(w_bar,q_bar,-alpha,q_bar_old,-beta);
				set(q_bar_old,q_bar);
				beta=sqrt(beta_new_squared);
				setScaled(q_bar,w_bar,T(1)/beta);//q_bar_it
				multiplyA(w_tilde,q_bar);
				alpha=T(dotProduct(q_bar,w_tilde));
				T mu=beta/d;
				t=-mu*t;
				d=alpha-d*mu*mu;
				scaleAndAdd(p_bar,-mu,q_bar);
				addScaled(x,p_bar,t/d);
				
			}
			std::cout << "PCG used max iterations = " << max_it << std::endl;
			std::cout << "PCG residual = " << r << std::endl;
			exitCG(setScaled,std::string("Exiting CG"));
		}

		template<typename Func>
		void exitCG(Func setScaled, const std::string output_message){
			std::cout << output_message << std::endl;
			setScaled(w_bar,w_bar,T(0));
			setScaled(w_tilde,w_tilde,T(0));
			setScaled(q_bar,q_bar,T(0));
			setScaled(q_bar_old,q_bar_old,T(0));
			setScaled(p_bar,p_bar,T(0));
		}
	};
}