#include "definitions.h"
#include <Eigen/Dense>

#ifndef MATHTOOLS_INCLUDED
#define MATHTOOLS_INCLUDED

namespace TGSL{

	inline void polarDecomposition(const T A00, const T A10, const T A01, const T A11, T& c, T& s,T& s00, T& s10, T& s11){
		T v0=A00+A11,v1=A01-A10;
		T denom=sqrt(v0*v0+v1*v1);
		if(denom==T(0)){
			c=T(1);s=T(0);
			s00=A00;s10=A10;s11=A11;
		}
		else{
			c=v0/denom;s=-v1/denom;
			s10=-A00*s + A10*c;
			s00=A00*c + A10*s;
			s11=-A01*s + A11*c;
		}
	}

	inline void symSchur2D(const T Aqq, const T App, const T Apq, T& c, T& s){
		if(Apq != 0){
			T tau = T(0.5)*(Aqq - App)/Apq;
			T t;
			if(tau >= 0)
				t = T(1)/(tau + std::sqrt(1 + tau*tau));
			else
				t = T(-1)/(-tau + std::sqrt(1 + tau*tau));
			c = T(1)/std::sqrt(1 + t*t);
			s = t*c;
		}
		else{
			c = 1;
			s = 0;
		}
	}

	inline void SVD(const T A00, const T A10, const T A01, const T A11, T& cu, T& su,T& cv, T& sv, T& s0, T& s1){
		T c,s,s00,s10,s11;
		polarDecomposition(A00,A10,A01,A11,c,s,s00,s10,s11);
		symSchur2D(s11,s00,s10,cv,sv);
		sv=-sv;
		s0=cv*(cv*s00+sv*s10)+sv*(cv*s10+sv*s11);s1=-sv*(-sv*s00+cv*s10)+cv*(-sv*s10+cv*s11);
		cu=c*cv-s*sv;su=s*cv+c*sv;
	}

	inline void Jacobi(const Eigen::Matrix<T, 2, 2>& B, Eigen::Matrix<T, 2, 2>& D, Eigen::Matrix<T, 2, 2>& V){
		T c, s;
		symSchur2D(B(1, 1), B(0, 0), B(0, 1), c, s);
		
		//Build Givens
		V(0, 0) = c;
		V(0, 1) = s;
		V(1, 0) = -s;
		V(1, 1) = c;

		D = V.transpose()*B*V;
	}

	template<typename F>
	T gaussianQuadrature2(F function, T a, T b) {
		std::array<T, 2> weights = { {T(1), T(1)} };
		std::array<T, 2> quad_points = { {T(1)/std::sqrt(T(3)), T(-1)/std::sqrt(T(3))} };

		T result = T(0);
		for (size_t i = 0; i < 2; ++i)
			result += weights[i] * function((b - a) / T(2)*quad_points[i] + (b + a) / 2);
		result *= (b - a) / T(2);

		return result;
	}

	template<typename F>
	T gaussianQuadrature4(F function, T a, T b) {
		std::array<T, 4> weights = { {(T(18) + std::sqrt(T(30))) / T(36), (T(18) + std::sqrt(T(30))) / T(36),
			(T(18) - std::sqrt(T(30))) / T(36), (T(18) - std::sqrt(T(30))) / T(36)} };
		std::array<T, 4> quad_points = { {std::sqrt(T(3) / T(7) - T(2) / T(7)*std::sqrt(T(6) / T(5))), T(-1)*std::sqrt(T(3) / T(7) - T(2) / T(7)*std::sqrt(T(6) / T(5))),
			std::sqrt(T(3) / T(7) + T(2) / T(7)*std::sqrt(T(6) / T(5))), T(-1)*std::sqrt(T(3) / T(7) + T(2) / T(7)*std::sqrt(T(6) / T(5)))} };

		T result = T(0);
		for (size_t i = 0; i < 4; ++i)
			result += weights[i] * function((b - a) / T(2)*quad_points[i] + (b + a) / 2);
		result *= (b - a) / T(2);

		return result;
	}

	template<typename F>
	T gaussianQuadrature2_2DSquare(F function, T a = T(0), T b = T(1)) {
		std::array<T, 2> weights = { {T(1), T(1)} };
		std::array<T, 2> quad_points = { {T(1) / std::sqrt(T(3)), T(-1) / std::sqrt(T(3))} };


		T result = T(0);
		for (size_t k = 0; k < 2; ++k)
			for (size_t l = 0; l < 2; ++l) {
				result += weights[k] * weights[l] * function(quad_points[k] * (b - a) / T(2) + (b + a) / T(2), quad_points[l] * (b - a) / T(2) + (b + a) / T(2));
			}
		result *= T(0.25);
		result *= (b - a)*(b - a);

		return result;
	}

	template<typename F>
	T gaussianQuadrature4_2DSquare(F function, T a = T(0), T b = T(1)){
		// for integrating over the square [a, b] x [a, b]
		std::array<T, 4> weights = { {(T(18) + std::sqrt(T(30))) / T(36), (T(18) + std::sqrt(T(30))) / T(36),
				(T(18) - std::sqrt(T(30))) / T(36), (T(18) - std::sqrt(T(30))) / T(36)} };
		std::array<T, 4> quad_points = { {std::sqrt(T(3) / T(7) - T(2) / T(7)*std::sqrt(T(6) / T(5))), T(-1)*std::sqrt(T(3) / T(7) - T(2) / T(7)*std::sqrt(T(6) / T(5))),
			std::sqrt(T(3) / T(7) + T(2) / T(7)*std::sqrt(T(6) / T(5))), T(-1)*std::sqrt(T(3) / T(7) + T(2) / T(7)*std::sqrt(T(6) / T(5)))} };

		T result = T(0);
		for(size_t k = 0; k < 4; ++k)
			for(size_t l = 0; l < 4; ++l)
				result += weights[k]*weights[l]*function(quad_points[k]*(b - a)/T(2) + (b + a)/T(2), quad_points[l]*(b - a)/T(2) + (b + a)/T(2));
		result *= T(0.25);
		result *= (b - a)*(b - a);

		return result;
	}

	template<typename F>
	T gaussianQuadrature5_2DSquare(F function, T a = T(0), T b = T(1)){
		// for integrating over the square [a, b] x [a, b]
		std::array<T, 5> weights = { {T(128)/T(225), (T(322) + T(13)*std::sqrt(T(70)))/T(900), (T(322) + T(13)*std::sqrt(T(70)))/T(900),
			(T(322) - T(13)*std::sqrt(T(70)))/T(900), (T(322) - T(13)*std::sqrt(T(70)))/T(900)} };
		std::array<T, 5> quad_points = { {T(0), T(1)/T(3)*std::sqrt(T(5) - T(2)*std::sqrt(T(10)/T(7))), T(-1)/T(3)*std::sqrt(T(5) - T(2)*std::sqrt(T(10)/T(7))),
			T(1)/T(3)*std::sqrt(T(5) + T(2)*std::sqrt(T(10)/T(7))), T(-1)/T(3)*std::sqrt(T(5) + T(2)*std::sqrt(T(10)/T(7)))} };

		T result = T(0);
		for(size_t k = 0; k < 5; ++k)
			for(size_t l = 0; l < 5; ++l)
				result += weights[k]*weights[l]*function(quad_points[k]*(b - a)/T(2) + (b + a)/T(2), quad_points[l]*(b - a)/T(2) + (b + a)/T(2));
		result *= T(0.25);
		result *= (b - a)*(b - a);

		return result;
	}

	inline int gcd(int a, int b){
		while(a != b){
			if(a > b)
				a = a - b;
			else
				b = b - a;
		}

		return a;
	}
}

#endif