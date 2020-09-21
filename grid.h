#pragma once

#include "definitions.h"

namespace TGSL{

struct Grid{
	Index gridN;
	T dx;
	Particle origin;
	enum interpType{
		linear=1,quadratic=2
	};
	interpType interp;
	size_t n_per_dir;	// number of interpolating points per side of triangle?

	Grid(){
		#ifdef LINEAR
			interp=linear;
		#elif defined QUADRATIC
			interp=quadratic;
		#endif
		if(interp==linear)
			n_per_dir = 2;
		else
			n_per_dir = 3;
	}

	Grid(const Index& in_gridN, const T in_dx, const Particle& in_origin) : gridN(in_gridN), dx(in_dx), origin(in_origin) {
		#ifdef LINEAR
			interp=linear;
		#elif defined QUADRATIC
			interp=quadratic;
		#endif
		if(interp==linear)
			n_per_dir = 2;
		else
			n_per_dir = 3;
	}

	void setInterp(interpType interp_in){
		interp=interp_in;
		if(interp==linear)
			n_per_dir = 2;
		else
			n_per_dir = 3;
	}

	T Nijk(T w, nm ii)const {
		#ifdef LINEAR
			return T(1) - w + T(ii)*(T(2)*w - T(1));
		#elif defined QUADRATIC
			return ((T(1.5)*(w*w - w) - T(0.25))*(T(ii) - 1) + (T(0.5)*w - T(0.25)))*(T(ii) - 1) - w*w + w + T(0.5);
		#endif
		TGSLAssert(false,"interpolation case not defined");
		return T(0);
	}

	T dNijk(T w, nm ii, T dx)const {
		#ifdef LINEAR
			return (T(2)*T(ii) - T(1))/dx;
		#elif defined QUADRATIC
			return ((w - 1)*(T(ii) - 1)*(T(ii) - 2)*T(0.5) + (2*w - 1)*T(ii)*(T(ii) - 2) + w*T(ii)*(T(ii) - 1)*T(0.5))/dx;
		#endif
		TGSLAssert(false,"interpolation case not defined");
		return T(0);
	}

	T NiProd(const Vector& Ni)const {
		return std::accumulate(Ni.begin(), Ni.end(), T(1), std::multiplies<T>());
	}

	void gradNi(const Vector& Ni, const Vector& dNi, Vector& result)const {
		for(size_t i = 0; i < d; ++i)
			result[i] = T(1);
		
		for(size_t i = 0; i < d; ++i)
			for(size_t j = 0; j < d; ++j)
				result[i] *= (i == j) ? dNi[j] : Ni[j];
	}

	Vector lengths()const{
		#ifdef ONE_D
		    return {dx*(gridN[0]-1)};
		#elif defined TWO_D
			return {dx*(gridN[0]-1),dx*(gridN[1]-1)};
		#else
			return {dx*(gridN[0]-1),dx*(gridN[1]-1),dx*(gridN[2]-1)};
		#endif

	}

	void baseNodeIndex(const Particle& x, Index& index, Vector& weights) const{
		for(size_t i = 0; i < d; ++i)
			index[i] = int(std::floor((x[i] - origin[i])/dx));
		for(size_t i = 0; i < d; ++i)
			weights[i] = (x[i] - origin[i])/dx - index[i];
	}

	bool outsideGridInterp(const Index& base) const{
		/*
		Use this for interpolating data defined over nodes of grid. Interpolation includes the case
		where only some of the grid nodes would be used (the remaining values would require out of bounds 
		indices).

		(9/3/19), added convnetion for quadratic.
		*/
		if(interp==linear){
			bool outside = false;
			for(size_t i = 0; i < d; ++i)
				if(base[i] < -1 || base[i] > gridN[i] - 1)
					outside = true;
			return outside;
		}
		else{
			bool outside = false;
			for(size_t i = 0; i < d; ++i)
				if(base[i] < 0 || base[i] > gridN[i] - 2)
					outside = true;
			return outside;
		}
	}

	bool validIndex(const Index& index) const{
		/*
		This assumes data will be stored at each index (i.e. numbered by flatIndex(index)). For linear interpolation over
		the grid we check that the index[i] \in [0,gridN[i]-1]. For quadratic interpolation, we check that the 
		index[i] \in [-1,gridN[i]-1].

		For linear interpolation, data is stored at the nodes of the grid, e.g. xi = i dx (i=0,1,...,gridN-1) and
		for quadratic interpolation, data is stored at the cell centers of the grid, e.g. xi = i dx + dx/2 (i=-1,0,1,...,gridN-1)
		*/

		if(interp==linear){
			bool result = true;
			for(size_t i = 0; i < d; ++i)
				if(index[i] < 0 || index[i] > (gridN[i]-1))
					result = false;
			return result;
		}
		else{
			bool result = true;
			for(size_t i = 0; i < d; ++i)
				if(index[i] < -1 || index[i] > (gridN[i]-1))
					result = false;
			return result;
		}
	}

	nm flatIndex(const Index& index) const{
		/*
		This assumes data will be stored at each index (i.e. numbered by flatIndex(index)). For linear interpolation over
		the grid we assume that the index[i] \in [0,gridN[i]-1]. For quadratic interpolation, we assume that the 
		index[i] \in [-1,gridN[i]-1].

		For linear interpolation, data is stored at the nodes of the grid, e.g. xi = i dx (i=0,1,...,gridN-1) and
		for quadratic interpolation, data is stored at the cell centers of the grid, e.g. xi = i dx + dx/2 (i=-1,0,1,...,gridN-1)
		*/

		#ifdef ONE_D
		    if(interp==linear){
				return index[0];
			}
			else{
				return index[0] + 1;
			}
		#elif defined TWO_D
			if(interp==linear){
				return index[0]*gridN[1] + index[1];
			}
			else{
				return (gridN[1] + 1)*(index[0] + 1) + index[1] + 1;
			}
		#else
			if(interp==linear){
				return index[0]*gridN[1]*gridN[2] + index[1]*gridN[2] + index[2];
			}
			else{
				return (gridN[1] + 1)*(gridN[2] + 1)*(index[0] + 1) + (gridN[2] + 1)*(index[1] + 1) + index[2] + 1;
			}
		#endif
		return -1;
	}

	Index loc2globIndex(const Index& index, const Index& locIndex) const{
		Index result;
		if(interp == linear){
			for(size_t i = 0; i < d; ++i)
				result[i] = index[i] + locIndex[i];
		}
		else{
			for(size_t i = 0; i < d; ++i)
				result[i] = index[i] + locIndex[i] - 1;
		}
		return result;
	}

	nm loc2globIndex(const nm index, const Index& locIndex) const{
		/*
		Converts the local multi-index to global flat index.
		index = flat index of base node (for linear this is lower left node, for quadratic this is the cell center)
		locIndex[i] = 0,1 for linear and 0,1,2 for quadratic.
		*/

		#ifdef ONE_D
		    if(interp==linear){
		    	return index + locIndex[0];
		    }
		    else{
		    	return index + locIndex[0] - 1;
		    }
		#elif defined TWO_D
			if(interp==linear){
				return index + gridN[1]*locIndex[0] + locIndex[1];
			}
			else{
				return index + (gridN[1] + 1)*(locIndex[0] - 1) + locIndex[1] - 1;
			}
		#else
			if(interp==linear){
				return index + gridN[1]*gridN[2]*locIndex[0] + gridN[2]*locIndex[1] + locIndex[2];
			}
			else{
				return index + (gridN[1] + 1)*(gridN[2] + 1)*(locIndex[0] - 1) + (gridN[2] + 1)*(locIndex[1] - 1) + locIndex[2] - 1;
			}
		#endif
		return -1;
	}

	void lin2multiIndex(const nm linIndex, Index& multiIndex) const{
		/* This inverts the conventions in flatIndex and gives the multiIndex from the flat index.
 		3*3 grid:
 		
		2 5 8
		1 4 7
		0 3 6
		*/

		#ifdef ONE_D
			if(interp==linear){
		    	multiIndex[0] = linIndex;
		    }
		    else{
		    	multiIndex[0] = linIndex - 1;
		    }
		#elif defined TWO_D
			if(interp==linear){
				multiIndex[0] = linIndex/gridN[1];
				multiIndex[1] = linIndex%gridN[1];
			}
			else{
				multiIndex[0] = linIndex/(gridN[1] + 1) - 1;
				multiIndex[1] = linIndex%(gridN[1] + 1) - 1;
			}
		#else
			if(interp==linear){
				multiIndex[0] = linIndex/(gridN[1]*gridN[2]);
				multiIndex[1] = (linIndex/gridN[2])%gridN[1];
				multiIndex[2] = linIndex%gridN[2];
			}
			else{
				multiIndex[0] = linIndex/((gridN[1] + 1)*(gridN[2] + 1)) - 1;
				multiIndex[1] = (linIndex/(gridN[2] + 1))%(gridN[1] + 1) - 1;
				multiIndex[2] = linIndex%(gridN[2] + 1) - 1;
			}
		#endif
	}

	Particle node(const Index& index) const{
		/* This gives the position of the grid node associated with multi index (index). These are nodes for linear 
		interpolation and cell centers for quadratic interpolation.
		*/

		Particle result;

		if(interp==linear){
			for(size_t i = 0; i < d; ++i){
				result[i] = T(index[i])*dx + origin[i];
			}
		}
		else{
			for(size_t i = 0; i < d; ++i){
				result[i] = (T(index[i]) + T(0.5))*dx + origin[i];
			}
		}

		return result;
	}

	Particle node(nm flatIndex) const{
		Particle result;
		Index index = Index();
		lin2multiIndex(flatIndex, index);

		if(interp==linear){
			for(size_t i = 0; i < d; ++i){
				result[i] = T(index[i])*dx + origin[i];
			}
		}
		else{
			for(size_t i = 0; i < d; ++i){
				result[i] = (T(index[i]) + T(0.5))*dx + origin[i];
			}
		}

		return result;
	}

	size_t size() const{
		size_t result = 1;

		if(interp == linear)
			for(size_t i = 0; i < d; ++i)
				result *= gridN[i];
		else
			for(size_t i = 0; i < d; ++i)
				result *= (gridN[i] + 1);

		return result;
	}

	bool interiorNode(const Index& index)const{
		bool interior = true;
		for(size_t c=0;c<d;c++){
			if(index[c]==-1 || index[c]==gridN[c]-1){
				interior=false;
			}
		}
		return interior;
	}

	template<int D, typename F1>
	std::array<T,D> interpolate(const Index& index, F1 grid_data) const{
		TGSLAssert(interp==quadratic,"Interpolation to grid nodes only makes sense for quadratic (and higher).");
		Vector3T weights={T(1)/T(8),T(3)/T(4),T(1)/T(8)};
		std::array<T,D> result;result.fill(T(0));
		TGSLAssert(interiorNode(index),"Interpolation to grid node only defined for interior nodes.");
		#ifdef ONE_D
			for(int ii=-1;ii<2;ii++){
				Index n={index[0]+ii};
				for(size_t c=0;c<D;c++)	
					result[c]+=weights[size_t(ii+1)]*grid_data(n)[c];
			}
		#elif TWO_D
			for(int ii=-1;ii<2;ii++){
				for(int jj=-1;jj<2;jj++){
					Index n={index[0]+ii,index[1]+jj};
					T value=weights[size_t(ii+1)]*weights[size_t(jj+1)];
					for(size_t c=0;c<D;c++)	
						result[c]+=value*grid_data(n)[c];
				}
			}
		#else
			for(int ii=-1;ii<2;ii++){
				for(int jj=-1;jj<2;jj++){
					for(int kk=-1;kk<2;kk++){
						Index n={index[0]+ii,index[1]+jj,index[2]+kk};
						T value=weights[size_t(ii+1)]*weights[size_t(jj+1)]*weights[size_t(kk+1)];
						for(size_t c=0;c<D;c++)	
							result[c]+=value*grid_data(n)[c];
					}
				}
			}
		#endif	
		return result;			
	}

	template<int D, typename F1, typename F2>
	std::array<T,D> interpolate(const Particle& xp, F2 grid_data, F1 bc, const T time) const{
		Vector local_weights;
		Index base_index;
		baseNodeIndex(xp, base_index, local_weights);

		if(outsideGridInterp(base_index))
			return bc(xp, time);
		else{
			using VectorDT = std::array<T, D>;
			VectorDT result;result.fill(T(0));

			nm flat=flatIndex(base_index);
			for(size_t ii = 0; ii < n_per_dir; ++ii){
				T Nii = Nijk(local_weights[0], ii);
				#ifdef ONE_D
					Vector Ni = {Nii};
					Index loc_index = {int(ii)};
					nm index = loc2globIndex(flat, loc_index);
					for(int c = 0; c < D; c++)
						result[c]+=NiProd(Ni)*grid_data(index)[c];
				#endif
				#ifndef ONE_D
					for(size_t jj = 0; jj < n_per_dir; ++jj){
						T Njj = Nijk(local_weights[1], jj);
						#ifdef TWO_D
							Vector Ni = {Nii, Njj};
							Index loc_index = {int(ii), int(jj)};
							nm index = loc2globIndex(flat, loc_index);
							for(int c=0;c<D;c++)
								result[c]+=NiProd(Ni)*grid_data(index)[c];
						#endif
						#ifndef TWO_D
							for(size_t kk = 0; kk < n_per_dir; ++kk){
								T Nkk = Nijk(local_weights[2], kk);
								Vector Ni = {Nii, Njj, Nkk};
								Index loc_index = {int(ii), int(jj), int(kk)};
								nm index = loc2globIndex(flat, loc_index);
								for(int c=0;c<D;c++)
									result[c]+=NiProd(Ni)*grid_data(index)[c];
							}
						#endif
					}
				#endif
			}
			return result;
		}
	}

	template<int D, typename F1, typename F2>
	std::array<T,D*d> diffInterpolant(const Particle& xp, F2 grid_data, F1 bc, const T time) const{
		Vector local_weights;
		Index base_index;
		baseNodeIndex(xp, base_index, local_weights);

		if(outsideGridInterp(base_index))
			return bc(xp, time);
		else{
			using VectorDdT = std::array<T, D*d>;
			VectorDdT result;result.fill(T(0));

			nm flat=flatIndex(base_index);
			for(size_t ii = 0; ii < n_per_dir; ++ii){
				T Nii = Nijk(local_weights[0],ii);
				#ifdef ONE_D
					T dNii=dNijk(local_weights[0],ii,dx);
					Index loc_index = {int(ii)};
					nm index = loc2globIndex(flat, loc_index);
					for(int c = 0; c < D; c++)
						result[c]+=dNii*grid_data(index)[c];
				#endif
				#ifndef ONE_D
					for(size_t jj = 0; jj < n_per_dir; ++jj){
						T Njj = Nijk(local_weights[1],jj);
						#ifdef TWO_D
							T dNii=dNijk(local_weights[0],ii,dx);
							T dNjj=dNijk(local_weights[1],jj,dx);
							Vector temp_gradNi;
							gradNi({Nii,Njj}, {dNii,dNjj},temp_gradNi);
							Index loc_index = {int(ii), int(jj)};
							nm index = loc2globIndex(flat, loc_index);
							for(size_t r=0;r<D;r++){
								for(size_t c=0;c<d;c++){
									result[d*r+c]+=grid_data(index)[r]*temp_gradNi[c];
								}
							}
						#endif
						#ifndef TWO_D
							for(size_t kk = 0; kk < n_per_dir; ++kk){
								T Nkk = Nijk(local_weights[2], kk);
								T dNii=dNijk(local_weights[0],ii,dx);
								T dNjj=dNijk(local_weights[1],jj,dx);
								T dNkk=dNijk(local_weights[2],kk,dx);
								Vector temp_gradNi;
								gradNi({Nii,Njj,Nkk}, {dNii,dNjj,dNkk},temp_gradNi);
								Index loc_index = {int(ii), int(jj), int(kk)};
								nm index = loc2globIndex(flat, loc_index);
								for(size_t r=0;r<D;r++){
									for(size_t c=0;c<d;c++){
										result[d*r+c]+=grid_data(index)[r]*temp_gradNi[c];
									}
								}
							}
						#endif
					}
				#endif
			}
			return result;
		}
	}
};
}
