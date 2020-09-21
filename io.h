#pragma once

#include <fstream>
#include <string>

#include "definitions.h"
using namespace TGSL;

template<typename VectorType, typename F1>
void writeXYZN(const std::vector<VectorType>& positions, const std::vector<VectorType>& normals, const size_t N, const std::string filename, F1 filter){
	std::ofstream outstream(filename);
	if(outstream.good()){
		for(size_t i = 0; i < N; ++i){
			if (filter(i)) {
				size_t j;
				for (j = 0; j < positions[i].size(); ++j)
					outstream << positions[i][j] << " ";
				for (; j < 3; ++j)
					outstream << "0 ";
				for (j = 0; j < normals[i].size(); ++j)
					outstream << normals[i][j] << " ";
				for (; j < 3; ++j)
					outstream << "0 ";
				outstream << std::endl;
			}
		}

		outstream.close();
	}
}

template<typename VectorType>
void writeXYZN(const std::vector<VectorType>& positions, const size_t N, const std::string filename){
	std::vector<VectorType> normals(0);
	writeXYZN(positions, normals, N, filename, [](const size_t i){return true;});
}

template<typename VectorType>
void writeXYZN(const std::vector<VectorType>& positions, const std::vector<VectorType>& normals, const size_t N, const std::string filename){
	writeXYZN(positions, normals, N, filename, [](const size_t i){return true;});
}

inline void writeCSV(const TV& x, const TV& y, const std::string filename){
	std::ofstream outstream(filename);
	if(outstream.good()){
		size_t N=x.size();
		TGSLAssert(N<=y.size(),"mismatched array sizes for write CSV");
		outstream << "x,y\n";
		for(size_t i = 0; i < N; ++i){
				outstream << x[i] << "," << y[i] << std::endl;
		}
		outstream.close();
	}
}

template<typename F1>
inline void writeOBJ(const TVP& positions, const std::string filename, F1 filter){
	sz N = positions.size();
	std::ofstream outstream(filename);
	if(outstream.good()){
		outstream << "# " << N << " points\n";
		for(sz i = 0; i < N; ++i){
			if (filter(i)) {
				outstream << "v";
				size_t j;
				for(j = 0; j < positions[i].size(); ++j)
					outstream << " " << positions[i][j];
				for(; j < 3; ++j)
					outstream << " 0";
				outstream << std::endl;
			}
		}
		outstream.close();
	}
}

inline void writeOBJ(const TVP& positions, const IV& mesh, const std::string filename){
	sz N = positions.size();
	std::ofstream outstream(filename);
	outstream << "H" << std::endl;
    if (not outstream.is_open())
    {
        std::cout << "Failed to open outputfile.\n";
    }
	if(outstream.good()){
		outstream << "# " << N << " points\n";
		for(sz i = 0; i < N; ++i){
			outstream << "v";
			sz j;
			for(j = 0; j < positions[i].size(); ++j)
				outstream << " " << positions[i][j];
			for(; j < 3; ++j)
				outstream << " 0";
			outstream << std::endl;
		}
		for(sz e = 0; e < mesh.size()/3; ++e){
			outstream << "f " << mesh[3*e]+1 << " " << mesh[3*e+1]+1 << " " << mesh[3*e+2]+1;
			outstream << std::endl;
		}
		outstream.close();
	}
}

template<typename VectorType, typename F1>
void writeOBJ(const std::vector<VectorType>& positions, const std::vector<VectorType>& normals, const size_t N, const std::string filename, F1 filter){
	std::ofstream outstream(filename);
	if(outstream.good()){
		outstream << "# " << N << " points\n" << "g\n";
		for(size_t i = 0; i < N; ++i){
			if (filter(i)) {
				outstream << "v";
				size_t j;
				for(j = 0; j < positions[i].size(); ++j)
					outstream << " " << positions[i][j];
				for(; j < 3; ++j)
					outstream << " 0";
				outstream << std::endl;
			}
		}
		if (normals.size()) {
			for(size_t i = 0; i < N; ++i){
				if (filter(i)) {
					outstream << "vn";
					size_t j;
					for(j = 0; j < normals[i].size(); ++j)
						outstream << " " << normals[i][j];
					for(; j < 3; ++j)
						outstream << " 0";
					outstream << std::endl;
				}
			}
		}
		outstream.close();
	}
}

template<typename VectorType>
void writeOBJ(const std::vector<VectorType>& positions, const size_t N, const std::string filename){
	std::vector<VectorType> normals(0);
	writeOBJ(positions, normals, N, filename, [](const size_t i){return true;});
}

template<typename VectorType>
void writeOBJ(const std::vector<VectorType>& positions, const std::vector<VectorType>& normals, const size_t N, const std::string filename){
	writeOBJ(positions, normals, N, filename, [](const size_t i){return true;});
}

template<typename VectorType>
void writePLY(const std::vector<VectorType>& positions, const TV& radii, const size_t N, const std::string filename)
{
	std::ofstream outstream(filename);
	if(outstream.good()){
		outstream << "ply" << std::endl;
		outstream << "format ascii 1.0" << std::endl;
		outstream << "element vertex " << N << std::endl;
		outstream << "property double x" << std::endl;
		outstream << "property double y" << std::endl;
		outstream << "property double z" << std::endl;
		outstream << "property double radii" << std::endl;
		outstream << "end_header" << std::endl;

		for(size_t i = 0; i < N; ++i){
			size_t j;
			for(j = 0; j < positions[i].size(); ++j)
				outstream << positions[i][j] << " ";
			for(; j < 3; ++j)
				outstream << "0 ";
			outstream << radii[i] << std::endl;
		}
		outstream.close();
	}
}

template<typename VectorType>
void writePLY(const std::vector<VectorType>& positions, const TV& radii, const std::vector<bool>& use_B, const size_t N, const std::string filename)
{
	std::ofstream outstream(filename);
	if(outstream.good()){
		outstream << "ply" << std::endl;
		outstream << "format ascii 1.0" << std::endl;
		outstream << "element vertex " << N << std::endl;
		outstream << "property double x" << std::endl;
		outstream << "property double y" << std::endl;
		outstream << "property double z" << std::endl;
		outstream << "property double radii" << std::endl;
		outstream << "property bool use_B" << std::endl;
		outstream << "end_header" << std::endl;

		for(size_t i = 0; i < N; ++i){
			size_t j;
			for(j = 0; j < positions[i].size(); ++j)
				outstream << positions[i][j] << " ";
			for(; j < 3; ++j)
				outstream << "0 ";
			outstream << radii[i] << std::endl;
			outstream << use_B[i] << std::endl;
		}
		outstream.close();
	}
}

template<typename VectorType>
void writePLY(const std::vector<VectorType>& positions, const TVV& velocity, const TV& curlv, const TV& radii, const std::vector<bool>& use_B, const size_t N, const std::string filename)
{
	std::ofstream outstream(filename);
	if(outstream.good()){
		outstream << "ply" << std::endl;
		outstream << "format ascii 1.0" << std::endl;
		outstream << "element vertex " << N << std::endl;
		outstream << "property double x" << std::endl;
		outstream << "property double y" << std::endl;
		outstream << "property double z" << std::endl;
		outstream << "property double vx" << std::endl;
		outstream << "property double vy" << std::endl;
		outstream << "property double vz" << std::endl;
		outstream << "property double curlv_x" << std::endl;
		outstream << "property double curlv_y" << std::endl;
		outstream << "property double curlv_z" << std::endl;
		outstream << "property double radii" << std::endl;
		outstream << "property bool use_B" << std::endl;
		outstream << "end_header" << std::endl;

		for(size_t i = 0; i < N; ++i){
			size_t j;
			for(j = 0; j < positions[i].size(); ++j)
				outstream << positions[i][j] << " ";
			for(; j < 3; ++j)
				outstream << "0 ";
			for(j = 0; j < positions[i].size(); ++j)
				outstream << velocity[i][j] << " ";
			for(; j < 3; ++j)
				outstream << "0 ";
			#ifdef TWO_D
				outstream << curlv[i] << " ";
			#else
				for(j = 0; j < positions[i].size(); ++j)
					outstream << curlv[i*d+j] << " ";
				for(; j < 3; ++j)
					outstream << "0 ";
			#endif
			outstream << radii[i] << " " << use_B[i] << std::endl;
		}
		outstream.close();
	}
}