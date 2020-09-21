#pragma once
#include <algorithm>
#include <iostream>

#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/task_scheduler_init.h>
#include <pstl/execution>
#include <pstl/algorithm>
#include <pstl/numeric>
#include <pstl/memory>

#include "definitions.h"

namespace TGSL{

/*

*/
void computeIncidentElements(const size_t element_size,const IV& mesh, std::vector<std::vector<int>>& incident_elements)
{
	// e.g. element_size = 3, mesh: 3*N, incident_element[i]: vector of indices of points that is adjacent to point i
	IV ordering(mesh.size());
	IV ranges(mesh.size()+1);
	sz num_nodes=sz(0);

	tbb::parallel_for(tbb::blocked_range<sz>(0, mesh.size()), [&](const tbb::blocked_range<sz>& range){
		for(sz p = range.begin(); p < range.end(); ++p){
			ordering[p] = p;
			ranges[p]=p;
		}
	});

	/*for(size_t n=0;n<mesh.size();n++)
		ordering[n]=n;*/

	std::sort(pstl::execution::par, std::begin(ordering), std::end(ordering), [&mesh](int a, int b) {
		return mesh[a]<mesh[b];});

	/*std::sort(std::begin(ordering), std::end(ordering), [&mesh](int& a, int& b) {
		return mesh[a]<mesh[b];
	});*/

	auto last = std::unique(pstl::execution::par, std::begin(ranges), std::begin(ranges) + mesh.size(), 
		[&](nm i, nm j){return mesh[ordering[i]] == mesh[ordering[j]];});

	num_nodes = last - ranges.begin();
	ranges[num_nodes] = nm(mesh.size());
	ranges.resize(num_nodes+1);

	/*IV ranges;
	ranges.emplace_back(0);
	for(size_t n=0;n<mesh.size()-1;n++){
		if(mesh[ordering[n]]!=mesh[ordering[n+1]]){
			ranges.emplace_back(n+1);
		}
	}
	ranges.emplace_back(mesh.size());*/

	incident_elements.resize(ranges.size()-1);
	tbb::parallel_for(tbb::blocked_range<sz>(0, sz(incident_elements.size())), [&](const tbb::blocked_range<sz>& range){
		for(sz p = range.begin(); p < range.end(); ++p){
			incident_elements[p].resize(ranges[p+1]-ranges[p]);
			for(int e=ranges[p];e<ranges[p+1];e++){
				incident_elements[p][e-ranges[p]]=ordering[e];
			}	
		}
	});

	/*incident_elements.resize(ranges.size()-1);
	for(size_t r=0;r<ranges.size()-1;r++){
		incident_elements[r].resize(ranges[r+1]-ranges[r]);
		for(int e=ranges[r];e<ranges[r+1];e++){
			incident_elements[r][e-ranges[r]]=ordering[e];
		}
	}*/
}

/*
Computes an array of boundary faces for a mesh of elements that can return element-wise faces and define ordering for faces.
*/
template<typename Func1, typename Func2, typename Func3, typename Face>
void computeBoundaryMesh(const size_t number_elements, Func1 addElementFaces, Func2 greaterThan, Func3 equal, std::vector<Face>& boundary_mesh)
{
	std::vector<Face> mesh_faces;
	for(size_t e=0;e<number_elements;e++)
		addElementFaces(e,mesh_faces);

	std::sort(std::begin(mesh_faces), std::end(mesh_faces), [&greaterThan](Face& a, Face& b) {
		return greaterThan(a,b); 
	});

	boundary_mesh.resize(0);
	int count=1;

	for(size_t f=0;f<mesh_faces.size();f++){
		int next_f=f+1;
		if(next_f == mesh_faces.size() && count == 1)
			boundary_mesh.emplace_back(mesh_faces[f]);
		else if(!equal(mesh_faces[next_f],mesh_faces[f])){
			if(count == 1)
				boundary_mesh.emplace_back(mesh_faces[f]);
			count = 1;
		}
		else
			count++;
	}	
}

void computeBoundaryMesh(const TVI3& tri_mesh,TVI2& boundary_mesh){
	auto addElementFaces = [&tri_mesh](const size_t e, TVI2& faces){
	Vector2I f0={tri_mesh[e][0],tri_mesh[e][1]};
	Vector2I f1={tri_mesh[e][1],tri_mesh[e][2]};
	Vector2I f2={tri_mesh[e][2],tri_mesh[e][0]};
	faces.emplace_back(f0);
	faces.emplace_back(f1);
	faces.emplace_back(f2);
};

auto greaterThan = [](Vector2I& i1, Vector2I& i2){
	bool greater_than = false;
	Vector2I is1,is2;
	if(i1[0]<i1[1])
		is1={i1[0],i1[1]};
	else
		is1={i1[1],i1[0]};
	if(i2[0]<i2[1])
		is2={i2[0],i2[1]};
	else
		is2={i2[1],i2[0]};

	if(is1[0]>is2[0])
		greater_than = true;
	else if(is1[0]==is2[0] && is1[1]>is2[1])
		greater_than = true;

	return greater_than;
};

auto equal = [](Vector2I& i1, Vector2I& i2){
	Vector2I is1,is2;
	if(i1[0]<i1[1])
		is1={i1[0],i1[1]};
	else
		is1={i1[1],i1[0]};
	if(i2[0]<i2[1])
		is2={i2[0],i2[1]};
	else
		is2={i2[1],i2[0]};

	return (is1[0]==is2[0] && is1[1]==is2[1]);
};

computeBoundaryMesh(tri_mesh.size(),addElementFaces,greaterThan, equal,boundary_mesh);
}


//Used in main.cpp
void computeBoundaryMesh(const IV& tri_mesh,TVI2& boundary_mesh){
	auto addElementFaces = [&tri_mesh](const size_t e, TVI2& faces){
		// store mesh triangle sides (i, j)
		Vector2I f0={tri_mesh[3*e],tri_mesh[3*e+1]};
		Vector2I f1={tri_mesh[3*e+1],tri_mesh[3*e+2]};
		Vector2I f2={tri_mesh[3*e+2],tri_mesh[3*e]};
		faces.emplace_back(f0);
		faces.emplace_back(f1);
		faces.emplace_back(f2);
	};

	auto greaterThan = [](Vector2I& i1, Vector2I& i2){
		bool greater_than = false;
		Vector2I is1,is2;
		if(i1[0]<i1[1])
			is1={i1[0],i1[1]};
		else
			is1={i1[1],i1[0]};
		if(i2[0]<i2[1])
			is2={i2[0],i2[1]};
		else
			is2={i2[1],i2[0]};

		if(is1[0]>is2[0])
			greater_than = true;
		else if(is1[0]==is2[0] && is1[1]>is2[1])
			greater_than = true;

		return greater_than;
	};

	auto equal = [](Vector2I& i1, Vector2I& i2){
		Vector2I is1,is2;
		if(i1[0]<i1[1])
			is1={i1[0],i1[1]};
		else
			is1={i1[1],i1[0]};
		if(i2[0]<i2[1])
			is2={i2[0],i2[1]};
		else
			is2={i2[1],i2[0]};

		return (is1[0]==is2[0] && is1[1]==is2[1]);
	};

	computeBoundaryMesh(tri_mesh.size()/3,addElementFaces,greaterThan, equal,boundary_mesh);
}

}