#include <vector>
#include <array>
#include <cmath>
#include <cstdint>
#include <string>
#include <iostream>
#include <random>

// Headers for error handling & debugging
#include <stdexcept>
#include <cassert>

#ifndef DEFINITIONS_INCLUDED
#define DEFINITIONS_INCLUDED

namespace TGSL{
	#ifdef BIGNUM
	typedef uint_fast64_t sz;
	typedef int_fast64_t nm;
	#else
	typedef size_t sz;
	typedef int nm;
	#endif

	#ifdef ONE_D
		constexpr size_t d = 1;
	#elif defined TWO_D
		constexpr size_t d = 2;
	#else
		constexpr size_t d = 3;
	#endif

	using T = double;

	using TV = std::vector<T>;
	using IV = std::vector<nm>;
	using IVV = std::vector<IV>;

	using Vector2T = std::array<T, 2>;
	using Vector3T = std::array<T, 3>;
	using Vector2I = std::array<nm, 2>;
	using Vector3I = std::array<nm, 3>;

	#ifdef ONE_D
		using Particle = std::array<T, 1>;
		using Vector = std::array<T, 1>;
		using Index = std::array<nm, 1>; 
	#elif defined TWO_D
		using Particle = std::array<T, 2>;
		using Vector = std::array<T, 2>;
		using Index = std::array<nm, 2>;
	#else
		using Particle = std::array<T, 3>;
		using Vector = std::array<T, 3>;
		using Index = std::array<nm, 3>;
	#endif

	using TVP3 = std::vector<Vector3T>;
	using TVP2 = std::vector<Vector2T>;
	using TVP = std::vector<Particle>;
	using TVV = std::vector<Vector>;
	using TVI = std::vector<Index>;
	using TVI3 = std::vector<Vector3I>;
	using TVI2 = std::vector<Vector2I>;

	constexpr double pi = 3.14159265358979323846;

	inline void TGSLAssert(const bool success, std::string flag){
		if(!success){
			std::cout << flag << std::endl;
			exit(1);
		}
	}
}
#endif
