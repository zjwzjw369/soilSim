#ifndef SETUP_FUNCTIONS_H
#define SETUP_FUNCTIONS_H

#include "Common.h"
#include "Parameters.h"
#include "Particle.h"
#include <time.h>

inline void createParticleGrid(std::vector<Particle>& particles, solverParams* sp, float3 lower, int3 dims, float radius, float mass, float3 velocity) {
	//srand(int(time(NULL)));
	srand(16);

	for (int x = 0; x < dims.x; x++) {
		for (int y = 0; y < dims.y; y++) {
			for (int z = 0; z < dims.z; z++) {
				float r1 = 0.001f + static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
				float r2 = 0.001f + static_cast <float>(rand()) / static_cast <float> (RAND_MAX);
				float r3 = 0.001f + static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
				float3 jitter = make_float3(r1, r2, r3) * radius;
				float3 pos = lower + make_float3(float(x), float(y), float(z)) * radius + jitter;
				particles.push_back(Particle(pos, velocity, mass));
			}
		}
	}
}

inline void createParticleSphereGrid(std::vector<Particle>& particles, solverParams* sp, float3 lower, int3 dims, float radius, float mass, float3 velocity, float SphereRadius) {
	//srand(int(time(NULL)));
	srand(16);

	for (int x = 0; x < dims.x; x++) {
		for (int y = 0; y < dims.y; y++) {
			for (int z = 0; z < dims.z; z++) {
				float r1 = 0.001f + static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
				float r2 = 0.001f + static_cast <float>(rand()) / static_cast <float> (RAND_MAX);
				float r3 = 0.001f + static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
				float3 jitter = make_float3(r1, r2, r3) * radius;
				float3 pos = lower + make_float3(float(x), float(y), float(z)) * radius + jitter;
				if (length(make_float2(pos.x, pos.z) - make_float2(lower.x + dims.x / 2.0* radius, lower.z + dims.z / 2.0* radius)) < SphereRadius) {
					particles.push_back(Particle(pos, velocity, mass));
				}
			}
		}
	}
}

//Some method for creating a snowball
inline void createSnowball(std::vector<Particle>& particles, float3 center, int3 dims, float radius, float mass, float3 velocity) {
	float sphereRadius = radius * (float)dims.x / 2.0f;
	for (int x = -dims.x / 2; x <= dims.x / 2; x++) {
		for (int y = -dims.y / 2; y <= dims.y / 2; y++) {
			for (int z = -dims.z / 2; z <= dims.z / 2; z++) {
				// generate a jittered point
				float r1 = 0.001f + static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
				float r2 = 0.001f + static_cast <float>(rand()) / static_cast <float> (RAND_MAX);
				float r3 = 0.001f + static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
				float3 jitter = make_float3(r1, r2, r3) * radius / 1.5;

				float3 pos = center + make_float3(float(x), float(y), float(z)) * radius + jitter;
				// see if pos is inside the sphere
				if (length(pos - center) < sphereRadius) {
					particles.push_back(Particle(make_float3(pos.x, pos.y, pos.z), velocity, mass));
				}
			}
		}
	}
}

inline void createTorus(std::vector<Particle>& particles, float3 center, int3 dims, float radius, float mass, float3 velocity) {
	float2 sphereRadius = make_float2(radius * (float)dims.x / 3.0f, radius * (float)dims.x / 6.0f);
	for (int x = -dims.x; x <= dims.x; x++) {
		for (int y = -dims.y; y <= dims.y; y++) {
			for (int z = -dims.z; z <= dims.z; z++) {
				// generate a jittered point
				float r1 = 0.001f + static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
				float r2 = 0.001f + static_cast <float>(rand()) / static_cast <float> (RAND_MAX);
				float r3 = 0.001f + static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
				float3 jitter = make_float3(r1, r2, r3) * radius / 1.5;

				float3 pos = center + make_float3(float(x), float(y), float(z)) * radius + jitter;
				// see if pos is inside the sphere
				float2 q = make_float2(length(make_float2((pos - center).x, (pos - center).z)) - sphereRadius.x, (pos - center).y);
				if (length(q) < sphereRadius.y) {
					particles.push_back(Particle(make_float3(pos.x, pos.y, pos.z), velocity, mass));
				}
			}
		}
	}
}
inline void createCustomModel(std::vector<Particle>& particles, std::string name, float mass, float3 scale, float3 center, float3 velocity) {
	std::ifstream file(name);
	if (!file.is_open())
	{
		std::cout << "error open file";
	}
	/*
	int x, y, z;
	file >> x >> y >> z;
	std::cout << x << " " << y << " " << z << std::endl;
	for (int i = 0; i < z; ++i) {
	for (int j = 0; j < y;++j) {
	for (int k = 0; k < x;++k) {
	float ccc;
	file >> ccc;
	if(ccc<=0)
	particles.push_back(Particle(center + scale * make_float3(i, j, k), velocity, mass));
	}
	}
	}*/

	while (!file.eof()) {
		int vCount;
		file >> vCount;
		float x, y, z;
		//std::cout << vCount << std::endl;
		for (int i = 0; i < vCount; ++i) {
			file >> x >> y >> z;
			particles.push_back(Particle(center + make_float3(x* scale.x, y * scale.y, z* scale.z), velocity, mass));
		}
	}
	file.close();
}

#endif