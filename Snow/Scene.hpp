#ifndef SCENE_H
#define SCENE_H

#include "Common.h"
#include "Parameters.h"
#include "SetupFunctions.hpp"
#include "Particle.h"
#include <iostream>
#include <sstream>
#include <fstream>
class Scene {
public:
	Scene(std::string name) : name(name) {}
	virtual void init(std::vector<Particle>& particles, solverParams* sp) {
		sp->deltaT = 5e-5f;
		sp->radius = 0.017f;
		sp->compression = 0.25f;
		sp->stretch = 0.010f;
		sp->hardening = 10.0f;
		sp->young = 3.5e5f;
		sp->poisson = 0.3f;
		sp->alpha = 0.95f;
		sp->density = 2200.0f;

		sp->lambda = getLambda(sp->poisson, sp->young);
		sp->mu = getMu(sp->poisson, sp->young);

		sp->h0 = 35; sp->h1 = 0; sp->h2 = 0.2; sp->h3 = 10;
		sp->gravity = make_float3(0, -9.8f, 0);

		sp->frictionCoeff = 1.0f;
		sp->stickyWalls = false;
	}

	float getLambda(float poisson, float young) {
		return (poisson * young) / ((1 + poisson) * (1 - 2 * poisson));
	}

	float getMu(float poisson, float young) {
		return young / (2 * (1 + poisson));
	}

	float getMass(float radius, float density) {
		return pow(radius, 3) * density / 4;
	}

	std::string name;
};

class GroundSmash : public Scene {
public:
	GroundSmash(std::string name) : Scene(name) {}

	virtual void init(std::vector<Particle>& particles, solverParams* sp) {
		Scene::init(particles, sp);

		const float restDistance = sp->radius * 1.f;

		int3 dims = make_int3(150);
		int3 snowDims = make_int3(40);

		sp->boxCorner1 = make_float3(0, 0.0f, 0);
		sp->boxCorner2 = make_float3((dims.x) * sp->radius, (dims.y) * sp->radius, (dims.z) * sp->radius);

		float3 lower = make_float3(dims.x / 2 * sp->radius, 0.5f, dims.z / 2 * sp->radius);
		createParticleGrid(particles, sp, lower, snowDims, restDistance, getMass(sp->radius, sp->density), make_float3(0, -5, 0));



		sp->numParticles = int(particles.size());
		sp->gridSize = dims.x * dims.y * dims.z;
		sp->gBounds = dims;
	}
};

class SnowballDrop : public Scene {
public:
	SnowballDrop(std::string name) : Scene(name) {}

	virtual void init(std::vector<Particle>& particles, solverParams* sp) {
		Scene::init(particles, sp);
		const float restDistance = sp->radius * 1.0f;

		int3 dims = make_int3(150);

		sp->boxCorner1 = make_float3(0, 0.0f, 0);
		sp->boxCorner2 = make_float3((dims.x) * sp->radius, (dims.y) * sp->radius, (dims.z) * sp->radius);

		int3 snowDims = make_int3(30);
		createSnowball(particles, make_float3(1.25f, 1.0f, 1.25f), snowDims, restDistance, getMass(sp->radius, sp->density), make_float3(0, 0.0f, 0));
		//createSnowball(particles, make_float3(1.25f, 2.0f, 1.25f), snowDims, restDistance, getMass(sp->radius, sp->density), make_float3(0, 0.0f, 0));

		sp->numParticles = int(particles.size());
		sp->gBounds = dims;
		sp->gridSize = dims.x * dims.y * dims.z;
	}
};

class CustomSceneLoad : public Scene {//name for path
public:
	CustomSceneLoad(std::string name) : Scene(name){}
	
	virtual void init(std::vector<Particle>& particles, solverParams* sp) {
		Scene::init(particles, sp);
		const float restDistance = sp->radius * 1.0f;

		int3 dims = make_int3(150);

		sp->boxCorner1 = make_float3(0, 0.0f, 0);
		sp->boxCorner2 = make_float3((dims.x) * sp->radius, (dims.y) * sp->radius, (dims.z) * sp->radius);
		
		std::ifstream file(name);
		if (!file.is_open())
		{
			std::cout << "error open file";
		}
		float3 transform = make_float3(1.0);
		float mass = getMass(sp->radius, sp->density);
		float3 velocity = make_float3(0.0f, -5.0f, 0.0f);
		while (!file.eof()) {
			int vCount;
			file >> vCount;
			float x, y, z;
			//std::cout << vCount << std::endl;
			for (int i = 0; i < vCount;++i) {
				file >> x >> y >> z;
				particles.push_back(Particle(transform +make_float3(x, y+0.5, z), velocity, mass));
			}
		}
		//int3 snowDims = make_int3(30);
		//createSnowball(particles, make_float3(1.25f, 1.0f, 1.25f), snowDims, restDistance, getMass(sp->radius, sp->density), make_float3(0, 0.0f, 0));
		//createSnowball(particles, make_float3(1.25f, 2.0f, 1.25f), snowDims, restDistance, getMass(sp->radius, sp->density), make_float3(0, 0.0f, 0));

		sp->numParticles = int(particles.size());
		sp->gBounds = dims;
		sp->gridSize = dims.x * dims.y * dims.z;
	}
};


class WallSmash : public Scene {
public:
	WallSmash(std::string name) : Scene(name) {}

	virtual void init(std::vector<Particle>& particles, solverParams* sp) {
		Scene::init(particles, sp);
		const float restDistance = sp->radius * 1.0f;

		int3 dims = make_int3(150);

		sp->boxCorner1 = make_float3(0.0f, 0.0f, 0);
		sp->boxCorner2 = make_float3((dims.x) * sp->radius, (dims.y) * sp->radius, (dims.z) * sp->radius);

		int3 snowDims = make_int3(30);
		createSnowball(particles, make_float3(1.0f, 1.0f, 1.25f), snowDims, restDistance, getMass(sp->radius, sp->density), make_float3(-10.0f, 0.0f, 0));

		sp->numParticles = int(particles.size());
		sp->gBounds = dims;
		sp->gridSize = dims.x * dims.y * dims.z;
		sp->stickyWalls = true;
	}
};

#endif