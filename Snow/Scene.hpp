#ifndef SCENE_H
#define SCENE_H

#include "Common.h"
#include "Parameters.h"
#include "SetupFunctions.hpp"
#include "Particle.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include "stb_image.h"
using namespace std;
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
		sp->gravity = make_float3(0, -9.80f, 0);

		sp->frictionCoeff = 0.5f;
		sp->stickyWalls = false;

		sp->terrainScale = make_float3(2.55f, 2.55f, 2.55f);
		sp->terrainTransform = make_float3(-1.0f, 0.0f, -1.0f);
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

		float3 lower = make_float3(0.2f, 1.5f, 0.2f);
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
		createSnowball(particles, make_float3(1.0f, 1.2f, 1.0f), snowDims, restDistance, getMass(sp->radius, sp->density), make_float3(0, -5.0f, 0));
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

class Landslide : public Scene {
public:
	Landslide(std::string terBeforePath,std::string terrAfterPath) : Scene(terBeforePath), terrainBeforePath(terBeforePath), terrainAfterPath(terrAfterPath){}

	virtual void init(std::vector<Particle>& particles, solverParams* sp) {
		Scene::init(particles, sp);
		const float restDistance = sp->radius * 1.0f;

		int3 dims = make_int3(150);

		sp->boxCorner1 = make_float3(0, 0.0f, 0);
		sp->boxCorner2 = make_float3((dims.x) * sp->radius, (dims.y) * sp->radius, (dims.z) * sp->radius);
		float mass = getMass(sp->radius, sp->density);
		float3 velocity = make_float3(0.0f, 0.0f, 0.0f);

		int width, height, nrChannels;
		stbi_us *rawImgAfter = stbi_load_16(terrainAfterPath.data(), &width, &height, &nrChannels, 0);
		if (!rawImgAfter)
		{
			cout << "open terrainAfterPath HeightMap file fail" << endl;
			return;
		}
		stbi_us *rawImgBefore = stbi_load_16(terrainBeforePath.data(), &width, &height, &nrChannels, 0);
		if (!rawImgBefore)
		{
			cout << "open terrainBeforePath HeightMap file fail" << endl;
			return;
		}
		float rrr = sp->radius ;
		int pacing = (int)(rrr / sp->boxCorner2.x * 2049);
		for (int i = 10; i < height -1; i+= pacing) {
			for (int j = 10; j < width -1; j+= pacing) {
				if (rawImgBefore[j*width * nrChannels + i * nrChannels]- rawImgAfter[j*width * nrChannels + i * nrChannels]>0) {
					float start = rawImgAfter[j*width * nrChannels + i * nrChannels] / (float)65535*sp->terrainScale.y;
					float end = rawImgBefore[j*width * nrChannels + i * nrChannels] / (float)65535* sp->terrainScale.y;
					for (float h = start+ sp->radius; h < end-0.0311;h+= rrr) {
						float r1 = 0.001f + static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
						float r2 = 0.001f + static_cast <float>(rand()) / static_cast <float> (RAND_MAX);
						float r3 = 0.001f + static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
						float3 jitter = make_float3(r1, r2, r3) * rrr;
						particles.push_back(Particle(make_float3(i / (float)height*sp->terrainScale.x, h, j / (float)width *sp->terrainScale.z)+jitter, velocity, mass));
					}
				}
			}
		}		
		sp->numParticles = int(particles.size());
		sp->gBounds = dims;
		sp->gridSize = dims.x * dims.y * dims.z;
	}
private:
	std::string terrainBeforePath, terrainAfterPath;
};






#endif