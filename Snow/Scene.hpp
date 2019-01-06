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
vector<int> particleTag;
class Scene {
public:
	Scene(std::string name) : name(name) {}
	virtual void init(std::vector<Particle>& particles, solverParams* sp) {
		sp->deltaT = 5e-5f;
		sp->radius = 0.017f;
		sp->compression = 0.25f;
		sp->stretch = 0.010f;
		sp->hardening = 10.0f;
		sp->young = 5e5;
		sp->poisson = 0.3f;
		sp->alpha = 0.95f;
		sp->density = 2000.0f;

		//³õÊ¼»¯¸ÕÌå
		sp->rig = new rigid();
		sp->rig->r = 0.2;
		sp->rig->vBall = make_float3(0.0f);
		sp->rig->wBall = make_float3(0.0f);
		sp->rig->force = make_float3(0.0f);
		sp->rig->center = make_float3(0.4, 0.56, 0.4);


		sp->lambda = getLambda(sp->poisson, sp->young);
		sp->mu = getMu(sp->poisson, sp->young);

		sp->h0 = 35; sp->h1 = 0; sp->h2 = 0.2; sp->h3 = 10;
		sp->gravity = make_float3(0, -9.80f, 0);

		//sp->frictionCoeff = 0.5f;
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
		const float restDistance = sp->radius;//16

		int3 dims = make_int3(150);

		sp->boxCorner1 = make_float3(0, 0.0f, 0);
		sp->boxCorner2 = make_float3((dims.x) * sp->radius, (dims.y) * sp->radius, (dims.z) * sp->radius);

		//int3 snowDims = make_int3(60.0);
		int3 snowDims = make_int3(30);
		//createTorus(particles, make_float3(1.275f, 1.8f, 1.275f), snowDims, restDistance, getMass(sp->radius, sp->density), make_float3(0, -5.0f, 0));
		//createTorus(particles, make_float3(1.25f, 1.4f, 1.25f), snowDims, restDistance, getMass(sp->radius, sp->density), make_float3(0, -5.0f, 0));
		//createTorus(particles, make_float3(1.225f, 1.0f, 1.225f), snowDims, restDistance, getMass(sp->radius, sp->density), make_float3(0, -5.0f, 0));
		//createSnowball(particles, make_float3(0.7f, 1.3f, 1.375f), snowDims, restDistance, getMass(sp->radius, sp->density), make_float3(0, -5.0f, 0));
		//createSnowball(particles, make_float3(1.0f, 0.6f, 1.475f), snowDims, restDistance, getMass(sp->radius, sp->density), make_float3(0, -5.0f, 0));
		//createSnowball(particles, make_float3(1.275f, 1.8f , 1.275f ), snowDims, restDistance, getMass(sp->radius, sp->density), make_float3(0, -10.0f, 0));
		//createParticleGrid(particles, sp, make_float3(0.0,0.02,0.0), make_int3(300,60,300), restDistance, getMass(sp->radius, sp->density), make_float3(0, 0, 0));
		//createSnowball(particles, make_float3(1.275f, 1.0f , 1.275f), snowDims, restDistance, getMass(sp->radius, sp->density), make_float3(0, 0.0f, 0));
		//createSnowball(particles, make_float3(1.25f, 2.0f, 1.25f), snowDims, restDistance, getMass(sp->radius, sp->density), make_float3(0, 0.0f, 0));
		//createParticleSphereGrid(particles, sp, make_float3(1.25,0.01,1.25), make_int3(100,2200,100), restDistance, getMass(sp->radius/8, sp->density), make_float3(0, 0, 0), sp->radius/1.5);//8


		//createSnowball(particles, make_float3(1.0f, 0.5f, 1.0f), snowDims, restDistance, getMass(sp->radius, sp->density), make_float3(0, -5.0f, 0));
		createParticleGrid(particles, sp, make_float3(0.0, 0.00, 0.0), make_int3(50, 20, 50), restDistance, getMass(restDistance, sp->density), make_float3(0, 0, 0));



		sp->numParticles = int(particles.size());
		sp->gBounds = dims;
		sp->gridSize = dims.x * dims.y * dims.z;
	}
};

class CustomSceneLoad : public Scene {//name for path
public:
	CustomSceneLoad(std::string name) : Scene(name) {}

	virtual void init(std::vector<Particle>& particles, solverParams* sp) {
		Scene::init(particles, sp);
		const float restDistance = sp->radius * 1.0f;

		int3 dims = make_int3(150);

		sp->boxCorner1 = make_float3(0, 0.0f, 0);
		sp->boxCorner2 = make_float3((dims.x) * sp->radius, (dims.y) * sp->radius, (dims.z) * sp->radius);

		//createCustomModel(particles, ".\\model\\dragon13w.txt", getMass(sp->radius, sp->density), make_float3(0.5, 0.5, 0.5), make_float3(0.6f, 0.5f, 0.6f), make_float3(0.0f, -5.0f, 0.0f));
		//createCustomModel(particles, ".\\model\\bunny27w.txt", getMass(sp->radius, sp->density), make_float3(3.0, 3.0, 3.0), make_float3(1.0f, 0.2f, 1.0f), make_float3(0.0f, -5.0f, 0.0f));
		//createCustomModel(particles, ".\\model\\bunny27w.txt", getMass(sp->radius, sp->density), make_float3(3.0, 3.0, 3.0), make_float3(1.0f, 0.7f, 1.0f), make_float3(0.0f, -5.0f, 0.0f));
		createCustomModel(particles, ".\\model\\bunny27w.txt", getMass(sp->radius, sp->density), make_float3(3.0, 3.0, 3.0), make_float3(1.0f, 1.2f, 1.0f), make_float3(0.0f, -5.0f, 0.0f));
		//createCustomModel(particles, ".\\model\\1.txt", getMass(sp->radius, sp->density), make_float3(0.005,0.005,0.005), make_float3(0.6f, 1.5f, 0.6f), make_float3(0.0f, -5.0f, 0.0f));
		//createCustomModel(particles, ".\\model\\bunny27w.txt", getMass(sp->radius, sp->density), make_float3(3.0, 3.0, 3.0), make_float3(0.6f, 1.5f, 1.2f), make_float3(0.0f, -5.0f, 0.0f));
		//createCustomModel(particles, ".\\model\\bunny27w.txt", getMass(sp->radius, sp->density), make_float3(3.0, 3.0, 3.0), make_float3(0.6f, 1.5f, 1.8f), make_float3(0.0f, -5.0f, 0.0f));
		//createCustomModel(particles, ".\\model\\bunny27w.txt", getMass(sp->radius, sp->density), make_float3(3.0, 3.0, 3.0), make_float3(1.2f, 1.5f, 0.6f), make_float3(0.0f, -5.0f, 0.0f));
		//createCustomModel(particles, ".\\model\\bunny27w.txt", getMass(sp->radius, sp->density), make_float3(3.0, 3.0, 3.0), make_float3(1.2f, 1.5f, 1.2f), make_float3(0.0f, -5.0f, 0.0f));
		//createCustomModel(particles, ".\\model\\bunny27w.txt", getMass(sp->radius, sp->density), make_float3(3.0, 3.0, 3.0), make_float3(1.2f, 1.5f, 1.8f), make_float3(0.0f, -5.0f, 0.0f));
		//createCustomModel(particles, ".\\model\\bunny27w.txt", getMass(sp->radius, sp->density), make_float3(3.0, 3.0, 3.0), make_float3(1.8f, 1.5f, 0.6f), make_float3(0.0f, -5.0f, 0.0f));
		//createCustomModel(particles, ".\\model\\bunny27w.txt", getMass(sp->radius, sp->density), make_float3(3.0, 3.0, 3.0), make_float3(1.8f, 1.5f, 1.2f), make_float3(0.0f, -5.0f, 0.0f));
		//createCustomModel(particles, ".\\model\\bunny27w.txt", getMass(sp->radius, sp->density), make_float3(3.0, 3.0, 3.0), make_float3(1.8f, 1.5f, 1.8f), make_float3(0.0f, -5.0f, 0.0f));


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
	Landslide(std::string terBeforePath, std::string terrAfterPath) : Scene(terBeforePath), terrainBeforePath(terBeforePath), terrainAfterPath(terrAfterPath) {}

	virtual void init(std::vector<Particle>& particles, solverParams* sp) {
		Scene::init(particles, sp);
		const float restDistance = sp->radius * 1.0f / 6.83;

		int3 dims = make_int3(150);

		sp->boxCorner1 = make_float3(0, 0.0f, 0);
		sp->boxCorner2 = make_float3((dims.x) * sp->radius, (dims.y) * sp->radius, (dims.z) * sp->radius);
		float3 velocity = make_float3(1.f, -1.0f, 1.f);

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
		int pacing = (int)(restDistance / sp->boxCorner2.x * 2049);
		for (int i = 0; i < height - 1; i += pacing) {
			for (int j = 0; j < width - 1; j += pacing) {
				if (rawImgBefore[j*width * nrChannels + i * nrChannels] - rawImgAfter[j*width * nrChannels + i * nrChannels] > 0) {
					float start = rawImgAfter[j*width * nrChannels + i * nrChannels] / (float)65535 * sp->terrainScale.y;
					float end = rawImgBefore[j*width * nrChannels + i * nrChannels] / (float)65535 * sp->terrainScale.y;
					for (float h = start + restDistance; h < end; h += restDistance) {
						float r1 = 0.001f + static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
						float r2 = 0.001f + static_cast <float>(rand()) / static_cast <float> (RAND_MAX);
						float r3 = 0.001f + static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
						float3 jitter = make_float3(r1, r2, r3) * restDistance;
						particles.push_back(Particle(make_float3(i / (float)height*sp->terrainScale.x, h, j / (float)width *sp->terrainScale.z) + jitter, velocity, getMass(restDistance, sp->density)));
						if (rand() / (double)RAND_MAX < 0.0005)
							particleTag.push_back(particles.size());
					}
				}
			}
		}

		//cout << "particleTag " <<particleTag.size()<<endl;
		sp->numParticles = int(particles.size());
		sp->gBounds = dims;
		sp->gridSize = dims.x * dims.y * dims.z;
	}
private:
	std::string terrainBeforePath, terrainAfterPath;
};






#endif