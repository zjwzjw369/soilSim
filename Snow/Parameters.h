#ifndef PARAMETERS_H
#define PARAMETERS_H

#include "Common.h"
struct rigid {
	//sphere struct
	float3 center;
	float r;
	float3 vBall;
	float3 force;
};
struct solverParams {
	float deltaT;
	float radius;
	float compression;
	float stretch;
	float hardening;
	float young;
	float poisson;
	float alpha;
	float density;

	rigid* rig;
	//lame parameters
	float lambda;
	float mu;

	float h0, h1, h2, h3;

	int numParticles;

	int gridSize;
	int3 gBounds;
	float3 gravity;

	float3 boxCorner1;
	float3 boxCorner2;
	float frictionCoeff;
	bool stickyWalls;
	float3 terrainScale;
	float3 terrainTransform;
};

#endif