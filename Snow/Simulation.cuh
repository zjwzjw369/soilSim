#ifndef SIMULATION_CUH
#define SIMULATION_CUH

#include "Parameters.h"
#include "Particle.h"
#include "Cell.h"

void update(Particle* particles, Cell* cells, int gridSize);
void getPositions(float* positionsPtr, Particle* particles);
void setParams(solverParams *tempParams);
void setTerrainTex(std::string terrainHeightPath, std::string terrainNormalpath);
void setTerrainTex16(std::string terrainHeightPath16, std::string terrainNormalpath);
#endif