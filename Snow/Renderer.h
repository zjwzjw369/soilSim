#ifndef RENDERER_H
#define RENDERER_H

#include "Common.h"
#include "Camera.hpp"
#include "Shader.h"
#include <cuda_gl_interop.h>
#include "Parameters.h"

struct snowBuffers {
	GLuint vao;
	GLuint positions;
	int numParticles;
};

struct planeBuffers {
	GLuint vao;
	GLuint vbo;
	GLuint ebo;
};

class Renderer {
public:
	cudaGraphicsResource *resource;

	Renderer(int width, int height, solverParams* sp);
	~Renderer();

	void setProjection(glm::mat4 &projection);
	void initSnowBuffers(int numParticles);
	void render(Camera& cam);
	void initTerrain(std::string rawFilename, std::string texFilename, int TerrainSize);//filename为terrain.raw所在的路径、TerrainSize为raw的大小 如大小为512*512，输入512即可

private:
	solverParams* sp;
	glm::mat4 mView, projection;
	int width, height;
	float aspectRatio;
	Shader plane;
	Shader snow;
	snowBuffers snowBuffers;
	planeBuffers wallBuffers;
	planeBuffers floorBuffers;

	void renderPlane(planeBuffers &buf);
	void renderSnow(Camera& cam);
};

#endif