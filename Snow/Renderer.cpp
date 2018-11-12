#include "Renderer.h"
#include <vector>
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
using namespace std;

static const float radius = 0.008f;
vector<GLfloat> terrainVertex;
vector<GLfloat> terrainTex;
Renderer::Renderer(int width, int height, solverParams* sp) :
	width(width),
	height(height),
	plane(Shader("plane.vert", "plane.frag")),
	snow(Shader("snow.vert", "snow.frag")),
	terrain(Shader("terrain.vert", "terrain.frag"))
{
	this->sp = sp;
	aspectRatio = float(width) / float(height);

	GLfloat floorVertices[] = {
		sp->boxCorner2.x, sp->boxCorner1.y, sp->boxCorner2.z,
		sp->boxCorner2.x, sp->boxCorner1.y, sp->boxCorner1.z,
		sp->boxCorner1.x, sp->boxCorner1.y, sp->boxCorner1.z,
		sp->boxCorner1.x, sp->boxCorner1.y, sp->boxCorner2.z
	};

	GLfloat wallVertices[] = {
		sp->boxCorner1.x, sp->boxCorner1.y, sp->boxCorner2.z,
		sp->boxCorner1.x, sp->boxCorner1.y, sp->boxCorner1.z,
		sp->boxCorner1.x, sp->boxCorner2.y, sp->boxCorner1.z,
		sp->boxCorner1.x, sp->boxCorner2.y, sp->boxCorner2.z
	};

	GLuint indices[] = {
		0, 1, 3,
		1, 2, 3
	};

	//Wall
	glGenVertexArrays(1, &wallBuffers.vao);

	glGenBuffers(1, &wallBuffers.vbo);
	glBindBuffer(GL_ARRAY_BUFFER, wallBuffers.vbo);
	glBufferData(GL_ARRAY_BUFFER, sizeof(wallVertices), wallVertices, GL_STATIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	glGenBuffers(1, &wallBuffers.ebo);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, wallBuffers.ebo);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

	//Floor
	glGenVertexArrays(1, &floorBuffers.vao);
	//cout << sp->boxCorner2.x << " " << sp->boxCorner2.y << " " << sp->boxCorner2.z << endl;
	glGenBuffers(1, &floorBuffers.vbo);
	glBindBuffer(GL_ARRAY_BUFFER, floorBuffers.vbo);
	glBufferData(GL_ARRAY_BUFFER, sizeof(floorVertices), floorVertices, GL_STATIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	glGenBuffers(1, &floorBuffers.ebo);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, floorBuffers.ebo);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

}

Renderer::~Renderer() {
	if (snowBuffers.vao != 0) {
		glDeleteVertexArrays(1, &snowBuffers.vao);
		glDeleteBuffers(1, &snowBuffers.positions);
	}

	if (wallBuffers.vao != 0) {
		glDeleteVertexArrays(1, &wallBuffers.vao);
		glDeleteBuffers(1, &wallBuffers.vbo);
		glDeleteBuffers(1, &wallBuffers.ebo);
	}

	if (floorBuffers.vao != 0) {
		glDeleteVertexArrays(1, &floorBuffers.vao);
		glDeleteBuffers(1, &floorBuffers.vbo);
		glDeleteBuffers(1, &floorBuffers.ebo);
	}
}

void Renderer::setProjection(glm::mat4 &projection) {
	this->projection = projection;
}

void Renderer::initTerrainBuffers() {

	glGenVertexArrays(1, &terBuffers.vao);
	glGenBuffers(1, &terBuffers.vbo);
	glBindBuffer(GL_ARRAY_BUFFER, terBuffers.vbo);
	glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*terrainVertex.size(), terrainVertex.data(), GL_STATIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glGenBuffers(1, &terBuffers.tbo);
	glBindBuffer(GL_ARRAY_BUFFER, terBuffers.tbo);
	glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*terrainTex.size(), terrainTex.data(), GL_STATIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	cout << terrainVertex.size() << endl;
	cout << terrainTex.size() << endl;

}

void Renderer::initSnowBuffers(int numParticles) {
	glGenVertexArrays(1, &snowBuffers.vao);

	glGenBuffers(1, &snowBuffers.positions);
	glBindBuffer(GL_ARRAY_BUFFER, snowBuffers.positions);
	glBufferData(GL_ARRAY_BUFFER, numParticles * 6 * sizeof(float), 0, GL_DYNAMIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	
	cudaGraphicsGLRegisterBuffer(&resource, snowBuffers.positions, cudaGraphicsRegisterFlagsWriteDiscard);

	snowBuffers.numParticles = numParticles;
}

void Renderer::render(Camera& cam) {
	//Set model view matrix
	mView = cam.getMView();

	//Clear buffer
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	//Plane
	renderPlane(wallBuffers);
	renderPlane(floorBuffers);

	//Snow
	renderSnow(cam);
	renderTerrain();

}

void Renderer::renderPlane(planeBuffers &buf) {
	glUseProgram(plane.program);
	glBindFramebuffer(GL_FRAMEBUFFER, 0);

	plane.setUniformmat4("mView", mView);
	plane.setUniformmat4("projection", projection);

	glBindVertexArray(buf.vao);
	glBindBuffer(GL_ARRAY_BUFFER, buf.vbo);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buf.ebo);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
	glEnableVertexAttribArray(0);
	//glBindTexture(GL_TEXTURE_2D, terBuffers.tId);

	glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);
}

void Renderer::initTerrain(std::string rawFilename, ::string texFilename) {
	int width, height, nrChannels;
	unsigned char *rawImg = stbi_load(rawFilename.data(), &width, &height, &nrChannels, 0);
	if (!rawImg)
	{
		cout << "open Terrain HeightMap file fail" << endl;
		return;
	}
	terrainVertex.clear();
	terrainTex.clear();
	terrainModel = glm::scale(terrainModel, glm::vec3(sp->terrainScale.x, sp->terrainScale.y, sp->terrainScale.z));

	//terrainModel = glm::translate(terrainModel, glm::vec3(sp->terrainTransform.x, sp->terrainTransform.y, sp->terrainTransform.z));
	for (int i = 0; i < height - 1; ++i) {
		for (int j = 0; j < width - 1; ++j) {
			float3 v1 = make_float3(i / (float)height, rawImg[j*width * 3 + i * 3] / (float)255, j / (float)width);
			float3 v2 = make_float3(i / (float)height, rawImg[j*width * 3 + i * 3] / (float)255, (j + 1) / (float)width);
			float3 v3 = make_float3((i + 1) / (float)height, rawImg[j*width * 3 + i * 3] / (float)255, (j + 1) / (float)width);
			float3 v4 = make_float3((i + 1) / (float)height, rawImg[j*width * 3 + i * 3] / (float)255, j / (float)width);
			terrainVertex.push_back(v1.x);
			terrainVertex.push_back(v1.y);
			terrainVertex.push_back(v1.z);
			terrainVertex.push_back(v2.x);
			terrainVertex.push_back(v2.y);
			terrainVertex.push_back(v2.z);
			terrainVertex.push_back(v3.x);
			terrainVertex.push_back(v3.y);
			terrainVertex.push_back(v3.z);
			terrainVertex.push_back(v1.x);
			terrainVertex.push_back(v1.y);
			terrainVertex.push_back(v1.z);
			terrainVertex.push_back(v3.x);
			terrainVertex.push_back(v3.y);
			terrainVertex.push_back(v3.z);
			terrainVertex.push_back(v4.x);
			terrainVertex.push_back(v4.y);
			terrainVertex.push_back(v4.z);
		}
	}
	for (int i = 0; i < terrainVertex.size() / 3; i++) {
		float2 tmp = make_float2(terrainVertex[i * 3], terrainVertex[i * 3 + 2]);
		terrainTex.push_back(tmp.x);
		terrainTex.push_back(tmp.y);
		//cout << tmp.x << "  " << tmp.y<<" "<< terrainVertex[i * 3 + 1] <<endl;
	}

	unsigned char *texData = stbi_load(texFilename.data(), &width, &height, &nrChannels, 0);
	if (!texData) {
		cout << "open Terrain Texture file fail" << endl;
		return;
	}
	glGenTextures(1, &terBuffers.tId);
	glBindTexture(GL_TEXTURE_2D, terBuffers.tId);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, texData);
	glGenerateMipmap(GL_TEXTURE_2D);

	stbi_image_free(rawImg);
	stbi_image_free(texData);
}

void Renderer::initTerrain16(std::string rFilename16, ::string texFilename) {
	int width, height, nrChannels;
	stbi_us *rawImg = stbi_load_16(rFilename16.data(), &width, &height, &nrChannels, 0);
	if (!rawImg)
	{
		cout << "open 16bit Terrain HeightMap file fail" << endl;
		return;
	}
	terrainVertex.clear();
	terrainTex.clear();
	terrainModel = glm::scale(terrainModel, glm::vec3(sp->terrainScale.x, sp->terrainScale.y, sp->terrainScale.z));
	for (int i = 0; i < height - 1; ++i) {
		for (int j = 0; j < width - 1; ++j) {
			float3 v1 = make_float3(i / (float)height, rawImg[j*width * nrChannels + i * nrChannels] / (float)65535, j / (float)width);
			float3 v2 = make_float3(i / (float)height, rawImg[j*width * nrChannels + i * nrChannels] / (float)65535, (j + 1) / (float)width);
			float3 v3 = make_float3((i + 1) / (float)height, rawImg[j*width * nrChannels + i * nrChannels] / (float)65535, (j + 1) / (float)width);
			float3 v4 = make_float3((i + 1) / (float)height, rawImg[j*width * nrChannels + i * nrChannels] / (float)65535, j / (float)width);
			terrainVertex.push_back(v1.x);
			terrainVertex.push_back(v1.y);
			terrainVertex.push_back(v1.z);
			terrainVertex.push_back(v2.x);
			terrainVertex.push_back(v2.y);
			terrainVertex.push_back(v2.z);
			terrainVertex.push_back(v3.x);
			terrainVertex.push_back(v3.y);
			terrainVertex.push_back(v3.z);
			terrainVertex.push_back(v1.x);
			terrainVertex.push_back(v1.y);
			terrainVertex.push_back(v1.z);
			terrainVertex.push_back(v3.x);
			terrainVertex.push_back(v3.y);
			terrainVertex.push_back(v3.z);
			terrainVertex.push_back(v4.x);
			terrainVertex.push_back(v4.y);
			terrainVertex.push_back(v4.z);
		}
	}
	for (int i = 0; i < terrainVertex.size() / 3; i++) {
		float2 tmp = make_float2(terrainVertex[i * 3], terrainVertex[i * 3 + 2]);
		terrainTex.push_back(tmp.x);
		terrainTex.push_back(tmp.y);
	}

	unsigned char *texData = stbi_load(texFilename.data(), &width, &height, &nrChannels, 0);
	if (!texData) {
		cout << "open Terrain Texture file fail" << endl;
		return;
	}
	glGenTextures(1, &terBuffers.tId);
	glBindTexture(GL_TEXTURE_2D, terBuffers.tId);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, texData);
	glGenerateMipmap(GL_TEXTURE_2D);

	stbi_image_free(rawImg);
	stbi_image_free(texData);
}

void Renderer::renderTerrain()
{

	glUseProgram(terrain.program);

	terrain.setUniformmat4("mView", mView);
	terrain.setUniformmat4("projection", projection);
	terrain.setUniformmat4("model", terrainModel);

	glBindVertexArray(terBuffers.vao);
	glBindBuffer(GL_ARRAY_BUFFER, terBuffers.vbo);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
	glEnableVertexAttribArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, terBuffers.tbo);
	glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 0, 0);
	glEnableVertexAttribArray(1);
	glBindTexture(GL_TEXTURE_2D, terBuffers.tId);
	//cout << terrainVertex.size() << endl;
	glDrawArrays(GL_TRIANGLES, 0, terrainVertex.size() / 3);

}

void Renderer::renderSnow(Camera& cam) {
	glUseProgram(snow.program);
	glBindFramebuffer(GL_FRAMEBUFFER, 0);

	snow.setUniformmat4("mView", mView);
	snow.setUniformmat4("projection", projection);
	snow.setUniformf("pointRadius", radius);
	snow.setUniformf("pointScale", width / aspectRatio * (1.0f / tanf(cam.zoom * 0.5f)));

	glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);
	//glEnable(GL_CULL_FACE);
	glEnable(GL_DEPTH_TEST);

	//Draw snow
	glBindVertexArray(snowBuffers.vao);
	glBindBuffer(GL_ARRAY_BUFFER, snowBuffers.positions);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat), 0);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat), (GLvoid*)(3 * sizeof(GLfloat)));
	glEnableVertexAttribArray(1);

	glDrawArrays(GL_POINTS, 0, GLsizei(snowBuffers.numParticles));
}
