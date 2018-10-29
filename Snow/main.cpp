#define GLEW_DYNAMIC
#include <GL/glew.h>
#include "Common.h"
#include <GLFW/glfw3.h>
#include "Camera.hpp"
#include "ParticleSystem.h"
#include "Renderer.h"
#include "Scene.hpp"
#include <cuda_profiler_api.h>

#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>

using namespace std;

static const int width = 1280;
static const int height = 720;
static const GLfloat lastX = (width / 2);
static const GLfloat lastY = (height / 2);
static float deltaTime = 0.0f;
static float lastFrame = 0.0f;
static int frameCounter = -1;
static int w = 0;
static bool video = false;
static bool startVideo = false;
static bool paused = true;
static bool spacePressed = false;

void handleInput(GLFWwindow* window, ParticleSystem &system, Camera &cam);
void mainUpdate(ParticleSystem& system, Renderer& renderer, Camera& cam, solverParams& params);

int main() {
	//Checks for memory leaks in debug mode
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);

	cudaGLSetGLDevice(0);

	if (!glfwInit()) exit(EXIT_FAILURE);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	glfwWindowHint(GLFW_RESIZABLE, GL_FALSE);
	GLFWwindow* window = glfwCreateWindow(width, height, "Snow Simulation", nullptr, nullptr);
	if (!window) {
		glfwTerminate();
		exit(EXIT_FAILURE);
	}

	glfwMakeContextCurrent(window);

	//Set callbacks for keyboard and mouse
	glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);

	glewExperimental = GL_TRUE;
	glewInit();
	glGetError();

	//Define the viewport dimensions
	glViewport(0, 0, width, height);

	//Create scenes
	vector<Particle> particles;

	//GroundSmash scene("GroundSmash");
	//SnowballDrop scene("SnowballDrop");
	//CustomSceneLoad scene("H:\\Research\\mpm\\project\\Snow2\\soilSim-master\\Snow\\model\\ZJU.txt");

	CustomSceneLoad scene(".\\model\\ZJU.txt");
	//WallSmash scene("Wallsmash");
	solverParams sp;

	scene.init(particles, &sp);

	Camera cam = Camera();
	cam.eye = glm::vec3(sp.boxCorner2.x, sp.boxCorner2.y, sp.boxCorner2.z);
	Renderer renderer = Renderer(width, height, &sp);
	renderer.setProjection(glm::infinitePerspective(cam.zoom, float(width) / float(height), 0.1f));

	ParticleSystem system = ParticleSystem(particles, sp);

	renderer.initTerrain(".\\model\\terrain\\terrain1.raw", ".\\model\\terrain\\terrain1.jpg", 513);
	//Initialize buffers for drawing snow
	renderer.initSnowBuffers(sp.numParticles);
	//Take 1 step for initialization
	system.updateWrapper(sp);

	//Tell ffmpeg to expect raw rgba 720p-60hz frames
	//-i - tells it to read frames from stdin
	string cmd = "output\\ffmpeg.exe -r " + std::to_string((1.0f / 60.0f) / sp.deltaT) + " -f rawvideo -pix_fmt rgba -s 1280x720 -i - "
		"-threads 0 -preset fast -y -pix_fmt yuv420p -crf 21 -vf vflip output\\snow.mp4";

	//open pipe to ffmpeg's stdin in binary write mode
	FILE* ffmpeg;
	int* buffer = new int[width*height];

	typedef BOOL(APIENTRY *PFNWGLSWAPINTERVALFARPROC)(int);
	PFNWGLSWAPINTERVALFARPROC wglSwapIntervalEXT = 0;
	wglSwapIntervalEXT = (PFNWGLSWAPINTERVALFARPROC)wglGetProcAddress("wglSwapIntervalEXT");
	//wglSwapIntervalEXT(1);//�򿪴�ֱͬ��������֡��
	wglSwapIntervalEXT(0);//�رմ�ֱͬ������ַ����Կ�����Ⱦ����


	float minX = particles[0].pos.x, maxX = particles[0].pos.x, minY = particles[0].pos.y, maxY = particles[0].pos.y, minZ = particles[0].pos.z, maxZ = particles[0].pos.z;

	for (int i = 0; i < particles.size(); ++i)
	{
		auto pos = particles[i].pos;
		if (pos.x < minX)
			minX = pos.x;
		if (pos.x > maxX)
			maxX = pos.x;
		if (pos.y < minY)
			minY = pos.y;
		if (pos.y > maxY)
			maxY = pos.y;
		if (pos.z < minZ)
			minZ = pos.z;
		if (pos.z > maxZ)
			maxZ = pos.z;
	}
	float XX = (minX + maxX) / 2;
	float YY = (minY + maxY) / 2;
	float ZZ = (minZ + maxZ) / 2;
	//bool isfffff = true;
	std::ofstream openfile("..\\VDBDatas\\vdbpoint.txt", std::ios::trunc);
	int countFrame = 0;
	int countT = 0;

	while (!glfwWindowShouldClose(window)) {
		//Set frame times
		float currentFrame = float(glfwGetTime());
		deltaTime = currentFrame - lastFrame;
		lastFrame = currentFrame;

		char fps[256];
		float ifps = 1.f / deltaTime;
		sprintf(fps, "mpm: %3.1f fps", ifps);
		glfwSetWindowTitle(window, fps);

		//Check and call events
		glfwPollEvents();
		handleInput(window, system, cam);

		//Step physics and render
		mainUpdate(system, renderer, cam, sp);


		if (countFrame++ == 800)
		{
			countFrame = 0;

			for (int i = 0; i < particles.size(); ++i)
			{
				openfile << particles[i].pos.x * 5000 - XX * 5000 << " " << particles[i].pos.y * 5000 - YY * 5000 << " " << particles[i].pos.z * 5000 - ZZ * 5000 << endl;
			}
			openfile << 99999.9f << " " << 99999.9f << " " << 99999.9f << endl;
			countT++;
		}

		if (countT == 3)
		{
			openfile.close();
			cout << "out put vdbpoint.txt finish";
		}



		//Swap the buffers
		glfwSwapBuffers(window);

		//Save video if turned on
		if (startVideo == true) {
			ffmpeg = _popen(cmd.c_str(), "wb");


			startVideo = false;
		}
		if (!paused) {
			if (frameCounter % (int)(1 / (sp.deltaT * 30 * 3)) == 0) {
				cout << lastFrame << endl;

			}
		}
		if (video == true) {
			glReadPixels(0, 0, width, height, GL_RGBA, GL_UNSIGNED_BYTE, buffer);
			fwrite(buffer, sizeof(int)* width * height, 1, ffmpeg);
		}


		//glfwSetCursorPos(window, lastX, lastY);
	}

	_pclose(ffmpeg);

	glfwTerminate();

	return 0;
}

void handleInput(GLFWwindow* window, ParticleSystem &system, Camera &cam) {
	if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
		glfwSetWindowShouldClose(window, GL_TRUE);

	if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
		cam.wasdMovement(FORWARD, deltaTime);

	if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
		cam.wasdMovement(BACKWARD, deltaTime);

	if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
		cam.wasdMovement(RIGHT, deltaTime);

	if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
		cam.wasdMovement(LEFT, deltaTime);

	if (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS)
		cam.wasdMovement(UP, deltaTime);

	if (glfwGetKey(window, GLFW_KEY_LEFT_CONTROL) == GLFW_PRESS)
		cam.wasdMovement(DOWN, deltaTime);

	if (glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS)
		spacePressed = true;

	if (glfwGetKey(window, GLFW_KEY_COMMA) == GLFW_PRESS) {
		startVideo = true;
		video = true;
	}

	if (glfwGetKey(window, GLFW_KEY_PERIOD) == GLFW_PRESS)
		video = false;

	if (glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_RELEASE) {
		if (spacePressed) {
			spacePressed = false;
			if (paused) {
				cout << "Running Simulation..." << endl;
			}
			else {
				cout << "Pausing Simulation..." << endl;
			}
			paused = !paused;
		}
	}

	double xpos, ypos;
	glfwGetCursorPos(window, &xpos, &ypos);
	if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS)
		cam.mouseMovement((float(xpos) - lastX), (lastY - float(ypos)), deltaTime);
}

void mainUpdate(ParticleSystem& system, Renderer& renderer, Camera& cam, solverParams& params) {
	//Step physics
	if (!paused) {
		system.updateWrapper(params);
		frameCounter++;
	}
	//Update the VBO
	void* positionsPtr;
	cudaCheck(cudaGraphicsMapResources(1, &renderer.resource));
	size_t size;
	cudaGraphicsResourceGetMappedPointer(&positionsPtr, &size, renderer.resource);
	system.getPositionsWrapper((float*)positionsPtr);
	cudaGraphicsUnmapResources(1, &renderer.resource, 0);

	//Render
	renderer.render(cam);
}