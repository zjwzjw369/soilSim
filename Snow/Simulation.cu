#ifndef SIMULATION_CU
#define SIMULATION_CU

#include "Common.h"
#include "Parameters.h"
#include "Particle.h"
#include "Cell.h"
#include "matrix.h"
#include "decomposition.h"
#include "stb_image.h"
#define cudaCheck(x) { cudaError_t err = x; if (err != cudaSuccess) { printf("Cuda error: %f in %s at %s:%f\n", err, #x, __FILE__, __LINE__); assert(0); } }

using namespace std;

static dim3 particleDims;
static dim3 gridDims;
static const int blockSize = 128;

rp3d::Vector3 gravity(0.0, -9.8, 0.0);
rp3d::DynamicsWorld *world;
rp3d::ProxyShape* proxyShape;
rp3d::RigidBody * body;
rp3d::RigidBody * boundingBox;


bool firstTime = true;
texture<float4, 2, cudaReadModeElementType> TerrainHeightTex;
texture<float4, 2, cudaReadModeElementType> TerrainNormalTex;

__device__ rigid r;
__constant__ solverParams sp;

__device__ float NX(const float &x) {
	if (x < 1.0f) {
		return 0.5f * (x * x * x) - (x * x) + (2.0f / 3.0f);
	}
	else if (x < 2.0f) {
		return (-1.0f / 6.0f) * (x * x * x) + (x * x) - (2.0f * x) + (4.0f / 3.0f);
	}
	else {
		return 0.0f;
	}
}

__device__ float dNX(const float &x) {
	float absx = fabs(x);
	if (absx < 1.0f) {
		return (1.5f * absx * x) - (2.0f * x);
	}
	else if (absx < 2.0f) {
		return -0.5f * (absx * x) + (2.0f * x) - (2.0f * x / absx);
	}
	else {
		return 0.0f;
	}
}

__device__ float weight(const float3 &xpgpDiff) {
	return NX(xpgpDiff.x) * NX(xpgpDiff.y) * NX(xpgpDiff.z);
}

__device__ float3 gradWeight(const float3 &xpgpDiff) {
	//TODO: verify this is correct, especially the 1/h scaling factor due to chain rule
	return (1.0f / sp.radius) * make_float3(dNX(xpgpDiff.x) * NX(fabs(xpgpDiff.y)) * NX(fabs(xpgpDiff.z)),
		NX(fabs(xpgpDiff.x)) * dNX(xpgpDiff.y) * NX(fabs(xpgpDiff.z)),
		NX(fabs(xpgpDiff.x)) * NX(fabs(xpgpDiff.y)) * dNX(xpgpDiff.z));
}

__device__ int getGridIndex(const int3 &pos) {
	return (pos.z * sp.gBounds.y * sp.gBounds.x) + (pos.y * sp.gBounds.x) + pos.x;
}


__device__ void applyCellBoundaryCollisions(Cell* cells, int index, float3& position, float3& velocity) {
	float vn;
	float3 vt;
	float3 normal;

	bool collision;

	//handle the 3 boundary dimensions separately
	for (int i = 0; i < 1; i++) {
		collision = false;

		if (i == 0) {
			if (position.x <= sp.boxCorner1.x) {
				collision = true;
				normal = make_float3(0);
				normal.x = 1.0f;
			}
			else if (position.x >= sp.boxCorner2.x) {
				collision = true;
				normal = make_float3(0);
				normal.x = -1.0f;
			}
		}
		if (i == 1) {
			if (position.y <= sp.boxCorner1.y) {
				collision = true;
				normal = make_float3(0);
				normal.y = 1.0f;
			}
			else if (position.y >= sp.boxCorner2.y) {
				collision = true;
				normal = make_float3(0);
				normal.y = -1.0f;
			}
		}
		if (i == 2) {
			if (position.z <= sp.boxCorner1.z) {
				collision = true;
				normal = make_float3(0);
				normal.z = 1.0f;
			}
			else if (position.z >= sp.boxCorner2.z) {
				collision = true;
				normal = make_float3(0);
				normal.z = -1.0f;
			}
		}

		//球体的碰撞
		float3 centerStar = r.center + sp.deltaT*r.vBall;
		if (length(position - centerStar) < r.r) {
			collision = true;
			normal = normalize(position - centerStar);
			float3 V = r.vBall + cross(r.wBall, position - centerStar);
			float3 vRela = V - velocity;
			//vRela *= 0.8;
			float3 vRelaN = dot(vRela, normal)*normal;
			float3 vRelaT = vRela - vRelaN;
			vRelaT= cross(position - centerStar, vRelaT);
			atomicAdd(&(r.force.x), -vRelaN.x*cells[index].mass);//f=mv 动量守恒
			atomicAdd(&(r.force.y), -vRelaN.y*cells[index].mass);//f=mv 动量守恒
			atomicAdd(&(r.force.z), -vRelaN.z*cells[index].mass);//f=mv 动量守恒
			atomicAdd(&(r.torque.x), -vRelaT.x*cells[index].mass);//动量矩守恒
			atomicAdd(&(r.torque.y), -vRelaT.y*cells[index].mass);
			atomicAdd(&(r.torque.z), -vRelaT.z*cells[index].mass);
			velocity += vRela;//1200.0是球的密度
							  //printf("%lf\n",velocity);
		}

		//结束球的碰撞

		//地形碰撞
		//float3 tmpPos = position - sp.deltaT*velocity/2.0f;
		/*
		float terrainHeight = tex2D(TerrainHeightTex, position.x / sp.terrainScale.x, position.z / sp.terrainScale.z).x;
		if (position.y<terrainHeight*sp.terrainScale.y) {
		collision = true;


		float strength = 8.0f;
		float tl = abs(tex2D(TerrainHeightTex, position.x-0.002f / sp.terrainScale.x, position.z - 0.002f / sp.terrainScale.z).x);
		float l = abs(tex2D(TerrainHeightTex, position.x - 0.002f / sp.terrainScale.x, position.z/ sp.terrainScale.z).x);
		float bl = abs(tex2D(TerrainHeightTex, position.x - 0.002f / sp.terrainScale.x, position.z - 0.002f / sp.terrainScale.z).x);
		float b = abs(tex2D(TerrainHeightTex, position.x / sp.terrainScale.x, position.z + 0.002f / sp.terrainScale.z).x);
		float br = abs(tex2D(TerrainHeightTex, position.x + 0.002f / sp.terrainScale.x, position.z + 0.002f / sp.terrainScale.z).x);
		float r = abs(tex2D(TerrainHeightTex, position.x + 0.002f / sp.terrainScale.x, position.z / sp.terrainScale.z).x);
		float tr = abs(tex2D(TerrainHeightTex, position.x + 0.002f / sp.terrainScale.x, position.z - 0.002f / sp.terrainScale.z).x);
		float t = abs(tex2D(TerrainHeightTex, position.x  / sp.terrainScale.x, position.z - 0.002f / sp.terrainScale.z).x);
		float dX = tr + 2 * r + br - tl - 2 * l - bl;
		float dY = bl + 2 * b + br - tl - 2 * t - tr;
		normal =  make_float3 (dX, 1.0f / strength, dY);


		}
		*/



		if (collision) {
			/*
			vn = dot(velocity, normal);
			vt = velocity - vn*normal;

			if (vn >= 0) {
			continue;
			}
			velocity = vt - vn  *normal;
			*/

			//if separating, do nothing	
			vn = dot(velocity, normal);

			if (vn >= 0) {
				continue;
			}

			//if get here, need to respond to collision in some way
			//if sticky, unconditionally set collision to 0
			if (sp.stickyWalls) {
				velocity = make_float3(0);
				return;
			}

			//get tangential component of velocity
			vt = velocity - vn*normal;

			//see if sticking impulse required
			//TODO: could have separate static and sliding friction
			if (length(vt) <= -sp.frictionCoeff * vn) {
				velocity = make_float3(0);
				return;
			}

			//apply dynamic friction
			velocity = vt + sp.frictionCoeff * vn * normalize(vt);
		}
	}
}
__device__ void applyBoundaryCollisions(float3& position, float3& velocity) {
	float vn;
	float3 vt;
	float3 normal;

	bool collision;

	//handle the 3 boundary dimensions separately
	for (int i = 0; i < 1; i++) {
		collision = false;

		if (i == 0) {
			if (position.x <= sp.boxCorner1.x) {
				collision = true;
				normal = make_float3(0);
				normal.x = 1.0f;
			}
			else if (position.x >= sp.boxCorner2.x) {
				collision = true;
				normal = make_float3(0);
				normal.x = -1.0f;
			}
		}
		if (i == 1) {
			if (position.y <= sp.boxCorner1.y) {
				collision = true;
				normal = make_float3(0);
				normal.y = 1.0f;
			}
			else if (position.y >= sp.boxCorner2.y) {
				collision = true;
				normal = make_float3(0);
				normal.y = -1.0f;
			}
		}
		if (i == 2) {
			if (position.z <= sp.boxCorner1.z) {
				collision = true;
				normal = make_float3(0);
				normal.z = 1.0f;
			}
			else if (position.z >= sp.boxCorner2.z) {
				collision = true;
				normal = make_float3(0);
				normal.z = -1.0f;
			}
		}

		//球体的碰撞
		/*
		float3 centerStar = center;
		if (length(position - centerStar) < r) {
		//collision = true;
		normal = normalize(position - centerStar);
		float restitution = 1.0;

		float3 tmpV= normal*(r - length(position - centerStar))*(1+restitution)*vba;
		velocity += tmpV;
		//atomicAdd(&(force), -tmpV.y*0.00002057);//f=mv 动量守恒
		}
		*/


		//结束球的碰撞


		//地形碰撞
		//float3 tmpPos = position - sp.deltaT*velocity/2.0f;
		/*
		float terrainHeight = tex2D(TerrainHeightTex, position.x / sp.terrainScale.x, position.z / sp.terrainScale.z).x;
		if (position.y<terrainHeight*sp.terrainScale.y) {
		collision = true;


		float strength = 1.0f;
		float tl = abs(tex2D(TerrainHeightTex, position.x-0.002f / sp.terrainScale.x, position.z - 0.002f / sp.terrainScale.z).x);
		float l = abs(tex2D(TerrainHeightTex, position.x - 0.002f / sp.terrainScale.x, position.z/ sp.terrainScale.z).x);
		float bl = abs(tex2D(TerrainHeightTex, position.x - 0.002f / sp.terrainScale.x, position.z - 0.002f / sp.terrainScale.z).x);
		float b = abs(tex2D(TerrainHeightTex, position.x / sp.terrainScale.x, position.z + 0.002f / sp.terrainScale.z).x);
		float br = abs(tex2D(TerrainHeightTex, position.x + 0.002f / sp.terrainScale.x, position.z + 0.002f / sp.terrainScale.z).x);
		float r = abs(tex2D(TerrainHeightTex, position.x + 0.002f / sp.terrainScale.x, position.z / sp.terrainScale.z).x);
		float tr = abs(tex2D(TerrainHeightTex, position.x + 0.002f / sp.terrainScale.x, position.z - 0.002f / sp.terrainScale.z).x);
		float t = abs(tex2D(TerrainHeightTex, position.x  / sp.terrainScale.x, position.z - 0.002f / sp.terrainScale.z).x);
		float dX = tr + 2 * r + br - tl - 2 * l - bl;
		float dY = bl + 2 * b + br - tl - 2 * t - tr;
		normal =  make_float3 (dX, 1.0f / strength, dY);


		}*/
		//地形碰撞结束



		if (collision) {
			/*
			vn = dot(velocity, normal);
			vt = velocity - vn*normal;

			if (vn >= 0) {
			continue;
			}
			velocity = vt - vn  *normal;
			*/

			//if separating, do nothing	
			vn = dot(velocity, normal);

			if (vn >= 0) {
				continue;
			}

			//if get here, need to respond to collision in some way
			//if sticky, unconditionally set collision to 0
			if (sp.stickyWalls) {
				velocity = make_float3(0);
				return;
			}

			//get tangential component of velocity
			vt = velocity - vn*normal;

			//see if sticking impulse required
			//TODO: could have separate static and sliding friction
			if (length(vt) <= -sp.frictionCoeff * vn) {
				velocity = make_float3(0);
				return;
			}

			//apply dynamic friction
			velocity = vt + sp.frictionCoeff * vn * normalize(vt);
		}
	}
}

__device__ mat3 calcStress(const mat3& fe, const mat3& fp) {
	float je = mat3::determinant(fe);
	float jp = mat3::determinant(fp);


	//float expFactor = expf(sp.hardening * (1 - jp));//old
	//float lambda = sp.lambda * expFactor;//old
	//float mu = sp.mu * expFactor;//old
	float lambda = sp.lambda;
	float mu = sp.mu;
	//Polar decomposition of fe using SVD
	mat3 re;
	computePD(fe, re);
	//See MPM tech report
	//float k = 100.0f;
	//float y = 7.0f;
	//return mat3(-k*(1 / pow(je, y) - 1.0));
	//mu = 0;
	return (2.0f * mu * mat3::multiplyABt(fe - re, fe)) + mat3(lambda * (je - 1) * je);
}

__device__ mat3 getVelocityGradient(Particle& particle, Cell* cells) {
	//Grad of velocity of particle p = sum over cell i of (v_inext * grad weight i for particle p)
	float hInv = 1.0f / sp.radius;
	int3 pos = make_int3(particle.pos * hInv);
	mat3 velocityGrad(0.0f);

	for (int z = -2; z < 3; z++) {
		for (int y = -2; y < 3; y++) {
			for (int x = -2; x < 3; x++) {
				int3 n = make_int3(pos.x + x, pos.y + y, pos.z + z);

				//Make sure within bounds of grid
				if (n.x >= 0 && n.x < sp.gBounds.x && n.y >= 0 && n.y < sp.gBounds.y && n.z >= 0 && n.z < sp.gBounds.z) {
					float3 diff = (particle.pos - (make_float3(n) * sp.radius)) * hInv;
					float3 gw = gradWeight(diff);
					int gIndex = getGridIndex(n);

					velocityGrad += mat3::outerProduct(cells[gIndex].velocityStar, gw);
				}
			}
		}
	}

	return velocityGrad;
}



__global__ void transferData(Particle* particles, Cell* cells) {
	int index = threadIdx.x + (blockIdx.x * blockDim.x);
	if (index >= sp.numParticles) return;

	//Compute current volume and stress factor
	//printf("%lf %lf %lf\n", calcStress(particles[index].fe, particles[index].fp)[0], calcStress(particles[index].fe, particles[index].fp)[4], calcStress(particles[index].fe, particles[index].fp)[8]);
	mat3 volumeStress = -particles[index].volume * calcStress(particles[index].fe, particles[index].fp);
	//particles[index].D = mat3(0.0f);
	float hInv = 1.0f / sp.radius;
	int3 pos = make_int3(particles[index].pos * hInv);
	for (int z = -2; z < 3; z++) {
		for (int y = -2; y < 3; y++) {
			for (int x = -2; x < 3; x++) {
				int3 n = make_int3(pos.x + x, pos.y + y, pos.z + z);
				//Make sure within bounds of grid
				if (n.x >= 0 && n.x < sp.gBounds.x && n.y >= 0 && n.y < sp.gBounds.y && n.z >= 0 && n.z < sp.gBounds.z) {
					float3 diff = (particles[index].pos - (make_float3(n) * sp.radius)) * hInv;//APIC 中的 xp-xi
					float wip = weight(fabs(diff));
					mat3 xixp = mat3((make_float3(n) * sp.radius) - particles[index].pos, make_float3(0.0f), make_float3(0.0f));//APIC
																																//	particles[index].D += wip*xixp*mat3::transpose(xixp);//APIC
				}
			}
		}
	}
	//Particle adds mass and velocity to all cells within 2h of itself
	for (int z = -2; z < 3; z++) {
		for (int y = -2; y < 3; y++) {
			for (int x = -2; x < 3; x++) {
				int3 n = make_int3(pos.x + x, pos.y + y, pos.z + z);
				//Make sure within bounds of grid
				if (n.x >= 0 && n.x < sp.gBounds.x && n.y >= 0 && n.y < sp.gBounds.y && n.z >= 0 && n.z < sp.gBounds.z) {
					float3 diff = (particles[index].pos - (make_float3(n) * sp.radius)) * hInv;//APIC 中的 xp-xi
					float3 gw = gradWeight(diff);
					int gIndex = getGridIndex(n);
					float wip = weight(fabs(diff));

					//float3 new_v = wip*particles[index].mass*(particles[index].velocity + particles[index].B*mat3::inverse(particles[index].D)*((make_float3(n) * sp.radius) - particles[index].pos));
					//atomicAdd(&(cells[gIndex].velocity.x), new_v.x);//APIC
					//atomicAdd(&(cells[gIndex].velocity.y), new_v.y);//APIC
					//atomicAdd(&(cells[gIndex].velocity.z), new_v.z);//APIC
					//Accumulate force at this cell node (the contribution from the current particle)
					float3 f = volumeStress * gw;

					float mi = particles[index].mass * wip;
					atomicAdd(&(cells[gIndex].mass), mi);
					atomicAdd(&(cells[gIndex].velocity.x), particles[index].velocity.x * mi);//FLIP&PIC
					atomicAdd(&(cells[gIndex].velocity.y), particles[index].velocity.y * mi);//FLIP&PIC
					atomicAdd(&(cells[gIndex].velocity.z), particles[index].velocity.z * mi);//FLIP&PIC
					atomicAdd(&(cells[gIndex].force.x), f.x);
					atomicAdd(&(cells[gIndex].force.y), f.y);
					atomicAdd(&(cells[gIndex].force.z), f.z);
				}
			}
		}
	}
}

__global__ void initialTransfer(Particle* particles, Cell* cells) {
	int index = threadIdx.x + (blockIdx.x * blockDim.x);
	if (index >= sp.numParticles) return;

	float hInv = 1.0f / sp.radius;
	int3 pos = make_int3(particles[index].pos * hInv);

	//Particle adds mass and velocity to all cells within 2h of itself
	for (int z = -2; z < 3; z++) {
		for (int y = -2; y < 3; y++) {
			for (int x = -2; x < 3; x++) {
				int3 n = make_int3(pos.x + x, pos.y + y, pos.z + z);

				//Make sure within bounds of grid
				if (n.x >= 0 && n.x < sp.gBounds.x && n.y >= 0 && n.y < sp.gBounds.y && n.z >= 0 && n.z < sp.gBounds.z) {
					float3 diff = (particles[index].pos - (make_float3(n) * sp.radius)) * hInv;
					int gIndex = getGridIndex(n);

					float mi = particles[index].mass * weight(fabs(diff));
					atomicAdd(&(cells[gIndex].mass), mi);
				}
			}
		}
	}
}

__global__ void computeVolumes(Particle* particles, Cell* cells) {
	int index = threadIdx.x + (blockIdx.x * blockDim.x);
	if (index >= sp.numParticles) return;

	float hInv = 1.0f / sp.radius;
	int3 pos = make_int3(particles[index].pos * hInv);
	float pDensity = 0.0f;
	float invCellVolume = hInv * hInv * hInv;

	//Cells up to 2h away contribute towards volume
	for (int z = -2; z < 3; z++) {
		for (int y = -2; y < 3; y++) {
			for (int x = -2; x < 3; x++) {
				int3 n = make_int3(pos.x + x, pos.y + y, pos.z + z);

				//Make sure within bounds of grid
				if (n.x >= 0 && n.x < sp.gBounds.x && n.y >= 0 && n.y < sp.gBounds.y && n.z >= 0 && n.z < sp.gBounds.z) {
					float3 diff = (particles[index].pos - (make_float3(n) * sp.radius)) * hInv;
					int gIndex = getGridIndex(n);
					pDensity += cells[gIndex].mass * invCellVolume * weight(fabs(diff));
				}
			}
		}
	}

	particles[index].volume = particles[index].mass / pDensity;
}



__global__ void updateVelocities(Cell* cells) {
	int index = threadIdx.x + (blockIdx.x * blockDim.x);
	if (index >= sp.gridSize) return;

	//Simple forward euler step for predicted velocities (10)
	if (cells[index].mass > 0.0f) {
		float invMass = 1.0f / cells[index].mass;
		cells[index].force += cells[index].mass * sp.gravity;
		cells[index].velocity *= invMass;
		cells[index].velocityStar = cells[index].velocity + sp.deltaT * invMass * cells[index].force;



	}
}

__global__ void bodyCollisions(Cell* cells) {
	int index = threadIdx.x + (blockIdx.x * blockDim.x);
	if (index >= sp.gridSize) return;

	float z = index / (int(sp.gBounds.y) * int(sp.gBounds.x));
	float y = index % (int(sp.gBounds.y) * int(sp.gBounds.x)) / (int(sp.gBounds.x));
	float x = index % int(sp.gBounds.x);

	//TODO: check that we are consistent with cell edge vs. center... e.g. is position index * radius or
	//index * radius + .5 * radius?

	float3 pos = make_float3(x, y, z) * sp.radius;
	applyCellBoundaryCollisions(cells, index, pos + sp.deltaT * cells[index].velocityStar, cells[index].velocityStar);
	//applyBoundaryCollisions(pos + sp.deltaT * cells[index].velocityStar, cells[index].velocityStar);

}

__global__ void updateDeformationGradient(Particle* particles, Cell* cells) {
	int index = threadIdx.x + (blockIdx.x * blockDim.x);
	if (index >= sp.numParticles) return;

	//The deformation gradient is updated as Fp_n+1 = (I + dt * velocityGradp_n+1)
	//See part 7 of MPM paper for more details on elastic and plastic part updates
	//First, temporarily define fe_next = (I + dt * vp)fe
	//fp_next = fp

	mat3 velocityGrad = getVelocityGradient(particles[index], cells);
	mat3 newFe = (mat3(1.0f) + (sp.deltaT * velocityGrad)) * particles[index].fe;
	mat3 newF = newFe * particles[index].fp;
	particles[index].fe = newFe;//new
								/*
								//Take the SVD of the elastic part
								mat3 U, S, V, Sinv;
								computeSVD(newFe, U, S, V);

								S = mat3(clamp(S[0], 1 - sp.compression, 1 + sp.stretch), 0.0f, 0.0f,
								0.0f, clamp(S[4], 1 - sp.compression, 1 + sp.stretch), 0.0f,
								0.0f, 0.0f, clamp(S[8], 1 - sp.compression, 1 + sp.stretch));
								Sinv = mat3(1.0f / S[0], 0.0f, 0.0f,
								0.0f, 1.0f / S[4], 0.0f,
								0.0f, 0.0f, 1.0f / S[8]);

								particles[index].fe = mat3::multiplyADBt(U, S, V);
								particles[index].fp = mat3::multiplyADBt(V, Sinv, U) * newF;
								*/
}

__global__ void updateParticleVelocities(Particle* particles, Cell* cells) {
	int index = threadIdx.x + (blockIdx.x * blockDim.x);
	if (index >= sp.numParticles) return;

	//For each particle, compute the PIC and FLIP components
	float hInv = 1.0f / sp.radius;
	int3 pos = make_int3(particles[index].pos * hInv);

	float3 velocityPic = make_float3(0.0f);
	float3 velocityFlip = particles[index].velocity;
	//particles[index].B = mat3(0.0f);
	//Particles influenced by all cells within 2h of itself
	for (int z = -2; z < 3; z++) {
		for (int y = -2; y < 3; y++) {
			for (int x = -2; x < 3; x++) {
				int3 n = make_int3(pos.x + x, pos.y + y, pos.z + z);

				//Make sure within bounds of grid
				if (n.x >= 0 && n.x < sp.gBounds.x && n.y >= 0 && n.y < sp.gBounds.y && n.z >= 0 && n.z < sp.gBounds.z) {
					float3 diff = (particles[index].pos - (make_float3(n) * sp.radius)) * hInv;
					int gIndex = getGridIndex(n);

					float w = weight(fabs(diff));

					mat3 xixp = mat3((make_float3(n) * sp.radius) - particles[index].pos, make_float3(0.0f), make_float3(0.0f));//APIC
																																//		particles[index].B += w*mat3(cells[gIndex].velocityStar, make_float3(0.0f), make_float3(0.0f))*mat3::transpose(xixp);//APIC

					velocityPic += cells[gIndex].velocityStar * w;//also used for APIC
					velocityFlip += (cells[gIndex].velocityStar - cells[gIndex].velocity) * w;
				}
			}
		}
	}

	//linear interp betwen PIC and FLIP
	particles[index].velocity = ((1 - sp.alpha) * velocityPic) + (sp.alpha * velocityFlip);

	//particles[index].velocity = velocityPic;//also used for APIC
}

__global__ void particleBodyCollisions(Particle* particles) {
	int index = threadIdx.x + (blockIdx.x * blockDim.x);
	if (index >= sp.numParticles) return;
	applyBoundaryCollisions(particles[index].pos + sp.deltaT * particles[index].velocity, particles[index].velocity);

}


__global__ void updatePositions(Particle* particles) {
	int index = threadIdx.x + (blockIdx.x * blockDim.x);
	if (index >= sp.numParticles) return;
	float3 tmpPos = particles[index].pos;
	particles[index].pos += sp.deltaT * particles[index].velocity;

}


__device__ mat3 MCC_project(mat3 S, float qc, float *dq, Particle *particle) {
	mat3 e = mat3(log(S[0]), 0.0f, 0.0f,
		0.0f, log(S[4]), 0.0f,
		0.0f, 0.0f, log(S[8]));//自然对数e=lnS,因为S为对角矩阵，可以仅对对角元素进行操作，减少计算量
	float trE = (e[0] + e[4] + e[8]);

	mat3 newE = e - trE / 3.0f * mat3(1.0f);//newE=e-tr(e)/3*I;
	if ((abs(newE[0])<0.00000000001&&abs(newE[4])<0.00000000001&&abs(newE[8])<0.00000000001) || trE > 0.0f) {//在屈服面之外
		(*dq) = 0;
		particle->color = make_float3(0.8f, 0.3f, 0.3f);
		return mat3(1.0f);
	}

	float newE_F = sqrt(mat3::innerProduct(newE, newE));//newE的F范数


	float p = trE*(3.0 * sp.lambda + 2.0 * sp.mu) / 3.0f / 2.0f / sp.mu / sqrt(1.5f);
	float q = sqrt(1.5f)*newE_F * 2.0 * sp.mu / 2.0f / sp.mu;
	//float M = 0.210128f;
	//float M = 0.5f;
	float M = 0.5f;
	float dy = q*q + M*M*p*(p + qc / 2.0f / sp.mu / sqrt(1.5f));


	if (dy <= 0.00000000001) {//在屈服面之内
		(*dq) = 0;
		particle->color = make_float3(0.3f, 0.8f, 0.3f);
		//printf("%lf  \n", dy);
		return S;
	}
	mat3 H;
	if (dy / (newE_F*newE_F)>1.0f) {
		H = e - newE;
		//particle->color = make_float3(0.0f, 0.0f, 0.0f);
	}
	else {
		H = e - (1.0f - sqrt(1.0f - dy / (newE_F*newE_F))) * newE;
		particle->color = make_float3(0.3f, 0.3f, 0.8f);//蓝
	}
	(*dq) = dy;
	if (abs(p)<qc / 2.0f / sp.mu / sqrt(1.5f) / 2.0) {
		(*dq) = -dy;
		//particle->color = make_float3(0.8f, 0.3f, 0.3f);//红
	}


	return mat3(exp(H[0]), 0.0f, 0.0f,
		0.0f, exp(H[4]), 0.0f,
		0.0f, 0.0f, exp(H[8]));//因为H为对角矩阵，可以仅对对角元素进行操作，减少计算量
}

__device__ mat3 project(mat3 S, float a, float *dq) {
	mat3 e = mat3(log(S[0]), 0.0f, 0.0f,
		0.0f, log(S[4]), 0.0f,
		0.0f, 0.0f, log(S[8]));//自然对数e=lnS,因为S为对角矩阵，可以仅对对角元素进行操作，减少计算量
	float trE = (e[0] + e[4] + e[8]);


	mat3 newE = e - trE / 3.0f * mat3(1.0f);//newE=e-tr(e)/3*I;
	if ((abs(newE[0])<0.000001f&&abs(newE[4])<0.000001f&&abs(newE[8])<0.000001f) || trE > 0.0f) {//在屈服面之外
		(*dq) = 0;
		return mat3(1.0f);
	}

	float newE_F = sqrt(mat3::innerProduct(newE, newE));//newE的F范数

	float c = 0.0f;
	float dy = newE_F + (3 * sp.lambda + 2 * sp.mu) / 2.0f / sp.mu*trE*a - c;//dy=||newE||F+(d*lambda+2*mu)/(2*mu)*tr(E)*a;


	if (dy <= 0) {//在屈服面之内
		(*dq) = 0;
		return S;
	}

	mat3 H = e - dy*(newE / newE_F);
	(*dq) = dy;
	float maxV = 1000.01f;//cone max
	return mat3(::clamp(exp(H[0]), 0.0f, maxV), 0.0f, 0.0f,
		0.0f, ::clamp(exp(H[4]), 0.0f, maxV), 0.0f,
		0.0f, 0.0f, ::clamp(exp(H[8]), 0.0f, maxV));//因为H为对角矩阵，可以仅对对角元素进行操作，减少计算量
}

__global__ void MCC_plasticity_hardening(Particle* particles) {
	int index = threadIdx.x + (blockIdx.x * blockDim.x);
	if (index >= sp.numParticles) return;
	mat3 Fe = particles[index].fe;
	mat3 U, S, V;
	computeSVD(Fe, U, S, V);
	float dq = 0.0f;
	mat3 T = MCC_project(S, particles[index].qc, &dq, &particles[index]);
	particles[index].fe = U*T*mat3::transpose(V);
	particles[index].fp = V*mat3::inverse(T)*S*mat3::transpose(V)*particles[index].fp;//后面的fp是new fp，这里可能存在问题
	particles[index].q = particles[index].q + dq;
	//particles[index].qc = particles[index].qc*exp(dq/0.1f);
}
__global__ void plasticity_hardening(Particle* particles) {
	int index = threadIdx.x + (blockIdx.x * blockDim.x);
	if (index >= sp.numParticles) return;
	mat3 Fe = particles[index].fe;
	mat3 U, S, V;
	computeSVD(Fe, U, S, V);
	float dq = 0.0f;
	mat3 T = project(S, particles[index].a, &dq);
	particles[index].fe = U*T*mat3::transpose(V);
	particles[index].fp = V*mat3::inverse(T)*S*mat3::transpose(V)*particles[index].fp;//后面的fp是new fp，这里可能存在问题
	particles[index].q = particles[index].q + dq;
	float fi_F = 40.0f;
	//float fi_F = sp.h0 + (sp.h1*particles[index].q - sp.h3)*exp(-sp.h2*particles[index].q);
	particles[index].a = 0.81649658f*(2 * sinf(fi_F / 180.0f*3.1415926f)) / (3 - sinf(fi_F / 180.0f*3.1415926f));//particles[index].a=sqrt(2/3)*(2 * sinf(fi_F)) / (3 - sinf(fi_F));//前面的一串数字是根号2/3
}

__global__ void getPos(float* positionsPtr, Particle* particles) {
	int index = threadIdx.x + (blockIdx.x * blockDim.x);
	if (index >= sp.numParticles) return;

	positionsPtr[6 * index + 0] = particles[index].pos.x;
	positionsPtr[6 * index + 1] = particles[index].pos.y;
	positionsPtr[6 * index + 2] = particles[index].pos.z;
	positionsPtr[6 * index + 3] = particles[index].color.x;
	positionsPtr[6 * index + 4] = particles[index].color.y;
	positionsPtr[6 * index + 5] = particles[index].color.z;
}

void update(Particle* particles, Cell* cells, int gridSize) {
	//Clear cell data
	cudaCheck(cudaMemset(cells, 0, gridSize * sizeof(Cell)));

	if (!firstTime) transferData << <particleDims, blockSize >> >(particles, cells);
	if (firstTime) {
		initialTransfer << <particleDims, blockSize >> >(particles, cells);
		computeVolumes << <particleDims, blockSize >> >(particles, cells);
		firstTime = false;
	}
	updateVelocities << <gridDims, blockSize >> >(cells);
	bodyCollisions << <gridDims, blockSize >> >(cells);
	updateDeformationGradient << <particleDims, blockSize >> >(particles, cells);
	updateParticleVelocities << <particleDims, blockSize >> >(particles, cells);
	particleBodyCollisions << <particleDims, blockSize >> >(particles);
	updatePositions << <particleDims, blockSize >> >(particles);
	MCC_plasticity_hardening << <particleDims, blockSize >> >(particles);
	//plasticity_hardening << <particleDims, blockSize >> >(particles);
}

__host__ void setTerrainTex(string terrainHeightPath, string terrainNormalpath) {
	int width, height, nrChannels;
	unsigned char *texData = stbi_load(terrainHeightPath.data(), &width, &height, &nrChannels, 0);
	if (!texData) {
		cout << "open Terrain Height Texture file fail" << endl;
		return;
	}
	unsigned char *dst = (unsigned char *)malloc(width *height * 4);
	for (int i = 0; i < height; ++i) {
		for (int j = 0; j < width; ++j) {
			dst[(i*width + j) * 4] = texData[(i*width + j) * 3];
			dst[(i*width + j) * 4 + 1] = texData[(i*width + j) * 3 + 1];
			dst[(i*width + j) * 4 + 2] = texData[(i*width + j) * 3 + 2];
		}
	}
	unsigned int size = width * height * 4 * sizeof(unsigned char);
	cudaChannelFormatDesc channelDesc =
		cudaCreateChannelDesc(8, 8, 8, 8, cudaChannelFormatKindUnsigned);
	cudaArray *cuArray;
	cudaCheck(cudaMallocArray(&cuArray,
		&channelDesc,
		width,
		height));

	cudaCheck(cudaMemcpyToArray(cuArray,
		0,
		0,
		dst,
		size,
		cudaMemcpyHostToDevice));

	TerrainHeightTex.addressMode[0] = cudaAddressModeWrap;
	TerrainHeightTex.addressMode[1] = cudaAddressModeWrap;
	//TerrainHeightTex.filterMode = cudaFilterModeLinear;
	TerrainHeightTex.normalized = true;    // access with normalized texture coordinates
	cudaCheck(cudaBindTextureToArray(TerrainHeightTex, cuArray, channelDesc));
	unsigned char *tex2Data = stbi_load(terrainNormalpath.data(), &width, &height, &nrChannels, 0);
	if (!tex2Data) {
		cout << "open Terrain Normal Texture file fail" << endl;
		return;
	}
	size = width * height * 4 * sizeof(unsigned char);
	unsigned char *dst2 = (unsigned char *)malloc(width *height * 4);
	for (int i = 0; i < height; ++i) {
		for (int j = 0; j < width; ++j) {
			dst2[(i*width + j) * 4] = tex2Data[(i*width + j) * 3];
			dst2[(i*width + j) * 4 + 1] = tex2Data[(i*width + j) * 3 + 1];
			dst2[(i*width + j) * 4 + 2] = tex2Data[(i*width + j) * 3 + 2];
		}
	}
	cudaArray *cuArray2;
	cudaCheck(cudaMallocArray(&cuArray2,
		&channelDesc,
		width,
		height));
	cudaCheck(cudaMemcpyToArray(cuArray2,
		0,
		0,
		dst2,
		size,
		cudaMemcpyHostToDevice));
	TerrainNormalTex.addressMode[0] = cudaAddressModeWrap;
	TerrainNormalTex.addressMode[1] = cudaAddressModeWrap;
	//TerrainNormalTex.filterMode = cudaFilterModeLinear;
	TerrainNormalTex.normalized = true;    // access with normalized texture coordinates
	cudaCheck(cudaBindTextureToArray(TerrainNormalTex, cuArray2, channelDesc));

	stbi_image_free(texData);
	stbi_image_free(tex2Data);
}

__host__ void setTerrainTex16(string terrainHeightPath16, string terrainNormalpath) {
	int width, height, nrChannels;
	stbi_us  *texData = stbi_load_16(terrainHeightPath16.data(), &width, &height, &nrChannels, 0);
	if (!texData) {
		cout << "open Terrain Height Texture file fail" << endl;
		return;
	}
	float *dst = (float *)malloc(width *height * 4 * sizeof(float));
	//cout << "nrChannels " << nrChannels << endl;
	for (int i = 0; i < height; ++i) {
		for (int j = 0; j < width; ++j) {
			dst[(i*width + j) * 4] = texData[(i*width + j)] / 65535.0;
			dst[(i*width + j) * 4 + 1] = texData[(i*width + j)] / 65535.0;
			dst[(i*width + j) * 4 + 2] = texData[(i*width + j)] / 65535.0;
		}
	}
	unsigned int size = width * height * 4 * sizeof(float);
	cudaChannelFormatDesc channelDesc =
		cudaCreateChannelDesc(32, 32, 32, 32, cudaChannelFormatKindFloat);
	cudaArray *cuArray;
	cudaCheck(cudaMallocArray(&cuArray,
		&channelDesc,
		width,
		height));

	cudaCheck(cudaMemcpyToArray(cuArray,
		0,
		0,
		dst,
		size,
		cudaMemcpyHostToDevice));

	TerrainHeightTex.addressMode[0] = cudaAddressModeWrap;
	TerrainHeightTex.addressMode[1] = cudaAddressModeWrap;
	TerrainHeightTex.filterMode = cudaFilterModeLinear;
	//cout << "TerrainHeightTex.filterMode" << TerrainHeightTex.filterMode << endl;
	TerrainHeightTex.normalized = true;    // access with normalized texture coordinates
	cudaCheck(cudaBindTextureToArray(TerrainHeightTex, cuArray, channelDesc));
	unsigned char *tex2Data = stbi_load(terrainNormalpath.data(), &width, &height, &nrChannels, 0);
	if (!tex2Data) {
		cout << "open Terrain Normal Texture file fail" << endl;
		return;
	}
	size = width * height * 4 * sizeof(float);
	float *dst2 = (float *)malloc(width *height * 4 * sizeof(float));
	for (int i = 0; i < height; ++i) {
		for (int j = 0; j < width; ++j) {
			dst2[(i*width + j) * 4] = tex2Data[(i*width + j) * 3] / 255.0f;
			dst2[(i*width + j) * 4 + 1] = tex2Data[(i*width + j) * 3 + 1] / 255.0f;
			dst2[(i*width + j) * 4 + 2] = tex2Data[(i*width + j) * 3 + 2] / 255.0f;
		}
	}
	cudaArray *cuArray2;
	cudaCheck(cudaMallocArray(&cuArray2,
		&channelDesc,
		width,
		height));
	cudaCheck(cudaMemcpyToArray(cuArray2,
		0,
		0,
		dst2,
		size,
		cudaMemcpyHostToDevice));
	TerrainNormalTex.addressMode[0] = cudaAddressModeWrap;
	TerrainNormalTex.addressMode[1] = cudaAddressModeWrap;
	//TerrainNormalTex.filterMode = cudaFilterModeLinear;
	TerrainNormalTex.normalized = false;    // access with normalized texture coordinates
	cudaCheck(cudaBindTextureToArray(TerrainNormalTex, cuArray2, channelDesc));

	stbi_image_free(texData);
	stbi_image_free(tex2Data);
}


__host__ void setParams(solverParams *params) {
	particleDims = int(ceil(params->numParticles / blockSize + 0.5f));
	gridDims = int(ceil(params->gridSize / blockSize + 0.5f));
	cudaCheck(cudaMemcpyToSymbol(sp, params, sizeof(solverParams)));
	glEnable(GL_DEPTH_TEST);
}

__host__ void updateRigid(solverParams *params) {
	float densityBall = 2000.0f;
	if (firstTime) {//第一次初始化场景包围盒刚体等
		rp3d::SphereShape *sphereShape = new rp3d::SphereShape(rp3d::decimal(params->rig->r));
		world = new rp3d::DynamicsWorld(gravity);
		rp3d::Transform transform = rp3d::Transform::identity();
		rp3d::decimal mass = rp3d::decimal(densityBall * 4 * 3.1415926* params->rig->r *params->rig->r*params->rig->r / 3.0 );
		body = world->createRigidBody(rp3d::Transform(rp3d::Vector3(params->rig->center.x, params->rig->center.y, params->rig->center.z), rp3d::Quaternion::identity()));
		body->setType(rp3d::BodyType::DYNAMIC);
		body->addCollisionShape(sphereShape, transform, mass);
		rp3d::Material & material = body->getMaterial();
		// Change the bounciness of the body
		material.setBounciness(rp3d::decimal(0.0));
		// Change the friction coefficient of the body
		material.setFrictionCoefficient(rp3d::decimal(0.2));

		//包围盒
		rp3d::BoxShape* bottomFloor = new rp3d::BoxShape(rp3d::Vector3(1.275f, 0.5, 1.275f));
		rp3d::BoxShape* topFloor = new rp3d::BoxShape(rp3d::Vector3(1.275f, 0.5, 1.275f));
		rp3d::BoxShape* leftFloor = new rp3d::BoxShape(rp3d::Vector3(0.5, 1.275f, 1.275f));
		rp3d::BoxShape* rightFloor = new rp3d::BoxShape(rp3d::Vector3(0.5, 1.275f, 1.275f));
		rp3d::BoxShape* frontFloor = new rp3d::BoxShape(rp3d::Vector3(1.275f, 1.275f, 0.5f));
		rp3d::BoxShape* afterFloor = new rp3d::BoxShape(rp3d::Vector3(1.275f, 1.275f, 0.5f));
		boundingBox = world->createRigidBody(rp3d::Transform::identity());
		boundingBox->setType(rp3d::BodyType::STATIC);
		boundingBox->addCollisionShape(bottomFloor, rp3d::Transform(rp3d::Vector3(1.275f, -0.5f, 1.275f), rp3d::Quaternion::identity()), 100000.0f);
		boundingBox->addCollisionShape(topFloor, rp3d::Transform(rp3d::Vector3(1.275f, 2.55f+0.5f, 1.275f), rp3d::Quaternion::identity()), 100000.0f);
		boundingBox->addCollisionShape(leftFloor, rp3d::Transform(rp3d::Vector3(-0.5f, 1.275f, 1.275f), rp3d::Quaternion::identity()), 100000.0f);
		boundingBox->addCollisionShape(rightFloor, rp3d::Transform(rp3d::Vector3(2.55f + 0.5f, -1.275f, 1.275f), rp3d::Quaternion::identity()), 100000.0f);
		boundingBox->addCollisionShape(frontFloor, rp3d::Transform(rp3d::Vector3(1.275f, 1.275f, -0.5f), rp3d::Quaternion::identity()), 100000.0f);
		boundingBox->addCollisionShape(afterFloor, rp3d::Transform(rp3d::Vector3(1.275f, 1.275f, 2.55+0.5f), rp3d::Quaternion::identity()), 100000.0f);
	}
	if (!firstTime)
		cudaCheck(cudaMemcpyFromSymbol(params->rig, r, sizeof(rigid)));
	/*用于测试
	//params->rig->vBall = make_float3(0.0f, -9.8 * params->deltaT, 0.0f) + params->rig->force / (densityBall * 4 * 3.1415926* params->rig->r *params->rig->r*params->rig->r / 3) + params->rig->vBall;
	//params->rig->center += params->rig->vBall *params->deltaT;
	*/
	printf("%lf %lf %lf %lf\n", params->rig->vBall.y, params->rig->center.x, params->rig->center.y, params->rig->center.z);

	params->rig->force /= params->deltaT;//这里的force其实是动量，因此除以dt求得力
	params->rig->torque /= params->deltaT;//这里的torque其实是动量，因此除以dt求得力
	body->applyForceToCenterOfMass(rp3d::Vector3(params->rig->force.x, params->rig->force.y, params->rig->force.z));
	body->applyTorque(rp3d::Vector3(params->rig->torque.x, params->rig->torque.y, params->rig->torque.z));
	world->update(params->deltaT);
	params->rig->center = make_float3(body->getWorldPoint(rp3d::Vector3(0.0, 0.0, 0.)).x, body->getWorldPoint(rp3d::Vector3(0.0, 0.0, 0.)).y, body->getWorldPoint(rp3d::Vector3(0.0, 0.0, 0.)).z);
	params->rig->vBall = make_float3(body->getLinearVelocity().x, body->getLinearVelocity().y, body->getLinearVelocity().z);
	params->rig->wBall = make_float3(body->getAngularVelocity().x, body->getAngularVelocity().y, body->getAngularVelocity().z);

	params->rig->force = make_float3(0.0f);
	params->rig->torque = make_float3(0.0f);
	cudaCheck(cudaMemcpyToSymbol(r, params->rig, sizeof(rigid)));

}

__host__ void getPositions(float* positionsPtr, Particle* particles) {
	getPos << <particleDims, blockSize >> >(positionsPtr, particles);
}

#endif