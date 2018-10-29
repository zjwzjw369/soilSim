#ifndef SIMULATION_CU
#define SIMULATION_CU

#include "Common.h"
#include "Parameters.h"
#include "Particle.h"
#include "Cell.h"
#include "matrix.h"
#include "decomposition.h"

#define cudaCheck(x) { cudaError_t err = x; if (err != cudaSuccess) { printf("Cuda error: %f in %s at %s:%f\n", err, #x, __FILE__, __LINE__); assert(0); } }

using namespace std;

static dim3 particleDims;
static dim3 gridDims;
static const int blockSize = 128;
bool firstTime = true;

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

__device__ void applyBoundaryCollisions(const float3& position, float3& velocity) {
	float vn;
	float3 vt;
	float3 normal;

	bool collision;

	//handle the 3 boundary dimensions separately
	for (int i = 0; i < 3; i++) {
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

		if (collision) {
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
	particles[index].D = mat3(0.0f);
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
					particles[index].D += wip*xixp*mat3::transpose(xixp);//APIC
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

					float3 new_v = wip*particles[index].mass*(particles[index].velocity + particles[index].B*mat3::inverse(particles[index].D)*((make_float3(n) * sp.radius) - particles[index].pos));
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
	applyBoundaryCollisions(pos + sp.deltaT * cells[index].velocityStar, cells[index].velocityStar);
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
	particles[index].B = mat3(0.0f);
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
					particles[index].B += w*mat3(cells[gIndex].velocityStar, make_float3(0.0f), make_float3(0.0f))*mat3::transpose(xixp);//APIC

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

	particles[index].pos += sp.deltaT * particles[index].velocity;
}


__device__ mat3 MCC_project(mat3 S, float qc, float *dq,Particle *particle) {
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
	float M = 0.210128f;
	float dy = q*q + M*M*p*(p + qc / 2.0f / sp.mu / sqrt(1.5f));
	//printf("%lf\n", ((3 * sp.lambda + 2 * sp.mu) / 2.0f / sp.mu)*0.210128);0.210128f
	//float M = 1.2f;
	//float dy = newE_F*newE_F + 2.0f/3.0f / 2.0f / sp.mu*((3.0f * sp.lambda + 2.0f * sp.mu) / 2.0f / sp.mu)*M*M  * trE/3.0f*((3.0f * sp.lambda + 2.0f * sp.mu)*trE/3.0f +qc );//dy=||newE||F+(d*lambda+2*mu)/(2*mu)*tr(E)*a;
	//if(trE>-0.001f&&trE<0.0f)
	//printf("%lf  trE\n", trE);
	if (dy <= 0.00000000001) {//在屈服面之内
		(*dq) = 0;
		particle->color = make_float3(0.3f, 0.8f, 0.3f);
		//printf("%lf  \n", dy);
		return S;
	}
	mat3 H;
	if (dy / (newE_F*newE_F)>1.0f) {
		H = e - newE;
		particle->color = make_float3(0.0f, 0.0f, 0.0f);
	}
	else {
		H = e - (1.0f - sqrt(1.0f - dy / (newE_F*newE_F))) * newE;
		particle->color = make_float3(0.3f, 0.3f, 0.8f);
	}
	(*dq) = dy;
	if (trE / 3<qc / 2) {
		//	(*dq) = -dy;
	}
	//if(dy>0.001)
	//printf("%lf\n", trE*(3 * sp.lambda + 2 * sp.mu) / 3.0);

	return mat3(exp(H[0]), 0.0f, 0.0f,
		0.0f, exp(H[4]), 0.0f,
		0.0f, 0.0f, exp(H[8]));//因为H为对角矩阵，可以仅对对角元素进行操作，减少计算量
}

__device__ mat3 project(mat3 S, float a, float *dq) {
	mat3 e = mat3(log(S[0]), 0.0f, 0.0f,
		0.0f, log(S[4]), 0.0f,
		0.0f, 0.0f, log(S[8]));//自然对数e=lnS,因为S为对角矩阵，可以仅对对角元素进行操作，减少计算量
	float trE = (e[0] + e[4] + e[8]);

	/*
	float maxE = -0.3f;//cone max
	if (e[0]<maxE||e[4]<maxE||e[8]<maxE) {
		//printf("%lf  %lf  %lf\n", e[0], e[4], e[8]);

		return mat3(clamp(S[0], maxE, 1000.0f), 0.0f, 0.0f,
		0.0f, clamp(S[4], maxE, 1000.0f), 0.0f,
		0.0f, 0.0f, clamp(S[8], maxE, 1000.0f));//因为H为对角矩阵，可以仅对对角元素进行操作，减少计算量
	}
	*/
	mat3 newE = e - trE / 3.0f * mat3(1.0f);//newE=e-tr(e)/3*I;
	if ((abs(newE[0])<0.000001f&&abs(newE[4])<0.000001f&&abs(newE[8])<0.000001f) || trE > 0.0f) {//在屈服面之外
		(*dq) = 0;
		return mat3(1.0f);
	}
	//printf("%lf  %lf  %lf      %lf\n", newE[0], newE[4], newE[8], mat3::innerProduct(newE, newE));

	float newE_F = sqrt(mat3::innerProduct(newE, newE));//newE的F范数

	//float dy = newE_F + (3 * sp.lambda + 2 * sp.mu) / 2.0f / sp.mu*trE*a;//dy=||newE||F+(d*lambda+2*mu)/(2*mu)*tr(E)*a;
																		 //  printf("%lf\n", newE_F);

	float c = 0.0f;
	float dy = newE_F + (3 * sp.lambda + 2 * sp.mu) / 2.0f / sp.mu*trE*a- c;//dy=||newE||F+(d*lambda+2*mu)/(2*mu)*tr(E)*a;


	if (dy <= 0) {//在屈服面之内
		(*dq) = 0;
		return S;
	}


	mat3 H = e - dy*(newE / newE_F);
	(*dq) = dy;
	float maxV = 1000.01f;//cone max
	return mat3(clamp(exp(H[0]),0.0f, maxV), 0.0f, 0.0f,
		0.0f, clamp(exp(H[4]), 0.0f, maxV), 0.0f,
		0.0f, 0.0f, clamp(exp(H[8]), 0.0f, maxV));//因为H为对角矩阵，可以仅对对角元素进行操作，减少计算量
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
	
	
	
	particles[index].qc = particles[index].qc*exp(dq);
	//float fi_F = 20.0f;
	//float fi_F = sp.h0 + (sp.h1*particles[index].q - sp.h3)*exp(-sp.h2*particles[index].q);
	//particles[index].a = 0.81649658f*(2 * sinf(fi_F / 180.0f*3.1415926f)) / (3 - sinf(fi_F / 180.0f*3.1415926f));//particles[index].a=sqrt(2/3)*(2 * sinf(fi_F)) / (3 - sinf(fi_F));//前面的一串数字是根号2/3
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
	//float fi_F = 20.0f;
	float fi_F = sp.h0 + (sp.h1*particles[index].q - sp.h3)*exp(-sp.h2*particles[index].q);
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

	//Rasterize particle data to grid (1)
	if (!firstTime) transferData << <particleDims, blockSize >> >(particles, cells);

	//First timestep only: compute particle volumes (2)
	if (firstTime) {
		initialTransfer << <particleDims, blockSize >> >(particles, cells);
		computeVolumes << <particleDims, blockSize >> >(particles, cells);
		firstTime = false;
	}

	//Update grid velocities (4)
	updateVelocities << <gridDims, blockSize >> >(cells);

	//Grid collisions (5)
	bodyCollisions << <gridDims, blockSize >> >(cells);

	//Update deformation gradient (7)
	updateDeformationGradient << <particleDims, blockSize >> >(particles, cells);

	//Update particle velocities (8)
	updateParticleVelocities << <particleDims, blockSize >> >(particles, cells);

	//Particle collisions (9)
	particleBodyCollisions << <particleDims, blockSize >> >(particles);

	//Update particle positions (10)
	updatePositions << <particleDims, blockSize >> >(particles);
	MCC_plasticity_hardening << <particleDims, blockSize >> >(particles);
}

__host__ void setParams(solverParams *params) {
	particleDims = int(ceil(params->numParticles / blockSize + 0.5f));
	gridDims = int(ceil(params->gridSize / blockSize + 0.5f));
	cudaCheck(cudaMemcpyToSymbol(sp, params, sizeof(solverParams)));
}

__host__ void getPositions(float* positionsPtr, Particle* particles) {
	getPos << <particleDims, blockSize >> >(positionsPtr, particles);
}

#endif