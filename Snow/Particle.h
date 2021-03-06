#ifndef PARTICLE_H
#define PARTICLE_H

#include "Common.h"
#include "matrix.h"

struct Particle {
	float3 pos;
	float3 velocity;
	float3 color;
	float mass;
	float volume;
	mat3 fe;
	mat3 fp;

	//mat3 B;	//affine matrix  used for APIC
	//mat3 D;	//helper matrix  used for APIC

	float qc;//for mcc hardening and softening
	float a;//for hardening 
	float q;

	Particle(float3 pos, float3 velocity, float mass, float3 c = make_float3(0.3f, 0.8f, 0.3f)) :
		pos(pos), velocity(velocity), mass(mass), volume(0),
		fe(mat3(1.0f)), fp(mat3(1.0f)),/* B(mat3(0.0f)), D(mat3(0.0f)),*/ a(0.0f), qc(50000.0f), color(c)//200000.0f
	{}
};//20000000.0f

#endif