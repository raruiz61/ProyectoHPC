#version 120
#extension GL_EXT_gpu_shader4 : enable

uniform vec4 SMAA_RT_METRICS;
uniform sampler2D colorTex;
uniform sampler2D blendTex;
uniform float isOutput;

varying vec2 texcoordAttr;
varying vec4 offset;

#define texture texture2D
#define textureLod texture2DLod
#define textureLodOffset texture2DLodOffset

#include smaa_gen_fs

void main()
{
	gl_FragColor = SMAANeighborhoodBlendingPS(texcoordAttr, offset, colorTex, blendTex);
}
