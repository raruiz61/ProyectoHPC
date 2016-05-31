#version 120
#extension GL_EXT_gpu_shader4 : enable

uniform vec4 SMAA_RT_METRICS;

uniform sampler2D colorTex;
uniform float isOutput;
varying vec2 texcoordAttr;
varying vec4 offset[3];

#define texture texture2D
#define textureLod texture2DLod
#define textureLodOffset texture2DLodOffset

#include smaa_gen_fs

void main()
{
	vec2 edges = SMAAColorEdgeDetectionPS(texcoordAttr, offset, colorTex);
//	vec2 edges = SMAALumaEdgeDetectionPS(texcoordAttr, offset, colorTex);
	gl_FragColor.rg = edges;
}
