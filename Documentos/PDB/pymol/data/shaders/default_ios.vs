precision highp float;

attribute vec4 a_Vertex;
attribute vec3 a_Normal;
attribute vec4 a_Color;
attribute float a_Accessibility; 

varying vec4 COLOR;
varying vec3 NORMAL;
varying vec3 normalizedViewCoordinate;
varying float fog;

uniform mat3 g_NormalMatrix;
uniform mat4 g_ModelViewMatrix;
uniform mat4 g_ProjectionMatrix;
uniform float g_pointSize;

uniform float g_Fog_end;
uniform float g_Fog_start;
uniform float g_Fog_scale;

void main()
{
  NORMAL = (normalize(g_NormalMatrix * a_Normal) / 2.0) + 0.5;
  COLOR = a_Color;
  gl_PointSize = g_pointSize;
  vec3 eye_pos = vec3(g_ModelViewMatrix * a_Vertex);
  gl_Position = g_ProjectionMatrix * g_ModelViewMatrix * a_Vertex;
  normalizedViewCoordinate = (gl_Position.xyz/gl_Position.w) / 2.0 + 0.5;
  fog = (g_Fog_end - abs(eye_pos.z)) * g_Fog_scale;
}
