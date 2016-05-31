precision highp float;

attribute highp vec4 a_Vertex;
attribute mediump vec4 a_Color;

varying mediump vec4 COLOR;

uniform mat4 g_ModelViewMatrix;
uniform mat4 g_ProjectionMatrix;
uniform float g_pointSize;

void main()
{
  COLOR = a_Color;
  gl_PointSize = g_pointSize;
  gl_Position = g_ProjectionMatrix * g_ModelViewMatrix * a_Vertex;
}
