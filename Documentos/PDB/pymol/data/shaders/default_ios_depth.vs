/* ios recent */
precision mediump float;

attribute mediump vec4 a_Vertex;
attribute mediump vec3 a_Normal;
attribute mediump vec4 a_Color;
attribute float a_Accessibility; /* This is for ambient occlusion, 1.0 by default */

varying vec4 COLOR ;

uniform mat4 g_ModelViewMatrix;
uniform mat4 g_ProjectionMatrix;

void main()
{
  COLOR = a_Color;
  vec4 position = g_ProjectionMatrix * g_ModelViewMatrix * a_Vertex;
  position.z = position.z * step(1., a_Color.a);
  gl_Position = position;
}
