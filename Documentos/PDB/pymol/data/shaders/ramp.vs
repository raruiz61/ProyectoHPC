attribute vec4 a_Vertex;
attribute vec4 a_Color;
attribute vec3 a_Normal;

varying vec4 COLOR ;
varying vec3 NORMAL ;

uniform vec2 t2PixelSize;

uniform vec3 offsetPt;

#ifdef PURE_OPENGL_ES_2
uniform mat4 g_ModelViewMatrix;
uniform mat4 g_ProjectionMatrix;
#endif

void main()
{
  COLOR = a_Color;
  NORMAL = a_Normal;
  vec4 vertex = a_Vertex + vec4(offsetPt,0.);
#ifdef PURE_OPENGL_ES_2
  gl_Position = g_ProjectionMatrix * /* g_ModelViewMatrix */ vertex;
#else
  gl_Position = gl_ModelViewProjectionMatrix * vertex;
#endif
}
