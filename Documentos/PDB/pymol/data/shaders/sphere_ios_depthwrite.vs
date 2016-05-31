attribute vec4 a_vertex_radius;
attribute vec2 a_rightUpFlags;
attribute vec4 a_Color;

uniform mat4 g_ModelViewMatrix;
uniform mat4 g_ProjectionMatrix;
uniform mediump mat4 g_orthographicMatrix;
uniform float sphere_size_scale;
uniform float depth_range_size;

//varying mediump vec4 sphereColor;

void main()
{
    mediump vec2 insetCoordinate = a_rightUpFlags * 0.91017; // Octagon
    vec4 position = vec4(a_vertex_radius.xyz, 1.);
    mediump float sphereRadius = a_vertex_radius.w * sphere_size_scale;
    vec4 transformedPosition;
    transformedPosition = g_ModelViewMatrix * position;
    transformedPosition.xy = transformedPosition.xy + insetCoordinate * vec2(sphereRadius);
    transformedPosition = transformedPosition * g_orthographicMatrix;

    gl_Position = g_ProjectionMatrix * transformedPosition;
}
