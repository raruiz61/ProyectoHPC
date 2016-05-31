
attribute mediump vec4 a_vertex_radius;
attribute mediump vec2 a_rightUpFlags;
attribute vec4 a_Color;
//attribute mediump vec2 ambientOcclusionTextureOffset;

varying mediump vec2 impostorSpaceCoordinate;
varying mediump vec3 normalizedViewCoordinate;
varying mediump vec2 depthLookupCoordinate;
//varying mediump vec2 ambientOcclusionTextureBase;
varying mediump float adjustedSphereRadius;
varying mediump vec4 sphereColor;

uniform mat4 g_ModelViewMatrix;
uniform mat4 g_ProjectionMatrix;
//uniform mediump mat3 modelViewProjMatrix;
uniform mediump mat4 g_orthographicMatrix;
uniform float sphere_size_scale;
uniform float depth_range_size;

void main()
{
    vec4 position = vec4(a_vertex_radius.xyz, 1.);
    mediump float sphereRadius = a_vertex_radius.w * sphere_size_scale;
    sphereColor = a_Color;

//    ambientOcclusionTextureBase = ambientOcclusionTextureOffset;
    vec4 transformedPosition = g_ModelViewMatrix * position;
    impostorSpaceCoordinate = a_rightUpFlags.xy;
    depthLookupCoordinate = (a_rightUpFlags.xy / 2.0) + 0.5;

    transformedPosition.xy = transformedPosition.xy + a_rightUpFlags.xy * vec2(sphereRadius);
    transformedPosition = transformedPosition * g_orthographicMatrix;
    
    float depthAdjustmentForOrthographicProjection = (vec4(0.0, 0.0, 0.5, 0.0) * g_orthographicMatrix).z;
    adjustedSphereRadius = sphereRadius * depthAdjustmentForOrthographicProjection/depth_range_size;

    transformedPosition = g_ProjectionMatrix * transformedPosition;
    transformedPosition.z = transformedPosition.z - adjustedSphereRadius * transformedPosition.w;
    gl_Position = transformedPosition;
    normalizedViewCoordinate = (transformedPosition.xyz/transformedPosition.w) / 2.0 + 0.5;
}
