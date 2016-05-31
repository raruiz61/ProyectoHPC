attribute vec4 attr_origin;
attribute vec4 attr_axis;

uniform mat4 g_ModelViewMatrix;
uniform mat4 g_ProjectionMatrix;
uniform mediump mat4 g_orthographicMatrix;

uniform float ortho;
uniform float uni_radius;

varying mediump vec2 impostorSpaceCoordinate;
varying mediump float depthOffsetAlongCenterAxis;
varying mediump float normalizedDisplacementAtEndCaps;
varying mediump float normalizedDepth;
varying mediump float depthAdjustmentForOrthographicProjection;

uniform float depth_range_size;

void main()
{
    vec2 rotationFactor;
    vec4 transformedDirection, transformedPosition, transformedOtherPosition;
    vec3 viewDisplacementForVertex, displacementDirectionAtEndCap;
    float displacementAtEndCaps, lengthOfCylinder, lengthOfCylinderInView;
    float cylinderRadius;
    vec2 inputImpostorSpaceCoordinate;
    float has_uni_radius = step(0.001, uni_radius);
    cylinderRadius = has_uni_radius * uni_radius * attr_origin.w;
    cylinderRadius = cylinderRadius + (1.0 - has_uni_radius) * attr_origin.w;

    float packed_flags = attr_axis.w;
    vec4 flags = mod(vec4(packed_flags/262144.0, packed_flags/4096.0, 
                          packed_flags/64.0, packed_flags), 64.0);
    inputImpostorSpaceCoordinate.x = (2.0 * step(.5, flags.y)) - 1.0;
    inputImpostorSpaceCoordinate.y = (2.0 * step(.5, flags.z)) - 1.0;

    depthAdjustmentForOrthographicProjection = (vec4(0.0, 0.0, 0.5, 0.0) * g_orthographicMatrix).z;
    transformedPosition = g_ModelViewMatrix * vec4(attr_origin.xyz, 1.0);
    transformedOtherPosition = g_ModelViewMatrix * vec4((attr_origin.xyz + attr_axis.xyz), 1.0);
    transformedPosition = vec4(transformedPosition.xyz/transformedPosition.w, 1.);
    transformedOtherPosition = vec4(transformedOtherPosition.xyz/transformedOtherPosition.w, 1.);
    transformedDirection = transformedPosition - transformedOtherPosition;

    lengthOfCylinder = length(transformedDirection.xyz);
    lengthOfCylinderInView = length(transformedDirection.xy);
    rotationFactor = transformedDirection.xy / lengthOfCylinderInView;

    displacementAtEndCaps = cylinderRadius * (transformedDirection.z) / lengthOfCylinder;
    normalizedDisplacementAtEndCaps = displacementAtEndCaps / lengthOfCylinderInView;

    depthOffsetAlongCenterAxis = cylinderRadius * lengthOfCylinder * inversesqrt(lengthOfCylinder * lengthOfCylinder - (transformedDirection.z) * (transformedDirection.z));
    depthOffsetAlongCenterAxis = depthOffsetAlongCenterAxis/depth_range_size;
//    depthOffsetAlongCenterAxis = clamp(depthOffsetAlongCenterAxis, 0.0, cylinderRadius * 2.0) / depth_range_size;

    displacementDirectionAtEndCap.xy = displacementAtEndCaps * rotationFactor;
    displacementDirectionAtEndCap.z = transformedDirection.z * displacementAtEndCaps / lengthOfCylinder;
    
    transformedDirection.xy = normalize(transformedDirection.xy);
    
    if ((displacementAtEndCaps * inputImpostorSpaceCoordinate.t) > 0.0) {
        viewDisplacementForVertex.x = inputImpostorSpaceCoordinate.x * transformedDirection.y * cylinderRadius + displacementDirectionAtEndCap.x;
        viewDisplacementForVertex.y = -inputImpostorSpaceCoordinate.x * transformedDirection.x * cylinderRadius + displacementDirectionAtEndCap.y;
        viewDisplacementForVertex.z = displacementDirectionAtEndCap.z;
        impostorSpaceCoordinate = vec2(inputImpostorSpaceCoordinate.s, inputImpostorSpaceCoordinate.t + 1.0 * normalizedDisplacementAtEndCaps);
    } else {
        viewDisplacementForVertex.x = inputImpostorSpaceCoordinate.x * transformedDirection.y * cylinderRadius;
        viewDisplacementForVertex.y = -inputImpostorSpaceCoordinate.x * transformedDirection.x * cylinderRadius;    
        viewDisplacementForVertex.z = 0.0;
        impostorSpaceCoordinate = vec2(inputImpostorSpaceCoordinate.s, inputImpostorSpaceCoordinate.t);
    }

    transformedPosition.xyz = (transformedPosition.xyz - viewDisplacementForVertex);
    transformedPosition *= g_orthographicMatrix;
    gl_Position = g_ProjectionMatrix * transformedPosition;
    normalizedDepth = (gl_Position.z/gl_Position.w) / 2.0 + 0.5;
}



