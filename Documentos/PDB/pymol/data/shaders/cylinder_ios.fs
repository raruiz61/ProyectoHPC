precision mediump float;

uniform sampler2D depthTexture;
uniform lowp sampler2D sphereDepthMap;
//uniform sampler2D ambientOcclusionTexture;
uniform mat3 inverseModelViewProjMatrix;
//uniform mediump float ambientOcclusionTexturePatchWidth;

varying mediump vec2 impostorSpaceCoordinate;
varying mediump vec3 normalAlongCenterAxis;
varying mediump float depthOffsetAlongCenterAxis;
varying mediump float normalizedDepthOffsetAlongCenterAxis;
varying mediump float normalizedDisplacementAtEndCaps;
varying mediump float normalizedRadialDisplacementAtEndCaps;
varying mediump vec2 rotationFactor;
varying mediump vec3 normalizedViewCoordinate;
//varying mediump vec2 ambientOcclusionTextureBase;
varying mediump float depthAdjustmentForOrthographicProjection;
varying mediump float normalizedDistanceAlongZAxis;
varying mediump vec3 cylinder_color1;
varying mediump vec3 cylinder_color2;


struct g_LightModelParameters {
       mediump vec4 ambient;
};
uniform g_LightModelParameters g_LightModel;


const mediump float oneThird = 1.0 / 3.0;

mediump float depthFromEncodedColor(mediump vec4 encodedColor)
{
    return oneThird * (encodedColor.r + encodedColor.g + encodedColor.b);
}

mediump vec2 textureCoordinateForCylinderSurfacePosition(mediump vec3 cylinderSurfacePosition)
{
    vec2 halfAbsoluteXY = abs(cylinderSurfacePosition.xy / 2.0);
    
    if (cylinderSurfacePosition.x >= 0.0)
    {
        return vec2(cylinderSurfacePosition.y / (4.0 * (halfAbsoluteXY.x + halfAbsoluteXY.y)) - 0.5, cylinderSurfacePosition.z);
    }
    else
    {
        return vec2(-cylinderSurfacePosition.y / (4.0 * (halfAbsoluteXY.x + halfAbsoluteXY.y)) + 0.5, cylinderSurfacePosition.z);
    }
}

void main()
{
    float adjustmentFromCenterAxis = sqrt(1.0 - impostorSpaceCoordinate.s * impostorSpaceCoordinate.s);
    float displacementFromCurvature = normalizedDisplacementAtEndCaps * adjustmentFromCenterAxis;
    float depthOffset = depthOffsetAlongCenterAxis * adjustmentFromCenterAxis * depthAdjustmentForOrthographicProjection;

    vec3 normal = vec3(normalizedRadialDisplacementAtEndCaps * rotationFactor.x * adjustmentFromCenterAxis + impostorSpaceCoordinate.s * rotationFactor.y,
                       -(normalizedRadialDisplacementAtEndCaps * rotationFactor.y * adjustmentFromCenterAxis + impostorSpaceCoordinate.s * rotationFactor.x),
                       normalizedDepthOffsetAlongCenterAxis * adjustmentFromCenterAxis);
    normal = normalize(normal);
    if ( (impostorSpaceCoordinate.t <= (-1.0 + displacementFromCurvature)) || (impostorSpaceCoordinate.t >= (1.0 + displacementFromCurvature))) {
        gl_FragColor = vec4(0.0);
    }
    else {
        float currentDepthValue = normalizedViewCoordinate.z - depthOffset + 0.0025;
        float previousDepthValue = depthFromEncodedColor(texture2D(depthTexture, normalizedViewCoordinate.xy));

        if ( (currentDepthValue - 0.002) > (previousDepthValue) ) {
            gl_FragColor = vec4(0.0);
        }
        else {
	    vec2 lookupCoordinate ;
            lookupCoordinate.x = normalizedRadialDisplacementAtEndCaps * rotationFactor.y * adjustmentFromCenterAxis + impostorSpaceCoordinate.s * rotationFactor.x;
            lookupCoordinate.y = (normalizedRadialDisplacementAtEndCaps * rotationFactor.x * adjustmentFromCenterAxis - impostorSpaceCoordinate.s * rotationFactor.y);
	    lowp vec4 precalculatedDepthAndLighting = texture2D(sphereDepthMap, (lookupCoordinate.xy / 2.0) + .5);
	    lowp vec3 lightingIntensity = g_LightModel.ambient.xyz + (0.5 * precalculatedDepthAndLighting.g);

	    vec3 cylinderColor = mix(cylinder_color1,cylinder_color2, step(0.0, normalizedDistanceAlongZAxis));
            lowp vec3 finalCylinderColor = cylinderColor * lightingIntensity;
	    finalCylinderColor = finalCylinderColor + (precalculatedDepthAndLighting.b);

            gl_FragColor = vec4(finalCylinderColor, 1.0);
        }            
   }        
}
