precision mediump float;

uniform sampler2D depthTexture;
uniform lowp sampler2D sphereDepthMap;

varying mediump vec2 impostorSpaceCoordinate;
varying mediump vec2 depthLookupCoordinate;
varying mediump vec3 normalizedViewCoordinate;
varying mediump float adjustedSphereRadius;
varying mediump vec4 sphereColor;

struct g_LightModelParameters {
       mediump vec4 ambient;
};
uniform g_LightModelParameters g_LightModel;

const mediump float oneThird = 1.0 / 3.0;

mediump float depthFromEncodedColor(mediump vec3 encodedColor)
{
    return (encodedColor.r + encodedColor.g + encodedColor.b) * oneThird;
}

const lowp vec3 stepValues = vec3(2.0, 1.0, 0.0);

void main()
{

    lowp vec4 precalculatedDepthAndLighting = texture2D(sphereDepthMap, depthLookupCoordinate);
    lowp float alphaComponent = 1.0;
  
    alphaComponent = step(0.5, precalculatedDepthAndLighting.a);

    float currentDepthValue = normalizedViewCoordinate.z + adjustedSphereRadius - adjustedSphereRadius * precalculatedDepthAndLighting.r;
    vec3 encodedColor = texture2D(depthTexture, normalizedViewCoordinate.xy).rgb;

    float previousDepthValue = depthFromEncodedColor(encodedColor);

        // Check to see that this fragment is the frontmost one for this area
    alphaComponent = alphaComponent * step((currentDepthValue - 0.002), previousDepthValue);
    
    // Ambient lighting            
    lowp vec3 lightingIntensity = g_LightModel.ambient.xyz + (0.5 * precalculatedDepthAndLighting.g);
    lowp vec3 finalSphereColor = sphereColor.xyz * lightingIntensity;

    // Specular lighting
    finalSphereColor = finalSphereColor + (precalculatedDepthAndLighting.b);
    
    gl_FragColor = vec4(finalSphereColor * alphaComponent, alphaComponent);
}
