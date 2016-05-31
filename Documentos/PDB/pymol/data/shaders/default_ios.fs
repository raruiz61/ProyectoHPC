precision highp float;

uniform float lighting_enabled;
uniform float lighting_not_enabled;
uniform sampler2D depthTexture;
uniform lowp float doNotUseDepthMap;
uniform lowp sampler2D lightingDepthMap;

#include ANAGLYPH_HEADER

varying vec4 COLOR;
varying vec3 NORMAL;
varying vec3 normalizedViewCoordinate;

uniform float fog_enabled;
varying float fog;

uniform sampler2D bgTextureMap;
uniform vec3 fogSolidColor;
uniform float fogIsSolidColor;

const float oneThird = 1.0 / 3.0;

float depthFromEncodedColor(vec3 encodedColor)
{
    return (encodedColor.r + encodedColor.g + encodedColor.b) * oneThird;
}

const lowp vec3 stepValues = vec3(2.0, 1.0, 0.0);

struct g_LightModelParameters {
       vec4 ambient;
};
uniform g_LightModelParameters g_LightModel;

uniform float isStretched;
uniform float isCentered;
uniform float isCenteredOrRepeated;
uniform float isTiled;
uniform vec2 tileSize;
uniform vec2 tiledSize;
uniform vec2 viewImageSize;
uniform vec2 pixelSize;
uniform vec2 halfPixel;

#define bgTextureLookup normalizedViewCoordinate.xy

#include ComputeFogColor

void main()
{
  lowp float alphaComponent;
  lowp vec4 precalculatedDepthAndLighting = texture2D(lightingDepthMap, NORMAL.xy);
  alphaComponent = step(0.5, precalculatedDepthAndLighting.a);
  vec3 encodedColor = texture2D(depthTexture, bgTextureLookup).rgb;
  float previousDepthValue = max(depthFromEncodedColor(encodedColor), doNotUseDepthMap);
  alphaComponent = step((normalizedViewCoordinate.z - 0.002), previousDepthValue);
  lowp vec3 lightingIntensity = g_LightModel.ambient.xyz + (0.5 * precalculatedDepthAndLighting.g);
  lowp vec3 finalColor = COLOR.rgb*lightingIntensity*lighting_enabled + lighting_not_enabled * COLOR.rgb;
  finalColor = finalColor + lighting_enabled * (precalculatedDepthAndLighting.b);

  float cfog = mix(1.0, clamp(fog, 0.0, 1.0), fog_enabled);

  vec4 fogColor = ComputeFogColor();
  vec4 fColor = vec4(mix(fogColor.xyz, alphaComponent*finalColor, cfog), alphaComponent*COLOR.a);

#include ANAGLYPH_BODY
}
