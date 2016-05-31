
varying vec3 NORMAL ;
#ifdef use_geometry_shaders
varying vec4 COLOR ;
#else
varying vec4 COLORIn ;
#define COLOR COLORIn
#endif
varying float fog;
varying vec2 bgTextureLookup;
#ifdef use_geometry_shaders
varying float lineEdge;
varying float aaCutoff;
#endif
uniform float fog_enabled;
uniform sampler2D bgTextureMap;
uniform vec3 fogSolidColor;
uniform float fogIsSolidColor;

uniform bool lighting_enabled;
uniform bool two_sided_lighting_enabled;
uniform int light_count;
uniform vec4 interior_color;
uniform float interior_color_threshold;
uniform float shininess;
uniform float shininess_0;
uniform bool use_interior_color_threshold;
uniform int spec_count;
uniform float spec_value;
uniform float spec_value_0;

#include ANAGLYPH_HEADER

uniform float isStretched;
uniform float isCentered;
uniform float isCenteredOrRepeated;
uniform float isTiled;
uniform vec2 tileSize;
uniform vec2 tiledSize;
uniform vec2 viewImageSize;
uniform vec2 pixelSize;
uniform vec2 halfPixel;
uniform float antialiasedLines;

#include ComputeFogColor

//#include ComputeColorForLight

float mysmoothstep(float edge0, float edge1, float x){
  float rets = step(edge0, edge1) * step(edge1, edge0);
  return rets * step(edge0, x) + (1.-rets) * smoothstep(edge0, edge1, x);
}

void main()
{
  vec4 final_color = vec4(0.);
#ifdef use_geometry_shaders
  float alpha = antialiasedLines * mysmoothstep(0., aaCutoff, 1. - abs(lineEdge)) * COLOR.a + (1.-antialiasedLines) * COLOR.a;
#else
  float alpha = COLOR.a;
#endif
  final_color = COLOR;

  float cfog = mix(1.0, clamp(fog, 0.0, 1.0), fog_enabled);
  vec4 fogColor = ComputeFogColor();
  vec4 fColor = vec4(mix(vec3(fogColor), final_color.rgb, cfog), alpha);

#include ANAGLYPH_BODY
}

