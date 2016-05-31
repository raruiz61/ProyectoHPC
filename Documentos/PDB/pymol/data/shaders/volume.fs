uniform sampler3D volumeTex;
uniform sampler2D colorTex2D;
uniform sampler1D colorTex1D;
uniform float volumeScale;
uniform float volumeBias;
uniform float sliceDist;
uniform mat4 TexMatrix;
uniform vec3 eyeposM;
uniform vec3 vDirM;
uniform sampler3D carvemask;
uniform bool carvemaskFlag;

varying vec3 vertexM;

uniform float fog_enabled;
uniform sampler2D bgTextureMap;
uniform vec3 fogSolidColor;
uniform float fogIsSolidColor;
varying float fog;

uniform float g_Fog_end;
uniform float g_Fog_scale;
varying vec2 bgTextureLookup;

uniform float isStretched;
uniform float isCentered;
uniform float isCenteredOrRepeated;
uniform float isTiled;
uniform vec2 tileSize;
uniform vec2 tiledSize;
uniform vec2 viewImageSize;
uniform vec2 pixelSize;
uniform vec2 halfPixel;

#include ANAGLYPH_HEADER

#include ComputeFogColor

bool iscarvemasked(vec3 t) {
  return carvemaskFlag && texture3D(carvemask, t).r > 0.5;
}

void main()
{
#ifdef volume_mode
#ifdef ortho
  vec3 tex0 = gl_TexCoord[0].xyz;
  vec3 tex1 = gl_TexCoord[1].xyz;
#else // ortho

  // normalized eye to vertex vector in model space
  vec3 eyeToVert = normalize(vertexM - eyeposM);

  // vector from vertex to the back slab plane
  vec3 slabOffset = eyeToVert * (sliceDist / dot(vDirM, eyeToVert));

  // back slab and front slab points in model space
  vec3 sB = vertexM - slabOffset;
  vec3 sF = vertexM + slabOffset;

  // texture coordinates
  vec3 tex0 = vec3(TexMatrix * vec4(sB, 1.));
  vec3 tex1 = vec3(TexMatrix * vec4(sF, 1.));

#endif // ortho

  if (iscarvemasked(tex0))
    discard;

  // texture lookup of map values
  vec2 v = vec2(
    texture3D(volumeTex, tex0).r,
    texture3D(volumeTex, tex1).r
  );

  // color lookup in (preintegrated) table
  v = v * volumeScale + volumeBias;
  vec4 color = texture2D(colorTex2D, v);

#else // volume_mode

  if (iscarvemasked(gl_TexCoord[0].xyz))
    discard;

  float v = texture3D(volumeTex, gl_TexCoord[0].xyz).r;
  v = v * volumeScale + volumeBias;
  if (v < 0. || v > 1.) discard;
  vec4 color = texture1D(colorTex1D, v);
#endif // volume_mode

  if (color.a == 0.0) discard;
  float cfog = mix(1.0, clamp(fog, 0.0, 1.0), fog_enabled);
  vec4 fogColor = ComputeFogColor();
  vec4 fColor = vec4(mix(vec3(fogColor), color.rgb, cfog), color.a);

#include ANAGLYPH_BODY
}

