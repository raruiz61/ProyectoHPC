#ifdef PURE_OPENGL_ES_2
precision highp float;
#endif

uniform sampler2D bgTextureMap;
varying vec2 bgTextureLookup ;

uniform vec3 fogSolidColor;
uniform float fogIsSolidColor;

uniform float isStretched;
uniform float isCentered;
uniform float isCenteredOrRepeated;
uniform float isTiled;
uniform vec2 tileSize;
uniform vec2 tiledSize;
uniform vec2 viewImageSize;
uniform vec2 pixelSize;
uniform vec2 halfPixel;

void main()
{
#ifdef bg_image_mode_stretched
  vec2 bgLookup = pixelSize * floor(bgTextureLookup / pixelSize);
  vec4 bgColor = texture2D(bgTextureMap, bgLookup);
  gl_FragColor = vec4(bgColor.rgb*bgColor.a + (1.-bgColor.a) * fogSolidColor.rgb, 1.);
#endif
#ifdef bg_image_mode_centered
  vec2 centeredBgLookup = pixelSize * floor((bgTextureLookup - halfPixel) / pixelSize);
  centeredBgLookup = halfPixel + centeredBgLookup / viewImageSize;
  vec2 bgLookup = pixelSize * floor(centeredBgLookup / pixelSize);
  float isCenteredButOutsideOrSolid = step(centeredBgLookup.x, 0.) + step(centeredBgLookup.y, 0.) +
                                      step(1., centeredBgLookup.x) + step(1., centeredBgLookup.y);
  isCenteredButOutsideOrSolid = step(.5, isCentered * isCenteredButOutsideOrSolid + (1.-isCentered)*fogIsSolidColor);
  vec4 bgColor = texture2D(bgTextureMap, bgLookup);
  bgColor = vec4(bgColor.rgb*bgColor.a + (1.-bgColor.a) * fogSolidColor.rgb, 1.); // background should not write any alpha
  gl_FragColor = (1.- isCenteredButOutsideOrSolid) * bgColor +
  	         isCenteredButOutsideOrSolid * vec4(fogSolidColor, 1.);
#endif
#ifdef bg_image_mode_tiled
  vec2 bgLookupT = mod(bgTextureLookup, tiledSize) / tiledSize;
  vec2 bgLookup = tileSize * floor(bgLookupT / tileSize);
  vec4 bgColor = texture2D(bgTextureMap, bgLookup);
  gl_FragColor = vec4(bgColor.rgb*bgColor.a + (1.-bgColor.a) * fogSolidColor.rgb, 1.); // background should not write any alpha
#endif
#ifdef bg_image_mode_solid
  gl_FragColor = vec4(fogSolidColor.rgb, 1.);
#endif
}

