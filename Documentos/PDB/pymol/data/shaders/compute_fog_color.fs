vec4 ComputeFogColor(){
#ifdef bg_image_mode_stretched
  vec2 bgLookup = pixelSize * floor(bgTextureLookup / pixelSize);
  vec4 bgColor = texture2D(bgTextureMap, bgLookup);
  return vec4(bgColor.rgb*bgColor.a + (1.-bgColor.a) * fogSolidColor.rgb, 1.);
#endif
#ifdef bg_image_mode_centered
  vec2 centeredBgLookup = pixelSize * floor((bgTextureLookup - halfPixel) / pixelSize);
  centeredBgLookup = halfPixel + centeredBgLookup / viewImageSize;
  vec2 bgLookup = pixelSize * floor(centeredBgLookup / pixelSize);
  float isCenteredButOutsideOrSolid = step(centeredBgLookup.x, 0.) + step(centeredBgLookup.y, 0.) +
                                      step(1., centeredBgLookup.x) + step(1., centeredBgLookup.y);
  isCenteredButOutsideOrSolid = step(.5, isCentered * isCenteredButOutsideOrSolid + (1.-isCentered)*fogIsSolidColor);
  vec4 bgColor = texture2D(bgTextureMap, bgLookup);
  bgColor = vec4(bgColor.rgb*bgColor.a + (1.-bgColor.a) * fogSolidColor.rgb, bgColor.a);
  return (1.- isCenteredButOutsideOrSolid) * bgColor +
         isCenteredButOutsideOrSolid * vec4(fogSolidColor, 1.);
#endif
#ifdef bg_image_mode_tiled
  vec2 bgLookup = mod(bgTextureLookup, tiledSize) / tiledSize;
  bgLookup = tileSize * floor(bgLookup / tileSize);
  vec4 bgColor = texture2D(bgTextureMap, bgLookup);
  return vec4(bgColor.rgb*bgColor.a + (1.-bgColor.a) * fogSolidColor.rgb, bgColor.a);
#endif
#ifdef bg_image_mode_solid
  return vec4(fogSolidColor.rgb, 1.);	
#endif
}
