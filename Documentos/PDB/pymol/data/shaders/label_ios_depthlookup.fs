/* ios recent */
precision mediump float;

uniform sampler2D textureMap;
uniform sampler2D depthTexture;
varying vec2 textureLookup ;
varying mediump vec3 normalizedViewCoordinate;

const mediump float oneThird = 1.0 / 3.0;

mediump float depthFromEncodedColor(mediump vec3 encodedColor)
{
  return (encodedColor.r + encodedColor.g + encodedColor.b) * oneThird;
}

void main()
{
  vec4 textureColor = texture2D(textureMap, textureLookup);

  vec3 encodedColor = texture2D(depthTexture, normalizedViewCoordinate.xy).rgb;
  float previousDepthValue = depthFromEncodedColor(encodedColor);
  float alphaComponent = step((normalizedViewCoordinate.z - 0.002), previousDepthValue);
  gl_FragColor = vec4(textureColor.xyz * alphaComponent, textureColor.a * alphaComponent);
}
