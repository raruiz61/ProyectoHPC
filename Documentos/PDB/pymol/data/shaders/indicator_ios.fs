precision mediump float;

uniform sampler2D textureMap;
uniform vec2 textureLookup ;
uniform vec2 textureScale ;

uniform mat3 matL;
uniform float gamma;

void main()
{
  vec4 fColor = texture2D(textureMap, textureLookup + gl_PointCoord * textureScale);
  // if not stereo, matL and gamma are identity
  fColor.rgb = matL * pow(fColor.rgb, vec3(gamma));
  gl_FragColor = fColor;
}

