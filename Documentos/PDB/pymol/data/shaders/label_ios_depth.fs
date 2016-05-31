/* ios recent */
precision mediump float;

uniform sampler2D textureMap;
varying vec2 textureLookup ;

const lowp vec3 stepValues = vec3(2.0, 1.0, 0.0);

void main()
{
  vec4 fontColor = texture2D(textureMap, textureLookup);
  float returnDepth = step(0.002, fontColor.a);
  float returnWhite = 1. - returnDepth;
  gl_FragColor = returnWhite * vec4(1.) + returnDepth * vec4(clamp(vec3(3. * gl_FragCoord.z) - stepValues, 0., 1.), 1.);
}

