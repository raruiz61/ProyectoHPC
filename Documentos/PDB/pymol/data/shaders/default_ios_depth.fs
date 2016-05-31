/* ios recent */
precision mediump float;


const lowp vec3 stepValues = vec3(2.0, 1.0, 0.0);

void main()
{
  gl_FragColor = vec4(clamp(vec3(3. * gl_FragCoord.z) - stepValues, 0., 1.), 1.);
}

