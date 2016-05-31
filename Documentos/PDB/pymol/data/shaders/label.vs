#ifdef PURE_OPENGL_ES_2
precision highp float;
#endif

attribute vec4 attr_targetpos;
attribute vec4 attr_worldpos;
attribute vec3 attr_screenoffset;
attribute vec2 attr_texcoords;
attribute vec3 attr_screenworldoffset;
attribute vec4 attr_pickcolor;
attribute float attr_relative_mode;

varying vec2 textureLookup ;
varying vec3 normalizedViewCoordinate;
varying vec4 pickcolor ;

#ifdef PURE_OPENGL_ES_2
uniform mat3 g_NormalMatrix;
uniform mat4 g_ModelViewMatrix;
uniform mat4 g_ProjectionMatrix;
uniform mat4 g_orthographicMatrix;
uniform float g_Fog_end;
uniform float g_Fog_start;
uniform float g_Fog_scale;
#endif

uniform vec2 screenSize;
uniform vec2 pixelSize;

uniform float aspectRatioAdjustment;
uniform float contentScaleFactor;

uniform float screenOriginVertexScale;

uniform float scaleByVertexScale;
uniform float labelTextureSize;
varying float fog;
uniform float fog_enabled; // actually bool
uniform float textureSize;

uniform float front;
uniform float clipRange;

vec4 normalizeVec4(vec4 point){
     return vec4(point.xyz/point.w, 1.);
}

float convertNormalZToScreenZ(float normalz){
   float a_centerN = (normalz + 1.) / 2.;
   float ptInPreProjectionZ = -(front + clipRange * a_centerN);
   vec4 ptInPreProjection = vec4(0., 0., ptInPreProjectionZ, 1.);
#ifdef PURE_OPENGL_ES_2
   vec4 projVect = g_ProjectionMatrix * ptInPreProjection;
#else
   vec4 projVect = gl_ProjectionMatrix * ptInPreProjection;
#endif
   return projVect.z / projVect.w;
}

void main()
{
  float drawConnector, isProjected, isScreenCoord, isPixelCoord, zTarget;

  drawConnector = step(1., mod(attr_relative_mode, 2.));
  isScreenCoord = step(2., mod(attr_relative_mode, 4.));
  isPixelCoord = step(4., mod(attr_relative_mode, 8.));
  zTarget = step(8., mod(attr_relative_mode, 16.));
  isProjected = step(isPixelCoord + isScreenCoord, 0.5);

  float screenVertexScale = scaleByVertexScale * screenOriginVertexScale * labelTextureSize + (1. - scaleByVertexScale);
#ifdef PURE_OPENGL_ES_2
  vec3 viewVector = vec3(vec4(0.,0.,-1.,0.) * g_ModelViewMatrix);
  vec4 transformedPosition = g_ProjectionMatrix * g_ModelViewMatrix * attr_worldpos;
  vec4 targetPosition = normalizeVec4(g_ProjectionMatrix * g_ModelViewMatrix * attr_targetpos);
#else
  vec3 viewVector = vec3(gl_ModelViewMatrixTranspose * vec4(0.,0., -1.,0.));
  vec4 transformedPosition = gl_ModelViewProjectionMatrix * attr_worldpos;
  vec4 targetPosition = normalizeVec4(gl_ModelViewProjectionMatrix * attr_targetpos);
#endif
  transformedPosition.xyz = transformedPosition.xyz/transformedPosition.w;
  vec4 a_center = (attr_worldpos + attr_screenworldoffset.z * vec4(viewVector, 0.));
#ifdef PURE_OPENGL_ES_2
  vec4 transformedPositionZ = g_ProjectionMatrix * g_ModelViewMatrix * a_center;
#else
  vec4 transformedPositionZ = gl_ModelViewProjectionMatrix * a_center;
#endif
  transformedPositionZ.xyz = transformedPositionZ.xyz/transformedPositionZ.w;
  transformedPositionZ.w = 1.;
  vec2 pixOffset = ((2. * attr_worldpos.xy / screenSize) - 1.);
  transformedPosition = isProjected * transformedPosition + isScreenCoord * attr_worldpos + isPixelCoord * vec4(pixOffset.x, pixOffset.y, -0.5, 0.);
  transformedPosition.xy = transformedPosition.xy + attr_screenworldoffset.xy/(screenSize*screenOriginVertexScale);

  transformedPosition.z = (1.-zTarget) * ((isProjected * transformedPositionZ.z) + (1.-isProjected) * convertNormalZToScreenZ(attr_worldpos.z)) + zTarget * targetPosition.z;

  transformedPosition.x = transformedPosition.x + aspectRatioAdjustment * attr_screenoffset.x * contentScaleFactor * 2./(screenSize.x*screenVertexScale);
  transformedPosition.y = transformedPosition.y + attr_screenoffset.y * contentScaleFactor * 2./(screenSize.y*screenVertexScale);

  transformedPosition.xy = (pixelSize * floor((transformedPosition.xy + 1.) / pixelSize))-1.;  // rounding to nearest pixel
  transformedPosition.w = 1.;
  gl_Position = transformedPosition;
  textureLookup = (attr_texcoords * textureSize) + vec2(.5, .5);
  textureLookup = floor(textureLookup) / textureSize;
  normalizedViewCoordinate = (gl_Position.xyz/gl_Position.w) / 2.0 + 0.5;
  if (fog_enabled > 0.5) {
#ifdef PURE_OPENGL_ES_2
  vec3 eye_pos = vec3(g_ModelViewMatrix * attr_worldpos);
  fog = (g_Fog_end - g_Fog_start) * g_Fog_scale;
#else
  vec3 eye_pos = vec3(isProjected * (gl_ModelViewMatrix * attr_worldpos) + ((1.-isProjected) * attr_worldpos));
  fog = (gl_Fog.end - abs(eye_pos.z)) * gl_Fog.scale;
#endif
    fog = max(fog, 0.0);
  } else {
    fog = 1.1; // >= 1.0
  }
  pickcolor = attr_pickcolor;
}
