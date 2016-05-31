uniform float SMAA_THRESHOLD;
uniform int SMAA_MAX_SEARCH_STEPS;
uniform int SMAA_MAX_SEARCH_STEPS_DIAG;
uniform int SMAA_CORNER_ROUNDING;
vec3 SMAAGatherNeighbours(vec2 texcoord,
                            vec4 offset[3],
                            sampler2D tex) {



    float P = texture(tex, texcoord).r;
    float Pleft = texture(tex, offset[0].xy).r;
    float Ptop = texture(tex, offset[0].zw).r;
    return vec3(P, Pleft, Ptop);

}




vec2 SMAACalculatePredicatedThreshold(vec2 texcoord,
                                        vec4 offset[3],
                                        sampler2D predicationTex) {
    vec3 neighbours = SMAAGatherNeighbours(texcoord, offset, predicationTex);
    vec2 delta = abs(neighbours.xx - neighbours.yz);
    vec2 edges = step(0.01, delta);
    return 2.0 * SMAA_THRESHOLD * (1.0 - 0.4 * edges);
}




void SMAAMovc(bvec2 cond, inout vec2 variable, vec2 value) {
    if (cond.x) variable.x = value.x;
    if (cond.y) variable.y = value.y;
}

void SMAAMovc(bvec4 cond, inout vec4 variable, vec4 value) {
    SMAAMovc(cond.xy, variable.xy, value.xy);
    SMAAMovc(cond.zw, variable.zw, value.zw);
}
void SMAAEdgeDetectionVS(vec2 texcoord,
                         inout vec4 offset[3]) {
    offset[0] = (SMAA_RT_METRICS.xyxy * vec4(-1.0, 0.0, 0.0, -1.0) + texcoord.xyxy);
    offset[1] = (SMAA_RT_METRICS.xyxy * vec4( 1.0, 0.0, 0.0, 1.0) + texcoord.xyxy);
    offset[2] = (SMAA_RT_METRICS.xyxy * vec4(-2.0, 0.0, 0.0, -2.0) + texcoord.xyxy);
}




void SMAABlendingWeightCalculationVS(vec2 texcoord,
                                     inout vec2 pixcoord,
                                     inout vec4 offset[3]) {
    pixcoord = texcoord * SMAA_RT_METRICS.zw;


    offset[0] = (SMAA_RT_METRICS.xyxy * vec4(-0.25, -0.125, 1.25, -0.125) + texcoord.xyxy);
    offset[1] = (SMAA_RT_METRICS.xyxy * vec4(-0.125, -0.25, -0.125, 1.25) + texcoord.xyxy);


    offset[2] = (SMAA_RT_METRICS.xxyy * vec4(-2.0, 2.0, -2.0, 2.0) * float(SMAA_MAX_SEARCH_STEPS) + vec4(offset[0].xz, offset[1].yw))

                                                       ;
}




void SMAANeighborhoodBlendingVS(vec2 texcoord,
                                inout vec4 offset) {
    offset = (SMAA_RT_METRICS.xyxy * vec4( 1.0, 0.0, 0.0, 1.0) + texcoord.xyxy);
}