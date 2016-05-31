




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
vec2 SMAALumaEdgeDetectionPS(vec2 texcoord,
                               vec4 offset[3],
                               sampler2D colorTex



                               ) {




    vec2 threshold = vec2(SMAA_THRESHOLD, SMAA_THRESHOLD);



    vec3 weights = vec3(0.2126, 0.7152, 0.0722);
    float L = dot(texture(colorTex, texcoord).rgb, weights);

    float Lleft = dot(texture(colorTex, offset[0].xy).rgb, weights);
    float Ltop = dot(texture(colorTex, offset[0].zw).rgb, weights);


    vec4 delta;
    delta.xy = abs(L - vec2(Lleft, Ltop));
    vec2 edges = step(threshold, delta.xy);


    if (dot(edges, vec2(1.0, 1.0)) == 0.0)
        discard;


    float Lright = dot(texture(colorTex, offset[1].xy).rgb, weights);
    float Lbottom = dot(texture(colorTex, offset[1].zw).rgb, weights);
    delta.zw = abs(L - vec2(Lright, Lbottom));


    vec2 maxDelta = max(delta.xy, delta.zw);


    float Lleftleft = dot(texture(colorTex, offset[2].xy).rgb, weights);
    float Ltoptop = dot(texture(colorTex, offset[2].zw).rgb, weights);
    delta.zw = abs(vec2(Lleft, Ltop) - vec2(Lleftleft, Ltoptop));


    maxDelta = max(maxDelta.xy, delta.zw);
    float finalDelta = max(maxDelta.x, maxDelta.y);


    edges.xy *= step(finalDelta, 2.0 * delta.xy);

    return edges;
}







vec2 SMAAColorEdgeDetectionPS(vec2 texcoord,
                                vec4 offset[3],
                                sampler2D colorTex



                                ) {




    vec2 threshold = vec2(SMAA_THRESHOLD, SMAA_THRESHOLD);



    vec4 delta;
    vec3 C = texture(colorTex, texcoord).rgb;

    vec3 Cleft = texture(colorTex, offset[0].xy).rgb;
    vec3 t = abs(C - Cleft);
    delta.x = max(max(t.r, t.g), t.b);

    vec3 Ctop = texture(colorTex, offset[0].zw).rgb;
    t = abs(C - Ctop);
    delta.y = max(max(t.r, t.g), t.b);


    vec2 edges = step(threshold, delta.xy);


    if (dot(edges, vec2(1.0, 1.0)) == 0.0)
        discard;


    vec3 Cright = texture(colorTex, offset[1].xy).rgb;
    t = abs(C - Cright);
    delta.z = max(max(t.r, t.g), t.b);

    vec3 Cbottom = texture(colorTex, offset[1].zw).rgb;
    t = abs(C - Cbottom);
    delta.w = max(max(t.r, t.g), t.b);


    vec2 maxDelta = max(delta.xy, delta.zw);


    vec3 Cleftleft = texture(colorTex, offset[2].xy).rgb;
    t = abs(C - Cleftleft);
    delta.z = max(max(t.r, t.g), t.b);

    vec3 Ctoptop = texture(colorTex, offset[2].zw).rgb;
    t = abs(C - Ctoptop);
    delta.w = max(max(t.r, t.g), t.b);


    maxDelta = max(maxDelta.xy, delta.zw);
    float finalDelta = max(maxDelta.x, maxDelta.y);


    edges.xy *= step(finalDelta, 2.0 * delta.xy);

    return edges;
}




vec2 SMAADepthEdgeDetectionPS(vec2 texcoord,
                                vec4 offset[3],
                                sampler2D depthTex) {
    vec3 neighbours = SMAAGatherNeighbours(texcoord, offset, depthTex);
    vec2 delta = abs(neighbours.xx - vec2(neighbours.y, neighbours.z));
    vec2 edges = step((0.1 * SMAA_THRESHOLD), delta);

    if (dot(edges, vec2(1.0, 1.0)) == 0.0)
        discard;

    return edges;
}





#ifndef SMAA_DISABLE_DIAG_DETECTION




vec2 SMAADecodeDiagBilinearAccess(vec2 e) {
    e.r = e.r * abs(5.0 * e.r - 5.0 * 0.75);
    return round(e);
}

vec4 SMAADecodeDiagBilinearAccess(vec4 e) {
    e.rb = e.rb * abs(5.0 * e.rb - 5.0 * 0.75);
    return round(e);
}




vec2 SMAASearchDiag1(sampler2D edgesTex, vec2 texcoord, vec2 dir, inout vec2 e) {
    vec4 coord = vec4(texcoord, -1.0, 1.0);
    vec3 t = vec3(SMAA_RT_METRICS.xy, 1.0);
    while (coord.z < float(SMAA_MAX_SEARCH_STEPS_DIAG - 1) &&
           coord.w > 0.9) {
        coord.xyz = (t * vec3(dir, 1.0) + coord.xyz);
        e = textureLod(edgesTex, coord.xy, 0.0).rg;
        coord.w = dot(e, vec2(0.5, 0.5));
    }
    return coord.zw;
}

vec2 SMAASearchDiag2(sampler2D edgesTex, vec2 texcoord, vec2 dir, inout vec2 e) {
    vec4 coord = vec4(texcoord, -1.0, 1.0);
    coord.x += 0.25 * SMAA_RT_METRICS.x;
    vec3 t = vec3(SMAA_RT_METRICS.xy, 1.0);
    while (coord.z < float(SMAA_MAX_SEARCH_STEPS_DIAG - 1) &&
           coord.w > 0.9) {
        coord.xyz = (t * vec3(dir, 1.0) + coord.xyz);



        e = textureLod(edgesTex, coord.xy, 0.0).rg;
        e = SMAADecodeDiagBilinearAccess(e);





        coord.w = dot(e, vec2(0.5, 0.5));
    }
    return coord.zw;
}





vec2 SMAAAreaDiag(sampler2D areaTex, vec2 dist, vec2 e, float offset) {
    vec2 texcoord = (vec2(20, 20) * e + dist);


    texcoord = ((1.0 / vec2(160.0, 560.0)) * texcoord + 0.5 * (1.0 / vec2(160.0, 560.0)));


    texcoord.x += 0.5;


    texcoord.y += (1.0 / 7.0) * offset;


    return textureLod(areaTex, texcoord, 0.0).rg;
}




vec2 SMAACalculateDiagWeights(sampler2D edgesTex, sampler2D areaTex, vec2 texcoord, vec2 e, vec4 subsampleIndices) {
    vec2 weights = vec2(0.0, 0.0);


    vec4 d;
    vec2 end;
    if (e.r > 0.0) {
        d.xz = SMAASearchDiag1(edgesTex, texcoord, vec2(-1.0, 1.0), end);
        d.x += float(end.y > 0.9);
    } else
        d.xz = vec2(0.0, 0.0);
    d.yw = SMAASearchDiag1(edgesTex, texcoord, vec2(1.0, -1.0), end);

   
    if (d.x + d.y > 2.0) {

        vec4 coords = (vec4(-d.x + 0.25, d.x, d.y, -d.y - 0.25) * SMAA_RT_METRICS.xyxy + texcoord.xyxy);
        vec4 c;
        c.xy = textureLodOffset(edgesTex, coords.xy, 0.0, ivec2(-1, 0)).rg;
        c.zw = textureLodOffset(edgesTex, coords.zw, 0.0, ivec2( 1, 0)).rg;
        c.yxwz = SMAADecodeDiagBilinearAccess(c.xyzw);
        vec2 cc = (vec2(2.0, 2.0) * c.xz + c.yw);


        SMAAMovc(bvec2(step(0.9, d.zw)), cc, vec2(0.0, 0.0));


        weights += SMAAAreaDiag(areaTex, d.xy, cc, subsampleIndices.z);
    }


    d.xz = SMAASearchDiag2(edgesTex, texcoord, vec2(-1.0, -1.0), end);
    if (textureLodOffset(edgesTex, texcoord, 0.0, ivec2(1, 0)).r > 0.0) {
        d.yw = SMAASearchDiag2(edgesTex, texcoord, vec2(1.0, 1.0), end);
        d.y += float(end.y > 0.9);
    } else
        d.yw = vec2(0.0, 0.0);

   
    if (d.x + d.y > 2.0) {

        vec4 coords = (vec4(-d.x, -d.x, d.y, d.y) * SMAA_RT_METRICS.xyxy + texcoord.xyxy);
        vec4 c;
        c.x = textureLodOffset(edgesTex, coords.xy, 0.0, ivec2(-1, 0)).g;
        c.y = textureLodOffset(edgesTex, coords.xy, 0.0, ivec2( 0, -1)).r;
        c.zw = textureLodOffset(edgesTex, coords.zw, 0.0, ivec2( 1, 0)).gr;
        vec2 cc = (vec2(2.0, 2.0) * c.xz + c.yw);


        SMAAMovc(bvec2(step(0.9, d.zw)), cc, vec2(0.0, 0.0));


        weights += SMAAAreaDiag(areaTex, d.xy, cc, subsampleIndices.w).gr;
    }

    return weights;
}
#endif
float SMAASearchLength(sampler2D searchTex, vec2 e, float offset) {


    vec2 scale = vec2(66.0, 33.0) * vec2(0.5, -1.0);
    vec2 bias = vec2(66.0, 33.0) * vec2(offset, 1.0);


    scale += vec2(-1.0, 1.0);
    bias += vec2( 0.5, -0.5);



    scale *= 1.0 / vec2(64.0, 16.0);
    bias *= 1.0 / vec2(64.0, 16.0);


    return textureLod(searchTex, (scale * e + bias), 0.0).r;
}




float SMAASearchXLeft(sampler2D edgesTex, sampler2D searchTex, vec2 texcoord, float end) {







    vec2 e = vec2(0.0, 1.0);
    while (texcoord.x > end &&
           e.g > 0.8281 &&
           e.r == 0.0) {
        e = textureLod(edgesTex, texcoord, 0.0).rg;
        texcoord = (-vec2(2.0, 0.0) * SMAA_RT_METRICS.xy + texcoord);
    }

    float offset = (-(255.0 / 127.0) * SMAASearchLength(searchTex, e, 0.0) + 3.25);
    return (SMAA_RT_METRICS.x * offset + texcoord.x);
}

float SMAASearchXRight(sampler2D edgesTex, sampler2D searchTex, vec2 texcoord, float end) {
    vec2 e = vec2(0.0, 1.0);
    while (texcoord.x < end &&
           e.g > 0.8281 &&
           e.r == 0.0) {
        e = textureLod(edgesTex, texcoord, 0.0).rg;
        texcoord = (vec2(2.0, 0.0) * SMAA_RT_METRICS.xy + texcoord);
    }
    float offset = (-(255.0 / 127.0) * SMAASearchLength(searchTex, e, 0.5) + 3.25);
    return (-SMAA_RT_METRICS.x * offset + texcoord.x);
}

float SMAASearchYUp(sampler2D edgesTex, sampler2D searchTex, vec2 texcoord, float end) {
    vec2 e = vec2(1.0, 0.0);
    while (texcoord.y > end &&
           e.r > 0.8281 &&
           e.g == 0.0) {
        e = textureLod(edgesTex, texcoord, 0.0).rg;
        texcoord = (-vec2(0.0, 2.0) * SMAA_RT_METRICS.xy + texcoord);
    }
    float offset = (-(255.0 / 127.0) * SMAASearchLength(searchTex, e.gr, 0.0) + 3.25);
    return (SMAA_RT_METRICS.y * offset + texcoord.y);
}

float SMAASearchYDown(sampler2D edgesTex, sampler2D searchTex, vec2 texcoord, float end) {
    vec2 e = vec2(1.0, 0.0);
    while (texcoord.y < end &&
           e.r > 0.8281 &&
           e.g == 0.0) {
        e = textureLod(edgesTex, texcoord, 0.0).rg;
        texcoord = (vec2(0.0, 2.0) * SMAA_RT_METRICS.xy + texcoord);
    }
    float offset = (-(255.0 / 127.0) * SMAASearchLength(searchTex, e.gr, 0.5) + 3.25);
    return (-SMAA_RT_METRICS.y * offset + texcoord.y);
}





vec2 SMAAArea(sampler2D areaTex, vec2 dist, float e1, float e2, float offset) {

    vec2 texcoord = (vec2(16, 16) * round(4.0 * vec2(e1, e2)) + dist);


    texcoord = ((1.0 / vec2(160.0, 560.0)) * texcoord + 0.5 * (1.0 / vec2(160.0, 560.0)));


    texcoord.y = ((1.0 / 7.0) * offset + texcoord.y);


    return textureLod(areaTex, texcoord, 0.0).rg;
}




void SMAADetectHorizontalCornerPattern(sampler2D edgesTex, inout vec2 weights, vec4 texcoord, vec2 d) {

#ifndef SMAA_DISABLE_CORNER_DETECTION
    vec2 leftRight = step(d.xy, d.yx);
    vec2 rounding = (1.0 - (float(SMAA_CORNER_ROUNDING) / 100.0)) * leftRight;

    rounding /= leftRight.x + leftRight.y;

    vec2 factor = vec2(1.0, 1.0);
    factor.x -= rounding.x * textureLodOffset(edgesTex, texcoord.xy, 0.0, ivec2(0, 1)).r;
    factor.x -= rounding.y * textureLodOffset(edgesTex, texcoord.zw, 0.0, ivec2(1, 1)).r;
    factor.y -= rounding.x * textureLodOffset(edgesTex, texcoord.xy, 0.0, ivec2(0, -2)).r;
    factor.y -= rounding.y * textureLodOffset(edgesTex, texcoord.zw, 0.0, ivec2(1, -2)).r;

    weights *= clamp(factor, 0.0, 1.0);

#endif
}

void SMAADetectVerticalCornerPattern(sampler2D edgesTex, inout vec2 weights, vec4 texcoord, vec2 d) {

#ifndef SMAA_DISABLE_CORNER_DETECTION
    vec2 leftRight = step(d.xy, d.yx);
    vec2 rounding = (1.0 - (float(SMAA_CORNER_ROUNDING) / 100.0)) * leftRight;

    rounding /= leftRight.x + leftRight.y;

    vec2 factor = vec2(1.0, 1.0);
    factor.x -= rounding.x * textureLodOffset(edgesTex, texcoord.xy, 0.0, ivec2( 1, 0)).g;
    factor.x -= rounding.y * textureLodOffset(edgesTex, texcoord.zw, 0.0, ivec2( 1, 1)).g;
    factor.y -= rounding.x * textureLodOffset(edgesTex, texcoord.xy, 0.0, ivec2(-2, 0)).g;
    factor.y -= rounding.y * textureLodOffset(edgesTex, texcoord.zw, 0.0, ivec2(-2, 1)).g;

    weights *= clamp(factor, 0.0, 1.0);

#endif
}




vec4 SMAABlendingWeightCalculationPS(vec2 texcoord,
                                       vec2 pixcoord,
                                       vec4 offset[3],
                                       sampler2D edgesTex,
                                       sampler2D areaTex,
                                       sampler2D searchTex,
                                       vec4 subsampleIndices) {
    vec4 weights = vec4(0.0, 0.0, 0.0, 0.0);

    vec2 e = texture(edgesTex, texcoord).rg;

   
    if (e.g > 0.0) {

#ifndef SMAA_DISABLE_DIAG_DETECTION


        weights.rg = SMAACalculateDiagWeights(edgesTex, areaTex, texcoord, e, subsampleIndices);



       
        if (weights.r == -weights.g) {

#endif

        vec2 d;


        vec3 coords;
        coords.x = SMAASearchXLeft(edgesTex, searchTex, offset[0].xy, offset[2].x);
        coords.y = offset[1].y;
        d.x = coords.x;




        float e1 = textureLod(edgesTex, coords.xy, 0.0).r;


        coords.z = SMAASearchXRight(edgesTex, searchTex, offset[0].zw, offset[2].y);
        d.y = coords.z;



        d = abs(round((SMAA_RT_METRICS.zz * d + -pixcoord.xx)));



        vec2 sqrt_d = sqrt(d);


        float e2 = textureLodOffset(edgesTex, coords.zy, 0.0, ivec2(1, 0)).r;



        weights.rg = SMAAArea(areaTex, sqrt_d, e1, e2, subsampleIndices.y);


        coords.y = texcoord.y;
        SMAADetectHorizontalCornerPattern(edgesTex, weights.rg, coords.xyzy, d);


#ifndef SMAA_DISABLE_DIAG_DETECTION
        } else
            e.r = 0.0;

#endif
    }

   
    if (e.r > 0.0) {
        vec2 d;


        vec3 coords;
        coords.y = SMAASearchYUp(edgesTex, searchTex, offset[1].xy, offset[2].z);
        coords.x = offset[0].x;
        d.x = coords.y;


        float e1 = textureLod(edgesTex, coords.xy, 0.0).g;


        coords.z = SMAASearchYDown(edgesTex, searchTex, offset[1].zw, offset[2].w);
        d.y = coords.z;


        d = abs(round((SMAA_RT_METRICS.ww * d + -pixcoord.yy)));



        vec2 sqrt_d = sqrt(d);


        float e2 = textureLodOffset(edgesTex, coords.xz, 0.0, ivec2(0, 1)).g;


        weights.ba = SMAAArea(areaTex, sqrt_d, e1, e2, subsampleIndices.x);


        coords.x = texcoord.x;
        SMAADetectVerticalCornerPattern(edgesTex, weights.ba, coords.xyxz, d);
    }

    return weights;
}




vec4 SMAANeighborhoodBlendingPS(vec2 texcoord,
                                  vec4 offset,
                                  sampler2D colorTex,
                                  sampler2D blendTex



                                  ) {

    vec4 a;
    a.x = texture(blendTex, offset.xy).a;
    a.y = texture(blendTex, offset.zw).g;
    a.wz = texture(blendTex, texcoord).xz;


   
    if (dot(a, vec4(1.0, 1.0, 1.0, 1.0)) < 1e-5) {
        vec4 color = textureLod(colorTex, texcoord, 0.0);
        return color;
    } else {
        bool h = max(a.x, a.z) > max(a.y, a.w);


        vec4 blendingOffset = vec4(0.0, a.y, 0.0, a.w);
        vec2 blendingWeight = a.yw;
        SMAAMovc(bvec4(h, h, h, h), blendingOffset, vec4(a.x, 0.0, a.z, 0.0));
        SMAAMovc(bvec2(h, h), blendingWeight, a.xz);
        blendingWeight /= dot(blendingWeight, vec2(1.0, 1.0));


        vec4 blendingCoord = (blendingOffset * vec4(SMAA_RT_METRICS.xy, -SMAA_RT_METRICS.xy) + texcoord.xyxy);



        vec4 color = blendingWeight.x * textureLod(colorTex, blendingCoord.xy, 0.0);
        color += blendingWeight.y * textureLod(colorTex, blendingCoord.zw, 0.0);
        return color;
    }
}




vec4 SMAAResolvePS(vec2 texcoord,
                     sampler2D currentColorTex,
                     sampler2D previousColorTex



                     ) {
    vec4 current = texture(currentColorTex, texcoord);
    vec4 previous = texture(previousColorTex, texcoord);
    return mix(current, previous, 0.5);

}
