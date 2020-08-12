// by Mike Nolta, from Elliptic.jl
const int N = 10;
const float two_to_N = 1024.;
float am(float u, float m) {
    float [N] ambuf;
    float a = 1.;
    float b = sqrt(1.-m);
    float c = sqrt(m);
    for (int n = 0; n < N; n++) {
        a = 0.5*(a+b);
        b = sqrt(a*b);
        c = 0.5*(a-b);
        ambuf[n] = c/a;
    }
    
    float phi = a * u * two_to_N;
    for (int n = N-1; n >= 0; n--) {
        phi = 0.5*(phi + asin(ambuf[n] * sin(phi)));
    }
    return phi;
}

const float TAU = 6.283185307179586;

void mainImage(out vec4 fragColor, in vec2 fragCoord) {
    vec2 p = TAU * fragCoord / iResolution.xy;
    if (p.y < am(p.x, 1./sqrt(2.))) {
        fragColor = vec4(0.5, 0.2, 0.9, 1.);
    } else {
        fragColor = vec4(1.);
    }
}
