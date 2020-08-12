// --- elliptic integral of the first kind ---
//
//   William Press and Saul Teukolsky, "Elliptic Integrals"
//   Computers in Physics, vol. 4, pp. 92--96, 1990
//   <doi:10.1063/1.4822893>

const int N = 12;

const float C1 = 1./24.;
const float C2 = 0.1;
const float C3 = 3./44.;
const float C4 = 1./14.;

float RF(vec3 u) {
    for (int n = 0; n < N; n++) {
      vec3 sqrt_u = sqrt(u);
      vec3 pair = sqrt_u * sqrt_u.yzx;
      float lambda = pair.x + pair.y + pair.z;
      u = 0.25*(u + vec3(lambda));
    }
    float avg = (u.x + u.y + u.z)/3.;
    vec3 off = u - vec3(avg);
    float e2 = off.x * off.y - off.z * off.z;
    float e3 = off.x * off.y * off.z;
    return (1. + (C1*e2 - C2 - C3*e3)*e2 + C4*e3) / sqrt(avg);
}

float F(float phi, float k) {
    float c = cos(phi);
    float s = sin(phi);
    float ks = k*s;
    return s * RF(vec3(c*c, 1. - ks*ks, 1.));
}

const float PI_2 = 1.5707963267948966;

void mainImage(out vec4 fragColor, in vec2 fragCoord) {
    vec2 p = fragCoord / iResolution.xy;
    vec3 color = vec3(1.);
    
    for (int n = 5; n >= 0; n--) {
        if (3.*p.y < F(PI_2*p.x, sqrt(0.2*float(n)))) {
            color = mix(color, vec3(0.5, 0.2, 0.25*float(n)), 0.2);
        }
    }
    
    fragColor = vec4(color, 1.);
}
