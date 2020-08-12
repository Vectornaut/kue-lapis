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

float rf(vec3 u) {
    for (int n = 0; n < N; n++) {
      vec3 sqrt_u = sqrt(u);
      vec3 pair = sqrt_u * sqrt_u.yzx;
      float lambda = pair.x + pair.y + pair.z;
      u = 0.25*(u + lambda);
    }
    float avg = (u.x + u.y + u.z)/3.;
    vec3 off = u - vec3(avg);
    float e2 = off.x * off.y - off.z * off.z;
    float e3 = off.x * off.y * off.z;
    return (1. + (C1*e2 - C2 - C3*e3)*e2 + C4*e3) / sqrt(avg);
}

void mainImage(out vec4 fragColor, in vec2 fragCoord) {
    vec2 p = fragCoord / iResolution.xy;
    if (4.*p.y < rf(vec3(p.x, 0.5, 0.5))) {
        fragColor = vec4(0.5, 0.2, 0.9, 1.);
    } else {
        fragColor = vec4(1.);
    }
}
