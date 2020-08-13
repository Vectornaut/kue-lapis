// --- elliptic integral of the first kind ---
//
// William Press and Saul Teukolsky, "Elliptic Integrals"
// Computers in Physics, vol. 4, pp. 92--96, 1990
// <doi:10.1063/1.4822893>

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

// --- Peirce projection ---
//
// James Pierpont, "Note on C. S. Peirce's paper on
//     'A Quincuncial Projection of the sphere'
// American Journal of Mathematics, vol. 18, no. 2, pp. 145--152, 1896

vec2 cn_coords(vec2 zeta) {
    float x_sq = zeta.x * zeta.x;
    float y_sq = zeta.y * zeta.y;
    float r_sq = x_sq + y_sq;
    float base = 2.*sqrt(1.-r_sq) + x_sq - y_sq;
    float flip = 2.*sqrt(1.-r_sq + x_sq*y_sq);
    return vec2(
        (r_sq - flip) / base,
        base / (r_sq + flip)
    );
}

void mainImage(out vec4 fragColor, in vec2 fragCoord) {
    float small_dim = min(iResolution.x, iResolution.y);
    vec2 p = 2.2*(fragCoord - 0.5*iResolution.xy)/small_dim;
    
    if (dot(p, p) < 1.) {
        vec2 z = cn_coords(p);
        
        // set test and coordinate switches
        float t = mod(iTime, 8.);
        float s = mod(t, 4.);
        
        // pick coordinate
        float u = t < 4. ? z.x : z.y;
        
        // test coordinate
        if (s < 2.) {
            // test bounds
            if (abs(u) > 1.) {
                fragColor = vec4(1.);
            } else if (abs(u) > 0.1) {
                fragColor = vec4(0.4);
            } else {
                fragColor = vec4(0.2);
            }
        } else {
            // show value
            if (u < 0.) {
                fragColor = vec4(1.+u, 1., 0.7, 1.);
            } else {
                fragColor = vec4(1, 1.-u, 0.7, 1.);
            }
        }
    } else {
        fragColor = vec4(vec3(0.1), 1.);
    }
}
