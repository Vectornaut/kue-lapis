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

const float PI = 3.141592653589793;

float K(float k) { return RF(vec3(0., 1. - k*k, 1.)); }

float F(float phi, float k) {
    // phi = phi_tile*pi + phi_off, with phi_off in [-pi/2, pi/2]
    float phi_tile = round(phi/PI);
    float phi_off = phi - phi_tile*PI;
    
    // integrate from zero to phi_tile*pi
    float val_tile = (phi_tile == 0.) ? 0. : phi_tile * 2.*K(k);
    
    // integrate from phi_tile*pi to phi
    float c = cos(phi_off);
    float s = sin(phi_off);
    float ks = k*s;
    float val_off = s * RF(vec3(c*c, 1. - ks*ks, 1.));
    
    return val_tile + val_off;
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

const float SQRT_1_2 = 0.7071067811865475;

vec2 peirce_proj(vec2 zeta) {
    vec2 angles = acos(cn_coords(zeta));
    return 0.5*vec2(F(angles.x, SQRT_1_2), F(angles.y, SQRT_1_2));
    /*return angles;*/
}

void mainImage(out vec4 fragColor, in vec2 fragCoord) {
    vec2 p = 2.*(fragCoord - 0.5*iResolution.xy)/iResolution.xy;
    if (8.*p.y < F(1.5*PI*p.x, SQRT_1_2)) {
        fragColor = vec4(0.6, 0.2, 0.9, 1.);
    } else {
        fragColor = vec4(1.);
    }
}

/*void mainImage(out vec4 fragColor, in vec2 fragCoord) {
    float small_dim = min(iResolution.x, iResolution.y);
    vec2 p = 2.2*(fragCoord - 0.5*iResolution.xy)/small_dim;
    vec3 color = vec3(0.1);
    
    if (p.x > 0. && p.y > 0. && dot(p, p) < 1.) {
        vec2 z = (2./K(SQRT_1_2))*peirce_proj(p);
        float u = mod(iTime, 4.) < 2. ? z.x : z.y;
        if (u < 0. || 1. < u) {
            color = vec3(1.);
        } else {
            color = vec3(u, 0., 0.7);
        }
    }
    fragColor = vec4(color, 1.);
}*/
