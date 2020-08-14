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

float RF(vec3 r) {
    for (int n = 0; n < N; n++) {
      vec3 sqrt_r = sqrt(r);
      vec3 pair = sqrt_r * sqrt_r.yzx;
      float lambda = pair.x + pair.y + pair.z;
      r = 0.25*(r + vec3(lambda));
    }
    float avg = (r.x + r.y + r.z)/3.;
    vec3 off = r - vec3(avg);
    float e2 = off.x * off.y - off.z * off.z;
    float e3 = off.x * off.y * off.z;
    return (1. + (C1*e2 - C2 - C3*e3)*e2 + C4*e3) / sqrt(avg);
}

const float PI = 3.141592653589793;

float K(float m) { return RF(vec3(0., 1. - m, 1.)); }

float F(float phi, float m) {
    // phi = phi_tile*pi + phi_off, with phi_off in [-pi/2, pi/2]
    float phi_tile = round(phi/PI);
    float phi_off = phi - phi_tile*PI;
    
    // integrate from zero to phi_tile*pi
    float val_tile = (phi_tile == 0.) ? 0. : phi_tile * 2.*K(m);
    
    // integrate from phi_tile*pi to phi
    vec2 u = vec2(cos(phi_off), sin(phi_off));
    float val_off = u.y * RF(vec3(u.x*u.x, 1. - m*u.y*u.y, 1.));
    
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

vec2 conj(vec2 z) { return vec2(z.x, -z.y); }

vec2 peirce_proj(vec2 zeta) {
    vec2 raw_angles = acos(clamp(cn_coords(zeta), -1., 1.));
    vec2 angles = sign(conj(zeta)) * (vec2(PI, 0.) - raw_angles);
    return 0.5*vec2(F(angles.x, 0.5), F(angles.y, 0.5));
}

void mainImage(out vec4 fragColor, in vec2 fragCoord) {
    vec2 p = 2.*(fragCoord - 0.5*iResolution.xy)/iResolution.xy;
    vec2 zeta;
    if (mod(iTime, 4.) < 2.) {
        zeta = vec2(p.x, 0.);
    } else {
        zeta = vec2(0., p.x);
    }
    vec2 z = peirce_proj(zeta);
    vec3 color = vec3(1., 1., 1.);
    if (2.*p.y < z.x) {
        color.y *= 0.2;
        color.z *= 0.5;
    }
    if (2.*p.y < z.y) {
        color.x *= 0.2;
        color.y *= 0.5;
    }
    fragColor = vec4(color, 1.);
}

/*void mainImage(out vec4 fragColor, in vec2 fragCoord) {
    float small_dim = min(iResolution.x, iResolution.y);
    vec2 p = 2.2*(fragCoord - 0.5*iResolution.xy)/small_dim;
    vec3 color = vec3(0.1);
    
    if (mod(iTime, 4.) < 2.) {
        if (dot(p, p) < 1.) {
            vec2 z = 0.5*(peirce_proj(p)/K(0.5) + vec2(1.));
            color = vec3(z.x, z.y, 0.7);
        }
    } else {
        if (abs(p.x) + abs(p.y) < 1.) {
            vec2 z = 0.5*(p + vec2(1.));
            color = vec3(z.x, z.y, 0.7);
        }
    }
    fragColor = vec4(color, 1.);
}*/
