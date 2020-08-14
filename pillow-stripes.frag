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
    float base = 2.*sqrt(max(1.-r_sq, 0.)) + x_sq - y_sq;
    float flip = 2.*sqrt(max(1.-r_sq, 0.) + x_sq*y_sq);
    return vec2(
        (r_sq - flip) / base,
        base / (r_sq + flip)
    );
}

vec2 conj(vec2 z) { return vec2(z.x, -z.y); }

vec2 peirce_proj(vec3 u) {
    vec2 zeta = u.xy;
    vec2 raw_angles = acos(clamp(cn_coords(zeta), -1., 1.));
    vec2 angles = sign(conj(zeta)) * (vec2(PI, 0.) - raw_angles);
    return 0.5*vec2(F(angles.x, 0.5), F(angles.y, 0.5));
}

// --- euler angles ---

mat3 rot_xy(float t) {
    return mat3(
         cos(t), sin(t), 0.0,
        -sin(t), cos(t), 0.0,
            0.0,    0.0, 1.0
    );
}

mat3 rot_yz(float t) {
    return mat3(
        1.0,     0.0,    0.0,
        0.0,  cos(t), sin(t),
        0.0, -sin(t), cos(t)
    );
}

// attitude = vec3(precession, nutation spin)
mat3 euler_rot(vec3 attitude) {
    return rot_xy(attitude[0]) * rot_yz(attitude[1]) * rot_xy(attitude[2]);
}

// --- pillowcase pattern ---


const vec3 color_a = vec3(1., 0.4, 0.6);
const vec3 color_b = vec3(1., 0.6, 0.6);
const vec3 color_c = vec3(1.);
const vec3 color_d = vec3(1., 0.75, 0.);
const vec3 color_e = vec3(1., 0.6, 0.2);

vec3 stripe(vec2 z, float sheet) {
    float s;
    if (sheet > 0.) {
        s = 8.*(z.x + 1.);
        if (s < 1.) {
            return color_e;
        } else if (s < 3.) {
            return color_d;
        } else if (s < 5.) {
            return color_c;
        } else if (s < 7.) {
            return color_b;
        } else if (s < 9.) {
            return color_a;
        } else if (s < 11.) {
            return color_b;
        } else if (s < 13.) {
            return color_c;
        } else if (s < 15.) {
            return color_d;
        } else {
            return color_e;
        }
    } else {
        s = 8.*(z.y + 1.);
        if (s < 1.) {
            return color_a;
        } else if (s < 3.) {
            return color_b;
        } else if (s < 5.) {
            return color_c;
        } else if (s < 7.) {
            return color_d;
        } else if (s < 9.) {
            return color_e;
        } else if (s < 11.) {
            return color_d;
        } else if (s < 13.) {
            return color_c;
        } else if (s < 15.) {
            return color_b;
        } else {
            return color_a;
        }
    }
}

void mainImage(out vec4 fragColor, in vec2 fragCoord) {
    float small_dim = min(iResolution.x, iResolution.y);
    vec2 p = 2.2*(fragCoord - 0.5*iResolution.xy)/small_dim - vec2(0.8, 0.);
    vec3 color = vec3(0.1);
    
    float r_sq = dot(p, p);
    if (r_sq < 1.) {
        mat3 orient = euler_rot(vec3(0., iTime, 0.));
        vec3 u = orient * vec3(p, sqrt(1. - r_sq));
        color = stripe(peirce_proj(u)/K(0.5), u.z);
    } else {
        vec2 p_mini = 2.*(p - vec2(-1.65, 0.25));
        if (abs(p_mini.x) + abs(p_mini.y) < 1.) {
            color = stripe(p_mini, 1.);
        }
        
        p_mini -= vec2(-1.);
        mat2 unfold = mat2(0., -1., -1., 0.);
        if (abs(p_mini.x) + abs(p_mini.y) < 1.) {
            color = stripe(unfold * p_mini, -1.);
        }
    }
    fragColor = vec4(color, 1.);
}

/*void mainImage(out vec4 fragColor, in vec2 fragCoord) {
    vec2 p = 2.*(fragCoord - 0.5*iResolution.xy)/iResolution.xy;
    vec2 zeta = 0.99999*vec2(cos(p.x*PI), sin(p.x*PI));
    vec2 z = peirce_proj(vec3(zeta, 0.))/K(0.5);
    vec3 color = vec3(1., 1., 1.);
    if (p.y < z.x) {
        color.y *= 0.2;
        color.z *= 0.5;
    }
    if (p.y < z.y) {
        color.x *= 0.2;
        color.y *= 0.5;
    }
    fragColor = vec4(color, 1.);
}*/
