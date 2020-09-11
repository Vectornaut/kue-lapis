// --- complex arithmetic ---

const vec2 ZERO = vec2(0.);
const vec2 ONE  = vec2(1., 0.);
const vec2 I    = vec2(0., 1.); //// UNUSED

//  the complex conjugate of `z`
vec2 conj(vec2 z) {
    return vec2(z.x, -z.y);
}

// the product of `z` and `w`
vec2 mul(vec2 z, vec2 w) {
    return mat2(z, conj(z).yx) * w;
}

// the reciprocal of `z`
vec2 rcp(vec2 z) {
    // 1/z = z'/(z'*z) = z'/|z|^2
    return conj(z) / dot(z, z);
}

// the square root of `z`
//
// Stanley Rabinowitz, "How to Find the Square Root of a Complex Number"
// Mathematics and Informatics Quarterly, vol. 3, pp. 54--56, 1993
//
// Timm Ahrendt, "Fast high-precision computation of complex square roots"
// Proceedings of ISSAC '96, pp. 142--149
// <doi:10.1145/236869.236924>
//
vec2 csqrt(vec2 z) {
    float r = length(z);
    return vec2(sqrt(0.5*(r + z.x)), sign(z.y)*sqrt(0.5*(r - z.x)));
}

vec2 csin(vec2 z) {
    return vec2(sin(z.x) * cosh(z.y), cos(z.x) * sinh(z.y));
}

// --- elliptic integral of the first kind ---
//
// B. C. Carlson, "Numerical computation of real or complex elliptic integrals"
// Numerical Algorithms, vol. 10, pp. 13--26, 1995
// <doi:10.1007/BF02198293>
//
// William Press and Saul Teukolsky, "Elliptic Integrals"
// Computers in Physics, vol. 4, pp. 92--96, 1990
// <doi:10.1063/1.4822893>

const int N = 12;

const vec2 C1 = 1./24. * ONE;
const vec2 C2 = 0.1    * ONE;
const vec2 C3 = 3./44. * ONE;
const vec2 C4 = 1./14. * ONE;

vec2 RF(vec2 x, vec2 y, vec2 z) {
    for (int n = 0; n < N; n++) {
        vec2 sqrt_x = csqrt(x);
        vec2 sqrt_y = csqrt(y);
        vec2 sqrt_z = csqrt(z);
        vec2 lambda = mul(sqrt_x, sqrt_y) + mul(sqrt_y, sqrt_z) + mul(sqrt_z, sqrt_x);
        x = 0.25*(x + lambda);
        y = 0.25*(y + lambda);
        z = 0.25*(z + lambda);
    }
    vec2 avg = (x + y + z)/3.;
    vec2 off_x = x - avg;
    vec2 off_y = y - avg;
    vec2 off_z = z - avg;
    vec2 e2 = mul(off_x, off_y) - mul(off_z, off_z);
    vec2 e3 = mul(mul(off_x, off_y), off_z);
    return mul(ONE + mul(mul(C1, e2) - C2 - mul(C3, e3), e2) + mul(C4, e3), rcp(csqrt(avg)));
}

const float PI = 3.141592653589793;

vec2 K(vec2 m) { return RF(ZERO, ONE - m, ONE); }

// to save work, we precompute K(m) and pass it as K_val
vec2 F(vec2 phi, vec2 m, vec2 K_val) {
    // phi = phi_tile*pi + phi_off, where phi_tile is an integer and phi_off.x
    // is in [-pi/2, pi/2]
    float phi_tile = round(phi.x / PI);
    vec2 phi_off = phi - phi_tile * vec2(PI, 0.);
    
    // integrate from zero to phi_tile*pi
    vec2 val_tile = phi_tile * 2.*K_val;
    
    // integrate from phi_tile*pi to phi
    vec2 s = csin(phi_off);
    vec2 s_sq = mul(s, s);
    vec2 val_off = mul(s, RF(ONE - s_sq, ONE - mul(m, s_sq), ONE));
    
    return val_tile + val_off;
}

// --- Peirce projection ---
//
// James Pierpont, "Note on C. S. Peirce's paper on
//     'A Quincuncial Projection of the sphere'
// American Journal of Mathematics, vol. 18, no. 2, pp. 145--152, 1896

//// ~~~ old version

vec2 cn_coords(vec3 u) {
    float x_sq = u.x * u.x;
    float y_sq = u.y * u.y;
    float z_sq = u.z * u.z;
    float r_sq = x_sq + y_sq;
    float base = 2.*abs(u.z) + x_sq - y_sq;
    float flip = 2.*sqrt(z_sq + x_sq*y_sq);
    return vec2(
        (r_sq - flip) / base,
        base / (r_sq + flip)
    );
}

const float K_1_2 = 1.854074677301372; // the quarter-period K(1/2)

/*vec2 old_peirce_proj(vec3 u) {
    // map into the top-sheet diamond,
    //
    //           (1,  1)
    //            .   .
    //          .       .
    //     (0,  0)     (2,  0)
    //          .       .
    //            .   .
    //           (1, -1)
    //
    vec2 raw_angles = acos(clamp(cn_coords(u), -1., 1.));
    vec2 angles = sign(conj(u.xy)) * (vec2(PI, 0.) - raw_angles);
    vec2 z = (0.5/K_1_2)*vec2(old_F(angles.x, 0.5, K_1_2), old_F(angles.y, 0.5, K_1_2)) + vec2(1., 0.);
    
    // if we're on the bottom sheet, reflect across the southwest edge of the
    // top-sheet diamond
    if (u.z > 0.) return z; else return -z.yx;
}*/

//// ~~~ new version

// inverse cosine, valid in the right half of the unit disk (Carlson 1995,
// equation 4.20)
vec2 cacos_right(vec2 z) {
    vec2 z_sq = mul(z, z);
    return mul(csqrt(ONE - z_sq), RF(z_sq, ONE, ONE));
}

// inverse cosine
vec2 cacos(vec2 z) {
    if (z.x >= 0.) {
        return cacos_right(z);
    } else {
        return PI*ONE - cacos_right(-z);
    }
}

vec2 peirce_proj(vec3 u) {
    // project stereographically onto the equatorial disk
    vec2 zeta = u.xy / (1. - abs(u.z));
    
    // map into the top-sheet diamond,
    //
    //           (1,  1)
    //            .   .
    //          .       .
    //     (0,  0)     (2,  0)
    //          .       .
    //            .   .
    //           (1, -1)
    //
    vec2 z = F(cacos(-zeta), 0.5*ONE, K_1_2*ONE);
    
    // if we're on the bottom sheet, reflect across the southwest edge of the
    // top-sheet diamond
    if (u.z > 0.) return z; else return -z.yx;
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

const vec3 color_a = vec3(1., 0.2, 0.5);
const vec3 color_b = vec3(1., 0.45, 0.7);
const vec3 color_c = vec3(1., 0.9, 0.95);
const vec3 color_d = vec3(1., 0.64, 0.);
const vec3 color_e = vec3(1., 0.45, 0.);

vec3 stripe(vec2 z, vec2 charge) {
    float s = 0.25 * dot(conj(charge), z.yx); // the signed area (1/4) * D(charge, z)
    float s_off = 16.*abs(s - round(s)); // fold s into the fundamental domain [0, 8]
    if (s_off < 1.) {
        return color_e;
    } else if (s_off < 3.) {
        return color_d;
    } else if (s_off < 5.) {
        return color_c;
    } else if (s_off < 7.) {
        return color_b;
    } else {
        return color_a;
    }
}

const vec2 charge = vec2(1., 2.);

vec3 raw_image(vec2 fragCoord) {
    float small_dim = min(iResolution.x, iResolution.y);
    vec2 p = 2.2*(fragCoord - 0.5*iResolution.xy)/small_dim - vec2(0.8, 0.);
    vec3 color = vec3(0.1, 0.0, 0.2);
    
    float r_sq = dot(p, p);
    if (r_sq < 1.) {
        vec3 attitude = iTime * vec3(1./(2.+PI), 1./2., 1./PI);
        mat3 orient = euler_rot(attitude);
        vec3 u = orient * vec3(p, sqrt(1. - r_sq));
        color = stripe(peirce_proj(u), charge);
        if (u.x > 0.) {
          color *= 0.5;
        }
    } else {
        vec2 p_mini = 2.*(p - vec2(-2.15, 0.25));
        vec2 p_rect = mat2(1., 1., -1., 1.) * p_mini;
        if (0. < p_rect.x && p_rect.x < 2. && abs(p_rect.y) < 2.) {
            color = stripe(p_mini, charge);
        }
    }
    return color;
}

void mainImage(out vec4 fragColor, in vec2 fragCoord) {
    vec2 jiggle = vec2(0.25);
    vec3 color_sum = vec3(0.);
    for (int sgn_x = 0; sgn_x < 2; sgn_x++) {
        for (int sgn_y = 0; sgn_y < 2; sgn_y++) {
            color_sum += raw_image(fragCoord + jiggle);
            jiggle.y = -jiggle.y;
        }
        jiggle.x = -jiggle.x;
    }
    fragColor = vec4(0.25*color_sum, 1.);
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
