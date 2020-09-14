// --- complex arithmetic ---

const vec2 ZERO = vec2(0.);
const vec2 ONE  = vec2(1., 0.);
const vec2 I    = vec2(0., 1.);

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

// the square root of `z`, from the complex arithmetic code listing in
// Appendix C of _Numerical Recipes in C_
//
// William Press, Saul Teukolsky, William Vetterling, and Brian Flannery,
// _Numerical Recipes in C_, 2nd edition. Cambridge University Press, 1992
//
vec2 csqrt(vec2 z) {
    // sqrt(0) = 0
    if (z.x == 0. && z.y == 0.) {
        return vec2(0.);
    }
    
    // calculate w
    vec2 a = abs(z);
    float w;
    if (a.x >= a.y) {
        float sl = a.y / a.x;
        w = sqrt(a.x) * sqrt(0.5*(1. + sqrt(1. + sl*sl)));
    } else {
        float sl = a.x / a.y;
        w = sqrt(a.y) * sqrt(0.5*(sl + sqrt(1. + sl*sl)));
    }
    
    // construct output
    if (z.x >= 0.) {
        return vec2(w, z.y / (2.*w));
    } else if (z.y >= 0.) {
        return vec2(z.y/(2.*w), w);
    } else {
        return -vec2(z.y/(2.*w), w);
    }
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

vec2 K(vec2 m) {
    return RF(ZERO, ONE - m, ONE);
}

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

// inverse cosine, valid in the right half of the unit disk (Carlson 1995,
// equation 4.20)
vec2 cacos_right(vec2 z) {
    vec2 z_sq = mul(z, z);
    return mul(csqrt(ONE - z_sq), RF(z_sq, ONE, ONE));
}

// inverse cosine
vec2 cacos(vec2 z) {
    if (z.x > 0.) return cacos_right(z); else return PI*ONE - cacos_right(-z);
}

// the generalized Peirce projection. its output z is defined by the equation
// cn(z, m) = -zeta, where zeta is the stereographic projection of the unit
// vector u from the south pole onto the equatorial plane. the image of the
// northern hemisphere looks like
//
//           (1,  1)
//            .   .
//          .       .
//     (0,  0)     (2,  0)
//          .       .
//            .   .
//           (1, -1)
//
// in the K(m), iK(1-m) frame
vec2 peirce_proj(vec3 u, vec2 m, vec2 K_val) {
    vec2 zeta = u.xy / (1. + u.z); /* there should be tons of roundoff error near the south pole! why don't we see it? */
    return F(cacos(-zeta), m, K_val);
}

/*vec2 peirce_proj(vec3 u, vec2 m, vec2 K_val) {
    // project stereographically onto the equatorial disk
    vec2 zeta = u.xy / (1. + abs(u.z));
    
    // map into the top-sheet diamond, which looks like
    //
    //           (1,  1)
    //            .   .
    //          .       .
    //     (0,  0)     (2,  0)
    //          .       .
    //            .   .
    //           (1, -1)
    //
    // in the K(m), iK(1-m) frame
    vec2 z = F(cacos(-zeta), m, K_val);
    
    // if we're on the bottom sheet, reflect across the southwest edge of the
    // top-sheet diamond
    if (u.z > 0.) return z; else return -z.yx;
}*/

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
    float s = 0.5 * dot(conj(charge), z.yx); // the signed area (1/2) * D(charge, z)
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

/*vec3 debug_stripe(vec2 z, vec2 charge) {
    float s = 0.5 * dot(conj(charge), z.yx); // the signed area (1/2) * D(charge, z)
    float s_off = 16.*abs(s - round(s)); // fold s into the fundamental domain [0, 8]
    vec3 color;
    if (s_off < 1.) {
        color = color_e;
    } else if (s_off < 3.) {
        color = color_d;
    } else if (s_off < 5.) {
        color = color_c;
    } else if (s_off < 7.) {
        color = color_b;
    } else {
        color = color_a;
    }
    if (z.y > 0.975) {
        color = mix(color, vec3(1., 1., 0.), 0.5);
    } else if (z.y < -0.975) {
        color = mix(color, vec3(0., 0., 1.), 0.5);
    }
    if (z.x < 0.025) {
      color = mix(color, vec3(0.25, 1., 0.), 0.5);
    } else if (z.x > 0.975) {
      color = mix(color, vec3(0., 0.75, 1.), 0.5);
    }
    return color;
}*/

const vec2 charge = vec2(1., -2.);

vec3 raw_image(
    vec2 fragCoord,
    float small_dim,
    mat3 orient,
    vec2 m,
    vec2 K_val,
    mat2 rectify
) {
    vec2 p = 2.2*(fragCoord - 0.5*iResolution.xy)/small_dim - vec2(0.65, 0.);
    vec3 color = vec3(0.1, 0.0, 0.2);
    float r_sq = dot(p, p);
    if (r_sq < 1.) {
        vec3 u = orient * vec3(p, sqrt(1. - r_sq));
        color = stripe(rectify * peirce_proj(u, m, K_val), charge);
        /*float t = mod(iTime, 4.);
        if (t < 2.) {
            color = stripe(rectify * peirce_proj(u, m, K_val), charge);
        } else {
            color = debug_stripe(rectify * peirce_proj(u, m, K_val), charge);
        }*/
    } else {
        vec2 p_mini = 1.2*(p + 2.25*ONE);
        if (0. < p_mini.x && p_mini.x < 1. && abs(p_mini.y) < 1.) {
            color = stripe(p_mini, charge);
            /*color = debug_stripe(p_mini, charge);*/
        }
    }
    return color;
}

void mainImage(out vec4 fragColor, in vec2 fragCoord) {
    // find smal dimension
    float small_dim = min(iResolution.x, iResolution.y);
    
    // set attitude
    vec3 attitude = iTime * 0.2 * vec3(1./(2.+PI), 1./2., 1./PI);
    mat3 orient = euler_rot(attitude);
    /*mat3 orient = mat3(1.);*/
    
    // set modulus
    vec2 m = vec2(0.5 + 0.4*sin(iTime), 0.);
    mat2 quarter_frame = mat2(K(m), mul(I, K(ONE - m)));
    mat2 rectify = inverse(quarter_frame * mat2(1., -1., 1., 1.));
    
    // mix subpixels
    vec2 jiggle = vec2(0.25);
    vec3 color_sum = vec3(0.);
    for (int sgn_x = 0; sgn_x < 2; sgn_x++) {
        for (int sgn_y = 0; sgn_y < 2; sgn_y++) {
            color_sum += raw_image(
                fragCoord + jiggle,
                small_dim,
                orient,
                m,
                quarter_frame[0],
                rectify
            );
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
