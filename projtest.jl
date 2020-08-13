using Elliptic, Plots

# --- elliptic integral of the first kind ---

C1 = 1/24.;
C2 = 0.1;
C3 = 3/44.;
C4 = 1/14.;

function RF(u, N)
    for n in 1:N
      sqrt_u = sqrt.(u)
      pair = sqrt_u .* sqrt_u[[2,3,1]]
      lambda = pair[1] + pair[2] + pair[3]
      u = 0.25*(u + fill(lambda, 3))
    end
    avg = (u[1] + u[2] + u[3])/3
    off = u - fill(avg, 3)
    e2 = off[1] * off[2] - off[3] * off[3]
    e3 = off[1] * off[2] * off[3]
    return (1 + (C1*e2 - C2 - C3*e3)*e2 + C4*e3) / sqrt(avg)
end

function my_F(phi, k, N = 12)
    c = cos(phi)
    s = sin(phi)
    ks = k*s
    return s * RF([c*c, 1 - ks*ks, 1], N)
end

# --- Peirce projection

function cn_coords(zeta)
    x_sq = real(zeta) * real(zeta);
    y_sq = imag(zeta) * imag(zeta);
    r_sq = x_sq + y_sq;
    base = 2*sqrt(1-r_sq) + x_sq - y_sq;
    flip = 2*sqrt(1-r_sq + x_sq*y_sq);
    return (r_sq - flip) / base + im*(base / (r_sq + flip))
end

SQRT_1_2 = 1/sqrt(2);

function my_peirce_proj(zeta, N)
    angles = acos.(reim(cn_coords(zeta)))
    return 0.5*[my_F(angles[1], SQRT_1_2, N), my_F(angles[2], SQRT_1_2, N)]
end

function good_peirce_proj(zeta)
    angles = acos.(reim(cn_coords(zeta)))
    return 0.5*[F(angles[1], 1/2), F(angles[2], 1/2)]
end

# --- test ---

function test_peirce_proj(N = 12)
    mesh = LinRange(0, 1, 200)
    my_z = [my_peirce_proj(zeta, N) for zeta in mesh]
    good_z = good_peirce_proj.(mesh)
    plot(mesh, [first.(my_z), first.(good_z)])
end
