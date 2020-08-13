# --- elliptic integral of the first kind ---

N = 12;

C1 = 1/24.;
C2 = 0.1;
C3 = 3/44.;
C4 = 1/14.;

function RF(u)
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

function my_F(phi, k)
    c = cos(phi)
    s = sin(phi)
    ks = k*s
    return s * RF([c*c, 1 - ks*ks, 1])
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
