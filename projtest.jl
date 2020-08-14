using Elliptic, Plots

# --- elliptic integral of the first kind ---

C1 = 1/24.;
C2 = 0.1;
C3 = 3/44.;
C4 = 1/14.;

function RF(u, N, use_taylor)
    for n in 1:N
        sqrt_u = sqrt.(u)
        pair = sqrt_u .* sqrt_u[[2,3,1]]
        lambda = pair[1] + pair[2] + pair[3]
        u = 0.25*(u + fill(lambda, 3))
    end
    avg = (u[1] + u[2] + u[3])/3
    if use_taylor
        off = u - fill(avg, 3)
        e2 = off[1] * off[2] - off[3] * off[3]
        e3 = off[1] * off[2] * off[3]
        return (1 + (C1*e2 - C2 - C3*e3)*e2 + C4*e3) / sqrt(avg)
    else
        return 1/sqrt(avg)
    end
end

my_K(k, N, use_taylor) = RF([0, 1 - k*k, 1], N, use_taylor)

function my_F(phi, k, N = 12, use_taylor = true)
    # phi = phi_tile*pi + phi_off, with phi_off in [-pi/2, pi/2]
    phi_tile = round(phi/pi);
    phi_off = phi - phi_tile*pi;
    
    # integrate from zero to phi_tile*pi
    val_tile = (phi_tile == 0) ? 0 : phi_tile * 2my_K(k, N, use_taylor);
    
    # integrate from phi_tile*pi to phi
    c = cos(phi_off);
    s = sin(phi_off);
    ks = k*s;
    val_off = s * RF([c*c, 1 - ks*ks, 1], N, use_taylor);
    
    return val_tile + val_off;
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

function get_angles(zeta)
    angles = collect(acos.(reim(cn_coords(zeta))))
    if real(zeta) > 0
        angles[1] = pi - angles[1]
    else
        angles[1] = angles[1] - pi
    end
    if imag(zeta) < 0
        angles[2] = -angles[2]
    end
    angles
end

function my_peirce_proj(zeta, N = 12, use_taylor = true)
    angles = get_angles(zeta)
    return 0.5*[my_F(angles[1], SQRT_1_2, N, use_taylor), my_F(angles[2], SQRT_1_2, N, use_taylor)]
end

function good_peirce_proj(zeta)
    angles = get_angles(zeta)
    return 0.5*[F(angles[1], 1/2), F(angles[2], 1/2)]
end

# --- test ---

function test_peirce_proj(dir = 1, N = 12, use_taylor = true)
    mesh = LinRange(-1, 1, 200)
    my_z = [my_peirce_proj(dir*zeta, N, use_taylor) for zeta in mesh]
    good_z = good_peirce_proj.(dir*mesh)
    x_plot = plot(
      mesh,
      [first.(my_z) first.(good_z)],
      linecolor = [RGB(0.5, 0, 0.5) RGB(1, 0.5, 0.8)],
      linestyle = [:solid :dash],
      ylims = (-2, 2),
      legend = false
    )
    y_plot = plot(
      mesh,
      [last.(my_z), last.(good_z)],
      linecolor = [RGB(0, 0.5, 0) RGB(0.5, 1, 0)],
      linestyle = [:solid :dash],
      ylims = (-2, 2),
      legend = false
    )
    plot(x_plot, y_plot, layout = (2, 1))
end

function test_F(N = 12, use_taylor = true)
    mesh = LinRange(-3pi/2, 3pi/2, 600)
    my_z = [my_F(phi, 1/sqrt(2), N, use_taylor) for phi in mesh]
    good_z = [F(phi, 1/2) for phi in mesh]
    comparison = plot(mesh, [first.(my_z), real.(good_z)], legend = false)
    crossover = plot(mesh, [cos.(mesh).^2, [1-0.5*sin(phi)^2 for phi in mesh]], legend = false)
    plot(comparison, crossover, layout = (2, 1))
end
