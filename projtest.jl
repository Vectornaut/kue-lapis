using Elliptic, Plots

# --- square root ---

function my_sqrt(t)
    r = abs(t)
    return sqrt(0.5*(r + real(t))) + im*sign(imag(t))*sqrt(0.5*(r - real(t)))
end

# --- elliptic integral of the first kind ---

C1 = 1/24.
C2 = 0.1
C3 = 3/44.
C4 = 1/14.

function RF(r, N, use_taylor)
    for n in 1:N
        sqrt_r = my_sqrt.(r)
        pair = sqrt_r .* sqrt_r[[2,3,1]]
        lambda = pair[1] + pair[2] + pair[3]
        r = 0.25*(r + fill(lambda, 3))
    end
    avg = (r[1] + r[2] + r[3])/3
    if use_taylor
        off = r - fill(avg, 3)
        e2 = off[1] * off[2] - off[3] * off[3]
        e3 = off[1] * off[2] * off[3]
        return (1 + (C1*e2 - C2 - C3*e3)*e2 + C4*e3) / my_sqrt(avg)
    else
        return 1/my_sqrt(avg)
    end
end

my_K(m, N = 12, use_taylor = true) = RF([0, 1 - m, 1], N, use_taylor)

function my_F(phi, m, N = 12, use_taylor = true)
    # phi = phi_tile*pi + phi_off, where phi_tile is an integer and
    # real(phi_off) is in [-pi/2, pi/2]
    phi_tile = round(real(phi)/pi)
    phi_off = phi - phi_tile*pi
    
    # integrate from zero to phi_tile*pi
    val_tile = phi_tile * 2my_K(m, N, use_taylor)
    
    # integrate from phi_tile*pi to phi
    s = sin(phi_off)
    s_sq = s*s
    val_off = s * RF([1 - s_sq, 1 - m*s_sq, 1], N, use_taylor)
    
    ### this version might be a bit slower in GLSL because of the extra
    ### conditional in `sign`
    ##c = 1/sin(phi_off)
    ##c *= c
    ##val_off = sign(real(phi_off)) * RF([c-1, c-m, c], N, use_taylor)
    
    return val_tile + val_off
end

# --- Peirce projection

function cn_coords(u)
    x_sq = u[1] * u[1]
    y_sq = u[2] * u[2]
    z_sq = u[3] * u[3]
    r_sq = x_sq + y_sq
    base = base = 2*abs(u[3]) + x_sq - y_sq
    flip = 2*sqrt(z_sq + x_sq*y_sq)
    return [
        (r_sq - flip) / base,
        base / (r_sq + flip)
    ]
end

get_angles(u) = [sign(u[1]), -sign(u[2])] .* ([pi, 0] - acos.(cn_coords(u)))

function good_peirce_proj(u)
    angles = get_angles(u)
    return 0.5*complex(F(angles[1], 1/2), F(angles[2], 1/2))
end

function my_peirce_proj(u, N = 12, use_taylor = true)
    angles = get_angles(u)
    return 0.5*(my_F(angles[1], 0.5, N, use_taylor) + im*my_F(angles[2], 0.5, N, use_taylor))
end

# Carlson 1995, equation 4.20
function my_acos(z, N = 12, use_taylor = true)
    if real(z) >= 0
        z_sq = z*z
        return my_sqrt(1 - z_sq) * RF([z_sq, 1, 1], N, use_taylor)
    else
        return pi - my_acos(-z, N, use_taylor)
    end
end

function cpx_peirce_proj(u, N = 12, use_taylor = true)
    zeta = complex(u[[1, 2]]...) / (1 - u[3])
    return my_F(my_acos(-zeta), 0.5, N, use_taylor)
end

# --- complexified cn (for testing) ---

rot_cn(u::Real, m::Real) = Jacobi.nc(u, 1-m)
rot_sn(u::Real, m::Real) = Jacobi.sc(u, 1-m)*im
rot_dn(u::Real, m::Real) = Jacobi.dc(u, 1-m)

function cn(u::Complex, m::Real)
    x = real(u)
    y = imag(u)
    
    cn_x = Jacobi.cn(x, m);  cn_iy = rot_cn(y, m)
    sn_x = Jacobi.sn(x, m);  sn_iy = rot_sn(y, m)
    dn_x = Jacobi.dn(x, m);  dn_iy = rot_dn(y, m)
    
    ((cn_x * cn_iy) - (sn_x * dn_x) * (sn_iy * dn_iy)) / (1 - m * sn_x^2 * sn_iy^2)
end

# --- test ---

function test_sqrt()
    ax_mesh = LinRange(-2, 2, 6)
    [my_sqrt(x+im*y) - sqrt(x+im*y) for y in reverse(ax_mesh), x in ax_mesh]
end

function test_acos()
    ax_mesh = LinRange(-2, 2, 6)
    [my_acos(x+im*y) - acos(x+im*y) for y in reverse(ax_mesh), x in ax_mesh]
end

# check values from section 3 of
#
#   B. C. Carlson, "Numerical computation of real or complex elliptic integrals"
#   Numerical Algorithms, vol. 10, pp. 13--26, 1995
#   <doi:10.1007/BF02198293>
#
function test_RF(N = 12, use_taylor = true)
    return [
        RF([1, 2, 0], N, use_taylor) - 1.3110287771461,
        RF([im, -im, 0], N, use_taylor) - 1.8540746773014,
        RF([0.5, 1, 0], N, use_taylor) - 1.8540746773014,
        RF([im-1, im, 0], N, use_taylor) - (0.79612586584234 - 1.213856669865im),
        RF([2, 3, 4], N, use_taylor) - 0.58408284167715,
        RF([im, -im, 2], N, use_taylor) - 1.0441445654064,
        RF([im-1, im, 1-im], N, use_taylor) - (0.9391205218619 - 0.53296252018635im)
    ]
end

function test_F(y = 0, N = 12, use_taylor = true)
    mesh = LinRange(-3pi/2, 3pi/2, 600)
    my_z = [my_F(phi + im*y, 0.5, N, use_taylor) for phi in mesh]
    good_z = [F(phi, 1/2) for phi in mesh]
    comparison = plot(
        mesh,
        [real.(my_z), good_z],
        linecolor = [RGB(0.5, 0, 0.5) RGB(1, 0.5, 0.8)],
        linestyle = [:solid :dash],
        legend = false
    )
    plot!(
        comparison,
        mesh,
        [imag.(my_z), imag.(good_z)],
        linecolor = [RGB(0, 0.5, 0) RGB(0.5, 1, 0)],
        linestyle = [:solid :dash]
    )
    crossover = plot(mesh, [cos.(mesh).^2, [1-0.5*sin(phi)^2 for phi in mesh]], legend = false)
    plot(comparison, crossover, layout = (2, 1))
end

function test_inverse(N = 12, use_taylor = true)
    ax_mesh = K(0.5)*LinRange(0, 1, 5)
    mesh = [r*(1+im) + s*(1-im) for r in reverse(ax_mesh), s in ax_mesh]
    [z - my_F(acos(cn(z, 0.5)), 0.5, N, use_taylor) for z in mesh]
end

function test_peirce_proj(long = 0, N = 12, use_taylor = true; use_complex = true)
    mesh = LinRange(-pi, 0, 200)
    longline = [[cos(long)*cos(lat), sin(long)*cos(lat), sin(lat)] for lat in mesh]
    good_z = 1 .+ good_peirce_proj.(longline) / K(0.5)
    if use_complex
        my_z = [cpx_peirce_proj(u, N, use_taylor) for u in longline]  / real(my_K(0.5, N, use_taylor))
    else
        my_z = 1 .+ [my_peirce_proj(u, N, use_taylor) for u in longline] / real(my_K(0.5, N, use_taylor))
    end
    x_plot = plot(
        mesh,
        [real.(my_z) real.(good_z)],
        linecolor = [RGB(0.5, 0, 0.5) RGB(1, 0.5, 0.8)],
        linestyle = [:solid :dash],
        ylims = (0, 2),
        legend = false
    )
    y_plot = plot(
        mesh,
        [imag.(my_z), imag.(good_z)],
        linecolor = [RGB(0, 0.5, 0) RGB(0.5, 1, 0)],
        linestyle = [:solid :dash],
        ylims = (-1, 1),
        legend = false
    )
    plot(x_plot, y_plot, layout = (2, 1))
end

function test_equator(N = 12, use_taylor = true; use_complex = true)
    mesh = LinRange(0, 2pi, 400)
    equator = [[cos(t), sin(t), 0] for t in mesh]
    good_z = 1 .+ good_peirce_proj.(equator) / K(0.5)
    if use_complex
        my_z = [cpx_peirce_proj(u, N, use_taylor) for u in equator]  / real(my_K(0.5, N, use_taylor))
    else
        my_z = 1 .+ [complex(my_peirce_proj(u, N, use_taylor)...) for u in equator] / real(my_K(0.5, N, use_taylor))
    end
    x_plot = plot(
        mesh,
        [real.(my_z) real.(good_z)],
        linecolor = [RGB(0.5, 0, 0.5) RGB(1, 0.5, 0.8)],
        linestyle = [:solid :dash],
        ylims = (0, 2),
        legend = false
    )
    y_plot = plot(
        mesh,
        [imag.(my_z), imag.(good_z)],
        linecolor = [RGB(0, 0.5, 0) RGB(0.5, 1, 0)],
        linestyle = [:solid :dash],
        ylims = (-1, 1),
        legend = false
    )
    plot(x_plot, y_plot, layout = (2, 1))
end
