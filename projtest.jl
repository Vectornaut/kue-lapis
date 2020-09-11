using Elliptic, Plots

# --- square root ---

function my_sqrt(t)
    r = abs(t)
    return sqrt(0.5*(r + real(t))) + im*sign(imag(t))*sqrt(0.5*(r - real(t)))
end

# --- elliptic integral of the first kind ---

C1 = 1/24.;
C2 = 0.1;
C3 = 3/44.;
C4 = 1/14.;

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
    phi_tile = round(real(phi)/pi);
    phi_off = phi - phi_tile*pi;
    
    # integrate from zero to phi_tile*pi
    val_tile = (phi_tile == 0) ? 0 : phi_tile * 2my_K(m, N, use_taylor);
    
    # integrate from phi_tile*pi to phi
    s = sin(phi_off)
    s_sq = s*s
    val_off = s * RF([1 - s_sq, 1 - m*s_sq, 1], N, use_taylor);
    
    ### this version might be a bit slower in GLSL because of the extra
    ### conditional in `sign`
    ##c = 1/sin(phi_off)
    ##c *= c
    ##val_off = sign(real(phi_off)) * RF([c-1, c-m, c], N, use_taylor)
    
    return val_tile + val_off;
end

# --- Peirce projection

function cn_coords(zeta)
    x_sq = real(zeta) * real(zeta);
    y_sq = imag(zeta) * imag(zeta);
    r_sq = x_sq + y_sq;
    base = 2*sqrt(max(1-r_sq, 0)) + x_sq - y_sq;
    flip = 2*sqrt(max(1-r_sq, 0) + x_sq*y_sq);
    return (r_sq - flip) / base + im*(base / (r_sq + flip))
end

get_angles(zeta) =
  sign.(reim(conj(zeta))) .* ([pi, 0] - collect(acos.(reim(cn_coords(zeta)))))

function my_peirce_proj(zeta, N = 12, use_taylor = true)
    angles = get_angles(zeta)
    return 0.5*[my_F(angles[1], 0.5, N, use_taylor), my_F(angles[2], 0.5, N, use_taylor)]
end

function good_peirce_proj(zeta)
    angles = get_angles(zeta)
    return 0.5*[F(angles[1], 1/2), F(angles[2], 1/2)]
end

# --- test ---

function test_equator(N = 12, use_taylor = true)
    mesh = LinRange(0, 2pi, 400)
    equator = cis.(mesh)
    my_z = [my_peirce_proj(zeta, N, use_taylor) for zeta in equator] / my_K(0.5, N, use_taylor)
    good_z = good_peirce_proj.(equator) / K(0.5)
    x_plot = plot(
      mesh,
      [first.(my_z) first.(good_z)],
      linecolor = [RGB(0.5, 0, 0.5) RGB(1, 0.5, 0.8)],
      linestyle = [:solid :dash],
      ylims = (-1, 1),
      legend = false
    )
    y_plot = plot(
      mesh,
      [last.(my_z), last.(good_z)],
      linecolor = [RGB(0, 0.5, 0) RGB(0.5, 1, 0)],
      linestyle = [:solid :dash],
      ylims = (-1, 1),
      legend = false
    )
    plot(x_plot, y_plot, layout = (2, 1))
end

function test_peirce_proj(dir = 1, N = 12, use_taylor = true)
    mesh = LinRange(-1, 1, 200)
    my_z = [my_peirce_proj(dir*zeta, N, use_taylor) for zeta in mesh] / my_K(0.5, N, use_taylor)
    good_z = good_peirce_proj.(dir*mesh) / K(0.5)
    x_plot = plot(
      mesh,
      [first.(my_z) first.(good_z)],
      linecolor = [RGB(0.5, 0, 0.5) RGB(1, 0.5, 0.8)],
      linestyle = [:solid :dash],
      ylims = (-1, 1),
      legend = false
    )
    y_plot = plot(
      mesh,
      [last.(my_z), last.(good_z)],
      linecolor = [RGB(0, 0.5, 0) RGB(0.5, 1, 0)],
      linestyle = [:solid :dash],
      ylims = (-1, 1),
      legend = false
    )
    plot(x_plot, y_plot, layout = (2, 1))
end

function test_sqrt()
    mesh = LinRange(-2, 2, 6)
    [my_sqrt(x+im*y) - sqrt(x+im*y) for y in reverse(mesh), x in mesh]
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
