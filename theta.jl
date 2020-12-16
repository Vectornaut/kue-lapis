# should be accurate as long as |Im(z)| < Im(tau) by a wide margin
function theta_series(z, tau, N)
    # start accumulating powers of eta
    eta_f = exp(2pi*im*z)
    eta_b = 1 / eta_f
    eta_f_pow = 1
    eta_b_pow = 1
    
    # start accumulating powers of q
    q = exp(pi*im*tau)
    q_sq = q^2
    q_step = q
    q_pow = 1
    
    # calculate terms
    terms = Array{Complex{Float64}}(undef, N)
    terms[1] = 1
    for n in 2:N
        # update q and eta factors
        eta_f_pow *= eta_f
        eta_b_pow *= eta_b
        q_pow *= q_step
        
        # calculate term
        terms[n] = q_pow * (eta_f_pow + eta_b_pow)
        
        # update q step
        q_step *= q_sq
    end
    
    # sum backwards
    sum(reverse(terms))
end

function theta(z, tau, N)
    # z = z_tile*tau + z_off, with imag(z_off) in [-imag(tau)/2, imag(tau)/2]
    z_tile = round(imag(z) / imag(tau))
    tile_base = z_tile*tau
    z_off = z - tile_base
    
    twist = -z_tile*pi*im * (tile_base + 2z_off)
    exp(twist) * theta_series(z_off, tau, N)
end

parameter(tau, N) = 1 - (theta(1/2, tau, N) / theta(0, tau, N))^4
