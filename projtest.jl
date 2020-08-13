function cn_coords(zeta)
    x_sq = real(zeta) * real(zeta);
    y_sq = imag(zeta) * imag(zeta);
    r_sq = x_sq + y_sq;
    base = 2*sqrt(1-r_sq) + x_sq - y_sq;
    flip = 2*sqrt(1-r_sq + x_sq*y_sq);
    return (r_sq - flip) / base + im*(base / (r_sq + flip))
end
