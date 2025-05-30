#=function cart_to_spherical(r_cart)
    # input: matrix with each column as a cartesian vector
    R = norm.(r_cart, 2)
    [R, arccos.(r_cart[3]./R), arctan.(r_cart[2]./r_cart[1])]
end=#

function sph_reflect(o, rh, C, Rad)
    # C dim 1: vector elements; C dim 2: sphere number
    # R dim 1: radius of sphere n

    N = size(C, 2)

    r = copy(o) # these will be overwritten if there are any valid collisions
    rh_new = copy(rh)

    rM = Vector{Float64}(undef, N)
    for i = 1:N
        c = C[i] .- o # position of sphere relative to o
        R = Rad[i]

        cprM = dot(c, rh) # projection of sphere center onto ray direction

        D = R^2 - norm(c - cprM*rh)^2

        if D > 0 # Ray intersects sphere
            rM[i] = cprM - sqrt(D) # Minus sign --> closest to o
        else
            rM[i] = Inf
            continue
        end

        if rM[i] <= 0
            rM[i] = Inf
        end
    end

    (rM_min, rM_min_i) = findmin(rM)

    if rM_min != Inf
        c_min = C[rM_min_i] # relative to global origin
        R_min = Rad[rM_min_i]

        r = o + rM_min*rh # relative to global origin

        n = (r .- c_min)/R_min
        rh_new = rh - 2*dot(rh, n)*n
    else
        r = copy(o)
        rh_new = copy(rh)
    end

    return (r, rh_new)
    # (o_new, rh_new)
end

function raytrace(o, rh, C, Rad, reflect_max)
    # Om = (a, b) initial ray angle; a is ccw of xz plane, b is above xy plane
    # rh = (x, y, z)

    #(a, b) = Om
    #rh = [sin(a), sin(b), sqrt(1 - (sin(a)^2 + sin(b)^2))]

    #(x, y, z) = rh
    #rh = cart_to_spherical(rh)

    for i = 1:reflect_max
        (o_new, rh_new) = sph_reflect(o, rh, C, Rad)

        if all(o_new .== o)
            break
        else
            o = copy(o_new)
            rh = copy(rh_new)
        end
    end

    return rh
end

function raytrace_area(x_int, y_int, C, Rad, reflect_max = 500)
    x_L = length(x_int)
    y_L = length(y_int)

    Rh = Array{Float64}(undef, x_L, y_L, 3)
    Rh_RGB = Array{RGB{Float64}}(undef, x_L, y_L)
    for i = 1:x_L
        x = x_int[i]
        for j = 1:y_L
            y = y_int[j]

            rh = [sin(x), sin(y), sqrt(1 - sin(x)^2 - sin(y)^2)]
            Rh[i, j, :] = raytrace([0; 0; 0], rh, C, Rad, reflect_max)
            Rh_ij_tuple = ((Rh[i, j, 1], Rh[i, j, 2], Rh[i, j, 3]) .+ 1) ./ 2
            Rh_RGB[i, j] = RGB{Float64}(Rh_ij_tuple[1], Rh_ij_tuple[2], Rh_ij_tuple[3])
        end
    end

    return (Rh, Rh_RGB)
end

function rand_spheres_nonint(N, C_lim, R_lim)
    # C_lim is an array of 3 limit tuples

    x_lim = C_lim[1]
    x_rng = x_lim[2] - x_lim[1]
    y_lim = C_lim[2]
    y_rng = y_lim[2] - y_lim[1]
    z_lim = C_lim[3]
    z_rng = z_lim[2] - z_lim[1]

    C = [x_lim[1] .+ x_rng*rand(Float64, (1, N));
         y_lim[1] .+ y_rng*rand(Float64, (1, N));
         z_lim[1] .+ z_rng*rand(Float64, (1, N))]

    R_rng = R_lim[2] - R_lim[1]

    R = R_lim[1] .+ R_rng*rand(Float64, (1, N))

    for n = 1:N
        for m = 1:N
            if n == m
                continue
            end

            while spheres_intersecting(C[:, n], C[:, m], R[n], R[m])
                C[:, n] = [x_lim[1] .+ x_rng*rand(Float64);
                           y_lim[1] .+ y_rng*rand(Float64);
                           z_lim[1] .+ z_rng*rand(Float64)]

                R[n] = R_lim[1] + R_rng*rand(Float64)
            end
        end
    end

    return (C, R)
end

function spheres_intersecting(C1, C2, R1, R2)
    Cd = C2 - C1
    Rt = R1 + R2
    if any(Cd .< Rt)
        D = norm(C2 - C1, 2)
        if D .< R1 + R2
            return true
        else
            return false
        end
    else
        return false
    end
end
