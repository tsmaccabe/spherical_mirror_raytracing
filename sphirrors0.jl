## Packages

using Pkg

# Built-In

using Plots
using LinearAlgebra
using Plots
using Images

# Functions

function cart_to_spherical(r_cart)
    # input: matrix with each column as a cartesian vector
    R = norm.(r_cart, 2)
    [R, arccos.(r_cart[3]./R), arctan.(r_cart[2]./r_cart[1])]
end

function sph_reflect(o, rh, C, Rad)
    # C dim 1: vector elements; C dim 2: sphere number
    # R dim 1: radius of sphere n

    N = length(C)

    r = o # these will be overwritten if there are any valid collisions
    rh_new = rh

    rM = Vector{Float64}(undef, N)
    for i = 1:N
        c = C[i]
        omc = o .- c
        if dot(rh, omc) <= 0
            rM[i] = Inf
            continue
        end
        rh_dot_omc = dot(rh, omc)
        cprh = rh_dot_omc*rh # center of sphere vector projected onto rh

        #omc_sq = dot(omc, omc)

        #D = rh_dot_omc^2 - omc_sq + Rad[i]^2
        #if D <= 0
        #    rM[i] = Inf
        #    continue
        #end

        #rM_p = -rh_dot_omC + sqrt(D)
        #rM_m = -rh_dot_omC - sqrt(D)
        #if rM_m < 0
        #    rM[i] = Inf
        #elseif rM_m >= 0 && rM_p >= 0
        #    rM[i] = rM_m
        #end
    end

    (rM_min, rM_min_i) = findmin(rM)

    if rM_min != Inf && !isnan(rM_min)
        c = C[rM_min_i]
        rad = Rad[rM_min_i]

        r = rM_min .* rh

        n = (r .- c) ./ rad

        rh_new = rh .- 2 * dot(rh, n) .* n
        rh_new = rh_new/norm(rh_new, 2) # for numerical stability
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
        if i > 2
            print(i)
            print(",")
        end
    end

    return rh
end

function raytrace_area(x_int, y_int, C, Rad, reflect_max = 500)
    x_L = length(x_int)
    y_L = length(y_int)

    Rh = zeros(x_L, y_L, 3)

    for i = 1:x_L
        x = x_int[i]
        for j = 1:y_L
            y = y_int[j]

            rh = [x, y, sqrt(1 - x^2 - y^2)]
            Rh[i, j, :] = raytrace([0; 0; 0], rh, C, Rad, reflect_max)
        end
    end

    return Rh
end

# Parameters

res = 256

x_res = res
y_res = res
N = x_res*y_res

a_max = pi/4
b_max = pi/4

x_max = sin(a_max)
y_max = sin(b_max)

dx = 2*x_max/x_res
dy = 2*y_max/y_res

x_int = -x_max + dx/2:dx:x_max - dx/2
y_int = -y_max + dy/2:dy:y_max - dy/2

a_int = asin.(x_int)
b_int = asin.(y_int)

#a_L = length(a_int)
#b_L = length(b_int)

# Initialize

d = 2.0
R = 1
C = [vec([-R -R d]), vec([-R R d]), vec([R -R d]), vec([R R d]),
     vec([-R -R d+2*R]), vec([-R R d+2*R]), vec([R -R d+2*R]), vec([R R d+2*R])]
Rad = vec([R, R, R, R, R, R, R, R, R])

# Perform raytracing

#=
for i = 1:a_L
    a = a_int[i]
    for j = 1:b_L
        b = b_int[j]
        Rh[i, j, :] = raytrace([0; 0; 0], (a, b), C, 1, 500)
    end
end
=#

Rh = raytrace_area(x_int, y_int, C, Rad, 500)

heatmap(x_int, x_int, Rh[:, :, 1],
    color = cgrad([:blue, :white, :orange]), clim = (-1, 1))
