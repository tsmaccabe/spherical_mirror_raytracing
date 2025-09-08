using Plots
using LinearAlgebra
using Plots
using Images
using Colors


# Parameters

res = 2048

# x - y def
#=
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
=#
# a - b undef

a_res = res
b_res = res
N = a_res*b_res

(a0, b0) = (0, 0)
a_max = pi/4
b_max = pi/4

da = 2*a_max/a_res
db = 2*b_max/b_res

a_int = (-a_max + da/2:da:a_max - da/2) .+ a0
b_int = (-b_max + db/2:db:b_max - db/2) .+ b0

x_int = sin.(a_int)
y_int = sin.(b_int)

a_L = length(a_int)
b_L = length(b_int)

# Initialize

d = 1.5
R = 1.0
C = hcat(vec([-R, -R, d]),
         vec([-R,  R, d]),
         vec([ R, -R, d]),
         vec([ R,  R, d]),
         vec([0.0, 0.0, d + R*sqrt(2)]))
Rad = fill(R, size(C, 2)) |> collect  # Vector{Float64}

(Rh, Rh_RGB) = raytrace_area(x_int, y_int, C, Rad, 5000)
plot(Rh_RGB)

# Perform Raytracing

#(Rh, Rh_RGB) = raytrace_area(a_int, b_int, C, R, 5000)
#Rh_RGB = 0.5*(Rh .+ 1)


# Plot Result

#heatmap(x_int, x_int, Rh[:, :, 3],
#    color = cgrad([:blue, :white, :orange]), clim = (-1, 1))
