using Plots
using LinearAlgebra
using Plots
using Images
using Colors

res = 2048

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

d = 1.2
R = 1.0
C = hcat(vec([-R, -R, d]),
         vec([-R,  R, d]),
         vec([ R, -R, d]),
         vec([ R,  R, d]),
         vec([0.0, 0.0, d + R*sqrt(2)]))
Rad = fill(R, size(C, 2)) |> collect  # Vector{Float64}

(Rh, Rh_RGB) = raytrace_area(x_int, y_int, C, Rad, 5000)

plot(Rh_RGB)