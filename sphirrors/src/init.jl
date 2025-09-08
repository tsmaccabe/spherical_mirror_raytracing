using Pkg

deps = Pkg.dependencies()
packages = [dep.name for dep in values(deps)]

check = ["Images", "LinearAlgebra", "Plots", "Images", "Colors"]
for name in check
    if !(name in packages)
        Pkg.add(name)
    end
end

