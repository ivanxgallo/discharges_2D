#-------------------- 2D discharges Script -----------------------#

# This is a simulation of discharges based on the Jaime's work
# The idea is to expand the first 1D model of
# Nonlocal self-organization of a dissipative system
# Version: 1.0.0
# Author: Ivan Gallo

using YAML
using Logging
include("system.jl")

# defining logger
logger = SimpleLogger(stdout, Logging.Info)
global_logger(logger)

@info "Initializating simulation parameters..."

parameters_file = "parameters.yml"
params = YAML.load_file(parameters_file)

# setting parameters
T       = params["T"]
M       = params["M"]
N       = params["N"]
q       = params["q"]
Ecx     = Float64(params["Ecx"])
Ecy     = Float64(params["Ecy"])
cloud   = params["cloud"]
k       = params["k"]
save    = params["save"]
x       = params["x"]
y       = params["y"]

@info "parameters are: $params"

for (k, v) in params
    println("$(k): $(v)")
end

# creating fields
fields = Fields(M,N)
# width of cloud
W = Int(cloud*N)
xi = Int((1 - cloud)*N)
xf = N
yi = Int(M/2) - W
yf = Int(M/2) + W


for t in 1:T
    @info "Executing iteration: t = $t"
    # introducing a random amount of charge
    set_points_charge(fields, q, xi, xf, yi, yf, 5)

    # calculating fields and potential
    computing_fields(fields)

    # t = t+1
    update_fields(fields)

    save_fields(fields, "$(t)_0")

    # searching superated umbrals
    superated_umbrals(fields, Ecx, Ecy)

    # iterating until the system is equilibrated
    tau = 0
    check_instability = true
    while check_instability   #(fields)
        @info "Carga total en sistema: $(sum(fields.rho))"
        tau += 1
        @info "Instability cycle: t = $t, tau = $tau"
        equilibrate_charges2(fields, k)
        computing_fields(fields)
        check_instability = check_rho(fields)
        update_fields(fields)
        superated_umbrals(fields, Ecx, Ecy)
        if tau%save == 0
            @info "Saving..."
            save_fields(fields, "$(t)_$tau")
        end
    end
    @info "Saving iteration t = $t"
    save_fields(fields, "$(t)_$tau")
end
@info "Finishing..."





