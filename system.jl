using CSV, DelimitedFiles

# defining class Fields
mutable struct Fields
    M::Int64
    N::Int64
    rho::Matrix{Int64}
    phi::Matrix{Float64}
    E::Matrix{Vector{Float64}}
    rho_t::Matrix{Int64}
    phi_t::Matrix{Float64}
    E_t::Matrix{Vector{Float64}}
    umbrals::Matrix{Vector{Tuple{Int64, Int64}}}
end
# constructor
function Fields(M, N)
    rho     = zeros(Int64, M, N)
    phi     = zeros(Float64, M, N)
    E       = fill([0., 0.], M, N)
    rho_t   = copy(rho)
    phi_t   = copy(phi)
    E_t     = copy(E)
    umbrals = [Vector{Tuple{Int, Int}}() for _ in 1:M, _ in 1:N]

    return Fields(M, N, rho, phi, E, rho_t, phi_t, E_t, umbrals)
end

function set_initial_charge(f::Fields, q::Int64, W::Int64)
    rho_aux = copy(f.rho_t)
    rho_aux[:,f.N-W+1:f.N] = copy(rand(0:q,f.M,W))
    f.rho_t += rho_aux
end

function set_point_charge(f::Fields, q::Int64, x::Int64, y::Int64)
    f.rho_t[y,x] += q
end

function set_points_charge(f::Fields, q::Int64, xi::Int64, xf::Int64, yi::Int64, yf::Int64, N::Int64)
    for n in 1:N
        x = rand(xi:xf)
        y = rand(yi:yf)

        f.rho_t[y,x] += rand(0:q)
    end
end

# remember that coordinates x and y are mapped to i, j
# but in a matrix i is row and j is column indexes (typically)
# meaning that i operate in y and j in x axis
# we will change the above and will say that: matrix[i,j] -> matriz[j,i] (matrix[y,x])
# being y the vertical axis and x the horizontal axis, as usually is.
function computing_fields(f::Fields)
    f.phi_t = zeros(Float64, M, N)
    f.E_t   = fill([0., 0.], M, N)
    for i in 1:f.N
        for j in 1:f.M
            for k in 1:f.N
                for l in 1:f.M
                    if i!=k || j!=l
                        #@info "PASE POR AQUI"
                        f.phi_t[j,i]  += f.rho_t[l,k]*G(i,k,j,l)
                        f.E_t[j,i][1] += f.rho_t[l,k]*DxG(i,k,j,l)
                        f.E_t[j,i][2] += f.rho_t[l,k]*DyG(i,k,j,l)
                    else
                        #@info "i = $i, j = $j"
                        #@info "k = $k, l = $l"
                        f.phi_t[j,i]  += f.rho_t[l,k]*I(i,k,j,l)
                        f.E_t[j,i][1] += f.rho_t[l,k]*IxG(i,k,j,l)
                        f.E_t[j,i][2] += f.rho_t[l,k]*IyG(i,k,j,l)
                    end
                end
            end
        end
    end
end

function update_fields(f::Fields)
    f.rho = copy(f.rho_t)
    f.phi = copy(f.phi_t)
    f.E   = copy(f.E_t)
    f.umbrals = [Vector{Tuple{Int, Int}}() for _ in 1:M, _ in 1:N]
end

function superated_umbrals(f::Fields, Ecx::Float64, Ecy::Float64)
    for i in 1:f.N
        for j in 1:f.M
            if abs(f.E_t[j,i][1]) > Ecx
                if f.E_t[j,i][1] < 0 && i >= 1
                    push!(f.umbrals[j,i], (j,i-1))
                else
                    if i < f.N
                        push!(f.umbrals[j,i], (j,i+1))
                    end
                end
            end
            if abs(f.E_t[j,i][2]) > Ecy
                if f.E_t[j,i][2] < 0 && j > 1
                    push!(f.umbrals[j,i], (j-1,i))
                else
                    if j < f.M
                        push!(f.umbrals[j,i], (j+1,i))
                    end
                end
            end
        end
    end
end

function check_instability(f::Fields)
    for u in f.umbrals
        if !isempty(u)
            return true
        end
    end
    return false
end

function check_rho(f::Fields)
    validation = !isequal(f.rho, f.rho_t)
    return validation
end

function clean_umbrals(f::Fields)
    f.umbrals = fill([], M, N)
end

function equilibrate_charges(f::Fields, k::Float64)
    for i in 1:f.N
        for j in 1:f.M
            # geting adjacent values
            up    = (j < f.M && in((j,i), f.umbrals[j+1,i])) ? trunc(Int64, k*f.rho[j+1,i]) : 0
            down  = (j > 1 && in((j,i), f.umbrals[j-1,i])) ? trunc(Int64, k*f.rho[j-1,i]) : 0
            left  = (i > 1 && in((j,i), f.umbrals[j,i-1])) ? trunc(Int64, k*f.rho[j,i-1]) : 0
            right = (i < f.N && in((j,i), f.umbrals[j,i+1])) ? trunc(Int64, k*f.rho[j,i+1]) : 0

            discharge = 0
            for c in 1:length(f.umbrals[j,i])
                discharge += trunc(Int64, k*f.rho[j,i])
            end

            f.rho_t[j,i] += up + down + left + right - discharge
        end
    end
end

function save_fields(f::Fields, name::String)
    writedlm("density_$name.csv", f.rho)
end




# -------------------------- defining green functions ------------------------ #
f(x,x0,y,y0)   = 1/sqrt((x-x0)^2 + (y-y0)^2)
G(x,x0,y,y0)   = f(x,x0,y,y0) - f(x,-x0,y,y0)
DxG(x,x0,y,y0) = (x-x0)*f(x,x0,y,y0)^3-(x+x0)*f(x,-x0,y,y0)^3
DyG(x,x0,y,y0) = (y-y0)*(f(x,x0,y,y0)^3 + f(x,-x0,y,y0)^3)
I(x,x0,y,y0)   = -f(x,-x0,y,y0)
IxG(x,x0,y,y0) = -(x+x0)*f(x,-x0,y,y0)^3
IyG(x,x0,y,y0) = (y-y0)*f(x,-x0,y,y0)^3

# --------------------------- other auxiliar functions ---------------------------- #
compare_one(x) = max(x, 1)


# --------------------------------- DEAD CODE ------------------------------------ #
#=
function equilibrate_charges(f::Fields, k::Float64, Ecx::Float64, Ecy::Float64)
    for i in 1:f.N
        for j in 1:f.M
            # geting adjacent values
            up    = (j < f.M && in((j,i), f.umbrals[j+1,i])) ? trunc(Int64, k*f.rho[j+1,i]) : 0
            down  = (j > 1 && in((j,i), f.umbrals[j-1,i])) ? trunc(Int64, k*f.rho[j-1,i]) : 0
            left  = (i > 1 && in((j,i), f.umbrals[j,i-1])) ? trunc(Int64, k*f.rho[j,i-1]) : 0
            right = (i < f.N && in((j,i), f.umbrals[j,i+1])) ? trunc(Int64, k*f.rho[j,i+1]) : 0

            # defining boundering conditions
            #if abs(f.E_t[j,i][1]) > abs(f.E_t[j,i][2])
                # x bounder: we dont care that charge leave through x = 0, but we dont want that leave through x = N
            if i < f.N
                discharge_x = (abs(f.E_t[j,i][1]) > Ecx) ? trunc(Int64, k*f.rho[j,i]) : 0
            else
                if abs(f.E_t[j,i][1]) > Ecx
                    discharge_x = (f.E_t[j,i][1] < 0) ? trunc(Int64, k*f.rho[j,i]) : 0
                else
                    discharge_x = 0
                end
            end
                #discharge = discharge_x
            #else
                # y bounder: we set y = 0 and y = M as limits (or walls) in the system
            if j < f.M && j > 1
                discharge_y = (abs(f.E_t[j,i][1]) > Ecy) ? trunc(Int64, k*f.rho[j,i]) : 0
            else
                if abs(f.E_t[j,i][1]) > Ecy
                    if j == f.M
                        discharge_y = (f.E_t[j,i][1] < 0) ? trunc(Int64, k*f.rho[j,i]) : 0
                    else
                        discharge_y = (f.E_t[j,i][1] > 0) ? trunc(Int64, k*f.rho[j,i]) : 0
                    end
                else
                    discharge_y = 0
                end
            end
                #discharge = discharge_y
            #end

            discharge = discharge_x + discharge_y
            f.rho_t[j,i] = f.rho_t[j,i] - discharge + up + down + left + right
        end
    end
end

=#
