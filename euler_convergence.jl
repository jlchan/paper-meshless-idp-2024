using Plots
using Trixi
using StaticArrays
using OrdinaryDiffEq

include("operator_setup.jl")

filename = "operators/matrices_75_250.mat"
filename = "operators/matrices_150_500.mat"
filename = "operators/matrices_300_1000.mat"
filename = "operators/matrices_600_2000.mat"
# filename = "operators/matrices_1200_4000.mat"

(; x, y, M, Qxy_norm, Qxy_normalized, wf, Fmask, normals)  = StartUpMeshfree(filename)

equations = CompressibleEulerEquations2D(1.4)
function exact_solution(x, y, t, equations::CompressibleEulerEquations2D)
    v1, v2 = 0.1, 0.2
    rho = 1 + 0.5 * sinpi(1/3 * (x + y - t * (v1 + v2)))
    p = 2.5
    return prim2cons(SVector(rho, v1, v2, p), equations)
end

u0 = exact_solution.(x, y, 0.0, equations)

function rhs!(du, u, p, t)
    (; Qxy_norm, Qxy_normalized, wf, Fmask, normals, 
       invMdiag, uP, du_threaded, inflow, outflow, 
       numerical_flux, equations) = p

    # note: we call "rowvals" here to get the column indices 
    # because Qxy is skew symmetric, so the indices of nonzero 
    # rows are the same as indices of nonzero columns.
    cols = rowvals(Qxy_normalized)
    Qxy_vals = nonzeros(Qxy_normalized)
    Qxy_norm_vals = nonzeros(Qxy_norm)
    @inbounds begin 
        Threads.@threads :static for i in axes(Qxy_norm, 1)
            u_i = u[i]
            du_i = du_threaded[Threads.threadid()]
            for ii in nzrange(Qxy_norm, i) 
                j = cols[ii]
                nij = Qxy_vals[ii] 
                norm_nij = Qxy_norm_vals[ii]

                u_j = u[j]                
                du_i += numerical_flux(u_i, u_j, nij, equations) * norm_nij
            end
            du[i] = du_i
        end

        @. uP = u[Fmask]

        # enforce {u} = exact solution at inflow, do nothing at outflow
        Threads.@threads :static for i in inflow
            u_i = u[Fmask[i]]
            u_exact = exact_solution(x[Fmask[i]], y[Fmask[i]], t, equations)
            uP[i] = 2 * u_exact - u_i
        end
        
        Threads.@threads :static for i in eachindex(Fmask)
            du[Fmask[i]] += wf[i] * numerical_flux(u[Fmask[i]], uP[i], normals[i], equations)
        end    
        Threads.@threads :static for i in eachindex(du)
            du[i] *= -invMdiag[i]
        end
    end # inbounds
end

# problem-dependent parameters
advection_velocity = 0.1, 0.2 # hard-coded for initial_condition_density_wave
a_dot_n = map(n -> dot(advection_velocity, n), normals)
inflow = findall(@. a_dot_n <= 0.0)
outflow = findall(@. a_dot_n > 0.0)

# numerical_flux = flux_lax_friedrichs
numerical_flux = flux_hllc

parameters = (; Qxy_norm, Qxy_normalized, wf, Fmask, normals, 
                invMdiag = inv.(M.diag), uP = similar(u0[Fmask]), 
                du_threaded = [zero(eltype(u0)) for _ in 1:Threads.nthreads()], 
                inflow, outflow, numerical_flux, equations)

tspan = (0, .7)

ode = ODEProblem(rhs!, u0, tspan, parameters)
sol = solve(ode, SSPRK43(), dt = 1e-8, abstol=1e-6, reltol=1e-4, 
            save_everystep = false, 
            callback=AliveCallback(alive_interval=50))

w = diag(M)
u = sol.u[end]
L2_error = sqrt(sum(w .* norm.(u - exact_solution.(x, y, tspan[end], equations)).^2))

