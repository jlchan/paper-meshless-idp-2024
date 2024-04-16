using Plots
using Trixi
using StaticArrays
using OrdinaryDiffEq

include("operator_setup.jl")

filename = "operators/matrices_150_500.mat"
filename = "operators/matrices_300_1000.mat"
filename = "operators/matrices_600_2000.mat"
# filename = "operators/matrices_1200_4000.mat"

(; x, y, M, Qxy_norm, Qxy_normalized, wf, Fmask, normals)  = StartUpMeshfree(filename)

equations = CompressibleEulerEquations2D(1.4)
function exact_solution(x, y, t, equations::CompressibleEulerEquations2D)
    rho = 1.0 #+ ((abs(x) < 1) * (abs(y) < 1))
    v1, v2 = 0.0, 0.0
    p = 1e-4 + ((abs(x) < 1) * (abs(y) < 1))

    return prim2cons(SVector(rho, v1, v2, p), equations)
end

u0 = exact_solution.(x, y, 0.0, equations)

function rhs!(du, u, p, t)
    (; Qxy_norm, Qxy_normalized, wf, Fmask, normals, 
       invMdiag, uP, du_threaded, 
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

        # wall BCs 
        Threads.@threads :static for i in eachindex(Fmask)        
            rho, v_1, v_2, p = cons2prim(u[Fmask[i]], equations)
            v_n = dot(SVector(v_1, v_2), normals[i])
            q = SVector(rho, v_1 - 2 * v_n, v_2 - 2 * v_n, p)
            uP[i] = prim2cons(q, equations)

            du[Fmask[i]] += wf[i] * numerical_flux(u[Fmask[i]], uP[i], normals[i], equations)
        end

        Threads.@threads :static for i in eachindex(du)
            du[i] *= -invMdiag[i]
        end
    end # @inbounds
end

sound_speed(u, equations::CompressibleEulerEquations2D) = 
    sqrt(equations.gamma * Trixi.pressure(u, equations) / Trixi.density(u, equations))

# function wave_speed_toro(u, equations::CompressibleEulerEquations2D)
# end

parameters = (; Qxy_norm, Qxy_normalized, wf, Fmask, normals, 
                invMdiag = inv.(M.diag), uP = similar(u0[Fmask]), 
                du_threaded = [zero(eltype(u0)) for _ in 1:Threads.nthreads()], 
                numerical_flux=flux_hllc, equations)

tspan = (0, .3)

ode = ODEProblem(rhs!, u0, tspan, parameters)
sol = solve(ode, SSPRK43(), dt=1e-10, abstol=1e-9, reltol=1e-7, 
            saveat=LinRange(tspan..., 25), 
            callback=AliveCallback(alive_interval=100))

w = diag(M)
u = sol.u[end]

id = findall(@. abs(y) < .02)
scatter(x[id], getindex.(u[id], 1), xlims=(0,3), ms=1)
scatter(x[id], Trixi.pressure.(u[id], equations), xlims=(-2,2), ms=1)
# pstar = ((a0L + a0R - (gamma-1)/2 * (uR-uL)) / (a0L * pL^(-(gamma-1)/(2 * gamma)) + a0R * pR^(-(gamma-1)/(2 * gamma))))^(2 * gamma / (gamma - 1))