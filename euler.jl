using Plots
using Trixi
using StaticArrays
using OrdinaryDiffEq

include("operator_setup.jl")

# filename = "operators/matrices_150_500.mat"
filename = "operators/matrices_300_1000.mat"
# filename = "operators/matrices_sharingan_150_500.mat"
# filename = "operators/matrices_600_2000.mat"
# filename = "operators/matrices_1200_4000.mat"

(; x, y, M, Qxy_norm, Qxy_normalized, wf, Fmask, normals)  = StartUpMeshfree(filename)

equations = CompressibleEulerEquations2D(1.4)
function exact_solution(x, y, t, equations::CompressibleEulerEquations2D)
    rho = .01 + exp(-25 * ((x-2)^2+y^2))
    v1, v2 = 0.0, 0.0
    p = rho^equations.gamma

    return prim2cons(SVector(rho, v1, v2, p), equations)
end

u0 = exact_solution.(x, y, 0.0, equations)

function rhs!(du, u, p, t)
    (; Qxy_norm, Qxy_normalized, 
       wf, Fmask, normals, uP, 
       invMdiag, du_threaded, 
       numerical_flux, equations) = p

    # note: we call "rowvals" here to get the column indices 
    # because Qxy is skew symmetric, so the indices of nonzero 
    # rows are the same as indices of nonzero columns.
    cols = rowvals(Qxy_normalized)
    Qxy_vals = nonzeros(Qxy_normalized)
    Qxy_norm_vals = nonzeros(Qxy_norm)

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
    Threads.@threads for i in eachindex(Fmask)        
        boundary_flux = 
            Trixi.boundary_condition_slip_wall(u[Fmask[i]], normals[i],
                                               SVector(x[Fmask[i]], y[Fmask[i]]), t,
                                               numerical_flux, equations)
        du[Fmask[i]] += wf[i] * boundary_flux        
    end

    Threads.@threads for i in eachindex(du)
        du[i] *= -invMdiag[i]
    end
end

parameters = (; Qxy_norm, Qxy_normalized, wf, Fmask, normals, 
                invMdiag = inv.(M.diag), uP = similar(u0[Fmask]), 
                du_threaded = [zero(eltype(u0)) for _ in 1:Threads.nthreads()], 
                numerical_flux=flux_lax_friedrichs, equations)

tspan = (0, 1.2)

ode = ODEProblem(rhs!, u0, tspan, parameters)
sol = solve(ode, SSPRK43(), dt=1e-7, abstol=1e-6, reltol=1e-4,
            saveat=LinRange(tspan..., 25), 
            callback=AliveCallback(alive_interval=10))

w = diag(M)
rho = getindex.(sol.u[end], 1)
scatter(x, y, zcolor=rho, ratio=1, leg=false, msw=0, ms=1, colorbar=true)
