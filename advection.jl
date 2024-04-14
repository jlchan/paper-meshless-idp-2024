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

equations = LinearScalarAdvectionEquation2D(1., 0.5)
exact_solution(x, y, t, ::LinearScalarAdvectionEquation2D) = 
    sin(pi/6 * (x - t)) * sin(pi/6 * (y - 0.5 * t))

u0 = exact_solution.(x, y, 0.0, equations)

function rhs!(du, u, p, t)
    (; Qxy_norm, Qxy_normalized, wf, Fmask, normals, 
       invMdiag, uP, du_threaded, inflow, equations) = p

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
            du_i += flux_lax_friedrichs(u_i, u_j, nij, equations) * norm_nij
        end
        du[i] = du_i
    end

    @. uP = u[Fmask]
    uP[inflow] = @. 2 * exact_solution.(x[Fmask[inflow]], y[Fmask[inflow]], t, equations) - u[Fmask[inflow]]
    @. du[Fmask] += wf * flux_lax_friedrichs(u[Fmask], uP, normals, equations)
    @. du *= -invMdiag    
end

# problem-dependent parameters
a_dot_n = map(n -> dot(equations.advection_velocity, n), normals)
inflow = findall(@. a_dot_n <= 0.0)

parameters = (; Qxy_norm, Qxy_normalized, wf, Fmask, normals, 
                invMdiag = inv.(M.diag), uP = similar(u0[Fmask]), 
                du_threaded = [zero(eltype(u0)) for _ in 1:Threads.nthreads()], 
                inflow, equations)

tspan = (0, .7)
ode = ODEProblem(rhs!, u0, tspan, parameters)
sol = solve(ode, SSPRK43(), abstol=1e-6, reltol=1e-3, 
            saveat=LinRange(tspan..., 25), 
            callback=AliveCallback(alive_interval=50))

w = diag(M)
u = sol.u[end]
L2_error = sqrt(sum(w .* (u - exact_solution.(x, y, tspan[end], equations)).^2))
