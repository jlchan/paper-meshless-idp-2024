using MAT
using SparseArrays
using LinearAlgebra

function get_operators(filename = "operators/matrices_150_500.mat")

    file = matread(filename)
    points = file["points"]
    x, y = points[:,1], points[:,2]
    Fmask = findall(vec(Bool.(points[:,3]) .== false))
    nx, ny = vec(file["normal_x"]), vec(file["normal_y"])
    Ex, Ey = file["E_x"], file["E_y"]
    Qx, Qy = file["Q_x"], file["Q_y"]
    M = Diagonal(vec(file["H_min"]))
    wf = Vector(diag(@. sqrt(Ex^2 + Ey^2))[Fmask])

    return (; x, y, M, Qx, Qy, wf, Fmask, nx, ny)
end

function StartUpMeshfree(filename)
    (; x, y, M, Qx, Qy, wf, Fmask, nx, ny) = get_operators(filename)

    normals = SVector.(nx, ny)

    # nonzeros(Qxy) stores the cols, not the rows. 
    # Qxy is skew symmetric, so we flip the sign 
    # so that Qxy stores the transpose instead for 
    # easy access of rows via sparse CSC format.
    Qxy = -SVector.((Qx - Qx'), (Qy - Qy'))

    # precompute quantities for RHS
    Qxy_norm = norm.(Qxy)
    Qxy_normalized = copy(Qxy)
    @. Qxy_normalized.nzval ./= norm(Qxy_normalized.nzval)

    # check to make sure sparsity patterns are the same 
    @assert all(Qxy_normalized.rowval .== Qxy_norm.rowval)

    return (; x, y, # used for the initial condition
              M, Qxy_norm, Qxy_normalized, wf, Fmask, normals) # used for rhs!
end
