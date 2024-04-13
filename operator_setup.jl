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