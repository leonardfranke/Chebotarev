using Nemo

det_dict = Dict{Tuple{Set{Int}, Set{Int}}, ZZPolyRingElem}()
ZZPoly, z = polynomial_ring(ZZ, "z");
function determinant(indizesRow, indizesCol, z)
    key = (Set(indizesRow), Set(indizesCol))
    if haskey(det_dict, key)
        found = det_dict[key]
        # println("Determinant: ", indizesRow, " - ", indizesCol)
        # println("Found: ", found)
        return found
    end

    if size(indizesRow, 1) == 1 && size(indizesCol, 1) == 1
        exp = (indizesRow[1]*indizesCol[1])
        sub = z^exp
        det_dict[key] = sub
        # println("Determinant: ", indizesRow, " - ", indizesCol)
        # println("Calced: ", sub)
        return sub
    end

    rowExpansionIndex = indizesRow[1]
    rowRest = setdiff(indizesRow, rowExpansionIndex)

    det = 0
    factor = 1
    for indexCol in indizesCol
        det += factor * z^(rowExpansionIndex*indexCol) * determinant(rowRest, setdiff(indizesCol, indexCol), z)
        factor = -factor
        # println()
    end
    # println(det)
    det_dict[key] = det
    # println("Determinant: ", indizesRow, " - ", indizesCol)
    # println("Calced: ", det)
    return det
end

q = 2
p = 5
   
K, _ = finite_field(q, 1, "y")
Kx, y = K["y"];
field, a = finite_field(sum(y^i for i in 0:p-1), "a")
PR, x = polynomial_ring(field, "x")

det1 = PR(determinant([0,1,2], [0,1,2], z))
det2 = PR(determinant([0,1,3], [0,1,3], z))
det3 = det2 / det1
println(det1)
println(det2)
println(det3)