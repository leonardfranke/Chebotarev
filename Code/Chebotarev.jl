using Nemo
using IterTools
using Combinatorics
using Primes
using DataFrames
using Printf
import DataStructures: SortedDict
using CSV
using StatsBase

function getPolyRing(p, n)
    field = GF(p, n)
    PR, x = polynomial_ring(field, "x")
    return field, PR, x
end

function getSubmatricies(A, max_size)
    n = size(A,1)
    Subs = []
    for size = 1:max_size
        combs = collect(combinations(2:n, size-1))
        len = length(combs)
        for i = 1:len
            for j = i:len
                indizesRow = combs[i]
                indizesCol = combs[j]
                prepend!(indizesRow, 1)
                prepend!(indizesCol, 1) 
                sub = A[indizesRow, indizesCol]
                push!(Subs, (sub, indizesRow.-1, indizesCol.-1))
            end
        end
    end
    return Subs
end

function determinant(A::AbstractMatrix{T}) where T
    n, m = size(A)
    if n != m
        error("Matrix muss quadratisch sein")
    end
    if n == 1
        return A[1, 1]
    elseif n == 2
        return A[1, 1] * A[2, 2] - A[1, 2] * A[2, 1]
    else
        det = 0
        for j in 1:n
            submatrix = Matrix{T}(undef, n-1, n-1)
            for i in 2:n
                subcol = 1
                for k in 1:n
                    if k != j
                        submatrix[i-1, subcol] = A[i, k]
                        subcol += 1
                    end
                end
            end
            cofactor = (-1)^(1 + j) * A[1, j] * determinant(submatrix)
            det += cofactor
        end
        return det
    end
end

function Fourier(x, n, group_exponent = nothing)
    F = Matrix{FqPolyRingElem}(undef, n, n)
    for i in 0:n-1
        for j in 0:n-1
            exp = group_exponent == nothing ? i*j : i*j % group_exponent
            F[i+1, j+1] = x^exp
        end
    end
    return F
end

ZR, _ = ZZ["z"]
function to_integer(elem, p, n)
    val = 0
    for i in 0:n-1
        coef = lift(ZZ, coeff(elem,i))
        val += coef * p^i
    end
    return val
end

function matrix_to_string(mat)
    str = "["
    for i in 1:size(mat, 1)
        for j in 1:size(mat, 2)
            str *= string(mat[i, j])
            if j < size(mat, 2)
                str *= " "  # Leerzeichen zwischen den Elementen einer Zeile
            end
        end
        str *= "; "  # Neue Zeile nach jeder Zeile der Matrix
    end
    return str * "]"
end

function vector_to_string(vec)
    str = "["
    for i in 1:length(vec)
        str *= string(vec[i]) * " "                    
    end
    return str * "]"
end

order_cache = Dict{FqFieldElem, Int}()
function ord(g)
    max_n = 400
    if haskey(order_cache, g)
        return order_cache[g]
    end

    if iszero(g)
        return Inf
    end
    
    n = 1
    h = g
    while !isone(h) && n < max_n  # Solange das Element nicht 1 ist
        h *= g
        n += 1
    end

    if n == max_n
        return Inf
    end
    
    order_cache[g] = n
    return n
end

function is_trivial(root, expRow, expCol)
    order = ord(root)
    if order == 1
        return true
    end
    if order == Inf
        return false
    end
    pairs = map(divisor -> (divisor, order/divisor), Primes.divisors(order))
    for pair in pairs
        k = pair[1]
        l = pair[2]
        remaindersRow = expRow .% k
        remaindersCol = expCol .% l

        for (remainders1, remainders2) in [(remaindersRow, remaindersCol), (remaindersCol, remaindersRow)]
            remaindersMap1 = countmap(remainders1)
            remaindersValues1 = values(remaindersMap1)
            remaindersKeys1 = keys(remaindersMap1)
            remaindersMap2 = countmap(remainders2)
            remaindersValues2 = values(remaindersMap2)
            remaindersKeys2 = keys(remaindersMap2)
    
            leftOver = length(remainders1) - maximum(remaindersValues1)
            surplus = sum(remaindersValues2) - length(remaindersKeys2) - 1
            hasEnoughPairs = surplus - leftOver >= 0
            if hasEnoughPairs
                return true
            end
        end        
    end
    
    return false
end

with_trivials = false
for p in [2]
    for n in primes(14,30)
        F, PR, x = getPolyRing(p, 1)
        
        # Erzeugt die Matrix mit den Einträgen x^{i*j}
        M = Fourier(x, n)
        
        # Berechnet alle Untermatrizen
        Subs = getSubmatricies(M, 4)

        # Sammelt alle Nullstellen der Minoren
        all_roots = []      
        for s in Subs
            minor = determinant(s[1])
            rots = roots(minor)
            ext_rots = product(rots, [s], [minor])
            push!(all_roots, ext_rots...)
        end
        
        dictionary = SortedDict{Int64, Dict{FqFieldElem, Vector{Tuple{Tuple{Matrix{FqPolyRingElem}, Vector{Int64}, Vector{Int64}},FqPolyRingElem}}}}()        
        for (root, matrix, minor) in all_roots

            # if root == 0
            #     @show root, matrix
            # end
            
            if !with_trivials && is_trivial(root, matrix[2], matrix[3])
                continue
            end
            
            mat_size = size(matrix[1], 1)
            if !haskey(dictionary, mat_size)
                dictionary[mat_size] = Dict{FqFieldElem, Vector{Tuple{Tuple{Matrix{FqPolyRingElem}, Vector{Int64}, Vector{Int64}},FqPolyRingElem}}}()
            end
            new_entry = (matrix, minor)
            if !haskey(dictionary[mat_size], root)
                dictionary[mat_size][root] = [new_entry]
            else                
                push!(dictionary[mat_size][root], new_entry)
            end
        end
        
        df = DataFrame()
        data = []
        max_rows = 0
        for (size, dict) in dictionary
            rows = []
            for (root, matricies) in dict
                order = ord(root)
                row_string = "$root ($order) => "
                for matrix in matricies
                    row_string *= "\n$(vector_to_string(matrix[1][2])) $(vector_to_string(matrix[1][3]))"
                end
                push!(rows, row_string)
            end
            if max_rows < length(rows)
                max_rows = length(rows)
            end
            push!(data, (size, rows))
        end
        for (size, rows) in data
            append!(rows, fill("", max_rows - length(rows)))
            df[!, string(size)] = rows
        end
        CSV.write("data/test-$(p)-$(n)-$(with_trivials ? "wt" : "nt").csv", df, delim="|")
        show(df, allcols=true)
    end
end