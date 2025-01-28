using Pkg
using Nemo
using Combinatorics
using Primes
using DataFrames
using CSV


function ord(a, field)
    if iszero(a)
        return nothing
    end
    group_order = Int(characteristic(field))^degree(field) - 1
    for d in Primes.divisors(group_order)
        # println(a, " ", d, " ", a^d)
        if a^d == field(1)
            return d
        end
    end
    return group_order  # Im schlimmsten Fall ist die Ordnung gleich der Gruppenordnung
end

function Fourier(x, n, ord)
    F = Matrix{FqPolyRingElem}(undef, n, n)
    for i in 0:n-1
        for j in 0:n-1
            exp = (i*j) % ord
            F[i+1, j+1] = x^exp
        end
    end
    return F
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

function normalize_polynomial(poly, x, order)
    normalized = sum(coeff(poly, i) * x^(i % order) for i in 0:degree(poly))
    return normalized
end



function find_zero_minors(q,p)
    p_root = nothing    
    K, _ = finite_field(q, 1, "x")
    Kx, x = K["x"];
    field, a = finite_field(sum(x^i for i in 0:p-1), "a")
    a_ord = ord(a, field)
    p_root = a
    println("My p-th root: ", p_root)
    
    PR, y = polynomial_ring(field, "y")
    M = Fourier(y, p, p)    
    println(M)

    max_size = round(Int, (p-1)/2)
    seenI = Set{Int}()
    for i in 3:(2^p - 1)
        mat_size = sum(digits(i, base=2))
        if mat_size > max_size
            continue
        end
        smallestI = minimum([((i >> k) | (i << (p - k))) & ((1 << p) - 1) for k in 0:(p - 1)])
        if smallestI in seenI
            continue
        end
        push!(seenI, smallestI)
        seenJ = Set{Int}()
        for j in i:(2^p-1)
            if sum(digits(j, base=2)) != mat_size
                continue
            end
            smallestJ = minimum([((j >> k) | (j << (p - k))) & ((1 << p) - 1) for k in 0:(p - 1)])
            if smallestJ in seenJ
                continue
            end
            push!(seenJ, smallestJ)

            indizesRow = Bool[digits(smallestI, base=2, pad=p)...]
            indizesCol = Bool[digits(smallestJ, base=2, pad=p)...]
            sub = M[indizesRow, indizesCol]            
            minor = determinant(sub)
            norm_minor = normalize_polynomial(minor, y, p)
            #println("Minor: ", minor)
            #println("Norm Minor: ", norm_minor)
            eval = norm_minor(p_root)
            if iszero(eval)
                indizesRow_int = findall(indizesRow).-1
                indizesCol_int = findall(indizesCol).-1
                println("Zero minor: ", indizesRow_int, " - ", indizesCol_int, " - ", minor, " - ", norm_minor)
                push!(results_df, (q, p, p_root, "[" * join(string.(indizesRow_int), " ") * "]", "[" * join(string.(indizesCol_int), " ") * "]"))
                return
            end
        end
    end 
    println("No zero minor")   
    push!(results_df, (q, p, p_root, "", ""))
end


results_df = DataFrame(q=Int[], p=Int[], p_root=Any[], row=String[], col=String[])
function writeCSV()
    println("Writing csv file")
    CSV.write("data.csv", results_df, delim="|")
end
atexit(writeCSV)

for q in primes(2,15)
    for p in primes(2,15)
        base_field = GF(p,1)
        ord_q = ord(base_field(q), base_field)    
        if ord_q != p-1
            println("Wrong order of q=",q," mod p=",p, ":  ", ord_q)
            push!(results_df, (q, p, "Wrong order of p mod q", "", ""))
            continue
        else
            println("Right order of q=",q," mod p=",p, ":  ", ord_q) 
        end
        
        find_zero_minors(q,p)
        println()
    end
end