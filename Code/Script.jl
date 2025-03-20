using Pkg
using Nemo
using Combinatorics
using Primes
using DataFrames
using CSV
using JLD2


function ord(a, field)
    if iszero(a)
        return nothing
    end
    group_order = BigInt(characteristic(field))^degree(field) - 1
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

function normalize_polynomial(poly, x, order)
    normalized = sum(coeff(poly, i) * x^(i % order) for i in 0:degree(poly))
    return normalized
end

function parse_polynomial(expr::String, z)
    # Entferne alle Leerzeichen
    expr = replace(expr, r"\s+" => "")
    expr = replace(expr, r"\*" => "")   # Entferne Multiplikationszeichen
    
    # Ersetze Minuszeichen, die als Pluszeichen mit negativem Vorzeichen interpretiert werden
    expr = replace(expr, r"(?<=\d)-" => "-")

    # Definiere ein Dictionary, um Exponenten zu speichern
    terms = Dict{Int, Int}()
    
    # Splitte den Ausdruck an "+" und "-" (um alle Terme zu trennen)
    split_terms = split(expr, r"(?=\+)|(?=-)")
    
    for term in split_terms
        if term == ""  # Überspringe leere Teile, die durch Split entstehen
            continue
        end
        
        if occursin("z^", term)
            # Extrahiere den Exponenten und den Koeffizienten
            parts = split(term, "z^")
            coeffi = parts[1] == "" ? 1 : (parts[1] == "-" ? -1 : (parts[1] == "+" ? 1 : parse(Int, parts[1])))  # Koeffizient vor "z^" (default 1, wenn leer)
            exponent = parse(Int, parts[2])  # Exponent nach "z^"
            terms[exponent] = coeffi
        elseif occursin("z", term)
            parts = split(term, "z")
            # Wenn es nur "z" ist, dann ist der Exponent 1
            coeffi = parts[1] == "" ? 1 : (parts[1] == "-" ? -1 : parse(Int, parts[1]))  # Koeffizient vor "z" (default 1, wenn leer)
            terms[1] = coeffi
        else
            # Konstante (z^0)
            terms[0] = parse(Int, term)
        end
    end

    poly = sum(coeffi * z^exponent for (exponent, coeffi) in terms)
    return poly
end


ZZPoly, z = polynomial_ring(ZZ, "z");
function find_zero_minors(q,p)
    p_root = nothing    
    K, _ = finite_field(q, 1, "y")
    Kx, y = K["y"];
    field, a = finite_field(sum(y^i for i in 0:p-1), "a")
    a_ord = ord(a, field)
    p_root = a
    
    PR, x = polynomial_ring(field, "x")

    max_size = round(Int, (p-1)/2)
    seenI = Set{Vector{Int64}}()
    for i in 3:(BigInt(2)^p - 1)
        mat_size = sum(digits(i, base=2))
        if mat_size > max_size
            continue
        end
        digits_i = digits(i, base=2, pad=p)
        smallestI = maximum([circshift(digits_i, k) for k in 0:(p - 1)])
        if smallestI in seenI
            continue
        end
        push!(seenI, smallestI)
        seenJ = Set{Vector{Int64}}()
        for j in i:(BigInt(2)^p-1)
            if sum(digits(j, base=2)) != mat_size
                continue
            end
            digits_j = digits(j, base=2, pad=p)
            smallestJ = maximum([circshift(digits_j, k) for k in 0:(p - 1)])
            if smallestJ in seenJ
                continue
            end
            push!(seenJ, smallestJ)

            indizesRow = Bool[smallestI...]
            indizesCol = Bool[smallestJ...]
            indizesRow_int = findall(indizesRow).-1
            indizesCol_int = findall(indizesCol).-1
            # println("Checking: ", indizesRow_int, " - ", indizesCol_int)
            minor = determinant(indizesRow_int, indizesCol_int, z)
            minor = PR(minor)
            # println(indizesRow_int, indizesCol_int, minor)
            norm_minor = normalize_polynomial(minor, x, p)
            #println("Minor: ", minor)
            #println("Norm Minor: ", norm_minor)
            eval = norm_minor(p_root)
            if iszero(eval)
                println("Zero minor: ", indizesRow_int, " - ", indizesCol_int, " - ", minor, " - ", norm_minor)
                @lock my_lock push!(results_df, (q, p, p_root, "[" * join(string.(indizesRow_int), " ") * "]", "[" * join(string.(indizesCol_int), " ") * "]"))
                writeCSV()
                return
            end
        end
    end 
    println("No zero minor")   
    @lock my_lock push!(results_df, (q, p, p_root, "", ""))
    writeCSV()
end

det_dict_lock = ReentrantLock();
filename = "det_dict.jls"
if isfile(filename)
    safe_dict = Dict{Tuple{Set{Int}, Set{Int}}, String}()
    @load filename safe_dict
    det_dict = Dict(k => parse_polynomial(v,z) for (k, v) in safe_dict)
else
    det_dict = Dict{Tuple{Set{Int}, Set{Int}}, ZZPolyRingElem}()  # Falls Datei fehlt, leeres Dict verwenden
end
println("Loaded dict: ", length(det_dict), " entries")
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
        @lock det_dict_lock det_dict[key] = sub
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
    @lock det_dict_lock det_dict[key] = det
    # println("Determinant: ", indizesRow, " - ", indizesCol)
    # println("Calced: ", det)
    return det
end


results_df = DataFrame(q=Int[], p=Int[], p_root=Any[], row=String[], col=String[])
function writeCSV()
    println("Writing csv file")
    @lock my_lock CSV.write("data3.csv", results_df, delim="|")
end

function exiting()
    safe_dict = Dict(k => string(v) for (k, v) in det_dict)
    @save "det_dict.jls" safe_dict
    writeCSV()
end
atexit(exiting)

println("Running with ", Threads.nthreads(), " threads!")
my_lock = ReentrantLock();
Threads.@threads for (q,p) in collect(Iterators.product(primes(3,3), primes(2,45)))
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