using Combinatorics
using Primes
using DataStructures

function generate_exp(p)
    norm_exps_list = []
    for J in Iterators.flatten([with_replacement_combinations(0:p-1, l) for l in 3:p-1])
        for k in 1:ceil(Int, length(J)/2)
            #println("k = $k")
            exps = []
            for _J in combinations(J, k)
                push!(exps, sum(_J))
            end
            #println("exps = $exps")
            norm_exps = mod.(exps.-minimum(exps), p)
            push!(norm_exps_list, norm_exps)
        end
    end
    return norm_exps_list
end

found_list = Dict{Tuple{Int, Int}, Bool}()
for p in primes(2,14)
    norm_exps_list = generate_exp(p)
    for norm_exps in norm_exps_list        
        counts = counter(norm_exps)
        size = length(counts)
        if(size != p)
            continue
        end 
        
        for q in primes(2,2)
            if(haskey(found_list, (p,q)))
                continue
            end
            complete = true
            for (exp, count) in counts
                if count % q != 1
                    complete = false
                    break
                end                
            end
            if complete
                push!(found_list, (p, q) => true)
                println("q = $q, p = $p, norm_exps = $norm_exps")
                continue
            end
        end
    end
end