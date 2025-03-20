using Combinatorics

p = 11
s = 2
set = [0, 1, 2, 6, 9]
subsets = combinations(set, s)
println(collect(subsets))
sums = map(subset -> sum(subset), subsets)
println(sums)
mod_sums = map(sum -> sum % p, sums)
println(mod_sums)
sorted_sums = sort(mod_sums)
println(sorted_sums)