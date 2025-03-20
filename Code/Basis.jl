using Combinatorics

J = [1,3,4]
k = 2
p = 11

subsets = combinations(J, k)
sums = map(sum, subsets)
println("sums = $sums")
norm_sums = mod.(sums, p)
sort_sums = sort(norm_sums)
println("sort_sums = $sort_sums")