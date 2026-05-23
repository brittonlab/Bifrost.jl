include(joinpath(@__DIR__, "juliacall-mcm.jl"))
using .MCMDemo

const OUT = joinpath(@__DIR__, "..", "..", "output", "juliacall-mcm.julia.csv")

# Accept T_sigma as an optional positional CLI arg so the Python driver can
# pass the exact value it used. Default matches MCMDemo.run's own default.
T_sigma = length(ARGS) >= 1 ? parse(Float64, ARGS[1]) : 10.0

rows = MCMDemo.run(threaded = Threads.nthreads() > 1, T_sigma = T_sigma)
MCMDemo.write_table(OUT, rows)
println("Wrote ", OUT)
