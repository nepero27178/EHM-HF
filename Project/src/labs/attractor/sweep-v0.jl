using DataFrames
using CSV

LAB_ROOT = @__DIR__
include(LAB_ROOT * "/../../modules/methods-simulating.jl")

Phase = "AF"
g = 0.05
Syms = Set{String}()
ModPars = DataFrame(Dict("t"=>1.0, "U"=>3.0, "V"=>6.0, "L"=>128, "β"=>100.0, "δ"=>0.0))
AlgPars = DataFrame(Dict("p"=>150, "Δv"=>DataFrame(Dict("uS"=>0.001, "m"=>0.001)), "Δn"=>1e-5, "g"=>g))
R0 = GetHFRun(Phase,Syms,ModPars,AlgPars;record=true,RBS=true,RBd=false)

vx = round(only(R0.v.m),digits=1)
vy = round(only(R0.v.uS),digits=1)

xx = round.([x for x in vx-0.2:0.1:vx+0.2],digits=1)
yy = round.([y for y in vy-0.2:0.1:vy+0.2],digits=1)
XY = [(x,y) for x in xx for y in yy]

TT = []
for xy in XY
   v0 = DataFrame("m"=>xy[1], "uS"=>xy[2])
   R = GetHFRun(Phase,Syms,ModPars,AlgPars;record=true,RBS=true,RBd=false,v0)
   push!(TT,R.Track)
   FilePathOut = "data-g=$(g)/vx=$(xy[1])_vy=$(xy[2]).csv"
   mkpath(dirname(FilePathOut))
   CSV.write(FilePathOut,R.Track)
end