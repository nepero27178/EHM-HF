using DataFrames
using CSV

include("/home/nepero27178/Thesis/EHM-HF/Project/src/modules/structs.jl")
include("/home/nepero27178/Thesis/EHM-HF/Project/src/modules/methods-IO.jl")
include("/home/nepero27178/Thesis/EHM-HF/Project/src/modules/methods-physics.jl")

FilePathIn = "/home/nepero27178/Thesis/EHM-HF/Project/data/refined/Mode=rs/Setup=B[128]/Phase=Normal/RB=S_Syms=.csv"
Setup,Phase,Syms,RB,_,_ = UnpackFilePath(FilePathIn)
DF = CSV.read(FilePathIn,DataFrame)
ff::Vector{Float64} = []
for row in eachrow(DF)
	μ = row.μ
	row = DataFrame(row)
	ModPars = select(row,[:t,:U,:V,:L,:β,:δ])
	v = select(row,[:uS])
	f = GetFreeEnergy(Phase,Syms,ModPars,v,μ;RBS=true,RBd=true)
	push!(ff,f)
end
DF.f .= ff
CSV.write(FilePathIn,DF)