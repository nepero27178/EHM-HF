#!/usr/bin/julia
using CSV
using Dates
using DataFrames
using Term

# Arguments handler
if length(ARGS) != 2
	println("How to use this program?
Type the following: \$ julia simulate.jl --mode --opt
Where:
· mode = hs / rs
· opt = up / mid / down")
	exit()
else
	UserInput = ARGS
	Mode = UserInput[1][3:end]
	Opt = UserInput[2][3:end]
end

# Includer
PROJECT_SRC_DIR = @__DIR__
if in(Mode, ["rs", "hs"])
	include(PROJECT_SRC_DIR * "/setup/" * Mode * "-setup.jl")
else
	@error "Invalid argument. Use: mode = hs / rs"
end

if !in(Opt, ["up", "mid", "down"])
	@error "Invalid argument. Use: opt = up / mid / down"
end

include(PROJECT_SRC_DIR * "/modules/structs.jl")
include(PROJECT_SRC_DIR * "/modules/methods-simulating.jl")
include(PROJECT_SRC_DIR * "/modules/methods-physics.jl")
include(PROJECT_SRC_DIR * "/modules/methods-optimizations.jl")
include(PROJECT_SRC_DIR * "/modules/methods-IO.jl")

function GetNewg(
	g::Float64;
	k::Float64=2.0,
	Opt::String="down"
)::Float64

	if g<=0.0 || g>=1.0
		@error "Invalid g @ GetNewg" g
		return NaN
	end

	if typeof(k) != Float64 || k <= 0
		@error "Invalid k @ GetNewg" k
		return NaN
	end

	if Opt=="down"
		return g/k

	elseif Opt=="mid"
		return 1/2*(1-1/k) + g/k

	elseif Opt=="up"
		return (1-1/k) + g/k

	end
end

function LayerHFScan(
	DF::DataFrame,
	Phase::String,
	Syms::Set{String},
	p::Int64,
	Δv::DataFrame,
	Δn::Float64,
	k::Float64;
	FilePathOut::String="",
	RBS::Bool=true,
	RBd::Bool=true,
	Opt::String="down",
	Maxp::Int64=50,
	OptBZ::Bool=true,
	record::Bool=false
)

	# Get Hartree Fock Parameters labels
	HFPs::Set{String} = GetHFPs(Phase,Syms,RBS,RBd)
	p*k < Maxp ? p = p*k : p=Maxp # Avoid exponential blowup

	# HF iterations
	i::Int64 = 1
	j::Int64 = 1
	c::Int64 = 0
	for (r,Row) in enumerate(eachrow(DF))

		g::Float64 = GetNewg(Row.g;k,Opt)
		Progress = @bold@yellow "[ LP: $(round(r/size(DF,1)*100,digits=1))% | CR: $(round(c/i*100,digits=1))% | Opt: $(Opt) | g: $(Row.g) ➔ $(g) ]"
		Setting = @default@white " Phase=$(Phase)  RB=$(RB...)  Syms=$(Syms...)  t=$(Row.t)  U=$(Row.U)  V=$(Row.V)  L=$(Row.L)  β=$(Row.β)  δ=$(Row.δ)"
		print(Panel(Progress * Setting;style="yellow",title="Row ($(r)/$(size(DF,1)))",title_justify=:right,fit=true))

		if !Row.Converged

			i += 1
			ModPars::DataFrame = DataFrame(Dict(
				"L" => Row.L,
				"U" => Row.U,
				"V" => Row.V,
				"t" => Row.t,
				"β" => Row.β,
				"δ" => Row.δ
			))

			AlgPars::DataFrame = DataFrame(Dict(
				"p" => p, # From setup
				"Δv" => Δv, # From setup
				"Δn" => Δn, # From setup
				"g" => g
			))

			# Initializers (use last converged value)
			LastRow = DataFrame(DF[j,:])
			v0 = select(LastRow, HFPs...) # DataFrame(LastRow[[HFP for HFP in HFPs]]) # Set{String} => Vector{String}

			# Main run
			R::HFRun = GetHFRun(Phase,Syms,ModPars,AlgPars;v0,RBS,RBd)
			Q::DataFrame = DataFrame(Dict(
				[ "Q"*x => first(R.Q[!,x]) for x in names(R.Q) ]
			))
			LayeredRow::DataFrame = hcat( # Horizontal concatenation of:
				ModPars,R.v,Q, # Already structured dfs
				DataFrame(Dict( # All the rest
					"ΔT" => R.ΔT,
					"I" => R.I,
					"μ" => R.μ,
					"g0" => first(Row.g),
					"g" => first(AlgPars.g),
					"f" => R.f,
					"Converged" => R.Cvd
				))
			)

			popat!(DF,r)
			insert!(DF,r,first(LayeredRow))

			R.Cvd ? c += 1 : false
			append::Bool = true
			if r == 1
				!R.Cvd ? c -= 1 : false
				append = false
			end

			if FilePathOut != "" #TODO Add no FilePathOut possibility
				# Write on file (respecting DF order!)
				CSV.write(FilePathOut,DataFrame(DF[r,:]);append)
			end

		elseif Row.Converged

			j = r # Advance to last converged row
			append = true
			if r == 1
				append = false
			end

			if FilePathOut != "" #TODO Add no FilePathOut possibility
				# Write on file
				CSV.write(FilePathOut,DataFrame(Row);append)
			end

		end

		cl::String = "\r\e[0K\e[1A"
		print(cl*cl*cl) # Carriage return and three lines up to overwrite

	end

	Completed = @bold@green "[ Completed | Convergence rate: $(round(c/i*100,digits=1))%   ] "
	Message = "The data have been saved at:"
	print(Panel(Completed * Message * "\n" * FilePathOut;style="green",fit=true))

end

# Main run
function main()

	# Create output directory
	# For data: Setup > Phase > Syms (to make comparable data in the same folder)
	# For plots: Phase > Setup > Syms (to make same-phase plots in the same folder)
	DirPathOut::String = dirname(PROJECT_SRC_DIR) * "/data/layered/Opt=$(Opt)/Mode=$(Mode)/Setup=$(Setup)/Phase=$(Phase)"

	# Read present layers and select the largest
	# l = try
	# 	Layers::Vector{Int64} = [U[end] for U in UnpackFilePath.(DirPathOut .* "/" .* readdir(DirPathOut); Layered=true)]
	# 	maximum(Layers)
	# catch
	# 	0
	# end

	OptsDirPathIn::String = dirname(PROJECT_SRC_DIR) * "/data/layered"
	OptKeys::Vector{String} = [S[2] for S in split.(readdir(OptsDirPathIn),'=')] # All possible layerings
	OptDict::Dict{String,Int64} = Dict([])
	for Key in OptKeys
		DirPathIn = OptsDirPathIn * "/Opt=$(Key)/Mode=$(Mode)/Setup=$(Setup)/Phase=$(Phase)"
		l = try
			Layers::Vector{Int64} = [U[end] for U in UnpackFilePath.(DirPathIn .* "/" .* readdir(DirPathIn); Layered=true)]
			maximum(Layers)
		catch
			0
		end
		OptDict[Key] = l
	end

	OptVals = [v for v in values(OptDict) if v > 0]
	if !allunique(OptVals)
		@error "Corrupted layering. Detected multiple MaxLayers" OptDict
	end
	l, OptIn = findmax(OptDict) # Get last layer

	# Generate FilePathIn and FilePathOut
	FilePathIn = nothing
	if l==0
		FilePathIn = dirname(PROJECT_SRC_DIR) * "/data/raw/Mode=$(Mode)/Setup=$(Setup)/Phase=$(Phase)/RB=$(RB...)_Syms=$(Syms...).csv"
	elseif l>0
		FilePathIn = OptsDirPathIn * "/Opt=$(OptIn)/Mode=$(Mode)/Setup=$(Setup)/Phase=$(Phase)/RB=$(RB...)_Syms=$(Syms...)_Layer=$(l).csv"
	end
	LogPathIn::String = replace(FilePathIn, ".csv" => ".log")
	FilePathOut::String = DirPathOut * "/RB=$(RB...)_Syms=$(Syms...)_Layer=$(l+1).csv"
	mkpath(dirname(FilePathOut))

	# Read previous p in order to increase
	LogIn::DataFrame = CSV.read(LogPathIn,DataFrame)
	p0::Int64 = first(LogIn.p)

	# Extract DataFrame
	DF::DataFrame = CSV.read(FilePathIn, DataFrame)

	# Warn user of memory-heavy simulations detection and give general info
	I::Int64 = length(UU) * length(VV) * length(LL) * length(ββ) * length(δδ)
	C::Int64 = length(Syms) + RBS + RBd + occursin("AF-",Phase)
	c::Float64 = round(sum(DF.Converged)/length(DF.Converged) * 100, digits=1)
	@info "Total iterations, algorithmic complexity (number of HFPs), source layer and source convergence rate" I C l c

	k::Float64 = 4.0
	RunStart::DateTime = now()
	TotalRunTime = @elapsed begin
		LayerHFScan(
			DF,
			Phase,Syms,
			p0,Δv,Δn,k;
			FilePathOut,
			RBS,RBd,
			OptBZ=false,Opt,record=false
		)
	end

	LogPathOut::String = replace(FilePathOut, ".csv" => ".log")
	Log::DataFrame = hcat(DataFrame(Dict(
		"FilePathIn" => FilePathIn,
		"SourceLayer" => l,
		"p0" => p0,
		"p" => k*p0
	)), Δv, DataFrame(Dict(
		"Δn" => Δn,
		"TotalRunTime" => TotalRunTime,
		"Machine" => gethostname(),
		"RunStart" => RunStart,
		"RunStop" => now()
	)))
	CSV.write(LogPathOut,Log)
end

if abspath(PROGRAM_FILE) == @__FILE__
	main()
end
