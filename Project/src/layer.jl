#!/usr/bin/julia
using CSV
using Dates
using DataFrames
using Term

# Arguments handler
if length(ARGS) != 1
	println("How to use this program?
Type the following: \$ julia simulate.jl --mode
Where:
· mode = hs / rs")
	exit()
else
	UserInput = ARGS
	Mode = UserInput[1][3:end]
end

# Includer
PROJECT_SRC_DIR = @__DIR__
if in(Mode, ["rs", "hs"])
	include(PROJECT_SRC_DIR * "/setup/" * Mode * "-setup.jl")
else
	@error "Invalid argument. Use: mode = hs / rs"
end
include(PROJECT_SRC_DIR * "/modules/structs.jl")
include(PROJECT_SRC_DIR * "/modules/methods-simulating.jl")
include(PROJECT_SRC_DIR * "/modules/methods-physics.jl")
include(PROJECT_SRC_DIR * "/modules/methods-optimizations.jl")
include(PROJECT_SRC_DIR * "/modules/methods-IO.jl")

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
	OptBZ::Bool=true,
	Optg::Bool=true,
	record::Bool=false
)

	# Get Hartree Fock Parameters labels
	HFPs::Set{String} = GetHFPs(Phase,Syms,RBS,RBd)

	# HF iterations
	i::Int64 = 1
	c::Int64 = 1
	for (r,Row) in enumerate(eachrow(DF))

		Progress = @bold@yellow "[ Layering progress: $(round(r/size(DF,1)*100,digits=1))% | Layering convergence rate: $(round(c/i*100,digits=1))% ]"
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

			k::Int64 = 2 # Half mixing, double max (does this make any sense?)
			AlgPars::DataFrame = DataFrame(Dict(
				"p" => p*k, # From setup
				"Δv" => Δv, # From setup
				"Δn" => Δn, # From setup
				"g" => Row.g/k
			))

			# Initializers
			v0 = DataFrame(Row[[HFP for HFP in HFPs]]) # Set{String} => Vector{String}

			# Main run
			R::HFRun = GetHFRun(Phase,Syms,ModPars,AlgPars;v0,RBS,RBd)
			Q::DataFrame = DataFrame(Dict(["Q"*x => first(R.Q[!,x]) for x in names(R.Q)]))
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
			insert!(DF,r,LayeredRow[1,:])

			R.Cvd ? c += 1 : false
			append::Bool = true
			if i == 1
				!R.Cvd ? c -= 1 : false
				append = false
			end

			if FilePathOut != "" #TODO Add no FilePathOut possibility
				# Write on file
				CSV.write(FilePathOut,LayeredRow;append)
			end
		end

		cl::String = "\r\e[0K\e[1A"
		print(cl*cl*cl) # Carriage return and three lines up to overwrite

	end

	Completed = @bold@green "[ Completed ] "
	Message = "The data have been saved at:"
	print(Panel(Completed * Message * "\n" * FilePathOut;style="green",fit=true))


end

# Main run
function main()

	# Create output directory
	# For data: Setup > Phase > Syms (to make comparable data in the same folder)
	# For plots: Phase > Setup > Syms (to make same-phase plots in the same folder)
	DirPathOut::String = dirname(PROJECT_SRC_DIR) * "/data/layered/Mode=$(Mode)/Setup=$(Setup)/Phase=$(Phase)"

	# Read present layers and select the largest
	l = try
		Layers::Vector{Int64} = [U[end] for U in UnpackFilePath.(DirPathOut .* "/" .* readdir(DirPathOut); Layered=true)]
		maximum(Layers)
	catch
		0
	end

	# Generate FilePathIn and FilePathOut
	FilePathIn = nothing
	if l==0
		FilePathIn = replace(DirPathOut, "layered" => "raw") * "/RB=$(RB...)_Syms=$(Syms...).csv"
	elseif l>0
		FilePathIn = DirPathOut * "/RB=$(RB...)_Syms=$(Syms...)_Layer=$(l).csv"
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

	k::Float64 = 2.0
	RunStart::DateTime = now()
	TotalRunTime = @elapsed begin
		LayerHFScan(
			DF,
			Phase,Syms,
			p0,Δv,Δn,k;
			FilePathOut,
			RBS,RBd,
			OptBZ=false,Optg=true,record=false
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
