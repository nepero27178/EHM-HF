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
include(PROJECT_SRC_DIR * "/modules/methods-simulating.jl")
include(PROJECT_SRC_DIR * "/modules/methods-physics.jl")
include(PROJECT_SRC_DIR * "/modules/methods-optimizations.jl")
include(PROJECT_SRC_DIR * "/modules/methods-IO.jl")



function RunHFScan(
	Phase::String,
	Syms::Set{String},
	tt::Vector{Float64},
	UU::Vector{Float64},
	VV::Vector{Float64},
	LL::Vector{Int64},
	δδ::Vector{Float64},
	ββ::Vector{Float64},
	p::Int64,
	Δv::DataFrame,
	Δn::Float64,
	g0::Float64;
	FilePathOut::String="",
	RBS::Bool=true,
	RBd::Bool=true,
	OptBZ::Bool=true,
	Optg::Bool=true,
	record::Bool=false
)

	# Warn user of memory-heavy simulations detection
	I::Int64 = length(UU) * length(VV) * length(LL) * length(ββ) * length(δδ)
	C::Int64 = length(Syms) + RBS + RBd + occursin("AF-",Phase)
	@info "Total iterations and algorithmic complexity" I C

	if FilePathOut == "" && I > 200
		@warn "Large simulation detected and no FilePathOut" I
	end

	# Get Hartree Fock Parameters labels
	HFPs::Set{String} = GetHFPs(Phase,Syms,RBS,RBd)

	# Initializers
	v0::DataFrame = DataFrame(Dict([
		HFP => 0.1 for HFP in HFPs
	]))

	# HF iterations
	i = 1
	for t in tt,
	L in LL,
	δ in δδ,
	β in ββ

		Uc::Float64 = 0.0
		if Optg && occursin("SU", Phase)
			Uc = GetUc(t,[L,L],δ,β)
		end

		g = g0
		for U in UU # TODO Evaluate possibility of computing dynamically g with the bare bands

			if U>(2/g0-1)*Uc && Optg
				Og = GetOptimalg(U,Uc) #TODO phase extension
				Og<g0 ? g=Og : 0
			end

			AlgPars::DataFrame = DataFrame(Dict(
				"p" => p,
				"Δv" => Δv,
				"Δn" => Δn,
				"g" => g
			))

			for V in VV

				ModPars::DataFrame = DataFrame(Dict(
					"t" => t,
					"U" => U,
					"V" => V,
					"L" => L,
					"β" => β,
					"δ" => δ
				))

				Progress = @bold@yellow "[ Progress: $(round(i/I*100,digits=3))% ]"
				Setting = @default@white " @$Phase  t=$t  U=$U  V=$V  L=$L  β=$β  δ=$δ"
				print(Panel(Progress * Setting;style="yellow",title="Run ($(i)/$(I))",title_justify=:right,width=200))
				print("\r\e[3A") # Carriage return and three lines up to overwrite

				# Main run
				R::HFRun = GetHFRun(Phase,Syms,ModPars,AlgPars;v0,RBS,RBd,OptBZ,record)
				Q::DataFrame = DataFrame(Dict(["Q"*x => first(R.Q[!,x]) for x in names(R.Q)]))
				Row::DataFrame = hcat(ModPars,R.v,Q,R.ΔT,R.I,R.μ,g0,g,R.f)
				append = i==1 ? false : true # Header only for first write
				CSV.write(FilePathOut,Row,append)
				i += 1
			end
		end
	end

	Completed = @bold@green "[ Completed ] "
	Location = "The data have been saved at:"
	print(Panel(Completed * Location * "\n" * FilePathOut;style="green",width=200))
end

# Main run
function main()

	# Create output directory
	# For simulations: Setup > Phase > Syms (to make comparable data in the same folder)
	# For plots: Phase > Setup > Syms (to make same-phase plots in the same folder)
	DirPathOut = PROJECT_SRC_DIR * "/../simulations/Mode=$(Mode)/Setup=$(Setup)/Phase=$(Phase)"
	FilePathOut = DirPathOut * "/Syms=$(Syms...).csv"
	mkpath(dirname(FilePathOut))

	# Filter out non half-filled simulations from AF phase
	occursin("AF-", Phase) ? filter!(==(0),δδ) : 0

	TotalRunTime = @elapsed begin
		RunHFScan(
			Phase,Syms,
			tt,UU,VV,
			LL,δδ,ββ,
			p,Δv,Δn,g;
			FilePathOut,
			RBS,RBd,
			OptBZ=false,Optg=false,record=false
		)
	end

	LogPathOut = DirPathOut * "/Syms=$(Syms...).log"
	Log::DataFrame = DataFrame(Dict(
		"tt" => tt,
		"UU" => UU,
		"VV" => VV,
		"LL" => LL,
		"ββ" => ββ,
		"δδ" => δδ,
		"p" => p,
		"Δv" => Δv,
		"Δn" => Δn,
		"TotalRunTime" => TotalRunTime,
		"Machine" => gethostname()
	))
	CSV.write(LogPathOut,Log)
end

if abspath(PROGRAM_FILE) == @__FILE__
	main()
end
