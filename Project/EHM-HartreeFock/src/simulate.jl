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
	g = g0

	# HF iterations
	i::Int64 = 1
	c::Int64 = 1
	for t in tt,
		L in LL,
		δ in δδ,
		β in ββ

		AlgPars::DataFrame = DataFrame(Dict(
			"p" => p,
			"Δv" => Δv,
			"Δn" => Δn,
			"g" => g
		))

		ModPars::DataFrame = DataFrame(Dict(
			"t" => t,
			"U" => 0.0,
			"V" => 0.0,
			"L" => L,
			"β" => β,
			"δ" => δ
		))

		Uc::Float64 = 0.0

		for U in UU
			ModPars.U .= U

			for 	V in VV
				ModPars.V .= V

				Progress = @bold@yellow "[ Progress: $(round(i/I*100,digits=1))% | Convergence rate: $(round(c/i*100,digits=1))% ]"
				Setting = @default@white " Phase=$(Phase)  RB=$(RB...)  Syms=$(Syms...)  t=$t  U=$U  V=$V  L=$L  β=$β  δ=$δ"
				print(Panel(Progress * Setting;style="yellow",title="Run ($(i)/$(I))",title_justify=:right,fit=true))

				# Main run
				R::HFRun = GetHFRun(Phase,Syms,ModPars,AlgPars;v0,RBS,RBd,OptBZ,record)
				Q::DataFrame = DataFrame(Dict(["Q"*x => first(R.Q[!,x]) for x in names(R.Q)]))
				Row::DataFrame = hcat( # Horizontal concatenation of:
					ModPars,R.v,Q, # Already structured dfs
					DataFrame(Dict( # All the rest
						"ΔT" => R.ΔT,
						"I" => R.I,
						"μ" => R.μ,
						"g0" => g0,
						"g" => AlgPars.g,
						"f" => R.f,
						"Converged" => R.Cvd
					))
				)

				R.Cvd ? c += 1 : false
				append::Bool = true
				if i == 1
					!R.Cvd ? c -= 1 : false
					append = false
				end

				# Write on file
				CSV.write(FilePathOut,Row;append)

				# Optimize g for next run
				if Optg && Phase=="SC-Singlet" && "s" in Syms
					Uc = GetUc(ModPars,R.v;how="DisSum",RBS,RBd,OptBZ)
					if U>(2/g0-1)*Uc
						Og = GetOptimalg(U,Uc)
						Og<g0 ? AlgPars.g .= Og : 0
					end
				end

				i += 1
				cl::String = "\r\e[0K\e[1A"
				print(cl*cl*cl) # Carriage return and three lines up to overwrite
			end
		end
	end

	Completed = @bold@green "[ Completed ] "
	Message = "The data have been saved at:"
	print(Panel(Completed * Message * "\n" * FilePathOut;style="green",fit=true))
end

# Main run
function main()

	# Create output directory
	# For simulations: Setup > Phase > Syms (to make comparable data in the same folder)
	# For plots: Phase > Setup > Syms (to make same-phase plots in the same folder)
	DirPathOut = replace(PROJECT_SRC_DIR,"/src"=>"/simulations") * "/Mode=$(Mode)/Setup=$(Setup)/Phase=$(Phase)"
	FilePathOut = DirPathOut * "/RB=$(RB...)_Syms=$(Syms...).csv"
	mkpath(dirname(FilePathOut))

	# Filter out non half-filled simulations from AF phase
	occursin("AF-", Phase) ? filter!(==(0.0),δδ) : false
	Phase=="Normal" ? filter!(==(0.0),UU) : false

	RunStart::DateTime = now()
	TotalRunTime = @elapsed begin
		RunHFScan(
			Phase,Syms,
			tt,UU,VV,
			LL,δδ,ββ,
			p,Δv,Δn,g;
			FilePathOut,
			RBS,RBd,
			OptBZ=false,Optg=true,record=false
		)
	end

	LogPathOut = DirPathOut * "/RB=$(RB...)_Syms=$(Syms...).log"
	Log::DataFrame = hcat(DataFrame(Dict(
		"tt" => [tt],
		"UU" => [UU],
		"VV" => [VV],
		"LL" => [LL],
		"ββ" => [ββ],
		"δδ" => [δδ],
		"p" => p
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
