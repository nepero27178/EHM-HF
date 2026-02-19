#!/usr/bin/julia
# NOTE: This script is standalone importable and imports all simulations methods.

using LinearAlgebra
using Roots
using Random
# using Integrals
# using Elliptic
using DataFrames
using CSV

LinearAlgebra.BLAS.set_num_threads(Threads.nthreads()) # Parallel optimization

PROJECT_METHODS_DIR = @__DIR__
include(PROJECT_METHODS_DIR * "/methods-physics.jl")
include(PROJECT_METHODS_DIR * "/methods-optimizations.jl")
include(PROJECT_METHODS_DIR * "/methods-IO.jl")
include(PROJECT_METHODS_DIR * "/structs.jl")

function FindRootμ(
	Phase::String,
	Syms::Set{String},
	Pars::DataFrame,
	v::DataFrame,
	n::Float64;
	Δn::Float64=1e-3,
	μ0::Float64=0.0,
	RBS::Bool=true,
	RBd::Bool=false,
	OptBZ::Bool=true,
	debug::Bool=false
)::Float64

	if n <= 0 || n >= 1
		@error "Invalid target density. Choose 0 < n < 1." n
		return
	end

	μ::Float64 = 0.0
	if n != 0.5 # Speed up at half-filling

		# Find root of function:
		δn(x::Float64) = GetDensity(Phase,Syms,Pars,v,x;RBS,RBd,OptBZ,debug) - n

		LowerBoundary::Float64 = μ0-0.5
		UpperBoundary::Float64 = μ0+0.5
		if δn(LowerBoundary) >= 0
			while δn(LowerBoundary) > 0
				if debug
					@warn "Moving down lower boundary" LowerBoundary
				end
				LowerBoundary -= 1.0
			end
			UpperBoundary = LowerBoundary + 1.0
		elseif δn(UpperBoundary) <= 0
			while δn(UpperBoundary) < 0
				if debug
					@warn "Moving up upper boundary" UpperBoundary
				end
				UpperBoundary += 1.0
			end
			LowerBoundary = UpperBoundary - 1.0
		end
		μ = find_zero(δn, (LowerBoundary, UpperBoundary))
	end

	if debug
		n = GetDensity(Phase,Syms,Pars,v,μ;RBS,RBd,OptBZ,debug)
		@info "Optimal chemical potential and density:" μ n
	end

	return μ
end

function GetHFStep(
	Phase::String,
	Syms::Set{String},
	Pars::DataFrame,
	v0::DataFrame;
	Δn::Float64=1e-3,
	μ0::Float64=0.0,
	RBS::Bool=true,
	RBd::Bool=false,
	OptBZ::Bool=true,
	debug::Bool=false
)::HFStep

	# Get pars
	t::Float64 = first(Pars.t)
	U::Float64 = first(Pars.U)
	V::Float64 = first(Pars.V)
	L::Int64 = first(Pars.L)
	β::Float64 = first(Pars.β)
	n::Float64 = 0.5+first(Pars.δ)

	# Initialize HFPs and find chemical potential
	v = copy(v0)
	μ = FindRootμ(Phase,Syms,Pars,v0,n;Δn,μ0,RBS,RBd,OptBZ,debug)

	# Get BZ
	K::Matrix{Vector{Float64}}, _, _ = GetK([L, L])
	LxLy::Int64 = L^2

	# Kinetics
	uS0::Float64 = Readv(v0,:uS;Cnd=RBS) # Read s*-wave
	ud0::Float64 = Readv(v0,:ud;Cnd=RBd) # Read d-wave
	εε::Dict{String,Float64} = Dict(
		"S" => -2*t + V*uS0,
		"d" => V*ud0
	)
	Getεk(k::Vector{Float64})::Float64 = GetFk(εε,k) # Function
	εK::Matrix{Float64} = Getεk.(K) # Get bare bands
	EK::Matrix{Float64} = zeros(size(εK)) # Initialize quasibands

	# Gap
	reΔΔ, imΔΔ = GetΔΔ(Phase,Syms,Pars,v0) # Get real and imaginary dicts
	GetreΔk(k::Vector{Float64})::Float64 = GetFk(reΔΔ,k) # Function
	GetimΔk(k::Vector{Float64})::Float64 = GetFk(imΔΔ,k) # Function
	reΔK::Matrix{Float64} = GetreΔk.(K) # Real gap function
	imΔK::Matrix{Float64} = GetimΔk.(K) # Imaginary gap function

	# Normal phase
	if Phase=="Normal"

		# Self-consistency equations
		RBS ? v.uS .= sum( StructureFactor.("S",K).*FermiDirac.(εK,μ,β) )/LxLy : false
		RBd ? v.ud .= sum( StructureFactor.("d",K).*FermiDirac.(εK,μ,β) )/LxLy : false

	# Symmetric AF phase
	elseif Phase=="AF-Symmetric"
		EK = sqrt.(εK.^2 .+ reΔK.^2 .+ imΔK.^2)
		m0::Float64 = try
			first(v0.m)
		catch
			@error "Magnetization not found in v0 @ HFStep" v0
			exit()
		end

		eK::Matrix{Float64} = εK./EK
		replace!(eK, NaN => 0.0) # Null field: <sz>=0

		rK::Matrix{Float64} = reΔK./EK
		replace!(rK, NaN => 0.0) # Null field: <sx>=0

		iK::Matrix{Float64} = imΔK./EK
		replace!(iK, NaN => 0.0) # Null field: <sy>=0

		# Self-consistency equations
		RBS ? v.uS .= -sum( StructureFactor.("S",K).*eK.*Th.("-",EK,μ,β) )/(2*LxLy) : false
		RBd ? v.ud .= -sum( StructureFactor.("d",K).*eK.*Th.("-",EK,μ,β) )/(2*LxLy) : false
		v.m .= sum( rK.*Th.("-",EK,μ,β) )/(2*LxLy)
		"S" in Syms ? v.vS .= sum( StructureFactor.("S",K).*iK.*Th.("-",EK,μ,β) )/(2*LxLy) : false
		"d" in Syms ? v.vd .= sum( StructureFactor.("d",K).*iK.*Th.("-",EK,μ,β) )/(2*LxLy) : false

	# Antisymmetric AF phase
	elseif Phase=="AF-Antisymmetric"
		EK = sqrt.(εK.^2 .+ reΔK.^2 .+ imΔK.^2)
		m0 = try # Pre-assigned data type
			first(v0.m)
		catch
			@error "Magnetization not found in v0 @ HFStep" v0
			exit()
		end

		eK = εK./EK # Pre-assigned data type
		replace!(eK, NaN => 0.0) # Null field: <sz>=0

		rK = reΔK./EK # Pre-assigned data type
		replace!(rK, NaN => 0.0) # Null field: <sx>=0

		iK = imΔK./EK # Pre-assigned data type
		replace!(iK, NaN => 0.0) # Null field: <sy>=0

		# Self-consistency equations
		RBS ? v.uS .= -sum( StructureFactor.("S",K).*eK.*Th.("-",EK,μ,β) )/LxLy : false
		RBd ? v.ud .= -sum( StructureFactor.("d",K).*eK.*Th.("-",EK,μ,β) )/LxLy : false
		v.m .= sum( rK.*Th.("-",EK,μ,β) )/(2*LxLy)
		"x" in Syms ? v.vx .= sum( StructureFactor.("x",K).*rK.*Th.("-",EK,μ,β) )/(2*LxLy) : false
		"y" in Syms ? v.vy .= sum( StructureFactor.("y",K).*rK.*Th.("-",EK,μ,β) )/(2*LxLy) : false

	# Singlet SC phase
	elseif Phase=="SC-Singlet"
		ξK::Matrix{Float64} = εK .- μ
		EK = sqrt.(ξK.^2 .+ reΔK.^2)

		tK::Matrix{Float64} = tanh.(EK .* β/2) ./ EK
		replace!(tK, NaN => β/2) # @ x~0 : tanh(x)/x~1

		# Self-consistency equations
		RBS ? v.uS .= sum( StructureFactor.("S",K).*(1 .- ξK.*tK) )/LxLy : false
		RBd ? v.ud .= sum( StructureFactor.("d",K).*(1 .- ξK.*tK) )/LxLy : false
		"s" in Syms ? v.ws .= sum( reΔK.*tK )/(2*LxLy) : false
		"S" in Syms ? v.wS .= sum( StructureFactor.("S",K).*reΔK.*tK )/(2*LxLy) : false
		"d" in Syms ? v.wd .= sum( StructureFactor.("d",K).*reΔK.*tK )/(2*LxLy) : false

	# Triplet SC phase
	elseif Phase=="SC-Triplet"
		ξK = εK .- μ # Pre-assigned data type
		EK = sqrt.(ξK.^2 .+ reΔK.^2)

		tK = tanh.(EK .* β/2) ./ EK # Pre-assigned data type
		replace!(tK, NaN => β/2) # @ x~0 : tanh(x)/x~1

		# Self-consistency equations
		RBS ? v.uS .= sum( StructureFactor.("S",K).*(1 .- ξK.*tK) )/LxLy : false
		RBd ? v.ud .= sum( StructureFactor.("d",K).*(1 .- ξK.*tK) )/LxLy : false
		"x" in Syms ? v.wx .= sum( StructureFactor.("x",K).*reΔK.*tK )/(2*LxLy) : false
		"y" in Syms ? v.wy .= sum( StructureFactor.("y",K).*reΔK.*tK )/(2*LxLy) : false
	end

	S::HFStep = HFStep(v,μ)
	return S
end



function GetHFRun(
	Phase::String,
	Syms::Set{String},
	ModPars::DataFrame,
	AlgPars::DataFrame;
	v0::DataFrame=DataFrame(),
	μ0::Float64=0.0,
	RBS::Bool=true,
	RBd::Bool=false,
	OptBZ::Bool=true,
	verbose::Bool=false,
	debug::Bool=false,
	record::Bool=false
)::HFRun

	if verbose
		@info "Running HF convergence algorithm" Phase Syms ModPars AlgPars
	end

	# Get model pars
	t::Float64 = first(ModPars.t)
	U::Float64 = first(ModPars.U)
	V::Float64 = first(ModPars.V)
	L::Int64 = first(ModPars.L)
	β::Float64 = first(ModPars.β)

	# Get algorithm pars
	p::Int64 = first(AlgPars.p)
	Δv::DataFrame = first(AlgPars.Δv)
	Δn::Float64 = first(AlgPars.Δn)
	g::Float64 = first(AlgPars.g)

	# Read or build v0
	Ns::Set{String} = Set(names(v0)) # Read names
	HFPs::Set{String} = GetHFPs(Phase,Syms,RBS,RBd) # Get phase HFPs
	if isempty(v0)
		v0::DataFrame = DataFrame(Dict(HFPs .=> 0.1)) # Full initialization
	elseif !isempty(v0)
		if issubset(Ns,HFPs)
			Ns = filter(!in(Ns),HFPs) # Filter missing entries
			v0 = hcat(v0,DataFrame(Dict(Ns .=> 0.1))) # Partial initialization
		elseif !issubset(Set(names(v0)),HFPs)
			@error "Invalid v0 columns" v0
			exit()
		end
	end
	select!(Δv,Cols(in(HFPs))) # Filter tolerances
	Q::DataFrame = DataFrame(Dict(names(v0) .=> 0.0)) # Initialize qualities
	Track::DataFrame =  copy(v0) # Initialize record track

	# Prepare dataframes outside of while loop
	v::DataFrame = copy(v0) # Otherwise chaos with pointers
	w::DataFrame = copy(v0) # Otherwise chaos with pointers

	# Recursive run
	μ::Float64 = 0.0
	i::Int64 = 1
	I::Int64 = p
	Cvd::Bool = false # Cvd="Converged" switch
	ΔT = @elapsed begin
		while true
			debug ? printstyled("\n[ Step $(i) ]\n", color=:yellow) : false
			S = GetHFStep(Phase,Syms,ModPars,v0;Δn,μ0=μ,RBS,RBd,OptBZ,debug) # Single step
			v = copy(S.v) # Otherwise chaos with pointers
			Q .= abs.(v.-v0)./Δv # Get qualities
			Cvd = all(first(Q.<=1)) # Compute converged switch
			record ? Track = vcat(Track,v) : false # Get record
			μ = S.μ

			if Cvd
				I = i
				verbose ? printstyled("\n[ Converged at step $I ]\n", color=:green) : false
				break  # Exit while

			elseif !Cvd
				w .= g.*v .+ (1-g).*v0
				if debug
					Tab::DataFrame = hcat(DataFrame(Dict("Object" => ["Initializer", "Current step", "Mix"])),vcat(v0,v,w))
					@info "Current step dataframes" Tab
				end
				i += 1
				i > p && break # Exit while
			end
			v0 = copy(w) # Otherwise chaos with pointers
		end
	end

	Tab = hcat(DataFrame(Dict("Object" => ["Final HFPs", "Qualities"])),vcat(v,Q)) # Pre-assigned data type
	f::Float64 = 0.0
	if Cvd
		if verbose
			@info "All converged" Tab
		end
		f = GetFreeEnergy(Phase,Syms,ModPars,v,μ;RBS,RBd,OptBZ,debug)

	elseif !Cvd
		if verbose
			@info "Not converged values saved as NaN" Tab
		end
		v = v .+ NaN.*(Q.>1.0) # Filter non-converged values and save as NaN
		f = NaN # Save free-energy as NaN
	end

	!record ? Track = v : false
	R::HFRun = HFRun(HFPs,v,Q,Track,Cvd,μ,f,ΔT,I)
	return R
end
