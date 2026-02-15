#!/usr/bin/julia
# NOTE: This script is standalone importable and imports all simulations methods.

using LinearAlgebra
using Roots
using Random
LinearAlgebra.BLAS.set_num_threads(Threads.nthreads()) # Parallel optimization
using Integrals
using Elliptic
using DataFrames
using DelimitedFiles

PROJECT_METHODS_DIR = @__DIR__
include(PROJECT_METHODS_DIR * "/methods-physics.jl")
include(PROJECT_METHODS_DIR * "/methods-optimizations.jl")
include(PROJECT_METHODS_DIR * "/methods-IO.jl")

function FindRootμ(
	Phase::String,
	Syms::Vector{String},
	Pars::Dict{String,Float64},
	v::Dict{String,Float64},
	n::Float64;
	Δn::Float64=1e-3,
	μ0::Float64=0.0,
	RBS::Bool=true,
	RBd::Bool=false,
	OptBZ::Bool=true
)::Float64

	if n <= 0 || n >= 1
		@error "Invalid target density. Choose 0 < n < 1." n
		return
	end

	μ::Float64 = 0.0
	if n != 0.5 # Speed up at half-filling

		# Find root of function:
		δn(x::Float64) = GetDensity(Phase,Syms,Pars,v,x,debug,RBS,RBd,OptBZ) - n

		LowerBoundary::Float64 = μ0-0.5
		UpperBoundary::Float64 = μ0+0.5
		if δn(LowerBoundary) > 0
			while δn(LowerBoundary) > 0
				if debug
					@warn "Moving down lower boundary" LowerBoundary
				end
				LowerBoundary -= 1.0
			end
			UpperBoundary = LowerBoundary + 1.0
			μ = find_zero(δn, (LowerBoundary, UpperBoundary))

		elseif δn(LowerBoundary) == 0.0
			μ = LowerBoundary

		elseif δn(UpperBoundary) < 0
			while δn(UpperBoundary) < 0
				if debug
					@warn "Moving up upper boundary" UpperBoundary
				end
				UpperBoundary += 1.0
			end
			LowerBoundary = UpperBoundary - 1.0
			μ = find_zero(δn, (LowerBoundary, UpperBoundary))

		elseif δn(UpperBoundary) == 0.0
			μ = UpperBoundary

		end

	end

	if debug
		n = GetDensity(Phase,Syms,Pars,v,μ,debug,RBS,RBd,OptBZ)
		@info "Optimal chemical potential and density:" μ n
	end

	return μ
end



function HFStep(
	Phase::String,
	Syms::Vector{String},
	Pars::DataFrame,
	v0::DataFrame,
	Δn::Float64=1e-3,
	μ0::Float64=0.0,
	RBS::Bool=true,
	RBd::Bool=false,
	OptBZ::Bool=true,
	debug::Bool=false
)::Tuple{DataFrame,Float64}

	v = copy(v0)	
	μ = FindRootμ(Phase,Syms,Pars,v,n;Δn,μ0,RBS,RBd,OptBZ)

	# Get pars
	t::Float64 = first(Pars.t)
	U::Float64 = first(Pars.U)
	V::Float64 = first(Pars.V)
	L::Int64 = first(Pars.L)
	β::Float64 = first(Pars.β)

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
	reΔΔ, imΔΔ = GetΔΔ(Phase,Syms,Pars,v) # Get real and imaginary dicts
	GetreΔk(k::Vector{Float64})::Float64 = GetFk(reΔΔ,k) # Function
	GetimΔk(k::Vector{Float64})::Float64 = GetFk(imΔΔ,k) # Function
	reΔK::Matrix{Float64} = GetreΔk.(K) # Real gap function
	imΔK::Matrix{Float64} = GetimΔk.(K) # Imaginary gap function

	# Normal phase
	if Phase=="Normal"

		# Self-consistency equations
		RBS ? v.uS .= sum( StructureFactor.("S",K).*FermiDirac.(εK,β,μ) )/LxLy : false
		RBd ? v.ud .= sum( StructureFactor.("d",K).*FermiDirac.(εK,β,μ) )/LxLy : false

	if Phase=="AF-Symmetric"
		EK = sqrt.(εK.^2 .+ reΔK.^2 .+ imΔK.^2)
		m0::Float64 = try
			first(v.m0)
		catch
			@error "Magnetization not found in v0 @ HFStep" v0
			exit()
		end

		eK::Matrix{Float64} = εK./Ek
		replace!(eK, NaN => 0.0) # Null field: <sz>=0

		rK::Matrix{Float64} = reΔK./Ek
		replace!(rK, NaN => 0.0) # Null field: <sy>=0

		iK::Matrix{Float64} = imΔK./Ek
		replace!(iK, NaN => 0.0) # Null field: <sy>=0

		# Self-consistency equations
		RBS ? v.uS .= sum( StructureFactor.("S",K).*eK.*Th("-",EK,β,μ) )/LxLy : false
		RBd ? v.ud .= sum( StructureFactor.("d",K).*eK.*Th("-",EK,β,μ) )/LxLy : false
		v.m .= sum( rK.*Th("-",EK,β,μ) )/(2*LxLy)
		"S" in Syms ? v.vS .= sum( StructureFactor.("S",K).*iK.*Th("-",EK,β,μ) )/(2*LxLy) : false
		"d" in Syms ? v.vd .= sum( StructureFactor.("d",K).*iK.*Th("-",EK,β,μ) )/(2*LxLy) : false

	elseif Phase=="AF-Antisymmetric"
		EK = sqrt.(εK.^2 .+ reΔK.^2 .+ imΔK.^2)
		m0::Float64 = try
			first(v.m0)
		catch
			@error "Magnetization not found in v0 @ HFStep" v0
			exit()
		end

		eK::Matrix{Float64} = εK./Ek
		replace!(eK, NaN => 0.0) # Null field: <sz>=0

		rK::Matrix{Float64} = reΔK./Ek
		replace!(rK, NaN => 0.0) # Null field: <sz>=0

		iK::Matrix{Float64} = imΔK./Ek
		replace!(iK, NaN => 0.0) # Null field: <sz>=0

		# Self-consistency equations
		RBS ? v.uS .= sum( StructureFactor.("S",K).*eK.*Th("-",EK,β,μ) )/LxLy : false
		RBd ? v.ud .= sum( StructureFactor.("d",K).*eK.*Th("-",EK,β,μ) )/LxLy : false
		v.m .= sum( rK.*Th("-",EK,β,μ) )/(2*LxLy)
		"px" in Syms ? v.vpx .= sum( StructureFactor.("px",K).*rK.*Th("-",EK,β,μ) )/(2*LxLy) : false
		"py" in Syms ? v.vpy .= sum( StructureFactor.("py",K).*rK.*Th("-",EK,β,μ) )/(2*LxLy) : false

	elseif Phase=="SC-Singlet"
		ξK::Matrix{Float64} = εK .- μ
		EK = sqrt.(ξK.^2 .+ reΔK.^2)

		tK::Matrix{Float64} = tanh(EK .* β/2) ./ EK
		replace!(tK, NaN => β/2) # @ x~0 : tanh(x)/x~1

		# Self-consistency equations
		RBS ? v.uS .= sum( StructureFactor.("S",K).*(1 .- ξK.*tK) )/LxLy : false
		RBd ? v.ud .= sum( StructureFactor.("d",K).*(1 .- ξK.*tK) )/LxLy : false
		"s" in Syms ? v.wS .= sum( reΔK.*tK )/(2*LxLy) : false
		"S" in Syms ? v.wS .= sum( StructureFactor.("S",K).*reΔK.*tK )/(2*LxLy) : false
		"d" in Syms ? v.wd .= sum( StructureFactor.("d",K).*reΔK.*tK )/(2*LxLy) : false

	elseif Phase=="SC-Triplet"
		ξK::Matrix{Float64} = εK .- μ
		EK = sqrt.(ξK.^2 .+ reΔK.^2)

		tK::Matrix{Float64} = tanh(EK .* β/2) ./ EK
		replace!(tK, NaN => β/2) # @ x~0 : tanh(x)/x~1

		# Self-consistency equations
		RBS ? v.uS .= sum( StructureFactor.("S",K).*(1 .- ξK.*tK) )/LxLy : false
		RBd ? v.ud .= sum( StructureFactor.("d",K).*(1 .- ξK.*tK) )/LxLy : false
		"px" in Syms ? v.px .= sum( StructureFactor.("px",K).*reΔK.*tK )/(2*LxLy) : false
		"py" in Syms ? v.py .= sum( StructureFactor.("py",K).*reΔK.*tK )/(2*LxLy) : false

	return v, μ
end



function HFRun(
	Phase::String,
	Syms::Vector{String},
	ModPars::DataFrame,
	AlgPars::DataFrame,
	v0i::DataFrame,
	Δn::Float64=1e-3,
	μ0::Float64=0.0,
	RBS::Bool=true,
	RBd::Bool=false,
	OptBZ::Bool=true,
	verbose::Bool=false,
	debug::Bool=false,
	record::Bool=false
)::Tuple{Dict{String,Any}, Dict{String,Any}}

	if verbose
		@info "Running HF convergence algorithm" Phase Syms Pars AlgPars
	end

	# Initialize HFPs
	HFPs::Vector{String} = GetHFPs(Phase;Syms)
	v0::DataFrame = DataFrame(Dict(HFPs .=> 0.1))
	# GO ON FROM HERE...

	# Initialize HF dictionaries
	v0::Dict{String,Float64} = Dict([])
	μ::Float64 = 0.0

	if v0i==Dict([])
		for HFP in HFPs
			v0[HFP] = rand()
		end
	elseif issubset(keys(v0i), HFPs)
		for HFP in HFPs
			v0[HFP] = copy(v0i[HFP])
		end
	end
	v = copy(v0) # Shallow copy of values
	Qs = copy(v0) # Copy NaN keys

	# Initialize record matrix
	Record::Dict{String,Vector{Float64}} = Dict([
		key => [ v0[key] ] for key in HFPs
	])

	# Recursive run
	i::Int64 = 1
	I::Int64 = p
	ΔT = @elapsed begin
		while i<=p

			if debug
				printstyled("\n---Step $i---\n", color=:yellow)
			end

			CurrentResults = PerformHFStep(
				Phase,
				Parameters,
				K,v0,n,β;
				Syms,
				RenormalizeBands,
				OptBZ,
				Δn,μ0=μ,
				debug,
			)
			v = copy(CurrentResults[1])
			μ = CurrentResults[2]
			for key in keys(v0)
				current = v[key]
				previous = v0[key]
				tolerance = Δv[key]
				Qs[key] = abs(current-previous) / tolerance
			end

			if all([Qs[key] for key in keys(v0)] .< 1)

				if verbose
					printstyled("\n---Converged at step $i---\n", color=:green)
				end
				I = i
				i = p+1

			elseif any([Qs[key] for key in keys(v0)] .>= 1)
				for key in keys(v0)
					current = v[key]
					previous = v0[key]
					if record
						Record[key] = vcat(Record[key],g*current + (1-g)*previous)
					end
					v[key] = g*current + (1-g)*previous
				end

				if debug
					@info "Initializer and current step after mixing" v0 v
				end
				i += 1

			end
			v0 = copy(v)
		end
	end

	fMFT::Float64 = 0.0
	if all([Qs[key] for key in keys(v0)] .<= 1)

		if verbose
			@info "Algorithm has converged." v Qs
		end
		fMFT = GetFreeEnergy(Phase,Parameters,K,v0,n,μ,β;RenormalizeBands,OptBZ)

	elseif any([Qs[key] for key in keys(v0)] .> 1)

		if verbose
			@info "Algorithm has not converged - v saved as NaN." v Qs Phase
		end

		# Substitute with NaN in order to plot blank points
		fMFT = NaN
		for key in keys(v0)
			v[key] = NaN
		end

	end

	Results::Dict{String,Any} = Dict([
		"HFPs" => v,
		"Record" => Record,
		"ChemicalPotential" => μ,
		"FreeEnergy" => fMFT
	])

	Performance::Dict{String,Any} = Dict([
		"Quality" => Qs,
		"Runtime" => ΔT,
		"Steps" => I
	])

	return Results, Performance
end
