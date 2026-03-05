using GLMakie
using CairoMakie
using DataFrames
using CSV

struct BrillouinZone
	K::Matrix{Vector{Float64}}
	MC::Vector{Vector{Float64}}
	NC::Vector{Vector{Float64}}
	MB::Vector{Vector{Float64}}
end

struct Simulation
	DF::DataFrame
	Setup::String
	Phase::String
	Syms::Set{String}
	RB::Set{String}
end

struct GroupedPlot
	H::Figure
	DF::SubDataFrame
	FileName::String
end

struct HFStep
	v::DataFrame			# Hartree-Fock vector
	μ::Float64			# Chemical potential
end

struct HFRun
	HFPs::Set{String} 	# Hartree-Fock parameters
	v::DataFrame			# Hartree-Fock vector
	Q::DataFrame			# Convergence quality
	Track::DataFrame		# Tracked evolution
	Cvd::Bool			# Converged flag
	μ::Float64			# Chemical potential
	f::Float64			# Free energy density
	ΔT::Float64			# Runtime
	I::Int64				# Total steps
end