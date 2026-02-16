using GLMakie
using CairoMakie
using DataFrames
using DelimitedFiles

struct Simulation
	DF::DataFrame
	Setup::String
	Phase::String
	Syms::Vector{String}
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