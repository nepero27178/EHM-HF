#!/usr/bin/julia

using CairoMakie
using LaTeXStrings
using ColorSchemes

# Includer
PROJECT_ROOT = @__DIR__
PROJECT_ROOT *= "/.."   # Up to the effective root
include(PROJECT_ROOT * "/modules/methods-simulating.jl")

CairoMakie.activate!()

MT = Makie.MathTeXEngine
MT_DIR = dirname(pathof(MT)) * "/../assets/fonts/NewComputerModern"

set_theme!(fonts = (
	regular = MT_DIR * "/NewCM10-Regular.otf",
	bold = MT_DIR * "/NewCM10-Bold.otf"
))


function PlotDensity(
		Phase::String,
		U::Float64,
		Δ::Float64,
		Lx::Int64,
		ββ::Vector{Float64},
		FilePathOut::String;
		cs::Symbol=:winter,
	)

	# Initialize plot
	Fig = Figure(size=(600,400),figure_padding = 1)
	ax = Axis(Fig[1, 1])

	ax.xlabel = L"$\mu/\Delta$"
	ax.ylabel = L"$N(\mu)/2L_xL_y$"
	ax.xticks = [-4,-2,-1,0,1,2,4]
	ax.title = L"%$(Phase) density ($t=1.0$, $U=%$(U)$, $L=%$(Lx)$)"

	Parameters::Dict{String,Float64} = Dict([
		"t" => 1.0,
		"U" => U
	])

	D = 2*Lx^2

	# Reciprocal space discretization (normalized to 1)
	Kx::Vector{Float64} = [kx for kx in -1:2/Lx:1]
	popfirst!(Kx)
	Ky::Vector{Float64} = [ky for ky in -1:2/Lx:1]
	popfirst!(Ky)
	K::Matrix{Vector{Float64}} = [ [kx,ky] for kx in Kx, ky in Ky ]

	v::Dict{String,Float64} = Dict("m" => Δ/U)
	vlines!(
		ax,
		[-Δ,Δ],
		linestyle = :dash,
		color = :gray,
		alpha = 0.5
	)
	text!(
		ax,
		[Point(-Δ-0.25,0.75),Point(Δ+0.25,0.25)],
		text = [L"\mu=-\Delta",L"\mu=\Delta"],
		fontsize = 12,
		color = :gray,
		align = (:center,:center),
		rotation = pi/2
	)

	xx = [x for x in -5:0.01:5]
	q = floor(Int64, length(colorschemes[cs]) / length(ββ) )
	for (j,β) in enumerate(ββ)
		n(μ::Float64) = sum( GetKPopulation(Phase,Parameters,K,v,μ,β) )/D
		yy = n.(xx)
		if β==Inf
			βlabel = "\\beta=\\infty"
		else
			βlabel = "\\beta=$(β)"
		end
		# scatter!(
		# 	ax, xx, yy,
		# 	marker = :circle,
		# 	color = colorschemes[cs][q*j],
		# 	markersize = 4,
		# 	label = L"%$(βlabel)",
		# )
		lines!(
			ax, xx, yy,
			color = colorschemes[cs][q*j],
			label = L"%$(βlabel)",
		)
	end

	xlims!(ax,-5,5)
	ylims!(ax,-0.05,1.05)
	Fig[1, 2] = Legend(Fig, ax, framevisible = false)

	# Save figure
	save(FilePathOut,Fig.scene)
end

if abspath(PROGRAM_FILE) == @__FILE__
	Phase = "AF"
	U = 10.0
	Δ = 1.0
	Lx = 256
	ββ = [Inf, 100.0, 10.0, 1.0, 0.1]
	FilePathOut = "n_$(Phase)_U=$(U)_Δ=$(Δ).pdf"
	PlotDensity(Phase,U,Δ,Lx,ββ,FilePathOut)
end