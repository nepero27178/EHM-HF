#!/usr/bin/julia

using CairoMakie
using LaTeXStrings
using ColorSchemes

# Includer
PROJECT_ROOT = @__DIR__
PROJECT_ROOT *= "/../.."   # Up to the effective root
include(PROJECT_ROOT * "/modules/methods-simulating.jl")
include(PROJECT_ROOT * "/setup/graphic-setup.jl")

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
		L::Int64,
		ββ::Vector{Float64},
		FilePathOut::String;
		cs::Symbol=:tabwarmcool,
	)

	# Initialize plot
	Fig = Figure(size=(600,400),figure_padding = 1)
	ax = Axis(Fig[1, 1])

	ax.xlabel = L"$\mu/\Delta$"
	ax.ylabel = L"$N(\mu)/2L_xL_y$"
	ax.xticks = [-4,-2,-1,0,1,2,4]
	ax.title = L"%$(Phase) density ($t=1.0$, $U=%$(U)$, $L=%$(Lx)$)"

	v::DataFrame = DataFrame(Dict("m" => Δ/U))
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
	RBS::Bool = false
	RBd::Bool = false
	Syms::Set{String} = Set{String}()
	q = floor(Int64, length(colorschemes[cs]) / length(ββ) )
	for (j,β) in enumerate(ββ)

		Pars::DataFrame = DataFrame(Dict(
			"t" => 1.0,
			"U" => U,
			"V" => 0.0,
			"L" => L,
			"β" => β,
			"δ" => 0.0
		))

		n(μ::Float64) = GetDensity(Phase,Syms,Pars,v,μ;RBS,RBd)
		yy = n.(xx)
		if β==Inf
			βlabel = "\\beta=\\infty"
		else
			βlabel = "\\beta=$(β)"
		end
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
	ββ = reverse([Inf, 100.0, 10.0, 1.0, 0.1])
	FilePathOut = "n_$(Phase)_U=$(U)_Δ=$(Δ).pdf"
	PlotDensity(Phase,U,Δ,Lx,ββ,FilePathOut)
end
