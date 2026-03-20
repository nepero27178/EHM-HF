#!/usr/bin/julia

using CairoMakie
using LaTeXStrings
using ColorSchemes

# Includer
PROJECT_ROOT = @__DIR__
PROJECT_ROOT *= "/../.."   # Up to the effective root
include(PROJECT_ROOT * "/setup/graphic-setup.jl")
include(PROJECT_ROOT * "/modules/structs.jl")
include(PROJECT_ROOT * "/modules/methods-physics.jl")
include(PROJECT_ROOT * "/modules/methods-IO.jl")

CairoMakie.activate!()

MT = Makie.MathTeXEngine
MT_DIR = dirname(pathof(MT)) * "/../assets/fonts/NewComputerModern"

set_theme!(fonts = (
	regular = MT_DIR * "/NewCM10-Regular.otf",
	bold = MT_DIR * "/NewCM10-Bold.otf"
))

Δ::Float64 = 1.0
μ::Float64 = 1.5
β::Float64 = 10.0
L::Vector{Int64} = [50,50]
K, Kx, Ky = GetK(L)
εε::Dict{String,Float64} = Dict(
	"S" => -2
)
Getεk(k::Vector{Float64})::Float64 = GetFk(εε,k) # Function
εK::Matrix{Float64} = Getεk.(K) # Get bare bands
EK::Matrix{Float64} = sqrt.( εK.^2 .+ Δ^2 )

H = Figure(size=(700,400),figure_padding=1)
ax = Axis3(H[1,1],
	azimuth=0.3*pi,
	xlabelrotation = 0, # Horizontal xlabel
	ylabelrotation = 0, # Horizontal ylabel
)
ax.xlabel = L"$k_x$"
ax.ylabel = L"$k_y$"
ax.zlabel = L"$E_\mathbf{k}/t$"
ax.title = L"Doped AF bands ($\mu/\Delta=%$(μ)$, $\beta=%$(β)$)"

ColorFunction = @. FermiDirac(-EK,μ,β)
surface!(ax,Kx,Ky,-EK,color=ColorFunction,colormap=colorschemes[:tabcoolrev],colorrange=(0,1))
ColorFunction = @. FermiDirac(EK,μ,β)
surface!(ax,Kx,Ky,EK,color=ColorFunction,colormap=colorschemes[:tabcoolrev],colorrange=(0,1))

Colorbar(H[1,2],colorrange=(0,1),colormap=colorschemes[:tabcoolrev],label=L"$f(E_\mathbf{k};\beta,\mu)$")
save("antiferromagnet-bands-3d.pdf",H)