using CairoMakie
using Colors
using ColorSchemes
using Contour

CairoMakie.activate!()
MT = Makie.MathTeXEngine
MT_DIR = dirname(pathof(MT)) * "/../assets/fonts/NewComputerModern"
set_theme!(fonts = (
	regular = MT_DIR * "/NewCM10-Regular.otf",
	bold = MT_DIR * "/NewCM10-Bold.otf"
))

LAB_ROOT = @__DIR__
include(LAB_ROOT * "/../../setup/graphic-setup.jl")
cmap = CoolWarm

FilePathOut = LAB_ROOT * "/normal-shrunk-bands.pdf"

function StructureFactor(
	Sym::String,
	k::Vector{Float64}
)::Float64

	AllSyms = ["s", "S", "x", "y", "d"]
	if !in(Sym, AllSyms)
		@error "Invalid Sym @ StructureFactor" Sym
		return
	end
	kx, ky = k

	if Sym=="s"
		return 1
	elseif Sym=="S"
		return cos(kx) + cos(ky)
	elseif Sym=="x"
		return sqrt(2) * sin(kx)
	elseif Sym=="y"
		return sqrt(2) * sin(ky)
	elseif Sym=="d"
		return cos(kx) - cos(ky)
	end
end

L = [150, 150]
Kx::Vector{Float64} = [kx for kx in -pi:2*pi/L[1]:pi]
Ky::Vector{Float64} = [ky for ky in -pi:2*pi/L[2]:pi]
K::Matrix{Vector{Float64}} = [ [kx,ky] for kx in Kx, ky in Ky ]
SS1::Matrix{Float64} = -2*StructureFactor.("S",K)
SS2::Matrix{Float64} = -1*StructureFactor.("S",K)
SS4::Matrix{Float64} = -1/2*StructureFactor.("S",K)

fig = Figure(size=(800,700),figure_padding = 1)
ax = Axis3(
	fig[1, 1],
	xlabelalign = (:center, :center),
	ylabelalign = (:center, :center),
	zlabelalign = (:center, :center),
	xlabelrotation = 0, # Horizontal xlabel
	ylabelrotation = 0, # Horizontal ylabel
	zlabelrotation = 0, # Horizontal zlabel
	azimuth=0.3*pi
)

surface!(ax,Kx,Ky,SS1,color=tabblue)
surface!(ax,Kx,Ky,SS2,color=tabred)
surface!(ax,Kx,Ky,SS4,color=tabgreen)

ax.xlabel = L"$k_x$"
ax.ylabel = L"$k_y$"
ax.xlabel = L"$\epsilon_\mathbf{k}$"
ax.xticks = ([-pi, -pi/2, 0, pi/2, pi], [L"$-\pi$", L"$-\pi/2$", L"$0$", L"$\pi/2$", L"$\pi$"])
ax.yticks = ([-pi, -pi/2, 0, pi/2, pi], [L"$-\pi$", L"$-\pi/2$", L"$0$", L"$\pi/2$", L"$\pi$"])

# fig[1, 2] = Legend(fig, ax, LegendLabel, framevisible = false)

save(FilePathOut, fig)
