using CairoMakie
using Colors
using ColorSchemes
using Contour

CairoMakie.activate!()
MT = Makie.MathTeXEngine
MT_DIR = dirname(pathof(MT)) * "/../assets/fonts/ComputerModern"
CairoMakie.set_theme!(fonts = (
		regular = MT_DIR * "/cmr10.ttf",
		bold = MT_DIR * "/cmb10.ttf",
		italic = MT_DIR * "/cmmi10.ttf",
		math = MT_DIR * "/cmsy10.ttf",
    ),
	latex_preamble = raw"""
	\usepackage{amsmath}
	\DeclareMathAlphabet{\mathcal}{OMS}{cmsy}{m}{n}
	"""
)

LAB_ROOT = @__DIR__
include(LAB_ROOT * "/../../setup/graphic-setup.jl")
FilePathOut = LAB_ROOT * "/normal-nematic-bands-pp.pdf"

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

fig = Figure(size=(500,400),figure_padding = 1)
ax = Axis(fig[1, 1], aspect=1)

# Create a colorbar axis on the left, spanning both rows
cbar = Colorbar(fig[1,2], limits=(-1, 1), width=15, colormap=cmap)

eS = 1.0
ed = 0.5
SS = eS.*StructureFactor.("S",K) + ed.*StructureFactor.("d",K)

heatmap!(ax,Kx,Ky,SS,colormap=cmap)
lines!(ax, Kx, abs.(Kx).-pi, color="white")
lines!(ax, Kx, -abs.(Kx).+pi, color="white")
lines!(ax, Kx./2, Kx./2, color="white")
lines!(ax, Kx./2, -Kx./2, color="white")

scatter!(ax, [-0.1*pi, -0.8*pi], [0.8*pi, -0.1*pi], color="white")
text!(ax,-0.1*pi,0.8*pi,text=L"$\mathbf{k}$",color="white",align=(:left,:top))
text!(ax,-0.8*pi,-0.1*pi,text=L"$\mathcal{R}\mathbf{k}$",color="white",align=(:left,:top))
text!(ax,0,0.5*pi,text=L"$\mathcal{A}$",color="white",align=(:center,:center))
text!(ax,0,-0.5*pi,text=L"$\mathcal{A}$",color="white",align=(:center,:center))
text!(ax,0.5*pi,0,text=L"$\mathcal{B}$",color="white",align=(:center,:center))
text!(ax,-0.5*pi,0,text=L"$\mathcal{B}$",color="white",align=(:center,:center))

ax.title = L"$\epsilon^{(s^*)} = 1$, $\epsilon^{(d)} = 1/2$"
ax.xlabel = L"$k_x$"
ax.ylabel = L"$k_y$"
ax.xticks = ([-pi, -pi/2, 0, pi/2, pi], [L"$-\pi$", L"$-\pi/2$", L"$0$", L"$\pi/2$", L"$\pi$"])
ax.yticks = ([-pi, -pi/2, 0, pi/2, pi], [L"$-\pi$", L"$-\pi/2$", L"$0$", L"$\pi/2$", L"$\pi$"])

save(FilePathOut, fig.scene)

# ---------

FilePathOut = LAB_ROOT * "/normal-nematic-bands-mp.pdf"

fig = Figure(size=(500,400),figure_padding = 1)
ax = Axis(fig[1, 1], aspect=1)

# Create a colorbar axis on the left, spanning both rows
cbar = Colorbar(fig[1,2], limits=(-1, 1), width=15, colormap=cmap)

eS = -1.0
ed = 0.5
SS = eS.*StructureFactor.("S",K) + ed.*StructureFactor.("d",K)

heatmap!(ax,Kx,Ky,SS,colormap=cmap)
lines!(ax, Kx, abs.(Kx).-pi, color="white")
lines!(ax, Kx, -abs.(Kx).+pi, color="white")
lines!(ax, Kx./2, Kx./2, color="white")
lines!(ax, Kx./2, -Kx./2, color="white")

scatter!(ax, [-0.1*pi, 0.8*pi], [0.8*pi, 0.1*pi], color="white")
text!(ax,-0.1*pi,0.8*pi,text=L"$\mathcal{R}\mathbf{k}$",color="white",align=(:left,:top))
text!(ax,0.8*pi,0.1*pi,text=L"$\mathbf{k}$",color="white",align=(:right,:top))
text!(ax,0,0.5*pi,text=L"$\mathcal{A}$",color="white",align=(:center,:center))
text!(ax,0,-0.5*pi,text=L"$\mathcal{A}$",color="white",align=(:center,:center))
text!(ax,0.5*pi,0,text=L"$\mathcal{B}$",color="white",align=(:center,:center))
text!(ax,-0.5*pi,0,text=L"$\mathcal{B}$",color="white",align=(:center,:center))

ax.title = L"$\epsilon^{(s^*)} = -1$, $\epsilon^{(d)} = 1/2$"
ax.xlabel = L"$k_x$"
ax.ylabel = L"$k_y$"
ax.xticks = ([-pi, -pi/2, 0, pi/2, pi], [L"$-\pi$", L"$-\pi/2$", L"$0$", L"$\pi/2$", L"$\pi$"])
ax.yticks = ([-pi, -pi/2, 0, pi/2, pi], [L"$-\pi$", L"$-\pi/2$", L"$0$", L"$\pi/2$", L"$\pi$"])

save(FilePathOut, fig.scene)