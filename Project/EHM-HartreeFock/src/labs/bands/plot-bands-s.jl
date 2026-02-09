using CairoMakie

LAB_ROOT = @__DIR__
include(LAB_ROOT * "/../../modules/methods-physics.jl")
FilePathOut = LAB_ROOT * "/plot-structure-factors.pdf"

# Activate backend
CairoMakie.activate!()

MT = Makie.MathTeXEngine
MT_DIR = dirname(pathof(MT)) * "/../assets/fonts/NewComputerModern"

set_theme!(fonts = (
	regular = MT_DIR * "/NewCM10-Regular.otf",
	bold = MT_DIR * "/NewCM10-Bold.otf"
))

L = [150, 150]
Kx::Vector{Float64} = [kx for kx in -pi:2*pi/L[1]:pi]
Ky::Vector{Float64} = [ky for ky in -pi:2*pi/L[2]:pi]
K::Matrix{Vector{Float64}} = [ [kx,ky] for kx in Kx, ky in Ky ]

fig = Figure(size=(800,700),figure_padding = 1)
axs = [Axis(fig[1, 1], aspect=1)]
linkxaxes!(axs...)
linkyaxes!(axs...)

# Create a colorbar axis on the left, spanning both rows
cbar = Colorbar(fig[1:2, 3], limits=(-1, 1), width=15)

ax = axs[1]
ss = -2*StructureFactor.(Sym,K)
surface!(ax,Kx,Ky,SS)
ax.title = L"%$(SymsStr[s])-wave structure factor"
ax.xlabel = L"$k_x$"
ax.ylabel = L"$k_y$"
ax.xticks = ([-pi, -pi/2, 0, pi/2, pi], [L"$-\pi$", L"$-\pi/2$", L"$0$", L"$\pi/2$", L"$\pi$"])
ax.yticks = ([-pi, -pi/2, 0, pi/2, pi], [L"$-\pi$", L"$-\pi/2$", L"$0$", L"$\pi/2$", L"$\pi$"])

# Hide x on top left
axs[1,1].xlabelvisible = false
axs[1,1].xticklabelsvisible = false
# Hide x,y on top right
axs[1,2].xlabelvisible = false
axs[1,2].xticklabelsvisible = false
axs[1,2].ylabelvisible = false
axs[1,2].yticklabelsvisible = false
# Hide y on bottom right
axs[2,2].ylabelvisible = false
axs[2,2].yticklabelsvisible = false

save(FilePathOut, fig.scene)
