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
FilePathOut = LAB_ROOT * "/zones.pdf"

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

fig = Figure(size=(800,700),figure_padding = 1)
axs = [Axis(fig[1, i], aspect=1) for i in 1:2]
linkxaxes!(axs...)
linkyaxes!(axs...)

# Create a colorbar axis on the left, spanning both rows
cbar = Colorbar(fig[1:2, 3], limits=(-1, 1), width=15, colormap=cmap)

ax = axs[1]

band!(ax, Kx, abs.(Kx).-pi, -abs.(Kx).+pi, color=(tabblue, 0.3))
lines!(ax, Kx, abs.(Kx).-pi, color=tabblue)
lines!(ax, Kx, -abs.(Kx).+pi, color=tabblue)

band!(ax, Kx, abs.(Kx), fill(pi,length(Kx)), color=(tabred, 0.3))
band!(ax, Kx, -abs.(Kx), fill(-pi,length(Kx)), color=(tabred, 0.3))
lines!(ax, Kx, Kx, color=tabred)
lines!(ax, Kx, -Kx, color=tabred)

text!(ax,pi*0.75,pi*0.25,text="MBZ",color=tabblue,align=(:right,:top))
text!(ax,pi*0.75,pi*0.75,text="NBZ",color=tabred,align=(:right,:bottom))

ax.xlabel = L"$k_x$"
ax.ylabel = L"$k_y$"
ax.xticks = ([-pi, -pi/2, 0, pi/2, pi], [L"$-\pi$", L"$-\pi/2$", L"$0$", L"$\pi/2$", L"$\pi$"])
ax.yticks = ([-pi, -pi/2, 0, pi/2, pi], [L"$-\pi$", L"$-\pi/2$", L"$0$", L"$\pi/2$", L"$\pi$"])

ax = axs[2]
eS = 1.0
ed = 0.5
SS = eS.*StructureFactor.("S",K) + ed.*StructureFactor.("d",K)

# Hide x,y on top right
axs[1,2].ylabelvisible = false
axs[1,2].yticklabelsvisible = false

save(FilePathOut, fig.scene)
