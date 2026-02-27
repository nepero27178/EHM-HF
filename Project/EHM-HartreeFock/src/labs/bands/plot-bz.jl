using CairoMakie
using Colors
using ColorSchemes
using Contour

# Get color scheme
TabColors = ColorSchemes.tab20 # .tab10
tabblue = TabColors[1]
tabgreen = TabColors[5] # TabColors[3] # tab10 -> 3; tab20 -> 5
tabred = TabColors[7] # TabColors[4] # tab10 -> 4; tab20 -> 7

# Convert to Lab color space
StartLab = Colors.convert(Lab{Float64}, tabblue)
StopLab = Colors.convert(Lab{Float64}, tabred)
MidLab = Colors.convert(Lab{Float64}, RGB(0.8,0.8,0.8))

# Interpolate in Lab space using the range syntax
n = 100  # number of colors in the colormap
Cool = [Colors.convert(RGB, c) for c in range(StartLab, stop=MidLab, length=n)]
Warm = [Colors.convert(RGB, c) for c in range(MidLab, stop=StopLab, length=n)]
cmap = vcat(Cool,Warm)
# cmap = :coolwarm

LAB_ROOT = @__DIR__
FilePathOut = LAB_ROOT * "/tmp.pdf"

# Activate backend
CairoMakie.activate!()

MT = Makie.MathTeXEngine
MT_DIR = dirname(pathof(MT)) * "/../assets/fonts/NewComputerModern"

set_theme!(fonts = (
	regular = MT_DIR * "/NewCM10-Regular.otf",
	bold = MT_DIR * "/NewCM10-Bold.otf"
))

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
axs = [Axis(fig[i, j], aspect=1) for i in 1:2, j in 1:2]
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

ax.title = "MBZ and NBZ"
ax.xlabel = L"$k_x$"
ax.ylabel = L"$k_y$"
ax.xticks = ([-pi, -pi/2, 0, pi/2, pi], [L"$-\pi$", L"$-\pi/2$", L"$0$", L"$\pi/2$", L"$\pi$"])
ax.yticks = ([-pi, -pi/2, 0, pi/2, pi], [L"$-\pi$", L"$-\pi/2$", L"$0$", L"$\pi/2$", L"$\pi$"])

ax = axs[2]
eS = 1.0
ed = 0.5
SS = eS.*StructureFactor.("S",K) + ed.*StructureFactor.("d",K)

heatmap!(ax,Kx,Ky,SS,colormap=cmap)
c = Contour.contours(Kx,Ky,SS)
for cl in Contour.levels(c)
   lvl = Contour.level(cl) # the z-value of this contour level
   for line in Contour.lines(cl)
       xs, ys = coordinates(line) # coordinates of this line segment
       lines!(ax, xs, ys, color="yellow")
   end
end

ax.title = L"$\epsilon^{(s^*)} = 1$, $\epsilon^{(d)} = 1/2$"
ax.xlabel = L"$k_x$"
ax.ylabel = L"$k_y$"
ax.xticks = ([-pi, -pi/2, 0, pi/2, pi], [L"$-\pi$", L"$-\pi/2$", L"$0$", L"$\pi/2$", L"$\pi$"])
ax.yticks = ([-pi, -pi/2, 0, pi/2, pi], [L"$-\pi$", L"$-\pi/2$", L"$0$", L"$\pi/2$", L"$\pi$"])

ax = axs[3]
eS = 1.0
ed = 1.0
SS = eS.*StructureFactor.("S",K) + ed.*StructureFactor.("d",K)

heatmap!(ax,Kx,Ky,SS,colormap=cmap)
c = Contour.contours(Kx,Ky,SS)
for cl in Contour.levels(c)
   lvl = Contour.level(cl) # the z-value of this contour level
   for line in Contour.lines(cl)
       xs, ys = coordinates(line) # coordinates of this line segment
       lines!(ax, xs, ys, color="yellow")
   end
end

ax.title = L"$\epsilon^{(s^*)} = 1$, $\epsilon^{(d)} = 1$"
ax.xlabel = L"$k_x$"
ax.ylabel = L"$k_y$"
ax.xticks = ([-pi, -pi/2, 0, pi/2, pi], [L"$-\pi$", L"$-\pi/2$", L"$0$", L"$\pi/2$", L"$\pi$"])
ax.yticks = ([-pi, -pi/2, 0, pi/2, pi], [L"$-\pi$", L"$-\pi/2$", L"$0$", L"$\pi/2$", L"$\pi$"])

ax = axs[4]
eS = 0.5
ed = 1.0
SS = eS.*StructureFactor.("S",K) + ed.*StructureFactor.("d",K)

heatmap!(ax,Kx,Ky,SS,colormap=cmap)
c = Contour.contours(Kx,Ky,SS)
for cl in Contour.levels(c)
   lvl = Contour.level(cl) # the z-value of this contour level
   for line in Contour.lines(cl)
       xs, ys = coordinates(line) # coordinates of this line segment
       lines!(ax, xs, ys, color="yellow")
   end
end

ax.title = L"$\epsilon^{(s^*)} = 1/2$, $\epsilon^{(d)} = 1$"
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
