using DataFrames
using CairoMakie
using CSV

CairoMakie.activate!()
MT = Makie.MathTeXEngine
MT_DIR = dirname(pathof(MT)) * "/../assets/fonts/NewComputerModern"
set_theme!(fonts = (
	regular = MT_DIR * "/NewCM10-Regular.otf",
	bold = MT_DIR * "/NewCM10-Bold.otf"
))

LAB_ROOT = @__DIR__
include(LAB_ROOT * "/../../setup/graphic-setup.jl")

F = Figure(size=(600,600), figure_padding=1)
ax = Axis(F[1,1])
ax.xlabel = L"$m$"
ax.ylabel = L"$u^{(s^*)}$"
ax.title = L"Algorithm paths ($t=1.0$, $U=3.0$, $V=6.0$, $L=128$, $\beta=100.0$, $\delta=0.0$)"

cs = colorschemes[:tabcoolwarm][1:8:end]

for (j,FilePathIn) in enumerate(readdir(LAB_ROOT * "/data"))
	DF = CSV.read(LAB_ROOT * "/data/" * FilePathIn,DataFrame)
	xx = DF.Im
	yy = DF.IuS
	lines!(
		ax,xx,yy,
		color=cs[j],
		linewidth=1
	)
	arrows2d!(
		ax, xx[1:end-1], yy[1:end-1],
		xx[2:end] .- xx[1:end-1],
		yy[2:end] .- yy[1:end-1],
		color=cs[j],
		transparency=0.8,
		tip=Point2f[(0, -0.25), (1.0, 0), (0, 0.25)],
		tipwidth=10,
		shaftwidth=1,
	)
end

vx = 0.449403
vy = 0.212453
scatter!(ax,[vx],[vy],color=:black)
text!(ax,vx,vy,text=L"$\mathbf{v}$",color=:black,align=(:right,:top))

FilePathOut = "attractor.pdf"
save(FilePathOut,F)