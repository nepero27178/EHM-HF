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
gg = [0.25, 0.1, 0.05]

F = Figure(size=(1300,500), figure_padding=1)
axs = [Axis(F[1,i],aspect=1) for i in 1:length(gg)]
linkyaxes!(axs...)

cs = colorschemes[:tabcoolwarm][1:8:end]
vx = 0.449403
vy = 0.212453
for (i,g) in enumerate(gg)
	ax = axs[i]
	ax.xlabel = L"$m$"
	ax.ylabel = L"$u^{(s^*)}$"
	ax.title = L"$g=%$(g)$"
	for (j,FilePathIn) in enumerate(readdir(LAB_ROOT * "/data-g=$(g)"))
		DF = CSV.read(LAB_ROOT * "/data-g=$(g)/" * FilePathIn,DataFrame)
		xx = DF.Im
		yy = DF.IuS
		scatter!(
			ax,xx[1:1],yy[1:1],
			color=cs[j]
		)
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

	scatter!(ax,[vx],[vy],color=:black)
	text!(ax,vx,vy,text=L"$\mathbf{v}$",color=:black,align=(:right,:top))
end

for i in [2,3]
	axs[i].ylabelvisible = false
	axs[i].yticklabelsvisible = false
end

Label(F[0,:], L"Algorithm paths ($t=1.0$, $U=3.0$, $V=6.0$, $L=128$, $\beta=100.0$, $\delta=0.0$)")
FilePathOut = "attractor.pdf"
save(FilePathOut,F)