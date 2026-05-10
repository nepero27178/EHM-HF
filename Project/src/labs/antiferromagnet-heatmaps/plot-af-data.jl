using CairoMakie
using Colors
using ColorSchemes

CairoMakie.activate!()
MT = Makie.MathTeXEngine
MT_DIR = dirname(pathof(MT)) * "/../assets/fonts/NewComputerModern"
set_theme!(fonts = (
	regular = MT_DIR * "/NewCM10-Regular.otf",
	bold = MT_DIR * "/NewCM10-Bold.otf"
))

LAB_ROOT = @__DIR__
include(LAB_ROOT * "/../../setup/graphic-setup.jl")
include(LAB_ROOT * "/../../modules/methods-3D-plotting.jl")

FilePathIn = "/home/nepero27178/Thesis/EHM-HF/Project/data/refined/Mode=rs/Setup=C[128]/Phase=AF/RB=S_Syms=_tailored2.csv"

function PlotSided(
	zVar::String;
	cs=colorschemes[:tabcoolerrev]
)
	FilePathOut = LAB_ROOT * "/af-$(zVar)-data-C128.pdf"

	H = Figure(size=(1300,400),figure_padding=1)
	P = Plot3D(FilePathIn;Print=true,xVar="U",yVar="V",zVar,colorbar=false)
	df = P[1].DF

	# FFS there is no simpler way to do this in Makie!?
	ax_P = P[1].H.content[1]
	ax_H = Axis(H[1,1])
	ax_H.xlabel = ax_P.xlabel.val
	ax_H.ylabel = ax_P.ylabel.val
	ax_H.xticks = ax_P.xticks.val
	ax_H.yticks = ax_P.yticks.val
	ylims!(ax_H,0.0,6.0)

	xx, yy, zz = ReshapeData(DataFrame(df);xVar="U",yVar="V",zVar)
	clims = (
		minimum( vcat(0.0,filter(!isnan,zz)) ),
		maximum( vcat(0.0,filter(!isnan,zz)) )
	)
	h = heatmap!(ax_H,xx,yy,zz,colormap=cs,colorrange=clims)
	text!(ax_H,3.9,5.75,text="High magnetisation",color=:white,align=(:right,:top))
	text!(ax_H,0.1,0.25,text="Low magnetisation",color=:black,align=(:left,:top),rotation=pi/2)

	ax_C = Axis(H[1,2])
	ax_C.xlabel = ax_P.xlabel.val
	ax_C.xticks = ax_P.xticks.val
	ax_C.yticks = ax_P.yticks.val
	ylims!(ax_C,0.0,6.0)

	heatmap!(ax_C,xx,yy,zz,colormap=cs,colorrange=clims)
	contour!(ax_C,xx,yy,zz,colormap=colorschemes[:blacklight],labels=true,levels=10,labelsize=10)
	ax_C.ylabelvisible = false
	ax_C.yticklabelsvisible = false

	linkyaxes!(ax_H,ax_C)
	Label(H[0,:], ax_P.title.val)
	Colorbar(H[1,3],h)
	save(FilePathOut, H)
end

PlotSided("m")
# PlotSided("uS")
