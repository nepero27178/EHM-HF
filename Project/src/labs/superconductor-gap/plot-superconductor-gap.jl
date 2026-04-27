using CairoMakie
using LaTeXStrings
using Colors
using ColorSchemes
using CSV
using DataFrames

LAB_ROOT = @__DIR__
include(LAB_ROOT * "/../../setup/graphic-setup.jl")
include(LAB_ROOT * "/../../modules/methods-3D-plotting.jl")

GLMakie.activate!()
MT = Makie.MathTeXEngine
MT_DIR = dirname(pathof(MT)) * "/../assets/fonts/NewComputerModern"
set_theme!(fonts = (
	regular = MT_DIR * "/NewCM10-Regular.otf",
	bold = MT_DIR * "/NewCM10-Bold.otf"
))

U = 0.0

FilePathIn = "/home/nepero27178/Thesis/EHM-HF/Project/data/refined/Mode=rs/Setup=B[128]/Phase=SC-Singlet/RB=Sd_Syms=Ssd.csv"
DF = CSV.read(FilePathIn,DataFrame)
df = filter(:U => x -> x==U, DF)

xx,yy,zzs = ReshapeData(df;xVar="δ",yVar="V",zVar="ws")
_,_,zzS = ReshapeData(df;xVar="δ",yVar="V",zVar="wS")
_,_,zzd = ReshapeData(df;xVar="δ",yVar="V",zVar="wd")
_,_,zzU = ReshapeData(df;xVar="δ",yVar="V",zVar="U")
_,_,zzV = ReshapeData(df;xVar="δ",yVar="V",zVar="V")

ZZs = -zzs .* zzU
ZZS = zzS .* zzV
ZZd = zzd .* zzV

F = Figure(size=(1300,800),figure_padding=1)
axs = Axis(
	F[1,1],
	xlabel = L"$\delta$",
	ylabel = L"$V$",
	title = L"$s$-wave",
	yaxisposition = :right
)
hs = heatmap!(axs,xx,yy,ZZs,colormap=colorschemes[:tabwarm])
Colorbar(F[1,0],hs,flipaxis=false)
text!(axs,0.45,0.0,text=L"$\tilde{\Delta}^{(s)}$",color=tabred,align=(:left,:bottom),fontsize=20,offset=(-20,10))

axS = Axis(
	F[2,1],
	xlabel = L"$\delta$",
	ylabel = L"$V$",
	title = L"$s^*$-wave",
	yaxisposition = :right
)
hS = heatmap!(axS,xx,yy,ZZS,colormap=colorschemes[:tabcoolrev])
Colorbar(F[2,0],hS,flipaxis=false)
text!(axS,0.45,0.0,text=L"$\tilde{\Delta}^{(s^*)}$",color=tabblue,align=(:left,:bottom),fontsize=20,offset=(-20,10))

axd = Axis(
	F[3,1],
	xlabel = L"$\delta$",
	ylabel = L"$V$",
	title = L"$d$-wave",
	yaxisposition = :right
)
axd.xlabel = L"$\delta$"
axd.ylabel = L"$V$"
axd.title = L"$d$-wave"
hd = heatmap!(axd,xx,yy,ZZd,colormap=colorschemes[:tabquietrev])
Colorbar(F[3,0],hd,flipaxis=false)
text!(axd,0.45,0.0,text=L"$\tilde{\Delta}^{(d)}$",color=tabgreen,align=(:left,:bottom),fontsize=20,offset=(-20,10))

linkxaxes!(axs,axS,axd)
axs.xlabelvisible = false
axs.xticklabelsvisible = false
axS.xlabelvisible = false
axS.xticklabelsvisible = false

ax3 = Axis3(
	F[1:3,2:4],
	aspect=(1,1,1),
	azimuth=0.32*pi,
	elevation=pi/10,
	xlabel=L"$\delta$",
	ylabel=L"$V$",
	zlabel=L"$\tilde{\Delta}^{(\gamma)}$",
	xlabelalign = (:center, :center),
	ylabelalign = (:center, :center),
	zlabelalign = (:center, :center),
	xlabelrotation = 0, # Horizontal xlabel
	ylabelrotation = 0, # Horizontal ylabel
	zlabelrotation = 0, # Horizontal zlabel
)
surface!(ax3,xx,yy,ZZs,colormap=colorschemes[:tabwarm],shading=false)
surface!(ax3,xx,yy,ZZS,colormap=colorschemes[:tabcoolrev],shading=false)
surface!(ax3,xx,yy,ZZd,colormap=colorschemes[:tabquietrev],shading=false)

Label(F[0,:], L"Singlet gap $\tilde{\Delta}_{\mathbf{k}}$ harmonics ($t=1.0$, $U=%$(U)$, $\beta=100.0$)",fontsize=20)

FilePathOut = "gap_U=$(U).png"
save(FilePathOut, F, px_per_unit=6)