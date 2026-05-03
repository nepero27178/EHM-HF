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

include("/home/nepero27178/Thesis/EHM-HF/Project/src/setup/graphic-setup.jl")
include("/home/nepero27178/Thesis/EHM-HF/Project/src/modules/methods-3D-plotting.jl")

FilePathIn = "/home/nepero27178/Thesis/EHM-HF/Project/data/refined/Mode=rs/Setup=B[128]/Phase=SC-Singlet/RB=Sd_Syms=Ssd.csv"
Setup,Phase,Syms,RB,Opt,Layer = UnpackFilePath(FilePathIn)
DF_SC = CSV.read(FilePathIn,DataFrame)
Sim = Simulation(DF_SC,Setup,Phase,Syms,RB)
EnlargeDF!(Sim)
U = 4.0
filter!(:U => x -> x==U, DF_SC)

FilePathIn = "/home/nepero27178/Thesis/EHM-HF/Project/data/refined/Mode=rs/Setup=B[128]/Phase=Normal/RB=S_Syms=.csv"
Setup,Phase,Syms,RB,Opt,Layer = UnpackFilePath(FilePathIn)
DF_N = CSV.read(FilePathIn,DataFrame)
Sim = Simulation(DF_N,Setup,Phase,Syms,RB)
EnlargeDF!(Sim)
filter!(:U => x -> x==0.0, DF_N)

H = Figure(size=(900,650),figure_padding=1)
ax = Axis3(
	H[1,1],
	xlabel = L"$\delta$",
	ylabel = L"$V$",
	zlabel = L"$\tilde{t}^{(s^*)}$",
	aspect=(1,1,1),
	azimuth=0.9*pi,
	elevation=pi/10,
	xlabelalign = (:center, :center),
	ylabelalign = (:center, :center),
	zlabelalign = (:center, :center),
	xlabelrotation = 0, # Horizontal xlabel
	ylabelrotation = 0, # Horizontal ylabel
	zlabelrotation = 0, # Horizontal zlabel
)

clims=(0,1)
xx,yy,tS_N = ReshapeData(DF_N;xVar="δ",yVar="V",zVar="tS")
_,_,tS_SC = ReshapeData(DF_SC;xVar="δ",yVar="V",zVar="tS")
ax.title = L"$\tilde{t}^{(s^*)}$ for SC and Normal phases ($t=1.0$, $U=%$(U)$, $L=128$, $\delta=0.0$, $\beta=100.0$)"

surface!(ax,xx,yy,tS_N,colormap=colorschemes[:tabcoolrev],shading=true,colorrange=clims)
# wireframe!(ax,xx,yy,tS_N,color=tabblue,linewidth=0.1)

surface!(ax,xx,yy,tS_SC,colormap=colorschemes[:tabquietrev],shading=true,colorrange=clims)
# wireframe!(ax,xx,yy,tS_SC,color=tabgreen,linewidth=0.1)

Colorbar(
	H[1,2],
	colormap=colorschemes[:tabcoolrev],
	colorrange=clims,
	label="Normal phase",
	labelcolor=tabblue
)
Colorbar(
	H[1,3],
	colormap=colorschemes[:tabquietrev],
	colorrange=clims,
	label="SC phase",
	labelcolor=tabgreen,
)

FilePathOut = "tS(N)+tS(SC).png"
save(FilePathOut, H, px_per_unit=6)