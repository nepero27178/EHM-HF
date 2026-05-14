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
CairoMakie.activate!()

FilePathIn = "/home/nepero27178/Thesis/EHM-HF/Project/data/refined/Mode=rs/Setup=B[128]/Phase=Normal/RB=S_Syms=.csv"
Setup,Phase,Syms,RB,Opt,Layer = UnpackFilePath(FilePathIn)
DF = CSV.read(FilePathIn,DataFrame)
Sim = Simulation(DF,Setup,Phase,Syms,RB)
EnlargeDF!(Sim)
U = 0.0
filter!(:U => x -> x==U, DF)

H = Figure(size=(450,400),figure_padding=1)
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
xx,yy,tS = ReshapeData(DF;xVar="δ",yVar="V",zVar="tS")
ax.title = L"$\tilde{t}^{(s^*)}$ ($t=1.0$, $U=%$(U)$, $L=128$, $\beta=100.0$)"
surface!(ax,xx,yy,tS,colormap=colorschemes[:tabcoolerrev],shading=true,colorrange=clims)

# Colorbar(
# 	H[1,2],
# 	colormap=colorschemes[:tabcoolrev],
# 	colorrange=clims,
# )

save("tS(N).pdf", H)