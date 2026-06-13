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
include("/home/nepero27178/Thesis/EHM-HF/Project/src/setup/graphic-setup.jl")
include("/home/nepero27178/Thesis/EHM-HF/Project/src/modules/structs.jl")
include("/home/nepero27178/Thesis/EHM-HF/Project/src/modules/methods-physics.jl")
include("/home/nepero27178/Thesis/EHM-HF/Project/src/modules/methods-IO.jl")
CairoMakie.activate!()

FilePathIn = LAB_ROOT * "/../../../data/refined/Mode=rs/Setup=A[128]/Phase=AF/RB=S_Syms=.csv"
Setup,Phase,Syms,RB,Opt,Layer = UnpackFilePath(FilePathIn)
DF_AF = CSV.read(FilePathIn,DataFrame)
Sim = Simulation(DF_AF,Setup,Phase,Syms,RB)
EnlargeDF!(Sim)

FilePathIn = LAB_ROOT * "/../../../data/refined/Mode=rs/Setup=B[128]/Phase=Normal/RB=S_Syms=.csv"
Setup,Phase,Syms,RB,Opt,Layer = UnpackFilePath(FilePathIn)
DF_N = CSV.read(FilePathIn,DataFrame)
filter!([:U,:δ] => (x,y) -> x==0.0 && y==0.0, DF_N)
df = copy(DF_N)
UU = unique(DF_AF.U)[2:end]
for U in UU
	df.U .= fill(U,size(df,1))
	global DF_N = vcat(DF_N,df)
end
Sim = Simulation(DF_N,Setup,Phase,Syms,RB)
EnlargeDF!(Sim)

H = Figure(size=(850,650),figure_padding=1)
ax = Axis3(
	H[1,1],
	xlabel = L"$U$",
	ylabel = L"$V$",
	zlabel = L"$\tilde{t}^{(s^*)}$",
	aspect=(1,1,1),
	azimuth=0.4*pi,
	elevation=pi/10,
	xlabelalign = (:center, :center),
	ylabelalign = (:center, :center),
	zlabelalign = (:center, :center),
	xlabelrotation = 0, # Horizontal xlabel
	ylabelrotation = 0, # Horizontal ylabel
	zlabelrotation = 0, # Horizontal zlabel
)

clims=(0,1)
xx,yy,tS_N = ReshapeData(DF_N;xVar="U",yVar="V",zVar="tS")
_,_,tS_AF = ReshapeData(DF_AF;xVar="U",yVar="V",zVar="tS")
ax.title = L"$\tilde{t}^{(s^*)}$ for AF and Normal phases ($t=1.0$, $L=128$, $\delta=0.0$, $\beta=100.0$)"

surface!(ax,xx,yy,tS_N,colormap=colorschemes[:tabwarm],shading=true,colorrange=clims)
# wireframe!(ax,xx,yy,tS_N,color=tabblue,linewidth=0.1)

surface!(ax,xx,yy,tS_AF,colormap=colorschemes[:tabcoolrev],shading=true,colorrange=clims)
# wireframe!(ax,xx,yy,tS_AF,color=tabgreen,linewidth=0.1)

Colorbar(
	H[1,2],
	colormap=colorschemes[:tabwarm],
	colorrange=clims,
	label="Normal phase",
	labelcolor=tabred
)
Colorbar(
	H[1,3],
	colormap=colorschemes[:tabcoolrev],
	colorrange=clims,
	label="AF phase",
	labelcolor=tabblue,
)

FilePathOut = LAB_ROOT * "/tS(N)+tS(AF).pdf"
save(FilePathOut, H)