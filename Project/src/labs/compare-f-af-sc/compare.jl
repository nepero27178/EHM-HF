using CairoMakie
using Contour

include("/home/nepero27178/Thesis/EHM-HF/Project/src/setup/graphic-setup.jl")
include("/home/nepero27178/Thesis/EHM-HF/Project/src/modules/structs.jl")
include("/home/nepero27178/Thesis/EHM-HF/Project/src/modules/methods-IO.jl")

CairoMakie.activate!()
MT = Makie.MathTeXEngine
MT_DIR = dirname(pathof(MT)) * "/../assets/fonts/NewComputerModern"
set_theme!(fonts = (
	regular = MT_DIR * "/NewCM10-Regular.otf",
	bold = MT_DIR * "/NewCM10-Bold.otf"
))

FilePathIn = "/home/nepero27178/Thesis/EHM-HF/Project/data/refined/Mode=rs/Setup=A[128]/Phase=SC-Singlet/RB=Sd_Syms=Ssd.csv"
Setup,Phase,Syms,RB,Opt,Layer = UnpackFilePath(FilePathIn)
DF_SC = CSV.read(FilePathIn,DataFrame)
Sim = Simulation(DF_SC,Setup,Phase,Syms,RB)
df_SC = filter(:δ => x -> x==0.0, DF_SC)

FilePathIn = "/home/nepero27178/Thesis/EHM-HF/Project/data/refined/Mode=rs/Setup=A[128]/Phase=AF/RB=S_Syms=.csv"
Setup,Phase,Syms,RB,Opt,Layer = UnpackFilePath(FilePathIn)
df_AF = CSV.read(FilePathIn,DataFrame)
Sim = Simulation(df_AF,Setup,Phase,Syms,RB)

H = Figure(size=(750,650),figure_padding=1)
ax = Axis3(
	H[1,1],
	xlabel = L"$U$",
	ylabel = L"$V$",
	zlabel = L"$[ f^{(\mathrm{AF})}-f^{(\mathrm{SC})} ]/t$",
	aspect=(1,1,1),
	azimuth=-0.1*pi,
	elevation=pi/9,
	xlabelalign = (:center, :center),
	ylabelalign = (:center, :center),
	zlabelalign = (:center, :center),
	xlabelrotation = 0, # Horizontal xlabel
	ylabelrotation = 0, # Horizontal ylabel
)

xx,yy,ff_AF = ReshapeData(df_AF;xVar="U",yVar="V",zVar="f")
_,_,ff_SC = ReshapeData(df_SC;xVar="U",yVar="V",zVar="f")
ax.title = L"$[ f^{(\mathrm{AF})}-f^{(\mathrm{SC})} ]/t$ ($t=1.0$, $L=128$, $\delta=0.0$, $\beta=100.0$)"

BlueSteps = 80
GreenSteps = 30
clims = (-4,1.5)
AsymCoolQuiet = vcat(
	[Colors.convert(RGB, c) for c in range(TabBlueLab, stop=WhiteLab, length=BlueSteps)],
	[Colors.convert(RGB, c) for c in range(WhiteLab, stop=TabGreenLab, length=GreenSteps)]
)
cs = AsymCoolQuiet
zz = ff_AF .- ff_SC
zlims!(ax,clims)
s = surface!(ax,xx,yy,zz,colormap=cs,shading=false,colorrange=clims)

lvs =  Contour.levels(Contour.contours(xx,yy,zz,[0]))
xc, yc = Contour.coordinates( Contour.lines(lvs[1])[1] )
zc = zeros(length(xc))

scatterlines!(ax,xc,yc,zc,color=:black,markersize=2,linewidth=0.1)
text!(ax,3.5,4,0.2,text="Phase boundary")
text!(ax,1,7.2,0.7,text="SC")
text!(ax,15,4,-2,text="AF")

Colorbar(H[1,2],s)
FilePathOut = "f(AF)-f(SC).pdf"
save(FilePathOut, H)