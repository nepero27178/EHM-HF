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
UU = unique(DF_SC.U)

FilePathIn = "/home/nepero27178/Thesis/EHM-HF/Project/data/refined/Mode=rs/Setup=B[128]/Phase=Normal/RB=S_Syms=.csv"
Setup,Phase,Syms,RB,Opt,Layer = UnpackFilePath(FilePathIn)
DF_N = CSV.read(FilePathIn,DataFrame)
filter!(:U => x -> x==0, DF_N)
Sim = Simulation(DF_N,Setup,Phase,Syms,RB)
EnlargeDF!(Sim)

H = Figure(size=(1300,500),figure_padding=1)
axs = [
	Axis3(
		H[1,j],
		xlabel = L"$\delta$",
		ylabel = L"$V$",
		zlabel = L"$\Delta f/t$",
		aspect=(1,1,1),
		azimuth=0.8*pi,
		elevation=pi/10,
		xlabelalign = (:center, :center),
		ylabelalign = (:center, :center),
		zlabelalign = (:center, :center),
		xlabelrotation = 0, # Horizontal xlabel
		ylabelrotation = 0, # Horizontal ylabel
		zlabelrotation = 0, # Horizontal zlabel
	) for j in 1:length(UU)
]

filter!(:U => x -> x==0.0, DF_N)
df_N = copy(DF_N)
for U in [4.0,12.0]
	df_N.U .= fill(U,size(df_N,1))
	global DF_N = vcat(DF_N,df_N)
end
clims=(
	minimum(DF_N.f .- DF_SC.f),
	maximum(DF_N.f .- DF_SC.f)
)

for (j,U) in enumerate(UU)
	df_SC = filter(:U => x -> x==U, DF_SC)
	xx,yy,ff_N = ReshapeData(df_N;xVar="δ",yVar="V",zVar="f")
	_,_,ff_SC = ReshapeData(df_SC;xVar="δ",yVar="V",zVar="f")
	ax = axs[j]
	ax.title = L"$\Delta f/t$ ($t=1.0$, $U=%$(U)$, $L=128$, $\delta=0.0$, $\beta=100.0$)"

	cs = colorschemes[:tabquiet]
	zz = ff_N .- ff_SC
	surface!(ax,xx,yy,zz,colormap=cs,shading=false,colorrange=clims)
end

Colorbar(H[1,0],colormap=colorschemes[:tabquiet],colorrange=clims)
FilePathOut = "f(N)-f(SC).png"
save(FilePathOut, H, px_per_unit=6)