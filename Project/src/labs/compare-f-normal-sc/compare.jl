using CairoMakie

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

FilePathIn = "/home/nepero27178/Thesis/EHM-HF/Project/data/refined/Mode=rs/Setup=B[128]/Phase=SC-Singlet/RB=Sd_Syms=Ssd.csv"
Setup,Phase,Syms,RB,Opt,Layer = UnpackFilePath(FilePathIn)
DF_SC = CSV.read(FilePathIn,DataFrame)
Sim = Simulation(DF_SC,Setup,Phase,Syms,RB)
df_SC = filter(:U => x -> x==4.0, DF_SC)
U = 4.0

FilePathIn = "/home/nepero27178/Thesis/EHM-HF/Project/data/refined/Mode=rs/Setup=B[128]/Phase=Normal/RB=S_Syms=.csv"
Setup,Phase,Syms,RB,Opt,Layer = UnpackFilePath(FilePathIn)
DF_N = CSV.read(FilePathIn,DataFrame)
df_N = filter(:U => x -> x==0, DF_N)
Sim = Simulation(DF_N,Setup,Phase,Syms,RB)

H = Figure(size=(750,650),figure_padding=1)
ax = Axis3(
	H[1,1],
	xlabel = L"$\delta$",
	ylabel = L"$V$",
	zlabel = L"$[ f^{(\mathrm{N})}-f^{(\mathrm{SC})} ]/t$",
	aspect=(1,1,1),
	azimuth=-0.2*pi,
	elevation=pi/9,
	xlabelalign = (:center, :center),
	ylabelalign = (:center, :center),
	zlabelalign = (:center, :center),
	xlabelrotation = 0, # Horizontal xlabel
	ylabelrotation = 0, # Horizontal ylabel
)

xx,yy,ff_N = ReshapeData(df_N;xVar="δ",yVar="V",zVar="f")
_,_,ff_SC = ReshapeData(df_SC;xVar="δ",yVar="V",zVar="f")
ax.title = L"$[ f^{(\mathrm{N})}-f^{(\mathrm{SC})} ]/t$ ($t=1.0$, $U=%$(U)$, $L=128$, $\beta=100.0$)"

cs = colorschemes[:tabquietrev]
zz = ff_N .- ff_SC
s = surface!(ax,xx,yy,zz,colormap=cs,shading=false)
# Colorbar(H[1,2],colormap=colorschemes[:tabquiet])

DF_n = CSV.read("/home/nepero27178/Thesis/EHM-HF/Project/src/labs/superconductor-region-boundaries/n_U=$(U).csv",DataFrame)
x_n = DF_n.x_n
y_n = DF_n.y_n
z_n = []
for (i,x) in enumerate(x_n)
	y = y_n[i]
	row_SC = filter([:δ,:V] => (u,v) -> u==x && v==y,df_SC)
	row_N = filter([:δ,:V] => (u,v) -> u==x && v==y,df_N)
	push!(z_n,only(row_N.f)-only(row_SC.f))
end
z_n = z_n .+ 0.005
scatterlines!(ax,x_n,y_n,z_n,color=:black,markersize=2,linewidth=0.1)
text!(ax,0.4,0.7,0.1,text="Normal",align=(:center,:center))

Colorbar(H[1,0],s)
FilePathOut = "f(N)-f(SC).pdf"
save(FilePathOut, H)