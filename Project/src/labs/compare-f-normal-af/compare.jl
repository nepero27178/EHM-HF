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

H = Figure(size=(800,650),figure_padding=1)
ax = Axis3(
	H[1,1],
	xlabel = L"$U$",
	ylabel = L"$V$",
	zlabel = L"$\Delta f/t$",
	title = L"$\Delta f/t$ ($t=1.0$, $L=128$, $\delta=0.0$, $\beta=100.0$)",
	aspect=(1,1,1),
	azimuth=0.8*pi,
	elevation=pi/10,
	xlabelalign = (:center, :center),
	ylabelalign = (:center, :center),
	zlabelalign = (:center, :center),
	xlabelrotation = 0, # Horizontal xlabel
	ylabelrotation = 0, # Horizontal ylabel
	zlabelrotation = 0, # Horizontal zlabel
)

cs = colorschemes[:tabquiet]
xx, yy, zz_AF = ReshapeData(DF_AF;xVar="U",yVar="V",zVar="f")
_, _, zz_N = ReshapeData(DF_N;xVar="U",yVar="V",zVar="f")
zz = zz_N .- zz_AF
clims=(
	minimum(0.0),
	maximum(zz)
)
h = surface!(ax,xx,yy,zz,colormap=cs,shading=false)
Colorbar(H[1,0], h)

FilePathOut = LAB_ROOT * "/f(N)-f(AF).png"
save(FilePathOut, H, px_per_unit=6)