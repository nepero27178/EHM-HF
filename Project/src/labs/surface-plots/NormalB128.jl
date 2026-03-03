using CairoMakie
using Colors
using ColorSchemes
using Contour

CairoMakie.activate!()
MT = Makie.MathTeXEngine
MT_DIR = dirname(pathof(MT)) * "/../assets/fonts/NewComputerModern"
set_theme!(fonts = (
	regular = MT_DIR * "/NewCM10-Regular.otf",
	bold = MT_DIR * "/NewCM10-Bold.otf"
))

LAB_DIR = @__DIR__
PROJECT_SRC_DIR = dirname(dirname(LAB_DIR))
include(PROJECT_SRC_DIR * "/setup/graphic-setup.jl")
include(PROJECT_SRC_DIR * "/modules/methods-3D-plotting.jl")

FilePathIn = dirname(PROJECT_SRC_DIR) * "/simulations/Mode=rs/Setup=B[128]/Phase=Normal/RB=S_Syms=.csv"
Setup, Phase, Syms, RB = UnpackFilePath(FilePathIn)
DF::DataFrame = CSV.read(FilePathIn,DataFrame)
filter!([:U,:V] => (x,y) -> x==0.0 && y!=4.0, DF)
Sim::Simulation = Simulation(DF,Setup,Phase,Syms,RB)
xx, yy, zz = ReshapeData(DF; xVar="δ", yVar="V", zVar="uS")

FilePathIn = dirname(PROJECT_SRC_DIR) * "/simulations/Mode=rs/Setup=B[128]a/Phase=Normal/RB=S_Syms=.csv"
Setup, Phase, Syms, RB = UnpackFilePath(FilePathIn)
DFa::DataFrame = CSV.read(FilePathIn,DataFrame)
Sim::Simulation = Simulation(DF,Setup,Phase,Syms,RB)
xxa, yya, zza = ReshapeData(DF; xVar="δ", yVar="V", zVar="uS")

yy = vcat(yy,yya)
zz = vcat(zz,zza)