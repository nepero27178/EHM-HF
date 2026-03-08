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

LAB_ROOT = @__DIR__
include(LAB_ROOT * "/../../setup/graphic-setup.jl")
include(LAB_ROOT * "/../../modules/methods-3D-plotting.jl")

FilePathIn = "/home/nepero27178/Thesis/EHM-HF/Project/data/refined/Mode=rs/Setup=B[128]/Phase=Normal/RB=S_Syms=.csv"
FilePathOut = LAB_ROOT * "/normal-mott-boundary.pdf"

uP = Plot3D(FilePathIn;Print=true,xVar="V",yVar="β",zVar="uS",cs=:tabcoolquiet)
tP = Plot3D(FilePathIn;Print=true,xVar="V",yVar="β",zVar="tS",cs=:tabcoolquiet)

save(FilePathOut, H)
