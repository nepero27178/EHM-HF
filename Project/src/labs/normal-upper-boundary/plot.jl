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

FilePathIn = "/home/nepero27178/Thesis/EHM-HF/Project/data/refined/Mode=hs/Setup=B[256]/Phase=Normal/RB=S_Syms=.csv"
FilePathOut = LAB_ROOT * "/normal-upper-boundary-uS.pdf"

P = Plot2D(FilePathIn;Print=true,xVar="V",pVar="δ",cVar="U",cs=:tabwarm)
H = P[1].H
df = filter(:U => x -> x==0.0, P[1].DF)
# xx = vcat(-0.05,df.V)
# yy = vcat(4.93, 2 ./ xx)

ax = H.content[1]
# text!(ax,0.1,7.0,text="Mott localisation",color=:black,align=(:center,:center))
# lines!(ax,xx,yy,color=:black,linestyle=(:dash,:dense),linewidth=1)
# xlims!(ax,0.0,0.49)
# ylims!(ax,0.0,8.0)

save(FilePathOut, H)