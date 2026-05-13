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
include(LAB_ROOT * "/../../modules/methods-2D-plotting.jl")

FilePathIn = "/home/nepero27178/Thesis/EHM-HF/Project/data/refined/Mode=hs/Setup=B[256]/Phase=Normal/RB=S_Syms=.csv"
FilePathOut = LAB_ROOT * "/normal-upper-boundary-uS.pdf"

P = Plot2D(FilePathIn;Print=true,xVar="V",yVar="uS",pVar="δ",cs=:tabwarm)
H = P[1].H
ax = H.content[1]
df = filter(:U => x -> x==0.0, P[1].DF)
xx = unique(df.V)
yy = 2 ./ xx
P1 = lines!(ax,xx,yy,color=:black,linestyle=:dash,linewidth=1)
yy = (2/pi)^2 .* ones(size(xx))
P2 = lines!(ax,xx,yy,color=tabgreen,linestyle=:dash,linewidth=1)
Legend(H[1,3],[P1,P2],[L"$2t/V$",L"$(2/\pi)^2$"],"Theoretical
boundaries:",framevisible=false,titlefont=:regular,valign=:center)

save(FilePathOut, H)