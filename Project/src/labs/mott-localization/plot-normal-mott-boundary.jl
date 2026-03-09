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
# FilePathOut = LAB_ROOT * "/normal-mott-boundary-uS.pdf"

# P = Plot3D(FilePathIn;Print=true,xVar="δ",yVar="V",zVar="tS",cs=:tabcoolerrev)
# H = P[1].H
# df = filter(:V => x -> x==0.0, P[1].DF)
# xx = df.δ
# yy = 2 ./ df.uS

# ax = H.content[1]
# text!(ax,0.1,7.0,text="Mott localization",color=tabblue,align=(:center,:center))
# scatterlines!(ax,xx,yy,color=tabblue,label=L"$V_{\mathrm{c}}[\beta=100,\delta]$")
# ylims!(ax,0.0,8.0)
# axislegend(ax,position=:rb)
# ax.yticks = ([2,4,4.93,6,8],["2","4","4.93","6","8"])

# save(FilePathOut, H)

FilePathOut = LAB_ROOT * "/normal-mott-boundary-f.pdf"

P = Plot3D(FilePathIn;Print=true,xVar="δ",yVar="V",zVar="f",cs=:tabcoolwarm)
H = P[1].H
df = filter(:V => x -> x==0.0, P[1].DF)
xx = vcat(-0.05,df.δ)
yy = vcat(4.93, 2 ./ df.uS)

ax = H.content[1]
text!(ax,0.1,7.0,text="Mott localization",color=:white,align=(:center,:center))
lines!(ax,xx,yy,color=:white,linestyle=(:dash,:dense),linewidth=1)
xlims!(ax,0.0,0.49)
ylims!(ax,0.0,8.0)

save(FilePathOut, H)