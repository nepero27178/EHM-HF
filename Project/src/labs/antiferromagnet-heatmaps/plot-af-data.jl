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

FilePathIn = "/home/nepero27178/Thesis/EHM-HF/Project/data/refined/Mode=rs/Setup=C[128]/Phase=AF/RB=S_Syms=_tailored2.csv"
FilePathOut = LAB_ROOT * "/af-data-C128.pdf"

H = Figure(size=(1300,400),figure_padding=1)

P = Plot3D(FilePathIn;Print=true,xVar="U",yVar="V",zVar="m",colorbar=false)
Hm = P[1].H
dfm = P[1].DF

# FFS there is no simpler way to do this in Makie!?
axm_P = Hm.content[1]
axm_H = Axis(H[1,1])
axm_H.xlabel = axm_P.xlabel.val
axm_H.ylabel = axm_P.ylabel.val
axm_H.title = axm_P.title.val
axm_H.xticks = axm_P.xticks.val
axm_H.yticks = axm_P.yticks.val

xx, yy, zz = ReshapeData(DataFrame(dfm);xVar="U",yVar="V",zVar="m")
hm = heatmap!(axm_H,xx,yy,zz,colormap=colorschemes[:tabcoolerrev])

P = Plot3D(FilePathIn;Print=true,xVar="U",yVar="V",zVar="uS",colorbar=false)
HuS = P[1].H
dfuS = P[1].DF

# FFS there is no simpler way to do this in Makie!?
axuS_P = HuS.content[1]
axuS_H = Axis(H[1,2])
axuS_H.xlabel = axuS_P.xlabel.val
axuS_H.ylabel = axuS_P.ylabel.val
axuS_H.title = axuS_P.title.val
axuS_H.xticks = axuS_P.xticks.val
axuS_H.yticks = axuS_P.yticks.val

xx, yy, zz = ReshapeData(DataFrame(dfuS);xVar="U",yVar="V",zVar="uS")
heatmap!(axuS_H,xx,yy,zz,colormap=colorschemes[:tabcoolerrev])

axuS_H.ylabelvisible = false
axuS_H.yticklabelsvisible = false
Colorbar(H[1,3],hm)

# text!(ax,0.1,7.0,text="Mott localization",color=:white,align=(:center,:center))
# lines!(ax,xx,yy,color=:white,linestyle=(:dash,:dense),linewidth=1)
# xlims!(ax,0.0,0.49)
# ylims!(ax,0.0,8.0)

save(FilePathOut, H)