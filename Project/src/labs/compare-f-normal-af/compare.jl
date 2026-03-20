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
P = Plot3D(FilePathIn;Print=true,xVar="U",yVar="V",zVar="f",colorbar=false)
ax_P = P[1].H.content[1]

H = Figure(size=(600,400),figure_padding=1)
ax = Axis(H[1,1])
ax.xlabel = ax_P.xlabel.val
ax.ylabel = ax_P.ylabel.val
ax.title = latexstring(replace(ax_P.title.val, L"$f_\mathrm{MFT}$" => L"$\Delta f_\mathrm{MFT}$"))

cs = colorschemes[:tabcoolwarm]

DFAF = P[1].DF
select!(DFAF, [:U,:V,:f])
xx, yy, zzAF = ReshapeData(DataFrame(DFAF);xVar="U",yVar="V",zVar="f")

FilePathIn = LAB_ROOT * "/../../../data/refined/Mode=rs/Setup=B[128]/Phase=Normal/RB=S_Syms=.csv"

# Read normal DF and filter out interface
DFNormal = CSV.read(FilePathIn,DataFrame)
DFNormal4 = filter(:V => x -> x==4.0, DFNormal)
filter(:g => x -> x!=0.5, DFNormal4)
filter!([:V,:g] => (x,y) -> x!=4.0, DFNormal)
DFNormal = vcat(DFNormal,DFNormal4)
filter!([:U,:δ] => (x,y) -> x==0.0 && y==0.0, DFNormal)
select!(DFNormal, [:U,:V,:f])

DF = copy(DFNormal)
UU = unique(DFAF.U)[2:end]
for U in UU
	DF.U .= fill(U,size(DF,1))
	global DFNormal = vcat(DFNormal,DF)
end

_, _, zzNormal = ReshapeData(DataFrame(DFNormal);xVar="U",yVar="V",zVar="f")
zz = zzNormal .- zzAF
clims=(
	minimum(0.0),
	maximum(zz)
)
h = heatmap!(ax,xx,yy,zz,colormap=cs)
Colorbar(H[1,2], h)

FilePathOut = LAB_ROOT * "/compare-f-normal-af.pdf"
save(FilePathOut, H)