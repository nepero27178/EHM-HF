using Contour
using DataFrames
using CSV

LAB_ROOT = @__DIR__
include(LAB_ROOT * "/../../modules/structs.jl")
include(LAB_ROOT * "/../../modules/methods-IO.jl")

function FindsMask(
	ss::Matrix{Float64},
	Δs::Float64,
	SS::Matrix{Float64},
	ΔS::Float64
)

	ss = abs.(ss)
	SS = abs.(SS)

	mask_s = ss .> Δs
	mask_S = SS .> ΔS
	mask_s = mask_s .|| mask_S # Either s or s*
	return mask_s
end

function FinddMask(
	dd::Matrix{Float64},
	Δd::Float64,
)

	dd = abs.(dd)
	mask_d = dd .> Δd
	return mask_d
end

FilePathIn = "/home/nepero27178/Thesis/EHM-HF/Project/data/refined/Mode=rs/Setup=B[128]/Phase=SC-Singlet/RB=Sd_Syms=Ssd.csv"
DF = CSV.read(FilePathIn,DataFrame)

U = 4.0
df = filter(:U => x -> x==U, DF)
VV = unique(df.V)
δδ = unique(df.δ)

xx,yy,ss = ReshapeData(df;xVar="δ",yVar="V",zVar="ws")
ss = -U .* Matrix(ss)
Δs = 0.01 # Change here

_,_,SS = ReshapeData(df;xVar="δ",yVar="V",zVar="wS")
SS = (ones(length(δδ)) * VV') .* Matrix(SS)
ΔS = 0.01 # Change here

_,_,dd = ReshapeData(df;xVar="δ",yVar="V",zVar="wd")
dd = (ones(length(δδ)) * VV') .* Matrix(dd)
Δd = 0.01 # Change here

mask_s = FindsMask(ss,Δs,SS,ΔS)
mask_d = FinddMask(dd,Δd)
mask_n = .!mask_s .* .!mask_d
mask_m = mask_s .* mask_d

lvs_s = Contour.levels(Contour.contours(xx,yy,mask_s,[0]))
xc_s, yc_s = Contour.coordinates( only(Contour.lines(lvs_s[1])) )

lvs_d = Contour.levels(Contour.contours(xx,yy,mask_d,[0]))
xc_d, yc_d = Contour.coordinates( only(Contour.lines(lvs_d[1])) )

lvs_n = Contour.levels(Contour.contours(xx,yy,mask_n,[0]))
xc_n, yc_n = Contour.coordinates( only(Contour.lines(lvs_n[1])) )

lvs_m = Contour.levels(Contour.contours(xx,yy,mask_m,[0]))
xc_m, yc_m = Contour.coordinates( only(Contour.lines(lvs_m[1])) )

FilePathOut = "s_U=$(U).csv"
DF_s = DataFrame(Dict(
	"x_s" => xc_s,
	"y_s" => yc_s,
))
CSV.write(FilePathOut,DF_s)

FilePathOut = "d_U=$(U).csv"
DF_d = DataFrame(Dict(
	"x_d" => xc_d,
	"y_d" => yc_d,
))
CSV.write(FilePathOut,DF_d)

FilePathOut = "n_U=$(U).csv"
DF_n = DataFrame(Dict(
	"x_n" => xc_n,
	"y_n" => yc_n,
))
CSV.write(FilePathOut,DF_n)

FilePathOut = "m_U=$(U).csv"
DF_m = DataFrame(Dict(
	"x_m" => xc_m,
	"y_m" => yc_m,
))
CSV.write(FilePathOut,DF_m)