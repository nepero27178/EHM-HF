using Contour
using Interpolations
using Dierckx

LAB_ROOT = @__DIR__
include(LAB_ROOT * "/../../setup/graphic-setup.jl")
include(LAB_ROOT * "/../../modules/methods-3D-plotting.jl")

CairoMakie.activate!()
MT = Makie.MathTeXEngine
MT_DIR = dirname(pathof(MT)) * "/../assets/fonts/NewComputerModern"
set_theme!(fonts = (
	regular = MT_DIR * "/NewCM10-Regular.otf",
	bold = MT_DIR * "/NewCM10-Bold.otf"
))

function FindGapMask(
	ss::Matrix{Float64},
	Δs::Float64,
	SS::Matrix{Float64},
	ΔS::Float64,
	dd::Matrix{Float64},
	Δd::Float64
)

	ss = abs.(ss)
	SS = abs.(SS)
	dd = abs.(dd)

	mask_s = ss .> Δs
	mask_S = SS .> ΔS
	mask_s = mask_s .|| mask_S # Either s or s*
	mask_d = dd .> Δd

	mask = 2*mask_s .+ mask_d

	return mask
end

function InterpolateBoundary(
	x_l::Vector{Float64},
	y_l::Vector{Float64};
	k::Int64=2,
	s::Float64=1.0,
	L::Int64=100
)

	t = range(0,stop=1,length=length(x_l))
	spline_x_l = Spline1D(t,x_l;k,s)
	spline_y_l = Spline1D(t,y_l;k,s)
	t_s = range(0, stop=1, length=L)
	x_s = spline_x_l(t_s)
	y_s = spline_y_l(t_s)
	return x_s, y_s

end

function InterpolateSurface(
	xx::Vector{Float64},
	yy::Vector{Float64},
	zz::Matrix{Float64};
	L::Int64=1000
)::Tuple{Vector{Float64},Vector{Float64},Matrix{Float64}}

	IntSurf = interpolate(
		(xx, yy), zz,
		( Gridded(Linear()), Gridded(Linear()) )
	)
	xx_I = range(extrema(xx)..., length=L)
	yy_I = range(extrema(yy)..., length=L)
	zz_I = [IntSurf(x,y) for x in xx_I, y in yy_I]

	return [x for x in xx_I], [y for y in yy_I], zz_I

end

FilePathIn = "/home/nepero27178/Thesis/EHM-HF/Project/data/refined/Mode=rs/Setup=B[128]/Phase=SC-Singlet/RB=Sd_Syms=Ssd.csv"
P = Plot3D(FilePathIn;Print=true,xVar="δ",yVar="V",zVar="td",cs=:tabcoolwarm)
J = 1
H = P[J].H
ax_H = H.content[1]
DF = DataFrame(P[J].DF)
U = only(unique(DF.U))

# Setup,Phase,Syms,RB,Layer = UnpackFilePath(FilePathIn)
# DF = CSV.read(FilePathIn,DataFrame)
# Sim = Simulation(DF,Setup,Phase,Syms,RB)
# EnlargeDF!(Sim)
# U = 4.0
# filter!(:U => x -> x==U, DF)

VV = unique(DF.V)
δδ = unique(DF.δ)

xx,yy,ss = ReshapeData(DF;xVar="δ",yVar="V",zVar="ws")
ss = -U .* Matrix(ss)

_,_,SS = ReshapeData(DF;xVar="δ",yVar="V",zVar="wS")
SS = (ones(length(δδ)) * VV') .* Matrix(SS)

_,_,dd = ReshapeData(DF;xVar="δ",yVar="V",zVar="wd")
dd = (ones(length(δδ)) * VV') .* Matrix(dd)

# _,_,ss = InterpolateSurface(xx,yy,ss,L=200)
Δs = 0.01 # Change here

# _,_,SS = InterpolateSurface(xx,yy,SS,L=200)
ΔS = 0.01 # Change here

# xx,yy,dd = InterpolateSurface(xx,yy,dd,L=200)
Δd = 0.01 # Change here

mask = FindGapMask(ss,Δs,SS,ΔS,dd,Δd)
lvs = Contour.levels(Contour.contours(xx,yy,mask,[2]))
xc, yc = Contour.coordinates( only(Contour.lines(lvs[1])) )
scatterlines!(ax_H,xc,yc,color=:white,markersize=5)
if J==2 || J==3
	text!(ax_H,0.08,3.2,text="Mixing boundary",align=(:center,:center),color=:white)
end

# x0_l, y0_l = Contour.coordinates( only(Contour.lines(lvs[1])) )
# x0_s, y0_s = InterpolateBoundary(x0_l,y0_l)
# lines!(ax_H,x0_s,y0_s,color=tabred)

# x1_l, y1_l = Contour.coordinates( only(Contour.lines(lvs[2])) )
# x1_s, y1_s = InterpolateBoundary(x1_l,y1_l)
# lines!(ax_H,x1_s,y1_s,color=tabgreen)

# try
# 	xNeg_l, yNeg_l = Contour.coordinates( only(Contour.lines(lvs[2])) )
# 	xNeg_s, yNeg_s = InterpolateBoundary(xNeg_l,yNeg_l)
# 	lines!(ax_H,xNeg_s,yNeg_s,color=tabblue)
# catch
# end

save("w-U=$(U)-curves.pdf",H)

F = Figure(size=(600,400),figure_padding=1)
ax_F = Axis(
	F[1,1],
	xlabel = L"$\delta$",
	ylabel = L"$V$",
	title = L"SC regions ($t=1.0$, $U=%$(U)$, $L=128$, $\beta=100.0$)"
)

if J==1
	heatmap!(ax_F,xx,yy,mask,colormap=CoolWarm)
	# text!(ax_F,0.1,6.0,text=L"$t^{(d)}\ge\Delta t$",color=:white,align=(:center,:center))
	# text!(ax_F,0.35,6.0,text=L"$t^{(d)}\le -\Delta t$",color=:white,align=(:center,:center))
	# text!(ax_F,0.25,1.0,text=L"$| t^{(d)} | < \Delta t$",color=:black,align=(:center,:center))
elseif J==2
	heatmap!(ax_F,xx,yy,mask,colormap=colorschemes[:tabcoolrev])
	text!(ax_F,0.35,0.75,text="Normal",align=(:center,:center))
	text!(ax_F,0.2,2.0,text=L"$d$-wave",align=(:center,:center))
	text!(ax_F,0.3,3.5,text="Mixed",align=(:center,:center),color=:white)
	text!(ax_F,0.42,2.25,text=L"$s \oplus s^*$-wave",align=(:center,:center),color=:white)
elseif J==3
	heatmap!(ax_F,xx,yy,mask,colormap=colorschemes[:tabcoolrev])
	text!(ax_F,0.375,1.0,text="Normal",align=(:center,:center))
	text!(ax_F,0.26,2.1,text=L"$d$-wave",align=(:center,:center))
	text!(ax_F,0.33,3.5,text="Mixed",align=(:center,:center),color=:white)
	text!(ax_F,0.43,2.6,text=L"$s \oplus s^*$-wave",align=(:center,:center),color=:white)
end

save("w-U=$(U)-solid.pdf",F)
