using Contour
using Interpolations
using Dierckx

LAB_ROOT = @__DIR__
include(LAB_ROOT * "/../../setup/graphic-setup.jl")
include(LAB_ROOT * "/../../modules/methods-3D-plotting.jl")

GLMakie.activate!()
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

	# mask = zeros(size(ss)) # INTRODUCES A BIAS
	# mask[ ss .>= Δs .&& SS .>= ΔS .&& dd .>= Δd ] .= 3.0
	# mask[ ss .>= Δs .&& SS .>= ΔS .&& dd .< Δd ] .= 2.0
	# mask[ ss .< Δs .&& SS .< ΔS .&& dd .>= Δd ] .= 1.0
	# mask[ ss .< Δs .&& SS .< ΔS .&& dd .< Δd ] .= 0.0

	mask_s = ss .> Δs
	mask_S = SS .> ΔS
	mask_s = mask_s .|| mask_S

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
J = 2
H = P[J].H
ax_H = H.content[1]
DF = DataFrame(P[J].DF)
U = only(unique(DF.U))

VV = unique(DF.V)
δδ = unique(DF.δ)

# Setup,Phase,Syms,RB,Layer = UnpackFilePath(FilePathIn)
# DF = CSV.read(FilePathIn,DataFrame)
# Sim = Simulation(DF,Setup,Phase,Syms,RB)
# EnlargeDF!(Sim)
# filter!(:U => x -> x==U, DF)

xx,yy,ss = ReshapeData(DF;xVar="δ",yVar="V",zVar="ws")
ss = -U .* Matrix(ss)

_,_,SS = ReshapeData(DF;xVar="δ",yVar="V",zVar="wS")
SS = (VV * ones(length(δδ))') .* Matrix(SS)

_,_,dd = ReshapeData(DF;xVar="δ",yVar="V",zVar="wd")
dd = (VV * ones(length(δδ))') .* Matrix(dd)

_,_,ss = InterpolateSurface(xx,yy,ss)
Δs = 0.01 # Change here

_,_,SS = InterpolateSurface(xx,yy,SS)
ΔS = 0.1 # Change here

xx,yy,dd = InterpolateSurface(xx,yy,dd)
Δd = 0.1 # Change here

mask = FindGapMask(ss,Δs,SS,ΔS,dd,Δd)
lvs = Contour.levels(Contour.contours(xx,yy,mask,[0,1,2,3]))

x0_l, y0_l = Contour.coordinates( only(Contour.lines(lvs[1])) )
x0_s, y0_s = InterpolateBoundary(x0_l,y0_l)
lines!(ax_H,x0_s,y0_s,color=tabred)

x1_l, y1_l = Contour.coordinates( only(Contour.lines(lvs[2])) )
x1_s, y1_s = InterpolateBoundary(x1_l,y1_l)
lines!(ax_H,x1_s,y1_s,color=tabgreen)

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
	xlabel = ax_H.xlabel,
	ylabel = ax_H.ylabel,
	title = ax_H.title
)

if J==1
	heatmap!(ax_F,xx,yy,mask,colormap=CoolWarm)
	# text!(ax_F,0.1,6.0,text=L"$t^{(d)}\ge\Delta t$",color=:white,align=(:center,:center))
	# text!(ax_F,0.35,6.0,text=L"$t^{(d)}\le -\Delta t$",color=:white,align=(:center,:center))
	# text!(ax_F,0.25,1.0,text=L"$| t^{(d)} | < \Delta t$",color=:black,align=(:center,:center))
elseif J==2 || J==3
	heatmap!(ax_F,xx,yy,mask,colormap=Warm)
	# text!(ax_F,0.25,3.7,text=L"$t^{(d)}\ge\Delta t$",color=:white,align=(:center,:center))
	# text!(ax_F,0.25,1.0,text=L"$| t^{(d)} | < \Delta t$",color=:black,align=(:center,:center))
end

save("w-U=$(U)-solid.pdf",F)
