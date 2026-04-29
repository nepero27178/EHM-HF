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

function FindNematicMask(
	tt::Matrix{Float64},
	Δt::Float64
)

	mask = zeros(size(tt))
	mask[ tt .>= Δt ] .= 1.0
	mask[ tt .<= -Δt ] .= -1.0
	mask[ tt .< Δt .&& tt .> -Δt ] .= 0.0

	return mask
end

function InterpolateBoundary(
	x_l::Vector{Float64},
	y_l::Vector{Float64};
	k::Int64=4,
	s::Float64=0.2,
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
J = 3
H = P[J].H
ax_H = H.content[1]
DF = DataFrame(P[J].DF)
U = only(unique(DF.U))

# Setup,Phase,Syms,RB,Layer = UnpackFilePath(FilePathIn)
# DF = CSV.read(FilePathIn,DataFrame)
# Sim = Simulation(DF,Setup,Phase,Syms,RB)
# EnlargeDF!(Sim)
# filter!(:U => x -> x==U, DF)

xx,yy,zz = ReshapeData(DF;xVar="δ",yVar="V",zVar="td")
zz = Matrix(zz)

xx,yy,zz = InterpolateSurface(xx,yy,zz;L=1000)
Δz = 0.002 # Change here
mask = FindNematicMask(zz,Δz)
lvs = Contour.levels(Contour.contours(xx,yy,mask,[0,-1]))

xPos_l, yPos_l = Contour.coordinates( only(Contour.lines(lvs[1])) )
xPos_s, yPos_s = InterpolateBoundary(xPos_l,yPos_l)
lines!(ax_H,xPos_s,yPos_s,color=tabred)

try
	xNeg_l, yNeg_l = Contour.coordinates( only(Contour.lines(lvs[2])) )
	xNeg_s, yNeg_s = InterpolateBoundary(xNeg_l,yNeg_l)
	lines!(ax_H,xNeg_s,yNeg_s,color=tabblue)
catch
end

save("t-U=$(U)-curves.pdf",H)

F = Figure(size=(600,400),figure_padding=1)
ax_F = Axis(
	F[1,1],
	xlabel = ax_H.xlabel,
	ylabel = ax_H.ylabel,
	title = ax_H.title
)

if J==1
	heatmap!(ax_F,xx,yy,mask,colormap=CoolWarm)
	text!(ax_F,0.1,6.0,text=L"$t^{(d)}\ge\Delta t$",color=:white,align=(:center,:center))
	text!(ax_F,0.35,6.0,text=L"$t^{(d)}\le -\Delta t$",color=:white,align=(:center,:center))
	text!(ax_F,0.25,1.0,text=L"$| t^{(d)} | < \Delta t$",color=:black,align=(:center,:center))
elseif J==2 || J==3
	heatmap!(ax_F,xx,yy,mask,colormap=Warm)
	text!(ax_F,0.25,3.7,text=L"$t^{(d)}\ge\Delta t$",color=:white,align=(:center,:center))
	text!(ax_F,0.25,1.0,text=L"$| t^{(d)} | < \Delta t$",color=:black,align=(:center,:center))
end

save("t-U=$(U)-solid.pdf",F)
