using Contour
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

	mask = zeros(size(ss))
	mask[ ss .>= Δs .&& SS .>= ΔS .&& dd .>= Δd ] .= 3.0
	mask[ ss .>= Δs .&& SS .>= ΔS .&& dd .< Δd ] .= 2.0
	mask[ ss .< Δs .&& SS .< ΔS .&& dd .>= Δd ] .= 1.0
	mask[ ss .< Δs .&& SS .< ΔS .&& dd .< Δd ] .= 0.0

	return mask
end

function InterpolateBoundary(
	x_l::Vector{Float64},
	y_l::Vector{Float64};
	k::Int64=3,
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

FilePathIn = "/home/nepero27178/Thesis/EHM-HF/Project/data/refined/Mode=rs/Setup=B[128]/Phase=SC-Singlet/RB=Sd_Syms=Ssd.csv"
P = Plot3D(FilePathIn;Print=true,xVar="δ",yVar="V",zVar="td",cs=:tabcoolwarm)
H = P[1].H
ax = H.content[1]
DF = DataFrame(P[1].DF)

# Setup,Phase,Syms,RB,Layer = UnpackFilePath(FilePathIn)
# DF = CSV.read(FilePathIn,DataFrame)
# Sim = Simulation(DF,Setup,Phase,Syms,RB)
# EnlargeDF!(Sim)
# filter!(:U => x -> x==0.0, DF)

xx,yy,zz = ReshapeData(DF;xVar="δ",yVar="V",zVar="td")
zz = Matrix(zz)
Δz = 0.002 # Change here
mask = FindNematicMask(zz,Δz)
lvs = Contour.levels(Contour.contours(xx,yy,mask,[0,-1]))

xPos_l, yPos_l = Contour.coordinates( only(Contour.lines(lvs[1])) )
xPos_s, yPos_s = InterpolateBoundary(xPos_l,yPos_l)

xNeg_l, yNeg_l = Contour.coordinates( only(Contour.lines(lvs[2])) )
xNeg_s, yNeg_s = InterpolateBoundary(xNeg_l,yNeg_l)

lines!(ax,xPos_s,yPos_s,color=tabred)
lines!(ax,xNeg_s,yNeg_s,color=tabblue)

save("test.pdf",H)