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
Setup,Phase,Syms,RB,Opt,Layer = UnpackFilePath(FilePathIn)
DF = CSV.read(FilePathIn,DataFrame)
Sim = Simulation(DF,Setup,Phase,Syms,RB)
EnlargeDF!(Sim)

# CHANGE WHEN DATA ARE FIXED
# filter!(:U => x -> x!=0.0, DF)
UU = unique(DF.U)
SymmetricRange = (-maximum(abs.(DF.td)),maximum(abs.(DF.td)))
AsymmetricRange = (-0.03,0.07) #extrema(DF.td)
BlueSteps = 30
RedSteps = 70
AsymCoolWarm = vcat(
	[Colors.convert(RGB, c) for c in range(TabBlueLab, stop=WhiteLab, length=BlueSteps)],
	[Colors.convert(RGB, c) for c in range(WhiteLab, stop=TabRedLab, length=RedSteps)]
)

H = Figure(size=(400*length(UU)+100,400),figure_padding=1)
axs_H = [Axis(
	H[1,i],
	xlabel = L"$\delta$",
	ylabel = L"$V$"
) for i in 1:length(UU)]
linkyaxes!(axs_H...)

F = Figure(size=(400*length(UU)+100,400),figure_padding=1)
axs_F = [Axis(
	F[1,i],
	xlabel = L"$\delta$",
	ylabel = L"$V$"
) for i in 1:length(UU)]
linkyaxes!(axs_F...)

for (i,U) in enumerate(UU)

	ax_H = axs_H[i]
	ax_H.title = L"$U=%$(U)$"
	ax_F = axs_F[i]
	ax_F.title = L"$U=%$(U)$"
	if i>1
		ax_H.ylabelvisible = false
		ax_H.yticklabelsvisible = false
		ax_F.ylabelvisible = false
		ax_F.yticklabelsvisible = false
	end
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

	_,_,zz = ReshapeData(df;xVar="δ",yVar="V",zVar="td")
	heatmap!(
		ax_H,
		xx,yy,zz,
		colormap=AsymCoolWarm,
		colorrange=AsymmetricRange,
	)

	mask = FindGapMask(ss,Δs,SS,ΔS,dd,Δd)
	lvs = Contour.levels(Contour.contours(xx,yy,mask,[2]))
	xc, yc = Contour.coordinates( only(Contour.lines(lvs[1])) )
	scatterlines!(ax_H,xc,yc,color=:black,markersize=2,linewidth=0.1)
	text!(ax_H,0.08,3.0,text="Mixing boundary",align=(:center,:top),color=:black)

	if U==0.0
		text!(ax_H,0.1,6.2,text=L"$\mathrm{N}_y$",align=(:center,:top),color=:black)
		text!(ax_H,0.35,6.7,text=L"$\mathrm{N}_x$",align=(:center,:top),color=:white)
		heatmap!(ax_F,xx,yy,mask,colormap=colorschemes[:tabcoolrev])
		text!(ax_F,0.3,0.6,text="Normal",align=(:center,:center))
		text!(ax_F,0.18,1.9,text=L"$d$-wave",align=(:center,:center))
		text!(ax_F,0.3,3.5,text="Mixed",align=(:center,:center),color=:white)
		text!(ax_F,0.42,2.25,text=L"$s^*$-wave",align=(:center,:center),color=:white)
	elseif U==4.0
		text!(ax_H,0.1,6.2,text=L"$\mathrm{N}_y$",align=(:center,:top),color=:black)
		text!(ax_H,0.37,6.7,text=L"$\mathrm{N}_x$",align=(:center,:top),color=:black)
		heatmap!(ax_F,xx,yy,mask,colormap=colorschemes[:tabcoolrev])
		text!(ax_F,0.35,0.75,text="Normal",align=(:center,:center))
		text!(ax_F,0.2,2.0,text=L"$d$-wave",align=(:center,:center))
		text!(ax_F,0.32,3.5,text="Mixed",align=(:center,:center),color=:white)
		text!(ax_F,0.43,2.35,text=L"$s \oplus s^*$-wave",align=(:center,:center),color=:white)
	elseif U==12.0
		text!(ax_H,0.18,6.2,text=L"$\mathrm{N}_y$",align=(:center,:top),color=:black)
		heatmap!(ax_F,xx,yy,mask,colormap=colorschemes[:tabcoolrev])
		text!(ax_F,0.375,1.0,text="Normal",align=(:center,:center))
		text!(ax_F,0.26,2.1,text=L"$d$-wave",align=(:center,:center))
		text!(ax_F,0.33,3.5,text="Mixed",align=(:center,:center),color=:white)
		text!(ax_F,0.43,2.52,text=L"$s \oplus s^*$-wave",align=(:center,:center),color=:white)
	end
end

Label(
	H[0,:],
	L"$\tilde{t}^{(d)}$ ($t=1.0$, $L=256$, $\beta=100.0$)",
	fontsize=15
)
Colorbar(
	H[1,length(UU)+1],
	colormap=AsymCoolWarm,
	colorrange=AsymmetricRange,
	ticks = [-0.03,0.0,0.03,0.06,0.07]
)
save("w-U=$(UU)-curves.pdf",H)

Label(
	F[0,:],
	L"SC and normal regions ($t=1.0$, $L=256$, $\beta=100.0$)",
	fontsize=15
)
save("w-U=$(UU)-solid.pdf",F)
