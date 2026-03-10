using DataFrames
using CairoMakie
using CSV
# using GLMakie
# GLMakie.activate!()

CairoMakie.activate!()
MT = Makie.MathTeXEngine
MT_DIR = dirname(pathof(MT)) * "/../assets/fonts/NewComputerModern"
set_theme!(fonts = (
	regular = MT_DIR * "/NewCM10-Regular.otf",
	bold = MT_DIR * "/NewCM10-Bold.otf"
))

LAB_ROOT = @__DIR__
include(LAB_ROOT * "/../../setup/graphic-setup.jl")

function ShowTrack(
	Track::DataFrame,
	pVar::String,
	ax::Axis
)::Axis

	xx = [x for x in 1:size(Track)[1]]
	yyI = Track[!,"I" * pVar]
	yyC = Track[!,"C" * pVar]
	yyM = Track[!,"M" * pVar]

	ax.xlabel = "Step"
	ax.ylabel = "Value"
	hlines!(ax,1/3,color=:gray,label=L"$2t/V$", linestyle=(:dash,:dense))

	slI = scatterlines!(ax,xx,yyI,color=tabred,label="Initializer")
	slC = scatterlines!(ax,xx,yyC,color=tabblue,label="Current step")
	slM = scatterlines!(ax,xx,yyM,color=tabgreen,label="Mixed value")

	return ax
end

function CompareTracks(
	Tracks::Vector{DataFrame},
	pVar::String
)

	H::Figure = Figure(size=(1500,400),figure_padding = 1)
	axs::Vector{Axis} = [Axis(H[1,i]) for i in 1:length(Tracks)]
	linkyaxes!(axs...)

	for (i,Track) in enumerate(Tracks)
		ax::Axis = axs[i]
		ax = ShowTrack(Track,pVar,ax)
		if i>1
			ax.ylabelvisible = false
			ax.yticklabelsvisible = false
		end
	end

	axs[1].title = L"$g=0.5$"
	axs[2].title = L"$g=0.25$"
	axs[3].title = L"$g=0.05$"

	Legend(H[1,4],axs[1],framevisible = false)
	Label(H[0,:],L"$u^{(s^*)}$ ($t=1.0$, $U=0.0$, $V=6.0$, $L=128$, $\beta=100.0$, $\delta=0.1$)")

	return H
end

FilePathIn05 = "data/g=0.5.track.csv"
FilePathIn025 = "data/g=0.25.track.csv"
FilePathIn005 = "data/g=0.05.track.csv"

Tracks = [
	CSV.read(FilePathIn05, DataFrame),
	CSV.read(FilePathIn025, DataFrame),
	CSV.read(FilePathIn005, DataFrame),
]
H = CompareTracks(Tracks,"uS")
save("tracks.pdf",H)