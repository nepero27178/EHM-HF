using DataFrames
using GLMakie
GLMakie.activate!()

function ShowTrack(
	Track::DataFrame,
	pVar::String
)::Tuple{Figure, Axis}

	xx = [x for x in 1:size(Track)[1]]
	yyI = Track[!,"I" * pVar]
	yyC = Track[!,"C" * pVar]
	yyM = Track[!,"M" * pVar]

	H = Figure(size=(700,500))
	ax = Axis(H[1,1])
	ax.xlabel = "Step"
	ax.ylabel = "Value"
	ax.title = pVar * " evolution"

	slI = scatterlines!(ax,xx,yyI,color="red")
	slC = scatterlines!(ax,xx,yyC,color="blue")
	slM = scatterlines!(ax,xx,yyM,color="cyan")

	Legend(H[1,2],
		[slI, slC, slM],
		["Initializers", "Current", "Mixed"]
	)

	return H, ax
end