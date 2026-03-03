#!/usr/bin/julia

using GLMakie
using CairoMakie
using LaTeXStrings
using ColorSchemes
using DataFrames
using CSV

PROJECT_MODULES_DIR = @__DIR__
include(PROJECT_MODULES_DIR * "/structs.jl")
include(PROJECT_MODULES_DIR * "/methods-IO.jl")

function Plot2D(
	FilePathIn::String;
	Print::Bool=false,
	xVar::String="V",
	yVar::String="m",
	pVar::String="δ",
	cs::Symbol=:imola50,
	Skip::Int64=0,
)::Vector{GroupedPlot}

	# List vars and pars
	xVars = ["t", "U", "V", "L", "δ", "β"]
	Pars = filter(!=(pVar), filter(!=(xVar), xVars))

	# Input safecheck
	!in(xVar, xVars) ? error("Invalid x variable, choose one of $(xVars)") : false
	!in(pVar, xVars) ? error("Invalid p variable, choose one of $(xVars)") : false
	xVar==pVar ? error("You have chosen xVar=pVar!") : false

	# Unpack filepath
	Setup, Phase, Syms, RB, _ = UnpackFilePath(FilePathIn)

	# Load data
	DF::DataFrame = CSV.read(FilePathIn,DataFrame)
	Sim::Simulation = Simulation(DF,Setup,Phase,Syms,RB)
	DF = EnlargeDF!(Sim) # Compute RMPs
	yVars = filter(!in(xVars), names(DF))
	!in(yVar,yVars) ? error("Invalid y variable, choose one of $(yVars)") : false

	if Print
		# Activate backend
		CairoMakie.activate!()

		MT = Makie.MathTeXEngine
		MT_DIR = dirname(pathof(MT)) * "/../assets/fonts/NewComputerModern"

		set_theme!(fonts = (
			regular = MT_DIR * "/NewCM10-Regular.otf",
			bold = MT_DIR * "/NewCM10-Bold.otf"
		))

		# Get LaTeX formatted labels
		VarLabels = GetLabels(Phase)
	elseif !Print
		# Activate backend
		GLMakie.activate!()
	end

	# Group data
	GroupedDF = groupby(DF,Pars)
	J = length(GroupedDF)
	C = floor(Int64, length(colorschemes[cs]) / J)
	PlotVec = GroupedPlot[]

	# Cycle over simulated points
	for (j,df) in enumerate(GroupedDF)

		# Select plot parameters and print on terminal
		PltPars = DataFrame(select(df, Symbol.(Pars))[1,:])
		InfoPars = copy(PltPars)
		InfoPars[!,"x"] .= xVar
		InfoPars[!,"y"] .= yVar
		InfoPars[!,"p"] .= pVar
		@info "\e[1;36mScan plot $(j)/$(J)\e[0m" InfoPars
		println()

		yLabel::String = ""
		if yVar=="m"
			yLabel = "Magnetization"
		elseif Print
			yLabel = "\$" * VarLabels[yVar] * "\$"
		elseif !Print
			yLabel = yVar
		end

		# Initialize plot
		Fig = Figure(size=(600,400),figure_padding = 1)
		ax = Axis(Fig[1, 1])
		if Print
			ax.xlabel = L"$%$(VarLabels[xVar])$"
			ax.ylabel = L"$%$(VarLabels[yVar])$"
		elseif !Print
			ax.xlabel = xVar
			ax.ylabel = yVar
		end

		# Create filename
		FileName = join( ["$(Pars[i])=$(df[!,Pars[i]][1])" for i in 1:length(Pars)], '_' )
		FileName = yVar * "_" * FileName

		# Create raw title string, either for printing or local plotting
		rawTitle::String = yLabel * " ("
		if Print
			ParTitle = [VarLabels[Par] * "=$(PltPars[!,Par][1])" for Par in Pars]
			ParTitle = ["\$" * Par * "\$" for Par in ParTitle]
		elseif !Print
			ParTitle = [Par * "=$(PltPars[!,Par][1])" for Par in Pars]
		end
		rawTitle *= join(ParTitle, ", ") * ")"

		# Include RB specifications
		if !Print
			r = split(rawTitle, ")")
			rawTitle = r[1] * ", RB=$(RB...))"
		end

		# Handle infinities
		Print ? rawTitle = replace(rawTitle, "Inf" => "\\infty") : false
		if Print
			ax.title = L"%$(rawTitle)"
		elseif !Print
			ax.title = rawTitle
		end

		Groupeddf = groupby(df,pVar)[1:(Skip+1):end]
		I = length(Groupeddf)
		C = floor(Int64, length(colorschemes[cs]) / I)
		C==0.0 ? error("Your cs (ColorScheme) is not large enough.") : false

		for (i,pdf) in enumerate(Groupeddf)
			p = pdf[!,pVar][1]
			xx = pdf[!,xVar]
			yy = pdf[!,yVar]

			if Print
				label = L"$%$(VarLabels[pVar])=%$(p)$"
			elseif !Print
				label = pVar * "=$(p)"
			end

			# Plot parametrically
			# in(yVar, GetHFPs(Phase,Syms,RBS,RBd)) ? yy = abs.(yy) : false
			scatterlines!(
				ax, xx, yy,
				marker = :circle,
				color = colorschemes[cs][C*i],
				markersize = 8,
				label = label
			)
		end

		if Print
			LegendLabel = L"Simulated $%$(pVar)$:"
		elseif !Print
			LegendLabel = "Simulated " * pVar * ":"
		end

		XX = filter(!isnan,df[!,xVar])
		YY = filter(!isnan,df[!,yVar])
		MinX = minimum(XX)
		MaxX = maximum(XX)
		MinY = minimum(YY)
		MaxY = maximum(YY)

		Pad = 5e-2
		MinXLim = MinX - Pad*(MaxX-MinX)
		MaxXLim = MaxX + Pad*(MaxX-MinX)
		MinYLim = MinY - Pad*(MaxY-MinY)
		MaxYLim = MaxY + Pad*(MaxY-MinY)

		xlims!(ax, MinXLim, MaxXLim)
		ylims!(ax, MinYLim, MaxYLim)
		Fig[1, 2] = Legend(Fig, ax, LegendLabel, framevisible = false)
		push!(PlotVec, GroupedPlot(Fig,df,FileName))
	end

	return PlotVec
end

function SavePlot2D(
	FilePathIn::String,
	DirPathOut::String;
	xVar::String="V",
	yVar::String="m",
	pVar::String="δ",
	cs::Symbol=:imola50,
	Skip::Int64=0,
	Extension::String="pdf"
)

	# Assert printing
	Print::Bool=true
	PlotVec = Plot2D(FilePathIn;Print,xVar,yVar,pVar,cs,Skip)

	# Initialize directory structure
	Setup, Phase, Syms, RB, _ = UnpackFilePath(FilePathIn)
	DirPathOut *= "/xVar=" * xVar * "_pVar=" * pVar * "/"
	mkpath(DirPathOut)

	# Save each plot
	for GP in PlotVec
		FilePathOut = DirPathOut * "/" * GP.FileName * "." * Extension
		with_theme(theme_latexfonts()) do #TODO Redundant?
			save(FilePathOut,GP.H.scene)
		end
	end

	@info "\e[1;36mPlots saved at:\e[0m" DirPathOut

end