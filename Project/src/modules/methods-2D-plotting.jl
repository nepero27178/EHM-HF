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
	xScale=identity,
	legend::Bool=true,
	compared::Bool=false,
	cVar::String=""
)::Vector{GroupedPlot}

	# List vars and pars
	xVars::Vector{String} = ["t", "U", "V", "L", "δ", "β"]
	Pars::Vector{String} = filter(!in([xVar,pVar]),xVars)

	# Input safecheck
	!in(xVar, xVars) ? error("Invalid x variable, choose one of $(xVars)") : false
	!in(pVar, xVars) ? error("Invalid p variable, choose one of $(xVars)") : false
	xVar==pVar ? error("You have chosen xVar=pVar!") : false

	# Unpack filepath and set graphics
	Setup, Phase, Syms, RB, _ = UnpackFilePath(FilePathIn)
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

	# Load data
	DF::DataFrame = CSV.read(FilePathIn,DataFrame)
	Sim::Simulation = Simulation(DF,Setup,Phase,Syms,RB)
	DF = EnlargeDF!(Sim) # Compute RMPs
	yVars = filter(!in(xVars), names(DF))
	!in(yVar,yVars) ? error("Invalid y variable, choose one of $(yVars)") : false

	# Assess if compared plot is possible (requires 3/6 variables to be constant)
	ParsDF::DataFrame = unique(select(DF,filter(!=(cVar),Pars)))
	if compared && nrow(ParsDF) != 1
		@error "Impossible to compare plot: bad data @ Plot2D" cVar
		return
	end

	# Initialize y label
	yLabel::String = ""
	if yVar=="m"
		yLabel = "Magnetization"
	elseif Print
		yLabel = "\$" * VarLabels[yVar] * "\$"
	elseif !Print
		yLabel = yVar
	end

	# Initialize legend label
	if legend
		LegendLabel = Print ? L"Simulated $%$(VarLabels[pVar])$:" : "Simulated " * pVar * ":"
	end

	# Group data
	GroupedDF::GroupedDataFrame = groupby(DF,Pars)
	J::Int64 = length(GroupedDF)
	PlotVec::Vector{GroupedPlot} = GroupedPlot[]

	# Initialize layout for comparison
	if compared

		if !in(cVar,Pars)
			@error "Invalid cVar @ Plot2D" cVar Pars
			return
		elseif J>4
			@warn "Many plots in a single row @ Plot2D" J
		end

		# If true, single legend; otherwise separated legends
		UniqueLegend = legend ? allequal( [Set(unique(GDF[!,pVar])) for GDF in GroupedDF] ) : false

		# Figure and axes
		cFig::Figure = Figure(size=(1500,400))
		caxs::Vector{Axis} = [Axis(cFig[1,j]) for j in 1:J]
		linkyaxes!(caxs...)

		# Create title (c=compared)
		if Print
			cParTitle::Vector{String} = [VarLabels[Par] * "=$(only(ParsDF[!,Par]))" for Par in filter(!=(cVar),Pars)]
			cParTitle = ["\$" * Par * "\$" for Par in cParTitle]
		end
		crawTitle::String = Print ? yLabel * " (" * join(cParTitle, ", ") * ")" : "Compared variable: " * cVar

		# Create filename (c=compared)
		cFileNameVec::Vector{String} = [Par * "=$(Par != cVar ? only(ParsDF[!,Par]) : "Compared")" for Par in Pars]
		cFileName::String = yVar * "_" * join(cFileNameVec,'_')

	elseif !compared && !isempty(cVar)
		@warn "Perhaps you meant compared=true @ Plot2D" cVar compared
		return
	end

	# Cycle over simulated points
	for (j,df) in enumerate(GroupedDF)

		# Select plot parameters and print on terminal
		PltPars::DataFrame = unique( select(df, Symbol.(Pars)) )
		InfoPars = copy(PltPars)
		InfoPars[!,"x"] .= xVar
		InfoPars[!,"y"] .= yVar
		InfoPars[!,"p"] .= pVar
		@info "\e[1;36mScan plot $(j)/$(J)\e[0m" InfoPars
		println()

		# Initialize plot
		Fig = Figure(size=(600,400),figure_padding = 1)
		ax = Axis(Fig[1, 1])
		AxisVec::Vector{Axis} = [ax]
		compared ? push!(AxisVec, caxs[j]) : false
		if Print
			[a.xlabel = L"$%$(VarLabels[xVar])$" for a in AxisVec]
			[a.ylabel = L"$%$(VarLabels[yVar])$" for a in AxisVec]
		elseif !Print
			[a.xlabel = xVar for a in AxisVec]
			[a.ylabel = yVar for a in AxisVec]
		end

		# Create filename
		FileNameVec::Vector{String} = [Par * "=$(first(df[!,Par]))" for Par in Pars]
		FileName::String = yVar * "_" * join(FileNameVec,'_')

		# Create raw title string, either for printing or local plotting
		if Print
			ParTitle = [VarLabels[Par] * "=$( only(PltPars[!,Par]) )" for Par in Pars]
			ParTitle = ["\$" * Par * "\$" for Par in ParTitle]
			if compared
				rawSubTitle = "\$$(VarLabels[cVar])=$(only(PltPars[!,cVar]))\$"
			end
		elseif !Print
			ParTitle = [Par * "=$(only(PltPars[!,Par]))" for Par in Pars]
			cParTitle = ParTitle
		end
		rawTitle::String = yLabel * " (" * join(ParTitle, ", ") * ")"

		# Include RB specifications
		if !Print
			r = split(rawTitle, ")")
			rawTitle = r[1] * ", RB=$(RB...))"
		end

		# Handle infinities
		Print ? rawTitle = replace(rawTitle, "Inf" => "\\infty") : false
		if Print
			ax.title = L"%$(rawTitle)"
			compared ? caxs[j].title = L"%$(rawSubTitle)" : false
		elseif !Print
			[a.title = rawTitle for a in AxisVec]
		end

		Groupeddf::GroupedDataFrame = groupby(df,pVar)[1:(Skip+1):end]
		I::Int64 = length(Groupeddf)
		C::Int64 = floor(Int64, length(colorschemes[cs]) / I)
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
			[scatterlines!(
				a, xx, yy,
				marker = :circle,
				color = colorschemes[cs][C*i],
				markersize = 8,
				label = label
			) for a in AxisVec]
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

		if occursin("log", string(xScale))
			MinXLim = Pad
			if MinXLim >= MaxXLim
				@error "Impossible to plot these data on log scale @ Plot2D" MinXLim MaxXLim
				return
			end
		end
		[xlims!(a, MinXLim, MaxXLim) for a in AxisVec] #TODO
		ylims!(ax, MinYLim, MaxYLim)
		legend ? Fig[1, 2] = Legend(Fig, ax, LegendLabel, framevisible = false) : false

		[a.xscale = xScale for a in AxisVec] #TODO
		push!(PlotVec, GroupedPlot(Fig,df,FileName))
	end

	if compared
		for (i,cax) in enumerate(caxs)
			if i>1
				cax.ylabelvisible = false
				cax.yticklabelsvisible = false
			end
		end
		cDF::DataFrame = unique(select(DF,cVar))
		Label(cFig[0,:], Print ? L"%$(crawTitle)" : crawTitle)
		if UniqueLegend
			cFig[1, J+1] = Legend(cFig, caxs[1], LegendLabel, framevisible = false)
		end

		push!(PlotVec, GroupedPlot(cFig,view(cDF,:,:),cFileName))
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
	xScale=identity,
	legend::Bool=true,
	compared::Bool=false,
	cVar::String="",
	Extension::String="pdf",
)

	# Assert printing
	Print::Bool=true
	PlotVec = Plot2D(FilePathIn;Print,xVar,yVar,pVar,cs,Skip,xScale,legend,compared,cVar)

	# Initialize directory structure
	Setup, Phase, Syms, RB, _ = UnpackFilePath(FilePathIn)
	DirPathOut *= "/xVar=" * xVar * "_pVar=" * pVar * "/"
	mkpath(DirPathOut)

	# Save each plot
	for GP in PlotVec
		FilePathOut = DirPathOut * "/" * GP.FileName * "." * Extension
		with_theme(theme_latexfonts()) do
			save(FilePathOut,GP.H.scene)
		end
	end

	@info "\e[1;36mPlots saved at:\e[0m" DirPathOut

end