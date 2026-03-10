#!/usr/bin/julia

using GLMakie
using CairoMakie
using LaTeXStrings
using ColorSchemes
using DataFrames
using CSV

PROJECT_METHODS_DIR = @__DIR__
include(PROJECT_METHODS_DIR * "/structs.jl")
include(PROJECT_METHODS_DIR * "/methods-IO.jl")
include(PROJECT_METHODS_DIR * "/methods-physics.jl")

function Plot3D(
	FilePathIn::String;
	Print::Bool=false,
	Mode::String="heatmap",
	xVar::String="U",
	yVar::String="V",
	zVar::String="f",
	cs::Symbol=:imola50,
	azm=false
)::Vector{GroupedPlot}

	if Mode=="heatmap" && isa(azm,Float64)
		@warn "Perhaps you wanted a surface plot?" Mode azm
	end

	# List vars and pars
	xyVars::Vector{String} = ["t", "U", "V", "L", "δ", "β"]
	Pars::Vector{String} = filter(!in([xVar,yVar]),xyVars)

	# Input safecheck
	!in(xVar, xyVars) ? error("Invalid x variable, choose one of $(xyVars)") : false
	!in(yVar, xyVars) ? error("Invalid y variable, choose one of $(xyVars)") : false
	xVar==yVar ? error("You have chosen xVar=yVar!") : false

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
	EnlargeDF!(Sim) # Compute RMPs
	zVars = filter(!in(xyVars), names(DF))
	!in(zVar,zVars) ? error("Invalid z variable, choose one of $(zVars)") : false

	# Initialize z label
	zLabel::String = ""
	if zVar=="m"
		zLabel = "Magnetization"
	elseif Print
		zLabel = "\$" * VarLabels[zVar] * "\$"
	elseif !Print
		zLabel = zVar
	end

	# Group data
	GroupedDF = groupby(DF,Pars)
	J::Int64 = length(GroupedDF)
	PlotVec::Vector{GroupedPlot} = GroupedPlot[]

	# Cycle over simulated points
	for (j,df) in enumerate(GroupedDF)

		# Select plot parameters and print on terminal
		PltPars::DataFrame = unique( select(df, Symbol.(Pars)) )
		InfoPars = copy(PltPars)
		InfoPars[!,"x"] .= xVar
		InfoPars[!,"y"] .= yVar
		InfoPars[!,"z"] .= zVar
		@info "\e[1;36m3D $(Mode) plot $(j)/$(J)\e[0m" InfoPars
		println()

		# Reshape local data
		xx, yy, zz = ReshapeData(DataFrame(df);xVar,yVar,zVar)

		# Initialize plot
		sz::Tuple{Float64,Float64} = (600,400)
		if Mode=="surface"
			sz = (600,500)
		end
		H::Figure = Figure(size = sz,figure_padding = 1)
		ax = nothing

		if Mode=="heatmap"
			ax = Axis(
				H[1, 1]
			)
		elseif Mode=="surface"
			ax = Axis3(
				H[1, 1],
				xlabelalign = (:center, :center),
				ylabelalign = (:center, :center),
				zlabelalign = (:center, :center),
				xlabelrotation = 0, # Horizontal xlabel
				ylabelrotation = 0, # Horizontal ylabel
				zlabelrotation = 0, # Horizontal zlabel
				azimuth=azm
			)
		else
			@error "Invalid Mode @ Plot3D" Mode
			return
		end

		if Print
			ax.xlabel = L"$%$(VarLabels[xVar])$"
			ax.ylabel = L"$%$(VarLabels[yVar])$"
			Mode=="surface" ? ax.zlabel = L"$%$(VarLabels[zVar])$" : false
		elseif !Print
			ax.xlabel = xVar
			ax.ylabel = yVar
			Mode=="surface" ? ax.zlabel = zVar : false
		end

		# Create filename
		FileName = zVar * "_" * join( ["$(Pars[i])=$(df[!,Pars[i]][1])" for i in 1:length(Pars)], '_' )
		Mode=="surface" ? FileName *= "_SurfacePlot" : false

		# Create raw title string, either for printing or local plotting
		rawTitle::String = zLabel * " ("
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

		# Plot parametrically
		clims::Tuple{Float64,Float64} = (
			minimum( vcat(0.0,filter(!isnan,zz)) ),
			maximum( vcat(0.0,filter(!isnan,zz)) )
		)
		if zVar=="tS"
			clims = (0.0,1.0)
		end

		h = nothing
		if Mode=="heatmap"
			h = CairoMakie.heatmap!(
				ax,
				xx, yy, zz,
				colormap=colorschemes[cs],
				colorrange=clims
			)
			Colorbar(H[1,2], h)
		elseif Mode=="surface"
			h = CairoMakie.surface!(
				ax,
				xx, yy, zz,
				colormap=colorschemes[cs],
				shading=false
			)
			w = CairoMakie.wireframe!(
				ax,
				xx, yy, zz,
				color=tabblue,
				linewidth=0.1
			)
			s = CairoMakie.scatter!(
				ax,
				xx, yy, zz,
				color=:black,
				markersize=1.5
			)
		else
			@error "Invalid Mode @ Plot3D" Mode
			return
		end

		push!(PlotVec, GroupedPlot(H,df,FileName))

	end

	return PlotVec
end

function SavePlot3D(
	FilePathIn::String,
	DirPathOut::String;
	Mode::String="heatmap",
	xVar::String="U",
	yVar::String="V",
	zVar::String="f",
	cs::Symbol=:imola50,
	azm=false,
	Extension::String="pdf"
)

	# Assert printing
	Print::Bool=true
	PlotVec = Plot3D(FilePathIn;Print,Mode,xVar,yVar,zVar,cs,azm)

	# Initialize directory structure
	Setup, Phase, Syms, RB, _ = UnpackFilePath(FilePathIn)
	DirPathOut *= "/xVar=" * xVar * "_yVar=" * yVar * "/"
	mkpath(DirPathOut)

	# Save each plot
	for GP in PlotVec
		FilePathOut = DirPathOut * "/" * GP.FileName * "." * Extension
		with_theme(theme_latexfonts()) do #TODO Redundant?
			save(FilePathOut,GP.H)
		end
	end

	@info "\e[1;36mPlots saved at:\e[0m" DirPathOut

end