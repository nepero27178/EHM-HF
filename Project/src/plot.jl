#!/usr/bin/julia

using GLMakie
using CairoMakie
using LaTeXStrings
using ColorSchemes
using DataFrames
using CSV

# Includer
PROJECT_SRC_DIR = @__DIR__
include(PROJECT_SRC_DIR * "/setup/graphic-setup.jl")
include(PROJECT_SRC_DIR * "/modules/structs.jl")
include(PROJECT_SRC_DIR * "/modules/methods-IO.jl")

# Arguments handler
if length(ARGS) != 2
	println("How to use this program?
Type the following: \$ julia plot.jl --mode --obj
Where:
· mode = hs / rs
· obj = HFPs / RMPs / Qs / phys")
	exit()
else
	UserInput = ARGS
	Mode::String = UserInput[1][3:end]
	Obj::String = UserInput[2][3:end]
end

# Process mode
if in(Mode, ["rs","hs"])
	include(PROJECT_SRC_DIR * "/setup/" * Mode * "-setup.jl")
	Mode=="hs" ? Dim = "2D" : Dim = "3D"
	include(PROJECT_SRC_DIR * "/modules/methods-" * Dim * "-plotting.jl")
else
	@error "Invalid mode. Use: mode = hs / rs"
	exit()
end

# Process obj
if Obj=="HFPs"
	objList = GetHFPs(Phase,Syms,RBS,RBd)
elseif Obj=="RMPs"
	objList = []
	RBS ? push!(objList,"tS") : false
	RBd ? push!(objList,"td") : false
elseif Obj=="Qs"
	QsList = ["Q$(HFP)" for HFP in GetHFPs(Phase,Syms,RBS,RBd)]
	RunList = ["ΔT", "I", "g0", "g"]
	objList = vcat(QsList, RunList)
elseif Obj=="phys"
	objList = ["μ", "f"]
else
	@error "Invalid obj. Use obj = HFPs / RMPs / Qs / phys"
	exit()
end

function main()
	Level::String = "raw"
	if Level=="refined"
		global Setup = split(Setup,'-')[1]
	end

	# Read files
	FilePathIn::String = dirname(PROJECT_SRC_DIR) * "/data/$(Level)/Mode=$(Mode)/Setup=$(Setup)/Phase=$(Phase)/RB=$(RB...)_Syms=$(Syms...).csv"

	# Create output directory
	# For simulations: Setup > Phase > RB+Syms (to make comparable data in the same folder)
	# For plots: Phase > Setup > Syms > RB (to make same-phase plots in the same folder)
	DirPathOut::String = dirname(PROJECT_SRC_DIR) * "/plots/$(Level)/Mode=$(Mode)/Phase=$(Phase)/Setup=$(Setup)/Syms=$(Syms...)/RB=$(RB...)/Obj=$(Obj)"
	mkpath(DirPathOut)

	JumpToData::String = "#!/usr/bin/bash\n\ncd $(dirname(FilePathIn))"
	open(dirname(DirPathOut) * "/jump-to-data.sh","w") do io
		write(io,JumpToData)
	end

	for obj in objList
		# Run plot modules for each HFP
		if Mode=="hs"
			SavePlot2D(
				FilePathIn,
				DirPathOut;
				xVar="V",
				yVar=obj,
				pVar="δ",
				cs=:tabwarm,
			)
			SavePlot2D(
				FilePathIn,
				DirPathOut;
				xVar="δ",
				yVar=obj,
				pVar="V",
				cs=:tabwarm,
				Skip=4
			)
		elseif Mode=="rs"
			SavePlot3D(
				FilePathIn,
				DirPathOut;
				xVar,
				yVar,
				zVar=obj,
				cs=:imola,
				# Extension="png"
			)
		end
	end
end

if abspath(PROGRAM_FILE) == @__FILE__
	main()
end
