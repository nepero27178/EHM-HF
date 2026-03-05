#!/usr/bin/julia
using CSV
using Dates
using DataFrames
using Term

# Arguments handler
if length(ARGS) != 1
	println("How to use this program?
Type the following: \$ julia simulate.jl --mode
Where:
· mode = hs / rs")
	exit()
else
	UserInput = ARGS
	Mode = UserInput[1][3:end]
end

# Includer
PROJECT_SRC_DIR = @__DIR__
if in(Mode, ["rs", "hs"])
	include(PROJECT_SRC_DIR * "/setup/" * Mode * "-setup.jl")
else
	@error "Invalid argument. Use: mode = hs / rs"
end
include(PROJECT_SRC_DIR * "/modules/structs.jl")
include(PROJECT_SRC_DIR * "/modules/methods-simulating.jl")
include(PROJECT_SRC_DIR * "/modules/methods-physics.jl")
include(PROJECT_SRC_DIR * "/modules/methods-optimizations.jl")
include(PROJECT_SRC_DIR * "/modules/methods-IO.jl")

# Main run
function main()

	# Create list of data to be imported
	FilePathsIn::Vector{String} = []

	# Get setup class (e.g. A[128]-a, A[128]-b, A[128] => A[128])
	SetupClass::String = split(Setup,'-')[1]

	# Create output directory
	# For data: Setup > Phase > Syms (to make comparable data in the same folder)
	# For plots: Phase > Setup > Syms (to make same-phase plots in the same folder)
	DirPathOut::String = dirname(PROJECT_SRC_DIR) * "/data/refined/Mode=$(Mode)/Setup=$(SetupClass)/Phase=$(Phase)"

	# Cycle over setups in the same class
	for Setup in filter(contains(SetupClass), AvailableSetups)

		# Look the layered data
		DirPathIn::String = dirname(PROJECT_SRC_DIR) * "/data/layered/Mode=$(Mode)/Setup=$(Setup)/Phase=$(Phase)"

		# Read present layers and select the largest
		l = try
			Layers::Vector{Int64} = [U[end] for U in UnpackFilePath.(
				DirPathIn .* "/" .* readdir(DirPathIn); Layered=true
			)]
			maximum(Layers)
		catch
			0
		end

		# Generate FilePathIn and FilePathOut
		FilePathIn = nothing
		if l==0
			FilePathIn = replace(DirPathIn, "layered" => "raw") * "/RB=$(RB...)_Syms=$(Syms...).csv"
		elseif l>0
			FilePathIn = DirPathIn * "/RB=$(RB...)_Syms=$(Syms...)_Layer=$(l).csv"
		end
		push!(FilePathsIn,FilePathIn)
	end

	# Merge data and write on file
	S::Simulation = MergeData(FilePathsIn)
	mkpath(DirPathOut)
	FilePathOut = DirPathOut * "/RB=$(S.RB...)_Syms=$(S.Syms...).csv"
	CSV.write(FilePathOut, S.DF)

	LogPathOut::String = replace(FilePathOut, ".csv" => ".log")
	Log::DataFrame = DataFrame(Dict("FilePathsIn" => FilePathsIn))
	CSV.write(LogPathOut,Log)
end

if abspath(PROGRAM_FILE) == @__FILE__
	main()
end
