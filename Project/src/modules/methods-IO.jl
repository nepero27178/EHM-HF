# General IO

function GetFk(
	FF::Dict{String,Float64},
	k::Vector{Float64};
)::Float64

	Fk::Float64 = 0.0
	for (Sym,FSym) in FF
		Fk += FSym * StructureFactor(Sym,k)
	end

	return Fk
end

function Readv(
	v::DataFrame,
	x::Symbol;
	Cnd::Bool=true
)::Float64

	x::Float64 = try v[!,x]
		Cnd ? first(v[!,x]) : 0.0
	catch
		0.0
	end

	return x

end

# Files IO

function UnpackFilePath(
	FilePathIn::String;
	Layered::Bool=false
)::Tuple{String,String,Set{String},Set{String},Int64}

	# Handle user forgetting to specify Layered
	if !Layered
		Layered = occursin("Layer=", FilePathIn) ? true : false
	end

	Str::String = replace(FilePathIn, ['.','_']=>'/')
	GetVal(x::String)::String = String(split(split(Str,"$(x)=")[2],"/")[1])

	Setup::String = GetVal("Setup")
	Phase::String = GetVal("Phase")
	Syms::Set{String} = Set([string(s) for s in GetVal("Syms")])
	RB::Set{String} = Set([string(s) for s in GetVal("RB")])
	Layer::Int64 = Layered ? parse(Int64,GetVal("Layer")) : 0

	return Setup,Phase,Syms,RB,Layer
end

function ReshapeData(
	DF::DataFrame;
	xVar::String="U",
	yVar::String="V",
	zVar::String="f"
)::Tuple{Any, Any, Any}

	# Sort DataFrame to avoid inconsistencies
	SDF::DataFrame = sort(DF, [xVar, yVar])

	xx = unique(SDF[!,xVar])
	yy = unique(SDF[!,yVar])
	zz = reshape(SDF[!,zVar],length(yy),length(xx))
	
	return xx, yy, zz'
end

function EnlargeDF!(
	Sim::Simulation
)::DataFrame

	DF::DataFrame = Sim.DF
	RBS::Bool = "S" in Sim.RB
	RBd::Bool = "d" in Sim.RB

	HFPs::Set{String} = GetHFPs(Sim.Phase,Sim.Syms,RBS,RBd)
	#TODO Consider the possibility of enlarging the DataFrame to left-out HFPs

	RBS ? DF.tS .= DF.t .- DF.V/2 .* DF.uS : false
	RBd ? DF.td .= -DF.V/2 .* DF.ud : false
	DF.s .= fill(0.0,size(DF,1))
	for Row in eachrow(DF)
		Pars::DataFrame = select(DataFrame(Row), [:L, :U, :V, :t, :β, :δ])
		v::DataFrame = select(DataFrame(Row), HFPs...)
		μ::Float64 = Row.μ
		Row.s = GetEntropy(Sim.Phase,Sim.Syms,Pars,v,μ;RBS,RBd)
	end

	return DF

end

function MergeData(
	FilePathsIn::Vector{String}
)::Simulation

	UDF::DataFrame = DataFrame(UnpackFilePath.(FilePathsIn), ["Setup", "Phase", "Syms", "RB", "Layer"])
	Uniform::Bool = allequal(UDF.Phase) && allequal(UDF.Syms) && allequal(UDF.RB) && allequal([s[1] for s in split.(UDF.Setup,'-')])
	if !Uniform
		@error "Non uniform FilePathsIn @ MergeData" UDF
		return

	elseif Uniform
		SetupClass = split(first(UDF.Setup),"-")[1]
		Phase = first(UDF.Phase)
		Syms = first(UDF.Syms)
		RB = first(UDF.RB)

		RBS::Bool = "S" in RB ? true : false
		RBd::Bool = "d" in RB ? true : false

		HFPs::Set{String} = GetHFPs(Phase,Syms,RBS,RBd) # Get phase HFPs
		QNames::Vector{String} = ["Q" * HFP for HFP in HFPs]

		DFs = CSV.read.(FilePathsIn, DataFrame)
		for i in 1:length(DFs) # Add layer information
			DFs[i].Layer .= fill(UDF.Layer[i], size(DFs[i],1))
		end

		DF = vcat(DFs...)

		# Find duplicated simulations => dDF = duplicates DataFrame
		ModPars::Vector{String} = ["L", "U", "V", "t", "β", "δ"]
		dDF = DF[ findall(nonunique(DF, ModPars)),: ][!,ModPars]

		# Check each duplicate
		for Row in eachrow(dDF)

			# Find repetitions and filter out non-converged values
			ii = findall( r==Row for r in eachrow(DF[!,ModPars]) )
			fDF = filter( :Converged => x -> x,DF[ii,:] )

			# Filter based on total quality (arbitrary)
			fDF.TotalQ .= sum([ fDF[!,Q] for Q in QNames ])
			_, j = findmin(fDF.TotalQ)
			fDF = select( DataFrame(fDF[j,:]), Not(:TotalQ) )

			# Map on the original dataframe and keep only the best
			i = findall( r==fDF[1,:] for r in eachrow(DF) )
			deleteat!(DF, filter( !=(only(i)),ii ))
		end

		return Simulation(DF,SetupClass,Phase,Syms,RB)
	end

end

# Physics IO

function GetHFPs(
	Phase::String,
	Syms::Set{String},
	RBS::Bool,
	RBd::Bool
)::Set{String}

	HFPs::Set{String} = Set()
	RBS ? push!(HFPs,"uS") : false
	RBd ? push!(HFPs,"ud") : false
	if Phase=="AF-Symmetric" || Phase=="AF-Antisymmetric"
		push!(HFPs,"m")

		if Phase=="AF-Symmetric" && !issubset(Syms,["S","d"])
			@error "Inconsistent Phase/Syms @ GetHFPs" Phase Syms
			exit()
		elseif Phase=="AF-Antisymmetric" && !issubset(Syms,["x","y"])
			@error "Inconsistent Phase/Syms @ GetHFPs" Phase Syms
			exit()
		end

		for Sym in Syms
			push!(HFPs,"v"*Sym)
		end

	elseif Phase=="SC-Singlet" || Phase=="SC-Triplet"

		if Phase=="SC-Singlet" && !issubset(Syms,["s","S","d"])
			@error "Inconsistent Phase/Syms @ GetHFPs" Phase Syms
			exit()
		elseif Phase=="SC-Triplet" && !issubset(Syms,["x","y"])
			@error "Inconsistent Phase/Syms @ GetHFPs" Phase Syms
			exit()
		end

		for Sym in Syms
			push!(HFPs,"w"*Sym)
		end

	end

	return HFPs

end


function GetLabels(
	Phase::String
)::Dict{String,String}

	VarLabels::Dict{String,String} = Dict([
		# Variables
		"t" => "t",
		"U" => "U",
		"V" => "V",
		"L" => "L",
		"δ" => "\\delta",
		"β" => "\\beta",
		# Other
		"ΔT" => "\\Delta T",
		"I" => "\\text{Total steps}",
		"μ" => "\\mu",
		"g0" => "g_0",
		"g" => "g",
		"f" => "f_\\mathrm{MFT}",
		"s" => "s_\\mathrm{MFT}",
		# Shared HFPs
		"uS" => "u^{(s^*)}",
		"ud" => "u^{(d)}",
		# Convergence
		"QuS" => "Q[u^{(s^*)}]",
		"Qud" => "Q[u^{(d)}]",
		# RMPs
		"tS" => "\\tilde{t}^{(s^*)}",
		"td" => "\\tilde{t}^{(d)}"
	])

	PhaseLabels::Dict{String,String} = Dict()
	if Phase=="AF-Symmetric"
		PhaseLabels = Dict([
			# HFPs
			"m" => "m",
			"vS" => "v^{(s^*)}",
			"vd" => "v^{(d)}",
			# Convergence
			"Qm" => "Q[m]",
			"QvS" => "Q[v^{(s^*)}]",
			"Qvd" => "Q[v^{(d)}]",
		])
	elseif Phase=="AF-Antisymmetric"
		PhaseLabels = Dict([
			# HFPs
			"m" => "m",
			"vx" => "v^{(p_x)}",
			"vy" => "v^{(p_y)}",
			# Convergence
			"Qm" => "Q[m]",
			"Qvx" => "Q[v^{(p_x)}]",
			"Qvy" => "Q[v^{(p_y)}]",
		])
	elseif Phase=="SC-Singlet"
		PhaseLabels = Dict([
			# HFPs
			"ws" => "w^{(s)}",
			"wS" => "w^{(s^*)}",
			"wd" => "w^{(d)}",
			# Convergence
			"Qws" => "Q[w^{(s)}]",
			"QwS" => "Q[w^{(s^*)}]",
			"Qwd" => "Q[w^{(d)}]",
		])
	elseif Phase=="SC-Triplet"
		PhaseLabels = Dict([
			# HFPs
			"wx" => "w^{(p_x)}",
			"wy" => "w^{(p_y)}",
			# Convergence
			"Qwx" => "Q[w^{(p_x)}]",
			"Qwy" => "Q[w^{(p_y)}]",
		])

	end

	return merge(VarLabels, PhaseLabels)

end