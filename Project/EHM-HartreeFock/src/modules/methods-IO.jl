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

function UnpackFilePath(
	FilePathIn::String
)::Tuple{String,String,Set{String},Set{String}}

	Str::String = replace(FilePathIn, ['.','_']=>'/')
	GetVal(x::String)::String = String(split(split(Str,"$(x)=")[2],"/")[1])

	Setup::String = GetVal("Setup")
	Phase::String = GetVal("Phase")
	Syms::Set{String} = Set([string(s) for s in GetVal("Syms")])
	RB::Set{String} = Set([string(s) for s in GetVal("RB")])

	return Setup,Phase,Syms,RB
end

function ReshapeData(
	DF::DataFrame;
	xVar::String="U",
	yVar::String="V",
	zVar::String="fMFT"
)::Tuple{Any, Any, Any}
	
	xx = unique(DF[!,xVar])
	yy = unique(DF[!,yVar])
	zz = reshape(DF[!,zVar],length(yy),length(xx))
	
	return xx, yy, zz
end

function EnlargeDF!(
	Sim::Simulation
)::DataFrame

	DF::DataFrame = Sim.DF
	# HFPs::Set{String} = GetHFPs(Sim.Phase,Sim.Syms,RBS,RBd)
	# Consider the possibility of enlarging the DataFrame to left-out HFPs

	if "S" in Sim.RB
		DF.tS .= DF.t .- DF.V/2 .* DF.uS
	end

	if "d" in Sim.RB
		DF.td .= -DF.V/2 .* DF.ud
	end

	return DF

end

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
		"fMFT" => "f_\\mathrm{MFT}",
		# Shared HFPs
		"uS" => "u^{(s^*)}",
		"ud" => "u^{(d)}",
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
			"Qm" => "Q_m",
			"QvS" => "Q_{v^{(s^*)}}",
			"Qvd" => "Q_{v^{(d)}}",
		])
	elseif Phase=="AF-Antisymmetric"
		PhaseLabels = Dict([
			# HFPs
			"m" => "m",
			"vx" => "v^{(p_x)}",
			"vy" => "v^{(p_y)}",
			# Convergence
			"Qm" => "Q_m",
			"Qvx" => "Q_{v^{(p_x)}}",
			"Qvy" => "Q_{v^{(p_y)}}",
		])
	elseif Phase=="SC-Singlet"
		PhaseLabels = Dict([
			# HFPs
			"ws" => "w^{(s)}",
			"wS" => "w^{(s^*)}",
			"wd" => "w^{(d)}",
			# Convergence
			"Qws" => "Q_{w^{(s)}}",
			"QwS" => "Q_{w^{(s^*)}}",
			"Qwd" => "Q_{w^{(d)}}",
		])
	elseif Phase=="SC-Triplet"
		PhaseLabels = Dict([
			# HFPs
			"wx" => "w^{(p_x)}",
			"wy" => "w^{(p_y)}",
			# Convergence
			"Qwx" => "Q_{w^{(p_x)}}",
			"Qwy" => "Q_{w^{(p_y)}}",
		])

	end

	return merge(VarLabels, PhaseLabels)

end