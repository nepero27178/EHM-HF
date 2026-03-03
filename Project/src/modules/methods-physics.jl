function GetK(
	L::Vector{Int64}
)::Tuple{Matrix{Vector{Float64}},Vector{Float64},Vector{Float64}}

	# Reciprocal space discretization (normalized to 1)
	Kx::Vector{Float64} = [kx for kx in -pi:2*pi/L[1]:pi]
	popfirst!(Kx)
	Ky::Vector{Float64} = [ky for ky in -pi:2*pi/L[2]:pi]
	popfirst!(Ky)

	return [ [kx,ky] for kx in Kx, ky in Ky ], Kx, Ky
end

function StructureFactor(
	Sym::String,
	k::Vector{Float64}
)::Float64

	AllSyms = ["s", "S", "x", "y", "d"]
	if !in(Sym, AllSyms)
		@error "Invalid Sym @ StructureFactor" Sym
		return
	end
	kx, ky = k

	if Sym=="s"
		return 1
	elseif Sym=="S"
		return cos(kx) + cos(ky)
	elseif Sym=="x"
		return sqrt(2) * sin(kx)
	elseif Sym=="y"
		return sqrt(2) * sin(ky)
	elseif Sym=="d"
		return cos(kx) - cos(ky)
	end
end

function FermiDirac(
	ε::Float64,							# Single-particle energy
	μ::Float64,							# Chemical potential
	β::Float64,							# Inverse temperature
)::Float64

	if β==Inf # Zero temperature distribution
		if ε<μ
			return 1
		elseif ε==μ
			return 0.5
		elseif ε>μ
			return 0
		end
	elseif β<Inf	 # Finite temperature distribution
		return 1 / ( exp(β * (ε-μ)) + 1 )
	end
end

function Th(
	Sign::String,
	ε::Float64,							# Single-particle energy
	μ::Float64,							# Chemical potential
	β::Float64,							# Inverse temperature
)::Float64

	if Sign=="+"
		return FermiDirac(-ε,μ,β) + FermiDirac(ε,μ,β)
	elseif Sign=="-"
		return FermiDirac(-ε,μ,β) - FermiDirac(ε,μ,β)
	else
		@error "Invalid Sign @ Th" Sign
		return
	end

end

function GetΔΔ(
	Phase::String,
	Syms::Set{String},
	Pars::DataFrame,
	v::DataFrame,
)::Tuple{Dict{String,Float64},Dict{String,Float64}}

	U = first(Pars.U)
	V = first(Pars.V)

	reΔΔ::Dict{String,Float64} = Dict()
	imΔΔ::Dict{String,Float64} = Dict()

	if Phase == "AF-Symmetric"

		# Real gap
		reΔΔ["s"] = U * first(v.m)

		# Imaginary gap
		vS::Float64 = Readv(v,:vS;Cnd="S" in Syms)
		vd::Float64 = Readv(v,:vd;Cnd="d" in Syms)

		imΔΔ["S"] = -V * vS
		imΔΔ["d"] = -V * vd

	elseif Phase == "AF-Antisymmetric"

		# Real gap
		reΔΔ["s"] = try
			U * first(v.m)
		catch
			@error "Magnetization not found in v @ GetΔΔ" v
			exit()
		end

		# Imaginary gap
		vx::Float64 = Readv(v,:vx;Cnd="x" in Syms)
		vy::Float64 = Readv(v,:vy;Cnd="y" in Syms)

		imΔΔ["x"] = -V * vx
		imΔΔ["y"] = -V * vy

	elseif Phase == "SC-Singlet"

		# Real gap
		ws::Float64 = Readv(v,:ws;Cnd="s" in Syms)
		wS::Float64 = Readv(v,:wS;Cnd="S" in Syms)
		wd::Float64 = Readv(v,:wd;Cnd="d" in Syms)

		reΔΔ["s"] = -U * ws
		reΔΔ["S"] = V * wS
		reΔΔ["d"] = V * wd

	elseif Phase == "SC-Triplet"

		# Imaginary gap
		wx::Float64 = Readv(v,:wx;Cnd="x" in Syms)
		wy::Float64 = Readv(v,:wy;Cnd="y" in Syms)

		imΔΔ["x"] = V * wx
		imΔΔ["y"] = V * wy

	elseif Phase != "Normal"

		@error "Invalid Phase @ GetΔΔ" Phase
		exit()

	end

	return reΔΔ, imΔΔ

end

function GetObj(
	Obj::String,
	Phase::String,
	Syms::Set{String},
	Pars::DataFrame,
	v::DataFrame,
	μ::Float64;
	RBS::Bool=true,
	RBd::Bool=false,
	OptBZ::Bool=true,
	debug::Bool=false
)::Float64

	# Get pars
	t::Float64 = first(Pars.t)
	U::Float64 = first(Pars.U)
	V::Float64 = first(Pars.V)
	L::Int64 = first(Pars.L)
	β::Float64 = first(Pars.β)

	# Get BZ
	K::Matrix{Vector{Float64}}, _, _ = GetK([L, L])
	LxLy::Int64 = L^2

	# Kinetics
	uS::Float64 = Readv(v,:uS;Cnd=RBS) # Read s*-wave
	ud::Float64 = Readv(v,:ud;Cnd=RBd) # Read d-wave
	εε::Dict{String,Float64} = Dict(
		"S" => -2*t + V*uS,
		"d" => V*ud
	)
	Getεk(k::Vector{Float64})::Float64 = GetFk(εε,k) # Function
	εK::Matrix{Float64} = Getεk.(K) # Get bare bands
	EK::Matrix{Float64} = zeros(size(εK)) # Initialize quasibands
	esK::Matrix{Float64} = zeros(size(εK))

	# Gap
	reΔΔ, imΔΔ = GetΔΔ(Phase,Syms,Pars,v) # Get real and imaginary dicts
	GetreΔk(k::Vector{Float64})::Float64 = GetFk(reΔΔ,k) # Function
	GetimΔk(k::Vector{Float64})::Float64 = GetFk(imΔΔ,k) # Function
	reΔK::Matrix{Float64} = GetreΔk.(K) # Real gap function
	imΔK::Matrix{Float64} = GetimΔk.(K) # Real gap function

	# Initializers
	f::Float64 = 0.0
	NK::Matrix{Float64} = zeros(size(εK))
	if Obj=="f"
		vu = RBS || RBd ? [ select( v, Cols(contains.("u")) )[1,:]... ] : false
		n::Float64 = 0.5 + first(Pars.δ)
		f += μ*n - V*sum(vu.^2)
	end

	# Custom xlogx function to handle numeric Inf
	xlogx(x::Float64) = x==0.0 ? 0.0 : x*log(x)

	# Normal phase
	if Phase=="Normal"

		# Density
		if Obj=="n"
			# Ignore spin degeneracy: it leads to *2 /2
			NK = FermiDirac.(εK,μ,β)
			return sum(NK)/LxLy

		# Free energy
		elseif Obj=="f" # Use explicit form to correct for numeric errors
			# Free energy from bands
			FK::Matrix{Float64} = FermiDirac.(εK,μ,β)
			f += sum((εK.-μ).*FK)/LxLy

			# Free energy from entropy
			esK = xlogx.(FK) + xlogx.(1 .- FK)
			f += sum(esK)/LxLy * 2/β
			return f

		end

	# Antiferromagnetic phases
	elseif Phase=="AF-Symmetric" || Phase=="AF-Antisymmetric"
		EK = sqrt.(εK.^2 .+ reΔK.^2 .+ imΔK.^2)
		m::Float64 = try
			first(v.m)
		catch
			@error "Magnetization not found in v @ GetObj" v
			exit()
		end

		# Density
		if Obj=="n"
			NK = Th.("+",EK,μ,β)
			return sum(NK)/(2*LxLy)

		# Free energy
		elseif Obj=="f" # Use explicit form to correct for numeric errors
			# Free energy from HFPs
			vv::Vector{Float64} = [ select(v, Cols(contains.("v")))[1,:]... ]
			f += U*m^2 - V*sum(vv.^2)

			# Free energy from bands
			FpK::Matrix{Float64} = FermiDirac.(EK,μ,β)
			FmK::Matrix{Float64} = FermiDirac.(-EK,μ,β)
			f += sum((EK.-μ).*FpK - (EK.+μ).*FmK)/LxLy # MBZ to BZ => /2, spin => *2

			# Free energy from entropy
			esK = xlogx.(FpK) + xlogx.(1 .- FpK) + xlogx.(FmK) + xlogx.(1 .- FmK)
			f += sum(esK)/LxLy * 1/β # MBZ to BZ => /2, spin => *2
			return f

		end

	# Superconducting phases
	elseif Phase=="SC-Singlet" || Phase=="SC-Triplet"
		ξK::Matrix{Float64} = εK .- μ
		EK = sqrt.(ξK.^2 .+ reΔK.^2)

		# Density
		if Obj=="n"
			tK::Matrix{Float64} = tanh.(EK .* β/2) ./ EK
			replace!(tK, NaN => β/2) # @ x~0 : tanh(x)/x~1
			NK = 1 .- ξK.*tK
			return sum(NK)/(2*LxLy)

		# Free energy
		elseif Obj=="f" # Use equivalent form since numeric errors are avoided by construction
			ws::Float64 = Readv(v,:ws;Cnd="s" in Syms)
			vw = [ select(select(v, Cols(contains("w"))), Cols(!contains("S")))[1,:]... ]
			f += -U*ws^2 + V*sum(vw.^2)
			esK = log.( 1 .- FermiDirac.(EK,0.0,β) ) # Entropic part
			ecK::Matrix{Float64} = ξK .- EK # Contraction part
			f += ( sum(esK)*2/β + sum(ecK) )/LxLy
			return f

		end

	else
		@error "Invalid Syms @ GetObj" Syms
		exit
	end

end

function GetDensity(
	Phase::String,
	Syms::Set{String},
	Pars::DataFrame,
	v::DataFrame,
	μ::Float64;
	RBS::Bool=true,
	RBd::Bool=false,
	OptBZ::Bool=true,
	debug::Bool=false
)::Float64

	return GetObj("n",Phase,Syms,Pars,v,μ;RBS,RBd,OptBZ,debug)
end

function GetFreeEnergy(
	Phase::String,
	Syms::Set{String},
	Pars::DataFrame,
	v::DataFrame,
	μ::Float64;
	RBS::Bool=true,
	RBd::Bool=false,
	OptBZ::Bool=true,
	debug::Bool=false
)::Float64

	return GetObj("f",Phase,Syms,Pars,v,μ;RBS,RBd,OptBZ,debug)
end