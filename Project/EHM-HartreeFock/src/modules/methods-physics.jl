function GetKGrid(
	L::Vector{Int64}
)::Tuple{Matrix{Vector{Float64}},Vector{Float64},Vector{Float64}}

	# Reciprocal space discretization (normalized to 1)
	Kx::Vector{Float64} = [kx for kx in -1:2/L[1]:1]
	popfirst!(Kx)
	Ky::Vector{Float64} = [ky for ky in -1:2/L[2]:1]
	popfirst!(Ky)

	return [ [kx,ky] for kx in Kx, ky in Ky ], Kx, Ky
end



function StructureFactor(
	Sym::String,                        # Symmetry
	k::Vector{Float64}                  # [kx, ky]
)::Float64

	AllSyms = ["s", "S", "px", "py", "d"]
	if !in(Sym, AllSyms)
		@error "Invalid symmetries."
		return
	end
	kx, ky = k

	if Sym=="s"
		return 1
	elseif Sym=="S"
		return cos(kx) + cos(ky)
	elseif Sym=="px"
		return sqrt(2) * sin(kx)
	elseif Sym=="py"
		return sqrt(2) * sin(ky)
	elseif Sym=="d"
		return cos(kx) - cos(ky)
	end
end



function Getεk(
	εε::Dict{String,Float64},
	k::Vector{Float64};
)::Float64

	εk::Float64 = 0.0
	for (Sym,εSym) in tt
		εk += εSym * StructureFactor(Sym,k)
	end

	return εk
end



function GetΔk(
	ΔΔ::Dict{String,Float64},
	k::Vector{Float64};
)::Float64

	Δk::Float64 = 0.0
	for (Sym,ΔSym) in ΔΔ
		Δk += ΔSym * StructureFactor(Sym,k)
	end

	return Δk
end



function FermiDirac(
	ε::Float64,							# Single-particle energy
	μ::Float64,							# Chemical potential
	β::Float64,							# Inverse temperature
)::Float64

	if β==Inf # Zero temperature distribution
		if ε<=μ
			return 1
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
)::Vector{Float64}

	if Sign=="+"
		return FermiDirac(ε,μ,β) + FermiDirac(-ε,μ,β)
	elseif Sign=="+"
		return FermiDirac(ε,μ,β) - FermiDirac(-ε,μ,β)
	else
		@error "Invalid sign."
		return
	end

end



function GetDensity(
	Phase::String,
	Syms::Vector{String},
	Pars::DataFrame,
	v::DataFrame,
	μ::Float64,
	debug::Bool=false,
	RBS::Bool=true,
	RBd::Bool=false,
	OptimizeBZ::Bool=true
)::Matrix{Float64}

	K, _, _ = GetKGrid(Pars.L)
	LxLy = prod(size(Pars.L))

	# Kinetics
	εS::Float64 = -2 * Pars.t
	try RBS && v.uS
		εS += V * v.uS
	catch
	end

	εd::Float64 = 0.0
	try RBd && v.uS
		εd += V * v.ud
	catch
	end

	εε::Dict{String,Float64} = Dict(
		"S" => εS,
		"d" => εd
	)

	if Phase=="Normal"

		N::Float64 = 0.0
		εK::Matrix{Float64} = Getεk.(εε,K)
		NK::Matrix{Float64} = 2 * FermiDirac.(εε,μ,Pars.β)
		return sum(NK)/(2*LxLy)

	elseif Phase=="AF-Symmetric"

		# Real gap
		reΔs::Float64 = U * v.m
		reΔΔ::Dict{String,Float64} = Dict(
			"s" => reΔs,
		)

		# Imaginary gap
		imΔS::Float64 = 0.0
		try "S" in Syms && v.vS
			imΔS += V * v.vS
		catch
		end

		imΔd::Float64 = 0.0
		try "d" in Syms && v.vd
			imΔd += V * v.vd
		catch
		end

		imΔΔ::Dict{String,Float64} = Dict(
			"S" => imΔS,
			"d" => imΔd
		)

		N::Float64 = 0.0
		εK::Matrix{Float64} = Getεk.(εε,K)
		reΔK::Matrix{Float64} = GetΔk.(reΔΔ,K)
		imΔK::Matrix{Float64} = GetΔk.(imΔΔ,K)
		EK::Matrix{Float64} = sqrt.(εK.^2 .+ reΔK.^2 .+ imΔK.^2)
		NK::Matrix{Float64} = Th.("+",Ek,μ,Pars.β)
		return sum(NK)/(2*LxLy)

	elseif Phase=="AF-Antisymmetric"

		# Real gap
		reΔs::Float64 = U * v.m
		reΔΔ::Dict{String,Float64} = Dict(
			"s" => reΔs,
		)

		# Imaginary gap
		imΔpx::Float64 = 0.0
		try "px" in Syms && v.vpx
			imΔpx += V * v.vpx
		catch
		end

		imΔpy::Float64 = 0.0
		try "py" in Syms && v.vpy
			imΔpy += V * v.vpy
		catch
		end

		imΔΔ::Dict{String,Float64} = Dict(
			"px" => imΔpx,
			"py" => imΔpy
		)

		N::Float64 = 0.0
		εK::Matrix{Float64} = Getεk.(εε,K)
		reΔK::Matrix{Float64} = GetΔk.(reΔΔ,K)
		imΔK::Matrix{Float64} = GetΔk.(imΔΔ,K)
		EK::Matrix{Float64} = sqrt.(εK.^2 .+ reΔK.^2 .+ imΔK.^2)
		NK::Matrix{Float64} = Th.("+",Ek,μ,Pars.β)
		return sum(NK)/(2*LxLy)

	elseif Phase=="SC-Singlet"

		# Real gap
		reΔs::Float64 = 0.0
		try "s" in Syms && v.ws
			reΔs += V * v.ws
		catch
		end

		reΔS::Float64 = 0.0
		try "S" in Syms && v.wS
			reΔS += V * v.wS
		catch
		end

		reΔd::Float64 = 0.0
		try "d" in Syms && v.wd
			reΔd += V * v.wd
		catch
		end

		reΔΔ::Dict{String,Float64} = Dict(
			"s" => reΔs,
			"S" => reΔS,
			"d" => reΔd
		)

		N::Float64 = 0.0
		εK::Matrix{Float64} = Getεk.(εε,K)
		reΔK::Matrix{Float64} = GetΔk.(reΔΔ,K)
		EK::Matrix{Float64} = sqrt.(εK.^2 .+ reΔK.^2 .+ imΔK.^2)
		NK::Matrix{Float64} = 1 .- εK./EK .* tanh(EK .* Pars.β ./ 2)
		return sum(NK)/(2*LxLy)

	elseif Phase=="SC-Triplet"

		# Imaginary gap
		imΔpx::Float64 = 0.0
		try "px" in Syms && v.wpx
			imΔpx += V * v.wpx
		catch
		end

		imΔpy::Float64 = 0.0
		try "py" in Syms && v.wpy
			imΔpy += V * v.wpy
		catch
		end

		imΔΔ::Dict{String,Float64} = Dict(
			"px" => imΔpx,
			"py" => imΔpy
		)

		N::Float64 = 0.0
		εK::Matrix{Float64} = Getεk.(εε,K)
		imΔK::Matrix{Float64} = GetΔk.(imΔΔ,K)
		EK::Matrix{Float64} = sqrt.(εK.^2 .+ imΔK.^2 .+ imΔK.^2)
		NK::Matrix{Float64} = 1 .- εK./EK .* tanh(EK .* Pars.β ./ 2)
		return sum(NK)/(2*LxLy)

	else
		@error "Invalid Syms."
		exit
	end

end

# @doc raw"""
# function GetKPopulation(
# 	Phase::String,
# 	Parameters::Dict{String,Float64},
# 	K::Matrix{Vector{Float64}},
# 	v::Dict{String,Float64},
# 	μ::Float64,
# 	β::Float64;
# 	debug::Bool=false,
# 	RenormalizeBands::Bool=true
# )::Matrix{Float64}
# """
# function GetKPopulation(
# 	Phase::String,						# Mean field phase
# 	Parameters::Dict{String,Float64},	# Model parameters t,U,V
# 	K::Matrix{Vector{Float64}},			# BZ grid
# 	v::Dict{String,Float64},				# HF parameters
# 	μ::Float64,							# Chemical potential
# 	β::Float64;							# Inverse temperature
# 	debug::Bool=false,
# 	RenormalizeBands::Bool=true,			# Conditional renormalization of t
# 	OptimizeBZ::Bool=true				# Conditional BZ optimization
# )::Matrix{Float64}

	# Sym = "S"
	# if in(Phase, ["AF", "FakeAF"])
	# 	Sym *= "-MBZ"
	# end

# 	Nk = zeros(size(K))
# 	wk::Int64 = 0
# 	Ek::Float64 = 0.0

	# Compute free energy hopping shift
	# LxLy::Int64 = prod(size(K))
	# cc::Matrix{Float64} = StructureFactor.("S",K.*pi)
	# εε::Matrix{Float64} = GetBareBandsrameters["t"],K.*pi)
	# ff::Matrix{Float64} = FermiDirac.(εε,μ,β)
	# w::Float64 = sum(cc.*ff)/LxLy # Bare bands correction

# 	for (i,q) in enumerate(K)

# 		wk = GetWeight(q; Sym, OptimizeBZ) # Avoid computational redundance
# 		k = q .* pi # Important: multiply k by pi

# 		if Phase=="Normal"

# 			t = Parameters["t"]
# 			if RenormalizeBands
# 				Conditional renormalization of bands
# 				t -= w * Parameters["V"]
# 			end
# 			εk = GetBareBands)
# 			Nk[i] = 2*FermiDirac(εk,μ,β)

# 		elseif in(Phase, ["AF", "FakeAF"]) && wk >= 1

# 			t = Parameters["t"]
# 			if RenormalizeBands
				# Conditional renormalization of bands
# 				t -= v["w0"] * Parameters["V"]
# 			end

			# Renormalized bands
# 			εk = GetBareBands)

			# Renormalized gap
# 			reΔk::Float64 = v["m"] * (Parameters["U"] + 8*Parameters["V"])
# 			imΔk::Float64 = 2*v["wp"]*Parameters["V"] * StructureFactor("S",k)

			# Renormalized gapped bands
# 			Ek = sqrt( εk^2 + reΔk^2 + imΔk^2 )

			# Local population
# 			Nk[i] = wk * (FermiDirac(-Ek,μ,β) + FermiDirac(Ek,μ,β))

# 		elseif in(Phase, ["SU-Singlet", "FakeSU-Singlet"]) && in(wk,[1,2,4])

# 			t = Parameters["t"]
# 			if RenormalizeBands && "gS" in keys(v)
				# Conditional renormalization of bands
# 				t -= v["gS"]/2 * Parameters["V"]
# 			end

# 			AllSyms = ["Δs", "ΔS", "Δd"]
# 			Fakev = copy(v)
# 			delete!(Fakev, "gS")
# 			delete!(Fakev, "gd")
# 			if !issubset(collect(keys(Fakev)), AllSyms)
# 				@error "Invalid set of symmetries. Please choose from $(AllSyms)."
# 			end

			# Free bands
# 			ξk::Float64 = GetBareBands) - μ
# 			if RenormalizeBands && "gd" in keys(v)
# 				ξk += Parameters["V"] * v["gd"] * StructureFactor("d",k)
# 			end

			# Gap
# 			Δk::Complex{Float64} = 0.0 + 1im * 0.0
# 			for (key, value) in v
# 				if !in(key, ["gS", "gd"])
# 					key = String(key)
# 					key = String(chop(key, head=1, tail=0))
# 					Δk += value * StructureFactor(key,k)
# 				end
# 			end

			# Renormalized gapped bands
# 			Ek = sqrt( ξk^2 + abs(Δk)^2 )

			# Local population
# 			ck::Float64 = 0.0
# 			if Ek!=0.0
# 				ck = ξk/Ek
# 			end
# 			Nk[i] = wk * (1 - ck * tanh(β*Ek/2))

# 		elseif Phase=="Su/Triplet"
# 			@error "Under construction"
# 			return
# 		end
# 	end

# 	return Nk

# end

function GetFreeEnergy(
	Phase::String,
	Syms::Vector{String},
	Pars::DataFrame,
	v::DataFrame,
	μ::Float64,
	debug::Bool=false,
	RenormalizeBands::Bool=true,
	OptimizeBZ::Bool=true
)::Matrix{Float64}

	K, _, _ = GetKGrid(Pars.L)
	LxLy = prod(size(Pars.L))

	# Kinetics
	uS::Float64 = 0.0
	try RenormalizeBands && v.uS
		uS += v.uS
	catch
	end

	ud::Float64 = 0.0
	try RenormalizeBands && v.ud
		ud += v.ud
	catch
	end

	f::Float64 = 0.0
	n::Float64 = 0.5 + Pars.δ
	μ -= n * (8*Pars.V - Pars.U)
	f += μ * n + n^2 * (8*Pars.V - Pars.U)

	ef::Float64 = 0.0

	if Phase=="Normal"

		N::Float64 = 0.0
		εK::Matrix{Float64} = Getεk.(εε,K)
		EF::Matrix{Float64} = log.(1-FermiDirac.(εε,μ,Pars.β))
		ef += sum(EF)/LxLy * 2/Pars.β

		try RBS && v.uS
			ef -= V * abs(v.uS)^2
		catch
		end

		try RBd && v.ud
			ef -= V * abs(v.ud)^2
		catch
		end

		#TODO

	elseif Phase=="AF-Symmetric"

		# Real gap
		reΔs::Float64 = U * v.m
		reΔΔ::Dict{String,Float64} = Dict(
			"s" => reΔs,
		)

		# Imaginary gap
		imΔS::Float64 = 0.0
		try "S" in Syms && v.vS
			imΔS += V * v.vS
		catch
		end

		imΔd::Float64 = 0.0
		try "d" in Syms && v.vd
			imΔd += V * v.vd
		catch
		end

		imΔΔ::Dict{String,Float64} = Dict(
			"S" => imΔS,
			"d" => imΔd
		)

		#TODO

	elseif Phase=="AF-Antisymmetric"

		# Real gap
		reΔs::Float64 = U * v.m
		reΔΔ::Dict{String,Float64} = Dict(
			"s" => reΔs,
		)

		# Imaginary gap
		imΔpx::Float64 = 0.0
		try "px" in Syms && v.vpx
			imΔpx += V * v.vpx
		catch
		end

		imΔpy::Float64 = 0.0
		try "py" in Syms && v.vpy
			imΔpy += V * v.vpy
		catch
		end

		imΔΔ::Dict{String,Float64} = Dict(
			"px" => imΔpx,
			"py" => imΔpy
		)

		#TODO

	elseif Phase=="SC-Singlet"

		# Real gap
		reΔs::Float64 = 0.0
		try "s" in Syms && v.ws
			reΔs += V * v.ws
		catch
		end

		reΔS::Float64 = 0.0
		try "S" in Syms && v.wS
			reΔS += V * v.wS
		catch
		end

		reΔd::Float64 = 0.0
		try "d" in Syms && v.wd
			reΔd += V * v.wd
		catch
		end

		reΔΔ::Dict{String,Float64} = Dict(
			"s" => reΔs,
			"S" => reΔS,
			"d" => reΔd
		)

		#TODO

	elseif Phase=="SC-Triplet"

		# Imaginary gap
		imΔpx::Float64 = 0.0
		try "px" in Syms && v.wpx
			imΔpx += V * v.wpx
		catch
		end

		imΔpy::Float64 = 0.0
		try "py" in Syms && v.wpy
			imΔpy += V * v.wpy
		catch
		end

		imΔΔ::Dict{String,Float64} = Dict(
			"px" => imΔpx,
			"py" => imΔpy
		)

		#TODO

	else
		@error "Invalid Syms."
		exit
	end

end

# @doc raw"""
# function GetFreeEnergy(
# 	Phase::String,
# 	Parameters::Dict{String,Float64},
# 	K::Matrix{Vector{Float64}},
# 	v::Dict{String,Float64},
# 	μ::Float64,
# 	β::Float64;
# 	debug::Bool=false,
# 	RenormalizeBands::Bool=true
# )Float64

# Returns: Free energy density for the given Phase.

# `GetFreeEnergy` takes as input `Phase` (string specifying the mean-field phase,
# the allowed are \"AF\", \"FakeAF\", \"SU-Singlet\", \"FakeSU-Singlet\",
# \"SU-Triplet\", \"FakeSU-Triplet\"), `Parameters`  (dictionary of model
# parameters containing `t`, `U`, `V`), `K` (k-points in the BZ), `v` (dictionary
# of real HF parameters), `μ` (chemical potential) and `β` (inverse temperature).
# It sums the free energy contribution per each couple of momenta. The boolean
# option `RenormalizeBands' allows for choosing to renormalize or not the hopping
# parameter.
# """
# function GetFreeEnergy(
# 	Phase::String,						# Mean field phase
# 	Parameters::Dict{String,Float64},	# Model parameters t,U,V
# 	K::Matrix{Vector{Float64}},			# BZ grid
# 	v::Dict{String,Float64},				# HF parameters
# 	n::Float64,							# Density
# 	μ::Float64,							# Chemical potential
# 	β::Float64;							# Inverse temperature
# 	RenormalizeBands::Bool=true,			# Conditional renormalization of t
# 	OptimizeBZ::Bool=true				# Conditional optimization of BZ
# )::Float64

# 	Sym = "S"
# 	if in(Phase, ["AF", "FakeAF"])
# 		Sym *= "-MBZ"
# 	end

	# f0MFT::Float64 = 0.0 # TODO Integrate free energy computation
# 	fMFT::Float64 = 0.0
# 	LxLy::Int64 = prod(size(K))

	# Compute free energy hopping shift
	# μ0::Float64 = FindRootμ("Free",Parameters,K,v,n,β;RenormalizeBands,OptimizeBZ)
	# cc::Matrix{Float64} = StructureFactor.("S",K.*pi)
	# εε::Matrix{Float64} = GetBareBandsrameters["t"],K.*pi)
	# ff::Matrix{Float64} = FermiDirac.(εε,μ0,β)
	# w::Float64 = sum(cc.*ff)/LxLy # Bare bands correction

# 	for (i,q) in enumerate(K)

# 		wk = GetWeight(q; Sym, OptimizeBZ) # Avoid computational redundance
# 		k = q .* pi # Important: multiply k by pi

# 		if Phase=="Free"

# 			@error "Under construction"
			# t = Parameters["t"]
			# if RenormalizeBands
				# Conditional renormalization of bands
			# 	t -= w * Parameters["V"]
			# end
			# εk = GetBareBands)
			# f0MFT += 2/β * log( 1-FermiDirac(εk,μ0,β) )

# 		elseif in(Phase, ["AF", "FakeAF"]) && wk >= 1
# 			@error "Under construction"

# 		elseif in(Phase, ["SU-Singlet", "FakeSU-Singlet"]) && in(wk,[1,2,4])

# 			t = Parameters["t"]
# 			if RenormalizeBands && "gS" in keys(v)
				# Conditional renormalization of bands
# 				t -= v["gS"]/2 * Parameters["V"]
# 			end

# 			AllSyms = ["Δs", "ΔS", "Δd"]
# 			Fakev = copy(v)
# 			delete!(Fakev, "gS")
# 			delete!(Fakev, "gd")
# 			if !issubset(collect(keys(Fakev)), AllSyms)
# 				@error "Invalid set of symmetries. Please choose from $(AllSyms)."
# 			end

			# Free bands
# 			ξk::Float64 = GetBareBands) - μ
# 			if RenormalizeBands && "gd" in keys(v)
# 				ξk += Parameters["V"] * v["gd"] * StructureFactor("d",k)
# 			end

			# Gap
# 			Δk::Complex{Float64} = 0.0 + 1im * 0.0
# 			for (key, value) in v
# 				if !in(key, ["gS", "gd"])
# 					key = String(key)
# 					key = String(chop(key, head=1, tail=0))
# 					Δk += value * StructureFactor(key,k)
# 				end
# 			end

			# Renormalized gapped bands
# 			Ek = sqrt( ξk^2 + abs(Δk)^2 )

			# Local population
# 			tk::Float64 = 0.0
# 			if Ek!=0.0
# 				tk = abs(Δk)^2/Ek * tanh(β*Ek/2)
# 			end
# 			fMFT += wk * (ξk - Ek + tk + 2/β * log( 1-FermiDirac(Ek,0.0,β) ))

# 		elseif Phase=="Su/Triplet"
# 			@error "Under construction"
# 		end
# 	end

# 	fMFT = fMFT/LxLy + 2*n*μ
# 	if Phase=="Free"
# 		fMFT = f0MFT

# 	elseif Phase=="SU-Singlet"
# 		Corr::Float64=0.0
# 		for gSym in ["gS", "gd"]
# 			try
# 				Corr += abs(v[gSym])^2
# 			catch
# 			end
# 		end

# 		fMFT -= 2 * Parameters["V"] * Corr # Lagrange correction
# 	end

# 	return fMFT
# end
