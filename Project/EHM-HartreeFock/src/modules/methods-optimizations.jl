function GetWeight(
	k::Vector{Float64};					# [kx, ky] in pi units
	Sym::String="S",						# Symmetry structure
	OptimizeBZ::Bool=true				# Option to disable optimization
)::Int64

	# Weights
	wk::Int64 = 0
	kx, ky = k

	if !OptimizeBZ
        wk = 1
	elseif Sym=="S-MBZ" # Half-sized Brillouin Zone

		if abs(kx)+abs(ky) <= 1
			if kx>=0 && ky>0 && kx+ky<1
				# Bulk (four times)
				wk = 4
			elseif kx+ky==1 && (kx!=1 && ky!=1)
				# Edge (two times due to nesting)
				wk = 2
			elseif kx==0 && (ky==0 || ky==1)
				# Special points
				wk = 1
			end
		end

	elseif Sym=="S"

        if (kx==0 && ky==0) || (kx==1 && ky==1)
            wk = 1
        elseif kx >= 0 && ky > 0
			if kx < 1 && ky < 1
				# Bulk (four times)
				wk = 4
			elseif kx==1 || ky==1
				# Edge (two times due to nesting)
				wk = 2
            end
        end
	end

	return wk

end

function GetUc(
	Pars::DataFrame,
	v::DataFrame;
	how::String="DisSum",
	Δn::Float64=1e-3,
	μ0::Float64=0.0,
	RBS::Bool=true,
	RBd::Bool=true,
	OptBZ::Bool=false
)::Float64

	# Get pars
	t::Float64 = first(Pars.t)
	L::Int64 = first(Pars.L)
	β::Float64 = first(Pars.β)

	# Find normal chemical potential
	μ = FindRootμ("Normal",Set{String}(),Pars,v;Δn,μ0,RBS,RBd,OptBZ,debug)

	if !in(Method, ["NumInt", "DisSum"])
		@error "Wrong method. Acceptable: \"NumInt\", \"DisSum\"."

	# Numerical integration
	elseif Method=="NumInt"

		F(x,m) = Elliptic.K( sqrt(1-x^2) ) * tanh( β/2 * (4*t*x - m) ) / (x -m/(4*t))
		DomainDown = (-1,0)
		ProblemDown = IntegralProblem(F, DomainDown, μ)
		SolutionDown = solve(ProblemDown, HCubatureJL())

		DomainUp = (0,1)
		ProblemUp = IntegralProblem(F, DomainUp, μ)
		SolutionUp = solve(ProblemUp, HCubatureJL())
		Uc = (2*pi)^2 / (SolutionUp.u + SolutionDown.u)

	# Discrete sum
	elseif Method=="DisSum"

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
		Getξk(k::Vector{Float64})::Float64 = GetFk(εε,k)-μ # Function
		ξK::Matrix{Float64} = Getξk.(K)
		tK::Matrix{Float64} = tanh.(ξK .* β/2) ./ ξK
		replace!(tK, NaN => β/2) # @ x~0 : tanh(x)/x~1
		Uc = 2*LxLy / sum(tK)

	end

    return Uc

end

function GetOptimalg(
    U::Float64,
    Uc::Float64;
    ΔU::Float64=1.0
)::Float64

    u::Float64 = U/Uc
    gc::Float64 = 2/(u+1)
    return gc * (1-gc/2 * ΔU/Uc)

end
