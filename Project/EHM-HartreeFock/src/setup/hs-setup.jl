#!/usr/bin/julia
const SetupFilePath::String = @__FILE__

# Phase
const AllPhases::Set{String} = Set(["Normal","AF-Symmetric","AF-Antisymmetric","SC-Singlet","SC-Triplet"])
const Phase::String = "AF-Symmetric" # Choose your phase
if !in(Phase, AllPhases)
	@error "Invalid phase, please modify at: " * SetupFilePath
	exit()
end

# Syms
const SymmetricStructures::Set{String} = Set(["s", "S", "d"])
const AntisymmetricStructures::Set{String} = Set(["x", "y"])
const Syms::Set{String} = Set(["s", "S", "d"]) # ← Change here

Err::Bool = false
(Phase=="Normal" && length(Syms)>0) ? Err = true : false
(Phase=="AF-Symmetric" && !issubset(Syms, SymmetricStructures)) ? Err = true : false
(Phase=="AF-Antisymmetric" && !issubset(Syms, AntisymmetricStructures)) ? Err = true : false
(Phase=="SC-Singlet" && !issubset(Syms, SymmetricStructures)) ? Err = true : false
(Phase=="SC-Triplet" && !issubset(Syms, AntisymmetricStructures)) ? Err = true : false
if Err
	@error "Invalid Syms=$(Syms...) for Phase=$(Phase). Please modify at: " SetupFilePath
end

# RB
const AllRB::Set{String} = Set(["S", "d"])
const RB::Set{String} = Set(["S", "d"]) # ← Change here
const RBS::Bool = "S" in RB ? true : false
const RBd::Bool = "d" in RB ? true : false

# Setup
const Setup::String = "Test[20]" # ← Change here
const AvailableSetups::Set{String} = Set([
	"Test[20]", # Test setup
	"A[256]", # UV plane
	"B[256]", # δV plane
])

const TestΔv::DataFrame = DataFrame(Dict([
	key => 5e-3 for key in [
		"uS","ud",
		"m","vS","vd","vx","vy",
		"ws","wS","wd","wx","wy"
	]
]))
const MainΔv::DataFrame = DataFrame(Dict([
	key => 5e-4 for key in [
		"uS","ud",
		"m","vS","vd","vx","vy",
		"ws","wS","wd","wx","wy"
	]
]))

if !in(Setup, AvailableSetups)
	@error "Invalid Setup=$(Setup), please modify at:" SetupFilePath
elseif Setup=="Test[20]"
	# TEST-SETUP
	tt = [1.0]
	UU = [0.0]
	VV = [1.0, 2.0, 3.0]
	LL = [20]
	δδ = [0.0, 0.1, 0.2, 0.3]
	ββ = [100.0]
	p = 20
	Δv = TestΔv
	Δn = 1e-2
	g = 0.1

# --- MAIN A RUN ---
elseif Setup=="A[256]"
	tt = [1.0]
	UU = [U for U in 0.0:0.25:5.0]
	VV = [V for V in 0.0:0.25:5.0]
	LL = [256]
	δδ = [0.0, 0.2, 0.4]
	ββ = [100.0]
	p = 200
	Δv = MainΔv
	Δn = 1e-2
	g = 0.2

# --- MAIN B RUN ---
elseif Setup=="B[256]"
	tt = [1.0]
	UU = [0.0, 5.0, 10.0]
	VV = [V for V in 0.0:0.1:5.0]
	LL = [256]
	δδ = [δ for δ in 0.0:0.05:0.45]
	ββ = [100.0]
	p = 200
	Δv = MainΔv
	Δn = 1e-2
	g = 0.2
end