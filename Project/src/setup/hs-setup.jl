#!/usr/bin/julia
SetupFilePath::String = @__FILE__

# Phase
AllPhases::Set{String} = Set(["Normal","AF-Symmetric","AF-Antisymmetric","SC-Singlet","SC-Triplet"])
Phase::String = "Normal" # Choose your phase
if !in(Phase, AllPhases)
	@error "Invalid phase, please modify at: " * SetupFilePath
	exit()
end

# Syms
SymmetricStructures::Set{String} = Set(["s", "S", "d"])
AntisymmetricStructures::Set{String} = Set(["x", "y"])
Syms::Set{String} = Set([]) # ← Change here

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
AllRB::Set{String} = Set(["S", "d"])
RB::Set{String} = Set(["S"]) # ← Change here
RBS::Bool = "S" in RB ? true : false
RBd::Bool = "d" in RB ? true : false

# Setup
Setup::String = "C[256]-a" # ← Change here
AvailableSetups::Set{String} = Set([
	"Test", # Test setup
	"A[256]-a", # UV plane
	"B[256]-a", # δV plane
	"B[256]-b", # δV plane
	"C[256]-a", # βV plane
])

TestΔv::DataFrame = DataFrame(Dict([
	key => 5e-3 for key in [
		"uS","ud",
		"m","vS","vd","vx","vy",
		"ws","wS","wd","wx","wy"
	]
]))
MainΔv::DataFrame = DataFrame(Dict([
	key => 5e-4 for key in [
		"uS","ud",
		"m","vS","vd","vx","vy",
		"ws","wS","wd","wx","wy"
	]
]))

xScale = identity
if !in(Setup, AvailableSetups)
	@error "Invalid Setup=$(Setup), please modify at:" SetupFilePath
elseif Setup=="Test"
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
	xVar = "V"
	pVar = "δ"

# --- MAIN UV plane RUN ---
elseif Setup=="A[256]-a"
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
	xVar = "U"
	pVar = "V"

# --- MAIN δV plane RUN ---
elseif Setup=="B[256]-a"
	tt = [1.0]
	UU = [0.0, 5.0, 10.0]
	VV = [V for V in 0.0:0.4:4.0]
	LL = [256]
	δδ = [δ for δ in 0.0:0.05:0.45]
	ββ = [100.0]
	p = 200
	Δv = MainΔv
	Δn = 1e-2
	g = 0.2
	xVar = "V"
	pVar = "δ"

# --- MAIN δV plane RUN ---
elseif Setup=="B[256]-b"
	tt = [1.0]
	UU = [0.0, 5.0, 10.0]
	VV = [V for V in 4.0:0.4:8.0]
	LL = [256]
	δδ = [δ for δ in 0.0:0.05:0.45]
	ββ = [100.0]
	p = 400
	Δv = MainΔv
	Δn = 1e-2
	g = 0.05
	xVar = "V"
	pVar = "δ"

# --- MAIN βV plane RUN ---
elseif Setup=="C[256]-a"
	tt = [1.0]
	UU = [0.0]
	VV = [V for V in 0.0:0.2:6.0]
	LL = [256]
	δδ = [0.0, 0.2, 0.4]
	ββ = [1.0, 2.0, 10.0, 20.0, 100.0]
	p = 100
	Δv = MainΔv
	Δn = 1e-2
	g = 0.1
	xVar = "V"
	pVar = "β"
	xScale = identity
end