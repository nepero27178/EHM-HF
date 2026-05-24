using CairoMakie
using CSV
using DataFrames

CairoMakie.activate!()
MT = Makie.MathTeXEngine
MT_DIR = dirname(pathof(MT)) * "/../assets/fonts/NewComputerModern"
set_theme!(fonts = (
	regular = MT_DIR * "/NewCM10-Regular.otf",
	bold = MT_DIR * "/NewCM10-Bold.otf"
))

LAB_ROOT = @__DIR__
include(LAB_ROOT * "/../../setup/graphic-setup.jl")

H = Figure(size=(1300,400),figure_padding = 1)
axs = [Axis(H[1,j]) for j in 1:2]
linkyaxes!(axs...)

ax = axs[1]
ax.xlabel = L"$U$"
ax.ylabel = L"$m$"
ax.title = "Half-filling hyperbolic approximation"
ax.titlefont = :regular

ax = axs[2]
ax.xlabel = L"$U$"
ax.ylabelvisible=false
ax.yticklabelsvisible=false
ax.title = "Half-filling exponential approximation"
ax.titlefont = :regular

ρρ = [0.12,0.14,0.16,0.18,0.2,0.22]
xx = [x for x in 0:0.01:16]
cs = CoolQuiet
for (i,ρ) in enumerate(ρρ)
	q = floor(Int64,length(cs)/length(ρρ)*i)
	ϵ = 1/(2*ρ)
	f(x) = ϵ/( x*sinh(1/(ρ*x)) )
	ff = f.(xx)
	lines!(axs[1],xx,ff,color=cs[q],label=L"$\rho_0=%$(ρ)$")

	g(x) = 2*ϵ/x * exp.( -1/(ρ*x) )
	gg = g.(xx)
	lines!(axs[2],xx,gg,color=cs[q])
end

FilePathIn = "/home/nepero27178/Thesis/EHM-HF/Project/data/refined/Mode=hs/Setup=D[256]/Phase=AF/RB=S_Syms=.csv"
DF = CSV.read(FilePathIn,DataFrame)

xx = unique(DF.U)
yy = filter(:δ => x -> x==0.0, DF).m
scatterlines!(axs[1],xx,yy,color=:black,label="Data")
scatterlines!(axs[2],xx,yy,color=:black)

H[1,3] = Legend(H,axs[1],framevisible=false)
FilePathOut = "m-analytic.pdf"
save(FilePathOut,H)