using CairoMakie

CairoMakie.activate!()
MT = Makie.MathTeXEngine
MT_DIR = dirname(pathof(MT)) * "/../assets/fonts/NewComputerModern"
set_theme!(fonts = (
	regular = MT_DIR * "/NewCM10-Regular.otf",
	bold = MT_DIR * "/NewCM10-Bold.otf"
))

LAB_ROOT = @__DIR__
include(LAB_ROOT * "/../../setup/graphic-setup.jl")

F = Figure(size=(800,400), figure_padding=1)
ax = Axis(F[1,1])
ax.xlabel = L"$x$"
ax.ylabel = "Distributions"
ax.title = L"Example $f(x)$ and $f\,'(x)$ distributions ($\beta=100.0$, $\Delta=0.2$)"
β = 100.0
Δ = 0.2
f(x,μ) = 1/( exp(β*(x-μ))+1 )
E(x) = sqrt(x^2+Δ^2)
g(x,μ) = -x/E(x) * (f(-E(x),μ)-f(E(x),μ))

ls = [:solid,:dash,:dot]
xx = [x for x in -1.0:0.01:1.0]

μμ = [0.0,0.3,0.6]
for (j,μ) in enumerate(μμ)
	ff = f.(xx,μ)
	lines!(ax,xx,ff,color=tabblue,linestyle=ls[j],label=L"$f(x)$, $\mu=%$(μ)$")
end

for (j,μ) in enumerate(μμ)
	gg = g.(xx,μ)
	lines!(ax,xx,gg,color=tabred,linestyle=ls[j],label=L"$f\,'(x)$, $\mu=%$(μ)$")
end

axislegend(ax,position=:lb)
FilePathOut = "pseudo-dirac.pdf"
save(FilePathOut,F)