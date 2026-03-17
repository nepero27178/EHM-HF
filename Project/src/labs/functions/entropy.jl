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
xlogx(x::Float64) = x==0.0 ? 0.0 : x*log(x)
xx = [x for x in 0.0:0.01:1.0]
yy = -xlogx.(xx) .- xlogx.(1 .- xx)

H = Figure(size=(500,400),figure_padding = 1)
ax = Axis(H[1,1])
ax.xlabel = L"$f$"
ax.ylabel = L"$s/k_\mathrm{B}$"
ax.title = L"Entropy density at occupation $f$: $s(f) = -k_\mathrm{B} [ f \log f + (1-f) \log (1-f) ]$"
lines!(ax,xx,yy,color=tabred)

FilePathOut = "entropy.pdf"
save(FilePathOut,H)