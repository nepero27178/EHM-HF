#!/usr/bin/julia

using CairoMakie
using LaTeXStrings
using ColorSchemes
using Elliptic

# Includer
PROJECT_ROOT = @__DIR__
PROJECT_ROOT *= "/../.."   # Up to the effective root
include(PROJECT_ROOT * "/setup/graphic-setup.jl")

MT = Makie.MathTeXEngine
MT_DIR = dirname(pathof(MT)) * "/../assets/fonts/NewComputerModern"

set_theme!(fonts = (
	regular = MT_DIR * "/NewCM10-Regular.otf",
	bold = MT_DIR * "/NewCM10-Bold.otf"
))

H = Figure(size=(600,400),figure_padding=1)
ax = Axis(H[1,1])
ax.xlabel = L"$x$"
ax.ylabel = L"$K(x)$"
ax.title = "Complete elliptic integral of the first kind"
ax.titlefont = :regular

xx = [x for x in 0.0:0.001:1.0]
KK = Elliptic.K.(xx)
lines!(ax,xx,KK,color=tabblue,label=L"$K(x)$")
xlims!(ax,0,1)
ylims!(ax,0,5)
axislegend(ax,position=:lt)
FilePathOut = "K(x).pdf"
save(FilePathOut,H)

H = Figure(size=(600,400),figure_padding=1)
ax = Axis(H[1,1])
ax.xlabel = L"$\epsilon/t$"
ax.ylabel = L"$\rho(\epsilon)$"
ax.title = "Square lattice DoS"
ax.titlefont = :regular

xx = [x for x in -4.0:0.001:4.0]
DD = Elliptic.K.( sqrt.(1 .- (xx/4).^2) ) ./ (2*pi^2)
lines!(ax,xx,DD,color=tabblue,label=L"$\rho(\epsilon)$")

Heaviside(x) = x >= 0 ? 1 : 0
x0 = 1.06
ρ0 = 0.47
HH = ρ0 * Heaviside.(x0 .- abs.(xx))
lines!(ax,xx,HH,color=tabred,linestyle=:dash,label=L"$\rho_0\mathrm{H}(\epsilon_0-|\epsilon|)$")

xlims!(ax,-4,4)
ylims!(ax,0,0.5)
ax.xticks = ([-4,-2,-x0,0,x0,2,4],["-4","-2",L"$-\epsilon_0$","0",L"$\epsilon_0$","2","4"])
ax.yticks = ([0.0,0.1,0.2,0.3,0.4,ρ0,0.5],["0.0","0.1","0.2","0.3","0.4",L"$\rho_0$","0.5"])
axislegend(ax,position=:rt)
FilePathOut = "D(x).pdf"
save(FilePathOut,H)