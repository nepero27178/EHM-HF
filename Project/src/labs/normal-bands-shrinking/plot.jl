using GLMakie
using CairoMakie
using LaTeXStrings

GLMakie.activate!()
MT = Makie.MathTeXEngine
MT_DIR = dirname(pathof(MT)) * "/../assets/fonts/NewComputerModern"
set_theme!(fonts = (
	regular = MT_DIR * "/NewCM10-Regular.otf",
	bold = MT_DIR * "/NewCM10-Bold.otf"
))

LAB_ROOT = @__DIR__
include(LAB_ROOT * "/../../setup/graphic-setup.jl")

H = Figure(size=(800,700),figure_padding=1)
ax = Axis3(
	H[1,1],
	xlabel=L"$k_x/\pi$",
	ylabel=L"$k_y/\pi$",
	zlabel=L"$\tilde{\epsilon}_{\mathbf{k}}/t$",
	xticks=[-1,0,1],
	yticks=[-1,0,1],
	zticks=[-4,-2,0,2,4],
	title="Bare bands shrinking",
	titlefont=:regular,
	aspect=(1,1,1),
	azimuth=0.3*pi,
	elevation=pi/8,
	xlabelalign = (:center, :center),
	ylabelalign = (:center, :center),
	zlabelalign = (:center, :center),
	xlabelrotation = 0, # Horizontal xlabel
	ylabelrotation = 0, # Horizontal ylabel
	zlabelrotation = 0, # Horizontal zlabel
)

xx = [x for x in -1:0.05:1]
yy = [y for y in -1:0.05:1]
zz = zeros(length(xx),length(yy))
for (i,x) in enumerate(xx)
	for (j,y) in enumerate(yy)
		zz[j,i] = -2*(cos(pi*x)+cos(pi*y))
	end
end
cs = CoolWarm
clims=(-4,4)
A = surface!(ax,xx,yy,zz,colormap=cs,colorrange=clims,alpha=0.9,shading=false)
B = surface!(ax,xx,yy,zz./2,colormap=cs,colorrange=clims,alpha=0.9,shading=false)
C = surface!(ax,xx,yy,zz./4,colormap=cs,colorrange=clims,alpha=0.9,shading=false)
Colorbar(H[1,0],colormap=cs,colorrange=clims)

FilePathOut = LAB_ROOT * "/normal-bands-shrinking.png"
save(FilePathOut, H, px_per_unit=6)