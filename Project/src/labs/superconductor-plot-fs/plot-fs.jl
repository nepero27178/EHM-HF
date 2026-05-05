using CairoMakie

include("/home/nepero27178/Thesis/EHM-HF/Project/src/setup/graphic-setup.jl")
include("/home/nepero27178/Thesis/EHM-HF/Project/src/modules/methods-3D-plotting.jl")

CairoMakie.activate!()
MT = Makie.MathTeXEngine
MT_DIR = dirname(pathof(MT)) * "/../assets/fonts/NewComputerModern"
set_theme!(fonts = (
	regular = MT_DIR * "/NewCM10-Regular.otf",
	bold = MT_DIR * "/NewCM10-Bold.otf"
))

FilePathIn = "/home/nepero27178/Thesis/EHM-HF/Project/data/refined/Mode=rs/Setup=B[128]/Phase=SC-Singlet/RB=Sd_Syms=Ssd.csv"
Setup,Phase,Syms,RB,Opt,Layer = UnpackFilePath(FilePathIn)
DF_SC = CSV.read(FilePathIn,DataFrame)
Sim = Simulation(DF_SC,Setup,Phase,Syms,RB)
EnlargeDF!(Sim)
df_SC = filter(:U => x -> x==4.0, DF_SC)
U = 4.0

H = Figure(size=(1100,400),figure_padding=1)
axs = [Axis(H[1,j],xlabel = L"$\delta$",ylabel = L"$V$") for j in 1:2]
linkyaxes!(axs...)

# Free energy
ax = axs[1]
xx,yy,ff = ReshapeData(df_SC;xVar="δ",yVar="V",zVar="f")
ax.title = L"$f/t$ ($t=1.0$, $U=%$(U)$, $L=128$, $\delta=0.0$, $\beta=100.0$)"

cs = colorschemes[:tabquieterrev]
f = heatmap!(ax,xx,yy,ff,colormap=cs)
Colorbar(H[1,0],f,flipaxis=false)

DF_s = CSV.read("/home/nepero27178/Thesis/EHM-HF/Project/src/labs/superconductor-region-boundaries/s_U=$(U).csv",DataFrame)
x_s = DF_s.x_s
y_s = DF_s.y_s
scatterlines!(ax,x_s,y_s,color=:black,markersize=2,linewidth=0.1)

DF_d = CSV.read("/home/nepero27178/Thesis/EHM-HF/Project/src/labs/superconductor-region-boundaries/d_U=$(U).csv",DataFrame)
x_d = DF_d.x_d
y_d = DF_d.y_d
scatterlines!(ax,x_d,y_d,color=:black,markersize=2,linewidth=0.1)

text!(ax,0.35,0.75,text="Normal",align=(:center,:center))
text!(ax,0.2,2.0,text=L"$d$-wave",align=(:center,:center),color=:white)
text!(ax,0.32,3.5,text="Mixed",align=(:center,:center),color=:white)
text!(ax,0.43,2.35,text=L"$s \oplus s^*$-wave",align=(:center,:center))

# Entropy
ax = axs[2]
ax.ylabelvisible = false
ax.yticklabelsvisible = false

xx,yy,ss = ReshapeData(df_SC;xVar="δ",yVar="V",zVar="s")
ax.title = L"$s/k_\mathrm{B}$ ($t=1.0$, $U=%$(U)$, $L=128$, $\delta=0.0$, $\beta=100.0$)"

cs = colorschemes[:tabwarm]
s = heatmap!(ax,xx,yy,ss,colormap=cs)
Colorbar(H[1,3],s)

scatterlines!(ax,x_s,y_s,color=:black,markersize=2,linewidth=0.1)
scatterlines!(ax,x_d,y_d,color=:black,markersize=2,linewidth=0.1)

text!(ax,0.35,0.75,text="Normal",align=(:center,:center))
text!(ax,0.2,2.0,text=L"$d$-wave",align=(:center,:center))
text!(ax,0.32,3.5,text="Mixed",align=(:center,:center))
text!(ax,0.43,2.35,text=L"$s \oplus s^*$-wave",align=(:center,:center))

save("fs_U=$(U).pdf", H)