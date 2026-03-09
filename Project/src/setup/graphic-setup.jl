using Colors
using ColorSchemes

# Get MatPlotLib base colors
const TabColors::ColorScheme{Vector{RGB{Float64}},String,String}  = ColorSchemes.tab10
const tabblue::RGB{Float64} = TabColors[1]
const tabgreen::RGB{Float64} = TabColors[3]
const tabred::RGB{Float64} = TabColors[4]

# Convert to Lab color space
const TabBlueLab::Lab{Float64} = Colors.convert(Lab{Float64}, tabblue)
const TabRedLab::Lab{Float64} = Colors.convert(Lab{Float64}, tabred)
const TabGreenLab::Lab{Float64} = Colors.convert(Lab{Float64}, tabgreen)
const WhiteLab::Lab{Float64} = Colors.convert(Lab{Float64}, RGB(0.85,0.85,0.85))
const WhiterLab::Lab{Float64} = Colors.convert(Lab{Float64}, RGB(0.98,0.98,0.98))

# Interpolate in Lab space to gain perceptual uniformity
const TotalLabSteps::Int64 = 100
const Cool::Vector{RGB{Float64}} = [Colors.convert(RGB, c) for c in range(TabBlueLab, stop=WhiteLab, length=TotalLabSteps)]
const Cooler::Vector{RGB{Float64}} = [Colors.convert(RGB, c) for c in range(TabBlueLab, stop=WhiterLab, length=TotalLabSteps)]
const Warm::Vector{RGB{Float64}} = [Colors.convert(RGB, c) for c in range(WhiteLab, stop=TabRedLab, length=TotalLabSteps)]
const CoolWarm::Vector{RGB{Float64}} = vcat(Cool,Warm)
const WarmCool::Vector{RGB{Float64}} = reverse(CoolWarm)
const Quiet::Vector{RGB{Float64}} = [Colors.convert(RGB, c) for c in range(TabGreenLab, stop=WhiteLab, length=TotalLabSteps)]
const CoolQuiet::Vector{RGB{Float64}} = [Colors.convert(RGB, c) for c in range(TabBlueLab, stop=TabGreenLab, length=TotalLabSteps)]

# Add colorschemes
colorschemes[:tabcool] = ColorScheme(Cool, "custom cool from matplotlib", "perceptually uniform sequential")
colorschemes[:tabcooler] = ColorScheme(Cooler, "custom cool from matplotlib", "perceptually uniform sequential")
colorschemes[:tabwarm] = ColorScheme(Warm, "custom warm from matplotlib", "perceptually uniform sequential")
colorschemes[:tabcoolwarm] = ColorScheme(CoolWarm, "custom coolwarm from matplotlib with grey midpoint", "two-tones perceptually uniform sequential")
colorschemes[:tabquiet] = ColorScheme(Quiet, "custom greens from matplotlib", "perceptually uniform sequential")
colorschemes[:tabcoolquiet] = ColorScheme(CoolQuiet, "custom blues-greens from matplotlib", "perceptually uniform sequential")

# Add reversed colorschemes
colorschemes[:tabcoolrev] = ColorScheme(reverse(Cool), "custom reversed cool from matplotlib", "perceptually uniform sequential")
colorschemes[:tabcoolerrev] = ColorScheme(reverse(Cooler), "custom reversed cool from matplotlib", "perceptually uniform sequential")
colorschemes[:tabwarmcool] = ColorScheme(WarmCool, "custom coolwarm from matplotlib with grey midpoint", "two-tones perceptually uniform sequential")