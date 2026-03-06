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
const WhiterLab::Lab{Float64} = Colors.convert(Lab{Float64}, RGB(1.0,1.0,1.0))

# Interpolate in Lab space to gain perceptual uniformity
const TotalLabSteps::Int64 = 100
const Cool::Vector{RGB{Float64}} = [Colors.convert(RGB, c) for c in range(TabBlueLab, stop=WhiteLab, length=TotalLabSteps)]
const Cooler::Vector{RGB{Float64}} = [Colors.convert(RGB, c) for c in range(TabBlueLab, stop=WhiterLab, length=TotalLabSteps)]
const Warm::Vector{RGB{Float64}} = [Colors.convert(RGB, c) for c in range(WhiteLab, stop=TabRedLab, length=TotalLabSteps)]
const CoolWarm::Vector{RGB{Float64}} = vcat(Cool,Warm)
const Quiet::Vector{RGB{Float64}} = [Colors.convert(RGB, c) for c in range(TabGreenLab, stop=WhiteLab, length=TotalLabSteps)]

# Add to full dictionary
colorschemes[:tabcool] = ColorScheme(Cool, "custom cool from matplotlib", "perceptually uniform sequential")
colorschemes[:tabcoolerrev] = ColorScheme(reverse(Cooler), "custom cool from matplotlib", "perceptually uniform sequential")
colorschemes[:tabwarm] = ColorScheme(Warm, "custom warm from matplotlib", "perceptually uniform sequential")
colorschemes[:tabcoolwarm] = ColorScheme(CoolWarm, "custom coolwarm from matplotlib with grey midpoint", "two-tones perceptually uniform sequential")
colorschemes[:tabquiet] = ColorScheme(Quiet, "custom cool from matplotlib", "perceptually uniform sequential")