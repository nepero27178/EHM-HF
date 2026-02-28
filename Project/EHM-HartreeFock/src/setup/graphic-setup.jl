using Colors
using ColorSchemes

# Get color scheme
const TabColors = ColorSchemes.tab10
const tabblue = TabColors[1]
const tabgreen = TabColors[3]
const tabred = TabColors[4]

# Convert to Lab color space
const StartLab = Colors.convert(Lab{Float64}, tabblue)
const StopLab = Colors.convert(Lab{Float64}, tabred)
const MidLab = Colors.convert(Lab{Float64}, RGB(0.85,0.85,0.85))

# Interpolate in Lab space using the range syntax
const TotalLabSteps = 100
const Cool = [Colors.convert(RGB, c) for c in range(StartLab, stop=MidLab, length=TotalLabSteps)]
const Warm = [Colors.convert(RGB, c) for c in range(MidLab, stop=StopLab, length=TotalLabSteps)]
const cmap = vcat(Cool,Warm)