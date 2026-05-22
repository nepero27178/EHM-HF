using DataFrames
using CSV

xx = [x for x in -4.5:0.01:4.5]
e(x,a) = exp(-(x/a)^2)
f(x,U) = (1-U)^2 * sqrt( maximum([0,1-x^2]) )
U = 2.5

f0(x,a0,A0) = A0 * f(x,U) * e(x,a0)
fL(x,a1,A1,s) = A1 * f(x+s,U) * e(x+s,a1)
fR(x,a1,A1,s) = A1 * f(x-s,U) * e(x-s,a1)
F(x,a0,A0,a1,A1,s) = fL(x,a1,A1,s) + fR(x,a1,A1,s) + f0(x,a0,A0)

yy1 = F.(xx,2,0.75,0,0,0)
yy2 = F.(xx,0.44,0.7,0.68,0.23,1.1) .+ 2
yy3 = F.(xx,0.15,0.7,0.5,0.4,1.4) .+ 4
yy4 = F.(xx,0,0,2,0.3,1.8) .+ 6

FilePathOut = "data.csv"
DF = DataFrame("x" => xx, "y1" => yy1, "y2" => yy2, "y3" => yy3, "y4" => yy4)
CSV.write(FilePathOut,DF)
