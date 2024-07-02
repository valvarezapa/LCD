include("LyapunovCycleDetector.jl")

using PyPlot
using Images
using Polynomials
using BenchmarkTools
using FileIO

# IMPORTANT: LCD benefits from multi-threading, so make sure you are running Julia
# with 4 threads! You can easily check the number of threads you are working with
# simply by executing the next line.
#Threads.nthreads()

# Sample polynomials
tradPol2=complex([-1.0,0.0,1.0]) # z^2-1
tradPol3=complex([-1.0,0.0,0.0,1.0]) # z^3-1
tradPol5=complex([-1.0,0.0,0.0,0.0,0.0,1.0]) # z^5-1
interestingFatou=complex([-1.0,-20.0,200.0,-550.0,600.0]) # 600z^4-550z^3+200z^2-20z-1
interestingJulia=complex([-4.0,0.0,-3.0,0.0,1.0]) # z^4-3z^2-4
#2-cycles
pol2cycle=complex([2.0,-2.0,0.0,1.0]) # z^3-2z+2
anotherPol2cycle=complex([-5.29,0.0,-4.29,0.0,1.0]) # (z-im)(z+im)(z-2.3)(z+2.3)
cubicPol2cycle=complex([-0.845-0.81*im,1.0,-0.845-0.81*im,1.0]) # (z-0.845-0.81*im)(z^2+1)
#3-cycles
pol3cycle=complex([-1.49175,1.0,0.0,-1.0]) # -z^3+z-1.49175
anotherPol3cycle=complex([-1.214-0.5895*im,1.0,-1.214-0.5895*im,1.0]) # (z-(1.214+0.5895*im))(z^2+1)
#4-cycles
pol4cycle=complex([1.422-0.407*im,1.0,1.422-0.407*im,1.0]) # (z-(-1.422+0.407*im))(z^2+1)
anotherPol4cycle=complex([-0.832-0.82*im,1.0,-0.832-0.82*im,1.0]) # (z-0.832-0.82*im)(z^2+1)

polCycle_chebyshev=complex([1.2866*im,-1.0,-1.2866*im,1.0]) # (z-1.2866*im)(z^2-1)
polExtra_chebyshev=complex([-2.5+1.0*im,-1.0,-1.0*(-2.5+1.0*im),1.0]) # (z-(-2.5+1.0*im))(z^2-1)
polAux_chebyshev=complex([1.157*im,-1.0,-1.0*(1.157*im),1.0]) # (z-(1.157*im))(z^2-1)

# Apply an iterative method  to the polynomial to obtain a rational function.
# Iterative methods currently supported: Newton, Chebyshev, Halley
pair=LCD.applyIterative_method(polAux_chebyshev,"chebyshev")
coefficientlistnum=pair[1]
coefficientlistden=pair[2]

# Computes the list of fixed points (not necessary for the following methods)
#fixedPoints=LCD.positionIteration_superAttractingFixedPoints(coefficientlistnum,coefficientlistden,(-1.5,1.5),(-1.5,1.5),200,5,8)[2]
#println(fixedPoints)

# Plot basins of attraction
#=
graphic=LCD.plotBasinsOfAttraction_Lyapunov(coefficientlistnum,coefficientlistden,(-1.5,1.5),(-1.5,1.5),200,5,8;coloringStrat="position")
img=graphic[1]
numColors=graphic[2]
matrix=graphic[3]
colors=graphic[4]
=#

# Plot ONLY the n-cycles
#LCD.plot_nCycles_Lyapunov(coefficientlistnum,coefficientlistden,(-1.5,1.5),(-1.5,1.5),100,2,8)

# Plot ONLY the n-cycles, with min<=n<=max
#LCD.plot_manyCycles_Lyapunov(coefficientlistnum,coefficientlistden,(-1.5,1.5),(-1.5,1.5),25,2,3,8)

# Plot basins of attraction with the n-cycles, given 2=<n<=max
#=
graphic=LCD.plotWithCycles_BasinsOfAttraction_Lyapunov_givenFixedPoints(coefficientlistnum,coefficientlistden,(-0.8,0.8),(-1.0,1.0),fixedPoints,200,5,8)
img=graphic[1]
numColors=graphic[2]
matrix=graphic[3]
colors=graphic[4]
=#

# Plot a neighborhood of 0 along with a neighborhood of infinity
#=
PyPlot.subplot(121)
graphic=LCD.plotBasinsOfAttraction_Lyapunov(coefficientlistnum,coefficientlistden,(-1.5,1.5),(-1.5,1.5),25,1,8;coloringStrat="position")
numColors=graphic[2]
ax1=graphic[1]
PyPlot.subplot(122)
ax2=LCD.plotBasinsOfAttraction_Lyapunov(coefficientlistnum,coefficientlistden,(-1.5,1.5),(-1.5,1.5),25,1,8;coloringStrat="position",inverted=true)[1]
=#

# Detect cyclic behaviour of cubic polynomials under iteration of Newton's method
#=
graphic=LCD.detectCathastrophicBehaviour_cubicPolynomials_Newton((0.81,0.87),(0.79,0.84),200,5,8)
img=graphic[1]
numColors=graphic[2]
matrix=graphic[3]
colors=graphic[4]
=#

# Detect cyclic behaviour of cubic polynomials under iteration of Newton's method, distinguishing cycles by its length
#=
graphic=LCD.detectCathastrophicBehaviour_distinguishCycles_cubicPolynomials_Newton((-1.4275,-1.4125),(0.4,0.4162),120,12,8)
img=graphic[1]
numColors=graphic[2]
matrix=graphic[3]
colors=graphic[4]
=#

# Detect cyclic behaviour of cubic polynomials under iteration of Chebyshev's method
#=
graphic=LCD.detectCathastrophicBehaviour_cubicPolynomials_Chebyshev((-2.5,2.5),(-2.5,2.5),300,5,8;behaviourOnly=true)
img=graphic[1]
numColors=graphic[2]
matrix=graphic[3]
colors=graphic[4]
=#

# Detect cyclic behaviour of cubic polynomials under iteration of Chebyshev's method, distinguishing cycles by its length
#=
graphic=LCD.detectCathastrophicBehaviour_distinguishCycles_cubicPolynomials_Chebyshev((-0.02,0.02),(1.14,1.18),300,8,8)
img=graphic[1]
numColors=graphic[2]
matrix=graphic[3]
colors=graphic[4]
=#

# Compute the Lyapunov constants of a region of the parameter space induced by the application of Newton's method to cubic polynomials
#=
graphic=LCD.plotLyapunovConstants_cubicPolynomials_Newton((0.82,0.865),(0.795,0.832),200,10,8)
img=graphic[1]
matrix=graphic[2]
=#

# Compute the Lyapunov constants of a region of the parameter space induced by the application of Chebyshev's method to cubic polynomials
#=
graphic=LCD.plotLyapunovConstants_cubicPolynomials_Chebyshev((-0.74,-0.68),(1.95,2.01),200,5,8;fcp="minus")
img=graphic[1]
matrix=graphic[2]
=#

# Add color bar
#PyPlot.colorbar(img,ticks=0:numColors)

# Print figure
#gcf()
clf(); # Can be useful to clear the current plot (usually when you plot n-cycles, subplots or color bars).

# Save figure
#=
savefig(fname="plot.png")
println("Figure saved succesfully!")
=#

# Save 4K resolution image ("img" is an integer matrix)
#=
FileIO.save("plot.jpg",img)
println("Figure saved succesfully!")
=#