module LCD # Lyapunov Cycle Detector for iteration of rational complex functions.


using Polynomials
using Colors
using PyPlot

# Bijection from C^2 (each representative (z,t) of each class [z:t] of P^1+(C)) to S^2+
# (the augmented 2-sphere).
function P1toSphere(point::Tuple{Complex{Float64},Complex{Float64}},precision::Int64=8)
  # Pre:
  # Post: Returns the coordinates in the augmented 2-sphere of the given complex tuple.
  z=point[1]
  t=point[2]
  den=real(conj(t)*t + conj(z)*z)
  if den<(1.0/10^precision)
      sphPoint=[0.0,0.0,0.0]
  else
      sphPoint=[real((conj(z)*t + conj(t)*z)/den),
      real((1*im*(conj(z)*t - conj(t)*z))/den),
      real((-conj(t)*t + conj(z)*z)/den)]
  end
  return sphPoint
end


#Example:
#P1toSphere((1.0+0.0*im, 0.0+1.0*im))


# Distance between two points of C^2. The pseudo-metric is induced by the chordal metric on
# P^1+(C), given by the previous bijection, considering on S^2+ the metric induced by the
# euclidean metric on R^3.
function chordalMetric(point1::Tuple{Complex{Float64},Complex{Float64}},
  point2::Tuple{Complex{Float64},Complex{Float64}},precision::Int64=8)
  # Pre:
  # Post: Returns the distance between the coordinates on the augmented 3-sphere of 2
  # points of C^2.
  norm=vector->sqrt(vector[1]^2+vector[2]^2+vector[3]^2)
  return norm((P1toSphere(point1,precision)-P1toSphere(point2,precision)))
end


#Example
#chordalMetric((1.0+0.0*im, 0.0+1.0*im),(0.0+0.0*im, 1.0+0.0*im))


# Builds a collection of points of a rectangle determined by the given intervals and precision,
# respecting the aspect ratio between the given intervals. The given precision is used to
# determine the number of points in the grid.
# Note that the collection contains P^1+(C) points [z:t] that, considering the P^1+(C) to C*
# (that is the Alexandroff compactification of C) bijection, correspond to z in C*.
function rectangle(xinterval::Tuple{Float64,Float64}=(-1.5,1.5),yinterval::Tuple{Float64,Float64}=(-1.5,1.5),
  precision=6)
  # Pre:
  # Post: Returns the collection of points.
  resolution=9.443329*(10.0^(precision)) # We take this 9-million-pixels-resolution as a reference.
  a=xinterval[1]
  b=xinterval[2]
  c=yinterval[1]
  d=yinterval[2]
  p1=floor(sqrt(resolution*((b-a)/(d-c)))) # Note that resolution=p1*p2
  p2=floor(p1*(d-c)/(b-a))
  collection=[(complex(r,i),complex(1.0)) for i=d:-(d-c)/p2:c, r=a:(b-a)/p1:b]
  return collection
end


#Example
#R2=rectangle((0.79,0.86),(0.79,0.86))


# Builds a collection of points of a rectangle determined by the given intervals and precision,
# respecting the aspect ratio between the given intervals. The given precision is used to
# determine the number of points in the grid.
# Note that the resulting collection is similar to the collection in the previous method,
# but in this case, each point [z:t] in the collection correspond to 1/t in C*; i.e. the 
# resulting collection of points is a (discretized) neighborhood of infinity in C*.
function invertedRectangle(xinterval::Tuple{Float64,Float64}=(-1.5,1.5),yinterval::Tuple{Float64,Float64}=(-1.5,1.5),
  precision=6)
  # Pre:
  # Post: Returns the collection of points.
  resolution=9.443329*(10.0^(precision)) # We take this 9-million-pixels-resolution as a reference.
  a=xinterval[1]
  b=xinterval[2]
  c=yinterval[1]
  d=yinterval[2]
  p1=floor(sqrt(resolution*((b-a)/(d-c))))
  p2=floor(p1*(d-c)/(b-a))
  collection=[(complex(1.0),complex(r,i)) for i=d:-(d-c)/p2:c, r=a:(b-a)/p1:b]
  return collection
end


#Example
#InvR2=invertedRectangle((0.79,0.86),(0.79,0.86))


# Builds a collection of points of a rectangle determined by the given intervals and precision.
# Note that the collection contains P^1+(C) points [z:t] that, considering the P^1+(C) to C*
# (that is the Alexandroff compactification of C) bijection, correspond to z in C*.
# The term rigid references that this method works fine with medium-sized rectangles (such that
# (-1.5,1.5)x(-1.5,1.5)), but doesn't work as intended with very small or very large rectangles. 
function rigidRectangle(xinterval::Tuple{Float64,Float64}=(-1.5,1.5),yinterval::Tuple{Float64,Float64}=(-1.5,1.5),
  precision=10)
  # Pre:
  # Post: Returns the collection of points.
  tol=1.0/(2^precision)
  a=xinterval[1]
  b=xinterval[2]
  c=yinterval[1]
  d=yinterval[2]
  collection=[(complex(r,i),complex(1.0)) for i=d:-tol:c, r=a:tol:b]
  return collection
end


#Example
#R1=rigidRectangle((-1.5,1.5),(-1.5,1.5))


# Builds a collection of points of a rectangle determined by the given intervals and precision.
# Note that the resulting collection is similar to the collection in the previous method,
# but in this case, each point [z:t] in the collection correspond to 1/t in C*; i.e. the 
# resulting collection of points is a (discretized) neighborhood of infinity in C*.
# The term rigid references that this method works fine with medium-sized rectangles (such that
# (-1.5,1.5)x(-1.5,1.5)), but doesn't work as intended with very small or very large rectangles. 
function rigidInvertedRectangle(xinterval::Tuple{Float64,Float64}=(-1.5,1.5),yinterval::Tuple{Float64,Float64}=(-1.5,1.5),
  precision=10)
  # Pre:
  # Post: Returns the collection of points.
  tol=1.0/(2^precision)
  a=xinterval[1]
  b=xinterval[2]
  c=yinterval[1]
  d=yinterval[2]
  invcollection=[(complex(1.0), complex(r,i)) for i=d:-tol:c, r=a:tol:b]
  return invcollection
end


#Example
#InvR1=rigidInvertedRectangle((-2.56,2.56),(-2.56,2.56),3)


# Applies an homogeneous normalization to a given point of C^2, given a nullity tolerance.
function homogeneousNormalization(point::Tuple{Complex{Float64},Complex{Float64}},
  tolerance::Int64=15)
  # Pre:
  # Post: Returns the homogeneous coordinates of a point of C^2.
  z=point[1]
  t=point[2]
  den=abs(z)+abs(t)
  if den<1.0/10^tolerance
      hpoint=(complex(0.0),complex(0.0))
  else
      hpoint=(z/den, t/den)
  end
  return hpoint
end


#Examples
#=
homogeneousNormalization((1.0+0.0*im, 0.0+1.0*im),9)

homogeneousNormalization((0.000000000000000000001+0.0*im,0.0+0.0*im),9)

homogeneousNormalization((1.0+0.0*im, 0.0+1.0*im),9)
=#


# Ensures that the two given lists have the same length, adding complex zeros if necessary.
function sameLength(list1,list2)
  # Pre: 
  # Post: The lists have been modified if necessary to have the same length.
  l1=length(list1)
  l2=length(list2)
  if l1!=l2
      if l1<l2
       while l1<l2
           push!(list1,complex(0.0))
             l1=l1+1
          end
      else
          while l2<l1
             push!(list2,complex(0.0))
             l2=l2+1
          end
      end
  end
  le=l1 # Note that at this point, it is ensured that l1=l2.
  while le>1 && list1[le]==complex(0.0) && list2[le]==complex(0.0)
      pop!(list1)
      pop!(list2)
      le=le-1
  end
end


#Example
#=
list1=[1.0,2.0]
list2=[1.0,2.0,3.0,4.0,5.0,0.0]
sameLength(list1,list2)
=#


# Normalizes a given homogeneous pair of polynomials, assuming that the norm of the pair
# is defined as the sum in absolute value of the coefficients of both polynomials.
function normalize_pairhpolynomials(numerator::Union{Vector{Float64},Vector{Complex{Float64}}},
  denominator::Union{Vector{Float64},Vector{Complex{Float64}}})
  # Pre:
  # Post: Retuns the normalized numerator and denominator of the given polynomials.
  norm=sum(abs,numerator)+sum(abs,denominator)
  newnum=numerator*(1.0/norm)
  newden=denominator*(1.0/norm)
  return [newnum,newden]
end


#Example
#=
tradNum=complex([1.0,0.0,0.0,2.0]) # f(z)=z^3-1
tradDen=complex([0.0,0.0,3.0])

num2cycle=complex([-2.0,0.0,0.0,2.0]) # f(z)=z^3-2z+2 (2-cycle in (complex(0.0),complex(1.0)))
den2cycle=complex([-2.0,0.0,3.0])

tradPair=normalize_pairhpolynomials(tradNum,tradDen)

pair2cycle=normalize_pairhpolynomials(num2cycle,den2cycle)
=#


# Evaluates in z and t a given homogeneous bivariate polynomial of degree d. 
function evaluate_hpoly(coefficientlist::Array{Complex{Float64},1}, 
  point::Tuple{Complex{Float64},Complex{Float64}}, d::Int)::Complex{Float64}
  # Pre: n=degree of the poly=l-1=length(coefficientlist)-1<=d
  # Post: Returns a complex value.
  l=length(coefficientlist)
  @assert l-1<=d
  z=point[1]
  t=point[2]
  complexcoefficients=complex(coefficientlist)
  p=complexcoefficients[1]*t^d
  for i in 2:l
    p=p+complexcoefficients[i]*z^(i-1)*t^(d-i+1)
  end
  return p
end


#Examples
#=
normTradNum=complex(tradPair[1])
normTradDen=complex(tradPair[2])
d=max(length(normTradNum),length(normTradDen))
evaluate_hpoly(normTradNum,(complex(3.0),complex(1.0)),d)

normNum2cycle=complex(pair2cycle[1])
normDen2cycle=complex(pair2cycle[2])
d=max(length(normNum2cycle),length(normDen2cycle))
evaluate_hpoly(normNum2cycle,(complex(3.0),complex(1.0)),d)
=#


# Evaluates in z and t a given homogeneous pair of bivariate polynomials
# that define a rational function, and returns the corresponding point [z:t]
# in P^{1+}(C). 
function nextpoint_rationalFunction(coefficientlistnum::Array{Complex{Float64},1},coefficientlistden::Array{Complex{Float64},1},
  point::Tuple{Complex{Float64},Complex{Float64}})::Tuple{Complex{Float64},Complex{Float64}}
  # Pre:
  # Post: Returns a complex tuple.
  d=max(length(coefficientlistnum),length(coefficientlistden))-1
  return homogeneousNormalization((evaluate_hpoly(coefficientlistnum,point,d),evaluate_hpoly(coefficientlistden,point,d)))
end


#Examples
#nextpoint_rationalFunction(normTradNum,normTradDen,(complex(3.0),complex(1.0)))


# Evaluates in z and t a given homogeneous pair of bivariate polynomials
# that define a rational function. 
function evaluate_rationalFunction(coefficientlistnum::Array{Complex{Float64},1},coefficientlistden::Array{Complex{Float64},1},
  point::Tuple{Complex{Float64},Complex{Float64}})::Complex{Float64}
  # Pre:
  # Post: Returns a complex value.
  d=max(length(coefficientlistnum),length(coefficientlistden))
  return evaluate_hpoly(coefficientlistnum,point,d)/evaluate_hpoly(coefficientlistden,point,d)
end


#Example
#evaluate_rationalFunction(normTradNum,normTradDen,(complex(3.0),complex(1.0)))


# Evaluates in (z,t) (representative of [z:t] from P^1(C)) the spherical derivative 
# of a given rational polynomical function (determined by the coefficients of its
# numerator and denominator).
# Note that, in the case that abs(t)~0, n (and m) are the degree of the numerator 
# (or denominator, respectively) plus 1, just for the shake of the simplicity of
# the code. This decision doesn't affect the result of the method, since the key
# for the calculation is the distance between n and m, not their specific values.
# This is because, in Julia, arrays begin in the position 1. 
# Also for the shake of simplicity, i'll denote the spherical derivative (associated to a 
# given rational polynomical function) in a point z in P^1(C) just by phi(z).
function sphericalDerivative(coefficientlistnum::Array{Complex{Float64},1},coefficientlistden::Array{Complex{Float64},1},
  point::Tuple{Complex{Float64},Complex{Float64}},tolerance::Int64=8)
  # Pre:
  # Post: Returns a non-negative real value.
  tol=1.0/10^tolerance
  z=point[1]
  t=point[2]
  n=length(coefficientlistnum)
  m=length(coefficientlistden)
  if abs(t)>tol
    d=max(n-1,m-1)
    u1=0.0
    @inbounds for i in 0:n-2
      @inbounds for j in 0:m-1
        u1=u1+(i+1)*coefficientlistnum[i+2]*coefficientlistden[j+1]*z^(i+j)*t^(n+m-3-i-j)
      end
    end
    u2=0.0
    @inbounds for i in 0:n-1
      @inbounds for j in 0:m-2
        u2=u2+(j+1)*coefficientlistnum[i+1]*coefficientlistden[j+2]*z^(i+j)*t^(n+m-3-i-j)
      end
    end
    u3=t*conj(t)+z*conj(z)
    u=abs(t)^(2*d-n-m+1)*(u1-u2)*u3
    v1=0.0
    @inbounds for i in 0:m-1
      @inbounds for j in 0:m-1
        v1=v1+coefficientlistden[i+1]*conj(coefficientlistden[j+1])*z^(i)*t^(d-i)*conj(z)^(j)*conj(t)^(d-j)
      end
    end
    v2=0.0
    @inbounds for i in 0:n-1
      @inbounds for j in 0:n-1
        v2=v2+coefficientlistnum[i+1]*conj(coefficientlistnum[j+1])*z^(i)*t^(d-i)*conj(z)^(j)*conj(t)^(d-j)
      end
    end
    v=v1+v2
    return abs(u/v)
  else
    if abs(n-m)>1
      return 0
    else
      if n==m
        pnum=Polynomial(coefficientlistnum)
        pden=Polynomial(coefficientlistden)
        pdernum=derivative(pnum)*pden-pnum*derivative(pden)
        coeffs_pdernum=complex(coeffs(pdernum))
        e_star=length(coeffs_pdernum)-1
        return abs(t^(2*n-4-e_star)*(coeffs_pdernum[e_star+1]*z^(e_star+1)*conj(z))/(z^(n-1)*conj(z^(n-1))*(coefficientlistnum[n]*conj(coefficientlistnum[n])+coefficientlistden[m]*conj(coefficientlistden[m]))))
      elseif n>m
        return abs(coefficientlistden[m]/conj(coefficientlistnum[n]))
      elseif m>n
        return abs(coefficientlistnum[n]/conj(coefficientlistden[m]))
      end
    end
  end
end


#Examples
#=
sphericalDerivative(normTradNum,normTradDen,(0.5+1.0*im,complex(1.0)))

sphericalDerivative(normTradNum,normTradDen,(complex(3.0),complex(1.0)))

sphericalDerivative(normTradNum,normTradDen,(complex(2.0),complex(0.0)))
=#


# Computes L^[0,numprod](phi) in a point of P^{1+}(C). It is important
# to note that, if the given point is equal to f^k(z) for some z in P^{1+}(C), 
# this method returns L^[k,k+numprod](phi)=(phi(f^k(z))phi(f^(k+1)(z))⋯phi(f^(k+numprod)(z)))^(1/numprod),
# i.e. a product of numprod factors (each factor is the spherical derivative in the point f^i(z),
# with i=k,...,k+numprod), to the power of 1/numprod.
function lyapunov_product(coefficientlistnum::Array{Complex{Float64},1},coefficientlistden::Array{Complex{Float64},1},
  point::Tuple{Complex{Float64},Complex{Float64}},numprod::Int64,tolerance::Int64=8)
  # Pre:
  # Post: Returns a positive real value.
  l=1.0
  p=point
  for i in 1:numprod
      l=l*sphericalDerivative(coefficientlistnum,coefficientlistden,p,tolerance)
      p=nextpoint_rationalFunction(coefficientlistnum,coefficientlistden,p)
  end
  return l^(1.0/numprod)
end


#Examples
#lyapunov_product(normTradNum,normTradDen,homogeneousNormalization((complex(3.0),complex(1.0))),7)


# Detects if the orbit of a given point (generated by iterating the given rational function) converges to 
# an n-cycle (given its length n), using the Lyapunov function L^[0,k](phi) and some of its properties.
# The resulting array contains TRUE in case an n-cycle is detected (and FALSE otherwise), 
# the n-cycle itself (if its the case), the number of iterations it took to detect the n-cycle, 
# and the constant associated with the n-cycle, given by the last computed value of the lyapunov function.
# Note that the resulting number of iterations will always be greater or equal to the length of the
# cycle we are trying to detect, because it will take at the very least n iterations to compute the 
# n-cycle itself.
# For example, to detect the 2-cycle (0,1), the minimum number of iterations required to detect it
# is 2.
function detect_nCycle_Lyapunov(coefficientlistnum::Array{Complex{Float64},1},coefficientlistden::Array{Complex{Float64},1},
  point::Tuple{Complex{Float64},Complex{Float64}},cycleLength::Int64,maxiterations::Int64=25,tolerance::Int64=8)
  # Pre: cycleLength>=1
  # Post: Returns an array containing a boolean, the n-cycle (a vector of P^{1+}(C) points),
  #       an integer and a real positive value.
  tol=1.0/10^tolerance
  p=point
  nextp=nextpoint_rationalFunction(coefficientlistnum,coefficientlistden,p)
  cycleDetected=false
  i=0
  while i<=maxiterations-cycleLength && !cycleDetected
    if abs(lyapunov_product(coefficientlistnum,coefficientlistden,p,cycleLength,tolerance)-lyapunov_product(coefficientlistnum,coefficientlistden,nextp,cycleLength,tolerance))<tol
      cycleDetected=true # In this case, the orbit of the given point has converged to an attracting cycleLength-cycle.
      cycle=[]
      pcycle=nextp
      for k in 1:cycleLength
        push!(cycle,pcycle)
        pcycle=nextpoint_rationalFunction(coefficientlistnum,coefficientlistden,pcycle)
      end
      result=[cycleDetected,cycle,i+cycleLength,lyapunov_product(coefficientlistnum,coefficientlistden,nextp,cycleLength,tolerance)]
    else # If the orbit has not converged yet, we moved onto the next point.
      p=nextp
      nextp=nextpoint_rationalFunction(coefficientlistnum,coefficientlistden,p)
    end
    i=i+1
  end
  if !cycleDetected # If the algorithm has ended without detecting any convergence:
    result=[cycleDetected,[],maxiterations,lyapunov_product(coefficientlistnum,coefficientlistden,p,cycleLength,tolerance)]
  end
  if result[4]<tol # Finally, if the computed constant is sufficiently small, up to some given tolerance, we take it as 0.
    result[4]=0.0
  end
  return result
end


#Examples
#=
detect_nCycle_Lyapunov(normTradNum,normTradDen,homogeneousNormalization((0.75-0.4*im,complex(1.0))),1)

detect_nCycle_Lyapunov(normTradNum,normTradDen,homogeneousNormalization((0.75-0.4*im,complex(1.0))),1)
# Note that, if the orbit of a point converges to a k-cycle, it will also converge to a (k+1)-cycle.

detect_nCycle_Lyapunov(normNum2cycle,normDen2cycle,(complex(0.0),complex(1.0)),1)

detect_nCycle_Lyapunov(normNum2cycle,normDen2cycle,(complex(0.0),complex(1.0)),2)
=#


# Given a point, a rational function and the maximum length of the n-cycles we
# want to detect, it searches the orbit of that point, described by the rational
# function, for n-cycles where n is between 1 and maxcycle.
# It is important to note that this method only 
# detects n-cycles for n<=maxcycle, and computes the n-cycle up to a given tolerance, 
# taking into account a given maximum number of iterations.
function cycleDetector_Lyapunov(coefficientlistnum::Array{Complex{Float64},1},coefficientlistden::Array{Complex{Float64},1},
  point::Tuple{Complex{Float64},Complex{Float64}},maxiterations::Int64=25,maxcycle::Int64=5,tolerance::Int64=8)
  # Pre: maxcycle>=1
  # Post: Returns a vector with the n-cycle (a vector of P^{1+}(C) points),
  #       an integer and a real positive value.
  p=point
  n=1
  cycleDetected=false
  ncycle=[false,[],maxiterations,-1]
  while n<=maxcycle && !cycleDetected
    ncycle=detect_nCycle_Lyapunov(coefficientlistnum,coefficientlistden,p,n,maxiterations,tolerance)
    if ncycle[1]==true  # If an n-cycle has been detected:
      cycleDetected=true
    else
      n=n+1
    end
  end
  result=[ncycle[2],ncycle[3],ncycle[4]]
  return result
end


#Examples
#cycleDetector_Lyapunov(normTradNum,normTradDen,homogeneousNormalization((complex(3.0),complex(1.0))))


# Checks if the given complex value is contained in the given list, up to
# some given tolerance.
function contains(list::Array{Any,1},element::Complex{Float64},tolerance::Int64=8)
  # Pre: list is an array of complex values.
  # Post: Returns TRUE if it is contained and FALSE in other case, and the position of the list
  #       where it was found.
  tol=1.0/10^tolerance
  l=length(list)
  contained=false
  i=1
  while i<=l && !contained
    if abs(list[i]-element)<tol
      contained=true
    else
      i=i+1
    end
  end
  return [contained,i]
end

# Given a rectangular grid of points in P^1(C) and a rational function,
# for each point of the grid, a vector with 2 components is computed.
# The first component is the "position" of the point; that is the 
# position of the super-attracting fixed point to which the point converges
# (under iteration of the rational function) in the computed list of fixed points,
# which is also returned as a result. The second component corresponds to the number 
# of iterations of the rational function that it took the point to converge.
# Note that if a point does not converge to any super-attracting fixed point 
# (and converges to a n-cycle with n>1), its position is considered to be 0.
# Also note that, since we are working on the Riemann sphere C*, the infinity
# point can also be considered as a fixed point with its own position and basin
# of attraction, if that's the case.
function positionIteration_superAttractingFixedPoints(coefficientlistnum::Array{Complex{Float64},1},coefficientlistden::Array{Complex{Float64},1},
  xinterval::Tuple{Float64,Float64},yinterval::Tuple{Float64,Float64},maxiterations::Int64=25,maxcycle::Int64=5,tolerance::Int64=8;inverted::Bool=false)
  # Pre: maxcycle>=1
  # Post: Returns a matrix of vectors with 2 components (position and iteration),
  #       the list of super-attracting fixed points and whether if infinity is
  #       a fixed point.
  if inverted==true
    rect=invertedRectangle(xinterval,yinterval)
  else
    rect=rectangle(xinterval,yinterval)
  end
  a=size(rect)[1]
  b=size(rect)[2]
  tol=1.0/10^tolerance
  d=max(length(coefficientlistnum),length(coefficientlistden))-1
  fixedPoints=[]
  if abs(evaluate_hpoly(coefficientlistnum,(complex(1.0),complex(0.0)),d))>tol && abs(evaluate_hpoly(coefficientlistden,(complex(1.0),complex(0.0)),d))<tol
    infFixedPoint=true # We check whether in this case the infinity is a fixed point of the given rational function.
    push!(fixedPoints,Inf)
  else
    infFixedPoint=false
  end
  positer=Array{Vector{Int64}}(undef,a,b)
  @inbounds Threads.@threads for j in 1:b
    @inbounds for i in 1:a
      cycle=cycleDetector_Lyapunov(coefficientlistnum,coefficientlistden,rect[i,j],maxiterations,maxcycle,tolerance)
      if infFixedPoint==true && length(cycle[1])==1 && abs(cycle[1][1][2])<tol # In case the orbit of the point converges to infinity.
        positer[i,j]=[1.0,cycle[2]] # Note that, if infinity is a fixed point, the position 1 on the list of fixed points is reserved to it.
      elseif cycle[3]<tol && length(cycle[1])==1 && abs(cycle[1][1][2])>tol # In case the orbit of the point converges to another super-attracting fixed point.
        fPoint=cycle[1][1][1]/cycle[1][1][2]
        if real(fPoint)<tol
          fPoint=imag(fPoint)*im
        elseif imag(fPoint)<tol
          fPoint=real(fPoint)
        end
        pos=contains(fixedPoints,complex(fPoint),tolerance)
        if pos[1]==false
          push!(fixedPoints,complex(fPoint))
          positer[i,j]=[length(fixedPoints),cycle[2]]
        else
          positer[i,j]=[pos[2],cycle[2]]
        end
      else # In case the orbit of the point does not converge to any point (1-cycle).
        positer[i,j]=[0.0,cycle[2]]
      end
    end
  end
  return [positer,fixedPoints,infFixedPoint]
end


#Examples
#positionIteration_superAttractingFixedPoints(normTradNum,normTradDen,(-1.5,1.5),(-1.5,1.5))


# This method works exactly in the same way that the previous one, but with the list of fixed
# points given beforehand. 
# Note that this collection of algorithms do not need any previous knowledge of the fixed points 
# in order to work, so the previous method should be used in almost every case instead of this one.
# This method was created to preserve the consistency of the color palette between different graphics.
# For example, if you generate a graphic using the previous method, and then one wants to generate a
# new one zooming in a particular area (or taking a neighborhood of infinity), the previous method
# will use the same colors that were used in the previous graphic, but it might assign them to a 
# different basin, which may lead to misunderstandings when it comes to interpret the graphics.
# On the contrary, using this method is considerably more time-consuming, but it ensures that every
# color in the graphics are consistent, no matter the grid that we choose.
function positionIteration_superAttractingFixedPoints_givenFixedPoints(coefficientlistnum::Array{Complex{Float64},1},coefficientlistden::Array{Complex{Float64},1},
  xinterval::Tuple{Float64,Float64},yinterval::Tuple{Float64,Float64},fixedPoints::Array{Any,1},maxiterations::Int64=25,maxcycle::Int64=5,
  tolerance::Int64=8;inverted::Bool=false)
  # Pre: maxcycle>=1
  # Post: Returns a matrix of vectors with 2 components (position and iteration),
  #       the list of super-attracting fixed points and whether if infinity is
  #       a fixed point.
  if inverted==true
    rect=invertedRectangle(xinterval,yinterval)
  else
    rect=rectangle(xinterval,yinterval)
  end
  a=size(rect)[1]
  b=size(rect)[2]
  tol=1.0/10^tolerance
  d=max(length(coefficientlistnum),length(coefficientlistden))-1
  if abs(evaluate_hpoly(coefficientlistnum,(complex(1.0),complex(0.0)),d))>tol && abs(evaluate_hpoly(coefficientlistden,(complex(1.0),complex(0.0)),d))<tol
    infFixedPoint=true # We check whether in this case the infinity is a fixed point of the given rational function.
  else
    infFixedPoint=false
  end
  positer=Array{Vector{Int64}}(undef,a,b)
  @inbounds Threads.@threads for j in 1:b
    @inbounds for i in 1:a
      cycle=cycleDetector_Lyapunov(coefficientlistnum,coefficientlistden,rect[i,j],maxiterations,maxcycle,tolerance)
      if infFixedPoint==true && length(cycle[1])==1 && abs(cycle[1][1][2])<tol # In case the orbit of the point converges to infinity.
        positer[i,j]=[1.0,cycle[2]] # Note that, if infinity is a fixed point, the position 1 on the list of fixed points is reserved to it.
      elseif cycle[3]<tol && length(cycle[1])==1 && abs(cycle[1][1][2])>tol # In case the orbit of the point converges to another super-attracting fixed point.
        fPoint=cycle[1][1][1]/cycle[1][1][2]
        if real(fPoint)<tol
          fPoint=imag(fPoint)*im
        elseif imag(fPoint)<tol
          fPoint=real(fPoint)
        end
        pos=contains(fixedPoints,complex(fPoint),tolerance)
        positer[i,j]=[pos[2],cycle[2]]
      else # In case the orbit of the point does not converge to any point (1-cycle).
        positer[i,j]=[0.0,cycle[2]]
      end
    end
  end
  return [positer,fixedPoints,infFixedPoint]
end


# Computes the position-iteration matrix described in the previous method,
# and transforms it into a matrix of real values, following a given coloring
# strategy. Such values will determine later the color of each corresponding 
# point in the given rectangular grid.
# Also, a color map is generated (taking into account the number of super-attracting
# fixed points found in the process or the maximum number of iterations, 
# when neccesary), and the resulting image is generated.
function colorRectangle_positionIteration(coefficientlistnum::Array{Complex{Float64},1},coefficientlistden::Array{Complex{Float64},1},
  xinterval::Tuple{Float64,Float64},yinterval::Tuple{Float64,Float64},maxiterations::Int64=25,maxcycle::Int64=5,tolerance::Int64=8;
  coloringStrat::AbstractString="positionIteration",inverted::Bool=false)
  # Pre: maxcycle>=1
  # Post: Returns the list of super-attracting fixed points, the seed used to
  #       generate the color map, the resulting image, and the number of colors
  #       used.
  positer=positionIteration_superAttractingFixedPoints(coefficientlistnum,coefficientlistden,xinterval,yinterval,maxiterations,maxcycle,tolerance;inverted)
  numFixedPoints=length(positer[2])
  a=size(positer[1])[1]
  b=size(positer[1])[2]
  if coloringStrat=="positionIteration"
    numColors=numFixedPoints+1
    @inbounds colorMatrix=[(positer[1][i,j][1]+(1.0-(positer[1][i,j][2])/maxiterations))+1.0 for i=1:a,j=1:b]
    colorMatrix[1,1]=1.0
    colorMatrix[a,b]=numColors-2
  elseif coloringStrat=="iteration"
    numColors=maxiterations+1
    @inbounds colorMatrix=[positer[1][i,j][2]+0.25 for i=1:a,j=1:b]
    colorMatrix[1,1]=maxiterations
  elseif coloringStrat=="position"
    numColors=numFixedPoints+1
    @inbounds colorMatrix=[positer[1][i,j][1]+0.25 for i=1:a,j=1:b]
    colorMatrix[1,1]=0.0
    colorMatrix[a,b]=numFixedPoints
  end
  if  positer[3]==false
    seed1=[Colors.RGB(0.0,0.0,0.0),Colors.RGB(0.5,0.5,0.5)]
  else
    seed1=[Colors.RGB(0.0,0.0,0.0),Colors.RGB(0.5,0.5,0.5),Colors.RGB(1.0,1.0,0.0)]
  end
  seed=union(seed1,[Colors.RGB(1.0,0.0,0.0),Colors.RGB(0.0,1.0,0.0),Colors.RGB(0.0,0.0,1.0)])
  gseed=Colors.distinguishable_colors(numColors,seed)
  fpcm=PyPlot.ColorMap(gseed)
  img=PyPlot.imshow(colorMatrix,cmap=fpcm,extent=[xinterval[1],xinterval[2],yinterval[1],yinterval[2]])
  return positer[2],gseed,img,numColors
end


#Examples
#colorRectangle_positionIteration(normTradNum,normTradDen,(-1.5,1.5),(-1.5,1.5);coloringStrat="position")


# The philosophy of this method is exactly the same than positionIteration_superAttractingFixedPoints_givenFixedPoints;
# it works exactly in the same way than the previous one, is considerably more time-consuming, and its purpose is to
# preserve the color consistency between different graphics.
function colorRectangle_positionIteration_givenFixedPoints(coefficientlistnum::Array{Complex{Float64},1},coefficientlistden::Array{Complex{Float64},1},
  xinterval::Tuple{Float64,Float64},yinterval::Tuple{Float64,Float64},fixedPoints::Array{Any,1},maxiterations::Int64=25,maxcycle::Int64=5,
  tolerance::Int64=8;coloringStrat::AbstractString="positionIteration",inverted::Bool=false)
  # Pre: maxcycle>=1
  # Post: Returns the list of super-attracting fixed points, the seed used to
  #       generate the color map, the resulting image, and the number of colors
  #       used.
  positer=positionIteration_superAttractingFixedPoints_givenFixedPoints(coefficientlistnum,coefficientlistden,xinterval,yinterval,fixedPoints,maxiterations,maxcycle,tolerance;inverted)
  numFixedPoints=length(positer[2])
  a=size(positer[1])[1]
  b=size(positer[1])[2]
  if coloringStrat=="positionIteration"
    numColors=numFixedPoints+1
    @inbounds colorMatrix=[(positer[1][i,j][1]+(1.0-(positer[1][i,j][2])/maxiterations))+1.0 for i=1:a,j=1:b]
    colorMatrix[1,1]=1.0
    colorMatrix[a,b]=numColors-2
  elseif coloringStrat=="iteration"
    numColors=maxiterations+1
    @inbounds colorMatrix=[positer[1][i,j][2]+0.25 for i=1:a,j=1:b]
    colorMatrix[1,1]=maxiterations
  elseif coloringStrat=="position"
    numColors=numFixedPoints+1
    @inbounds colorMatrix=[positer[1][i,j][1]+0.25 for i=1:a,j=1:b]
    colorMatrix[1,1]=0.0
    colorMatrix[a,b]=numFixedPoints
  end
  if  positer[3]==false
    seed1=[Colors.RGB(0.0,0.0,0.0),Colors.RGB(0.5,0.5,0.5)]
  else
    seed1=[Colors.RGB(0.0,0.0,0.0),Colors.RGB(0.5,0.5,0.5),Colors.RGB(1.0,1.0,0.0)]
  end
  seed=union(seed1,[Colors.RGB(1.0,0.0,0.0),Colors.RGB(0.0,1.0,0.0),Colors.RGB(0.0,0.0,1.0)])
  gseed=Colors.distinguishable_colors(numColors,seed)
  fpcm=PyPlot.ColorMap(gseed)
  img=PyPlot.imshow(colorMatrix,cmap=fpcm,extent=[xinterval[1],xinterval[2],yinterval[1],yinterval[2]])
  return positer[2],gseed,img,numColors
end


# Given the list of coefficients of a complex polynomial, this method applies the
# complex Newton-Raphson method to it and returns the obtained rational function.
function complexNewtonRaphson_method(coefficientlistpol::Array{Complex{Float64},1})
  # Pre:
  # Post: Retuns 2 lists of complex coefficients.
  p=Polynomial(coefficientlistpol)
  polX=Polynomial([0.0,1.0])
  pnum=polX*derivative(p)-p  # Applies the complex Newton-Raphson method.
  pden=derivative(p)
  coefficientlistnum=complex(coeffs(pnum))
  coefficientlistden=complex(coeffs(pden))
  return coefficientlistnum,coefficientlistden
end


#Example
#pair=complexNewtonRaphson_method(tradPol)


# Given the list of coefficients of a complex polynomial, this method applies the
# complex Chebyshev method to it and returns the obtained rational function.
function complexChebyshev_method(coefficientlistpol::Array{Complex{Float64},1})
  # Pre:
  # Post: Retuns 2 lists of complex coefficients.
  p=Polynomial(coefficientlistpol)
  derP=derivative(p)
  polX=Polynomial([0.0,1.0])
  pnum=2*polX*derP^3-2*p*derP^2-derivative(derP)*p^2  # Applies the complex Chebyshev method.
  pden=2*derP^3
  coefficientlistnum=complex(coeffs(pnum))
  coefficientlistden=complex(coeffs(pden))
  return coefficientlistnum,coefficientlistden
end


#Example
#pair=complexChebyshev_method(tradPol)


# Given the list of coefficients of a complex polynomial, this method applies the
# complex Halley method to it and returns the obtained rational function.
function complexHalley_method(coefficientlistpol::Array{Complex{Float64},1})
  # Pre:
  # Post: Retuns 2 lists of complex coefficients.
  p=Polynomial(coefficientlistpol)
  derP=derivative(p)
  polX=Polynomial([0.0,1.0])
  pnum=2*polX*derP^2-derivative(derP)*p*polX-2*p*derP  # Applies the complex Halley method.
  pden=2*derP^2-p*derivative(derP)
  coefficientlistnum=complex(coeffs(pnum))
  coefficientlistden=complex(coeffs(pden))
  return coefficientlistnum,coefficientlistden
end


#Example
#pair=complexHalley_method(tradPol)


# Given the list of coefficients of a complex polynomial, applies the
# chosen method to it and returns the obtained rational function.
function applyIterative_method(coefficientlistpol::Array{Complex{Float64},1},method::AbstractString="newton")
  # Pre: Only Newton's, Chebyshev's and Halley's method are currently supported.
  # Post: Retuns 2 lists of complex coefficients.
  if method=="newton"
    pair=complexNewtonRaphson_method(coefficientlistpol)
  elseif method=="chebyshev"
    pair=complexChebyshev_method(coefficientlistpol)
  elseif method=="halley"
    pair=complexHalley_method(coefficientlistpol)
  else
    error("The chosen method is not currently supported.")
  end
  return pair
end


#Example
##pair=applyIterative_method(tradPol)


# Plots the computed basins of attraction of super-attracting fixed
# points of the given rational function by calling the previous method.
function plotBasinsOfAttraction_Lyapunov(coefficientlistnum::Array{Complex{Float64},1},coefficientlistden::Array{Complex{Float64},1},
  xinterval::Tuple{Float64,Float64},yinterval::Tuple{Float64,Float64},maxiterations::Int64=25,maxcycle::Int64=5,tolerance::Int64=8;
  coloringStrat::AbstractString="positionIteration",inverted::Bool=false)
  # Pre: maxcycle>=1
  # Post: Returns an image and plots it by using PyPlot.
  #       Also returns the number of plotted colors.
  coloredRect=colorRectangle_positionIteration(coefficientlistnum,coefficientlistden,xinterval,yinterval,maxiterations,maxcycle,tolerance;coloringStrat,inverted)
  return coloredRect[3], coloredRect[4] # Returns the resulting image, and the number of colors used.
end


#Examples
#=
tradPol=complex([-1.0,0.0,0.0,1.0])
pol2cycle=complex([2.0,-2.0,0.0,1.0])
img=plotBasinsOfAttraction_Lyapunov(normTradNum,normTradDen,(-1.5,1.5),(-1.5,1.5),25,1,8;coloringStrat="position")[1]

gcf()
=#


# The philosophy of this method is exactly the same than positionIteration_superAttractingFixedPoints_givenFixedPoints;
# it works exactly in the same way than the previous one, is considerably more time-consuming, and its purpose is to
# preserve the color consistency between different graphics.
function plotBasinsOfAttraction_Lyapunov_givenFixedPoints(coefficientlistnum::Array{Complex{Float64},1},coefficientlistden::Array{Complex{Float64},1},
  xinterval::Tuple{Float64,Float64},yinterval::Tuple{Float64,Float64},fixedPoints::Array{Any,1},maxiterations::Int64=25,maxcycle::Int64=5,tolerance::Int64=8;
  coloringStrat::AbstractString="positionIteration",inverted::Bool=false)
  # Pre: maxcycle>=1
  # Post: Returns an image and plots it by using PyPlot.
  #       Also returns the number of plotted colors.
  coloredRect=colorRectangle_positionIteration_givenFixedPoints(coefficientlistnum,coefficientlistden,xinterval,yinterval,fixedPoints,maxiterations,maxcycle,tolerance;coloringStrat,inverted)
  return coloredRect[3], coloredRect[4] # Returns the resulting image, and the number of colors used.
end


#=
The previous lines of code are completely functional on itselves. You could use them 
to generate Newton fractals and to compute the basins of attraction of super-attracting
fixed points of any rational complex function.
The following part of the code is destined specifically to provide methods to compute
and plot the basins of attraction of attracting n-cycles of any rational complex function.
We begin this second part of the code describing a few auxiliary methods.
=#


# Checks if the infinity point is contained in a given list of P^1(C) points,
# up to some given tolerance.
function IsInfIn(list::Array{Any},tolerance::Int64=8)
  # Pre:
  # Post: Returns TRUE if it's the case and FALSE if not.
  tol=1.0/10^tolerance
  n=length(list)
  i=1
  inf=false
  while i<=n && !inf
    if abs(list[i][2])<tol
      inf=true
    else
      i=i+1
    end
  end
  return inf
end


# Auxiliar method that transforms a list of points of P^1(C) (which corresponds
# to an n-cycle) into a list of complex values.
function P1toComplex_cycle(cycle::Array{Any})
  # Pre: Forall i in 1:n, abs(cycle[i][2])>0; i.e. the infinity point is not contained in the given cycle.
  # Post: Returns a list of complex values.
  complexCycle=[]
  l=length(cycle)
  for i in 1:l
    push!(complexCycle,cycle[i][1]/cycle[i][2])
  end
  return complexCycle
end


# Checks if a given list of lists (nCycles) contains another given list (cycle).
function contains_list(nCycles::Array{Any},cycle::Array{Any},tolerance::Int64=6)
  # Pre: All the lists involved contain complex values.
  # Post: Returns TRUE if nCycles contains the list cycle, and FALSE in other case.
  numCycles=length(nCycles)
  if numCycles==0
    return false
  else
    contained=false
    l=length(cycle)
    i=1
    while i<=numCycles && !contained # Checks whether if one of the cycles in nCycles is equivalent to the given cycle.
      ok=true
      j=1
      while j<=l && ok
        if length(nCycles[i])!=l
          ok=false
        else
          containsElement=contains(nCycles[i],cycle[j],tolerance)
          if containsElement[1]==false
            ok=false
          else
            j=j+1
          end
        end
      end
      if ok
        contained=true
      else
        i=i+1
      end
    end
    return contained
  end
end


# When a cycle has an element that appears in the cycle 2 or more times we say that
# the cycle is degenerated, essentially because, if said cycle is has n elements and
# it is degenerated, it means that it represents an m-cycle with m<n.
# Our interest here is to compute only the "true" or "canonical" n-cycles, given
# a certain n.
# Note that if we have an n-cycle, we can naturally extend it to a degenerated (2n)-cycle. 
function degeneratedCycle(cycle::Array{Any},tolerance::Int64=8)
  # Pre:
  # Post: Returns TRUE if the given cycle is degenerated, and FALSE in other case.
  degenerated=false
  n=length(cycle)
  i=1
  while i<=n && !degenerated
    if contains(cycle,cycle[i],tolerance)[1]==true && contains(cycle,cycle[i],tolerance)[2]!=i # If there is an element of the cycle that appears 2 or more times in it.
      degenerated=true
    else
      i=i+1
    end
  end
  return degenerated
end


# Computes all the attracting n-cycles (for a given n) detected for a given rational 
# function in a given rectangle; that is, if there is some point in the given rectangle 
# whose orbit converges to an attracting n-cycle.
function getAttractingNCycles_Lyapunov(coefficientlistnum::Array{Complex{Float64},1},coefficientlistden::Array{Complex{Float64},1},
  xinterval::Tuple{Float64,Float64},yinterval::Tuple{Float64,Float64},maxiterations::Int64=25,n::Int64=2,tolerance::Int64=8)
  # Pre: 
  # Post: Returns a list with all the attracting n-cycles detected in the
  #       given rectangle.
  rect=rectangle(xinterval,yinterval)
  a=size(rect)[1]
  b=size(rect)[2]
  tol=1.0/10^tolerance
  nCycles=[]
  @inbounds Threads.@threads for j in 1:b
    @inbounds for i in 1:a
      cycle=detect_nCycle_Lyapunov(coefficientlistnum,coefficientlistden,rect[i,j],n,maxiterations,tolerance)
      if cycle[1]==true && length(cycle[2])==n && cycle[4]<1 && IsInfIn(cycle[2])==false # If an attracting n-cycle (not including the infinity point) is detected.
        complexCycle=P1toComplex_cycle(cycle[2])
        if !degeneratedCycle(complexCycle,tolerance) && !contains_list(nCycles,complexCycle,6) # Checks that the detected cycle isn't degenerated and hasn't been already computed.
          push!(nCycles,complexCycle)
        end
      end
    end
  end
  return nCycles
end


#Example
#=
pol2cycle=complex([2.0,-2.0,0.0,1.0])
getAttractingNCycles_Lyapunov(pol2cycle,(-1.5,1.5),(-1.5,1.5),25,2,8)
=#


# Computes the list of attracting n-cycles, and plots them into 
# the given rectangle, with a given color.
# Each cycle is represented in the rectangle as a closed polygonal curve
# connecting every element of the cycle. 
function plot_nCycles_Lyapunov(coefficientlistnum::Array{Complex{Float64},1},coefficientlistden::Array{Complex{Float64},1},
  xinterval::Tuple{Float64,Float64},yinterval::Tuple{Float64,Float64},maxiterations::Int64=25,n::Int64=2,tolerance::Int64=8,plotcolor::AbstractString="blue")
  # Pre:
  # Post: Plots the attracting n-cycles using PyPlot.
  nCycles=getAttractingNCycles_Lyapunov(coefficientlistnum,coefficientlistden,xinterval,yinterval,maxiterations,n,tolerance)
  l=length(nCycles)
  for j in 1:l
    cycle=nCycles[j]
    for i in 1:n
      if i==n
        x1=real(cycle[i])
        y1=imag(cycle[i])
        x2=real(cycle[1])
        y2=imag(cycle[1])
      else
        x1=real(cycle[i])
        y1=imag(cycle[i])
        x2=real(cycle[i+1])
        y2=imag(cycle[i+1])
      end
      x_values=[x1,x2]
      y_values=[y1,y2]
      PyPlot.plot(x_values,y_values,color=plotcolor,linestyle="-")
    end
  end
end


#Example
#=
tradPol=complex([-1.0,0.0,0.0,1.0])
pol2cycle=complex([2.0,-2.0,0.0,1.0])
pol3cycle=complex([-1.49175,1.0,0.0,-1.0])
plot_nCycles_Lyapunov(pol2cycle,(-1.5,1.5),(-1.5,1.5),25,2,8)

gcf()

# clf(); can be useful to clear the current plot.
=#


# Computes the list of attracting n-cycles, with n in a given range (from 
# min to max).
function getManyCycles_Lyapunov(coefficientlistnum::Array{Complex{Float64},1},coefficientlistden::Array{Complex{Float64},1},xinterval::Tuple{Float64,Float64},yinterval::Tuple{Float64,Float64},
  maxiterations::Int64=25,min::Int64=2,max::Int64=2,tolerance::Int64=8)
  # Pre:
  # Post: Returns a list with all the attracting n-cycles (min<=n<=max)
  #       detected in the given rectangle.
  rect=rectangle(xinterval,yinterval)
  a=size(rect)[1]
  b=size(rect)[2]
  tol=1.0/10^tolerance
  nCycles=[]
  @inbounds Threads.@threads for j in 1:b
    @inbounds for i in 1:a
      cycle=cycleDetector_Lyapunov(coefficientlistnum,coefficientlistden,rect[i,j],maxiterations,max,tolerance)
      if length(cycle[1])>=min && length(cycle[1])<=max && cycle[3]<1 && IsInfIn(cycle[1])==false # If an attracting n-cycle (not including the infinity point) is detected.
        complexCycle=P1toComplex_cycle(cycle[1])
        if !degeneratedCycle(complexCycle,tolerance) && !contains_list(nCycles,complexCycle,6) # Checks that the detected cycle isn't degenerated and hasn't been already computed.
          push!(nCycles,complexCycle)
        end
      end
    end
  end
  return nCycles
end


#Example
#=
pol2cycle=complex([2.0,-2.0,0.0,1.0])
getManyCycles_Lyapunov(pol2cycle,(-1.5,1.5),(-1.5,1.5),25,2,3,8)
=#


# Computes the list of attracting n-cycles, with n in a given range (from 
# min to max), and plots them into the given rectangle, assigning a different
# color for each group of cycles of the same length.
# Each cycle is represented in the rectangle as a closed polygonal curve
# connecting every element of the cycle.
# Note that there's only 7 possible colors so, if abs(min-max)>7, some colors
# will repeat. In general, the first 7 colors to appear are blue, yellow, magenta,
# cyan, green, black and red, in that specific order (for cycles of length n), and
# then the sequence will repeat for the following groups of n-cycles.
function plot_manyCycles_Lyapunov(coefficientlistnum::Array{Complex{Float64},1},coefficientlistden::Array{Complex{Float64},1},xinterval::Tuple{Float64,Float64},yinterval::Tuple{Float64,Float64},
  maxiterations::Int64=25,min::Int64=2,max::Int64=2,tolerance::Int64=8)
  # Pre:
  # Post: Plots the attracting n-cycles using PyPlot.
  nCycles=getManyCycles_Lyapunov(coefficientlistnum,coefficientlistden,xinterval,yinterval,maxiterations,min,max,tolerance)
  plotcolor=["black","red","blue","yellow","magenta","cyan","green"]
  l=length(nCycles)
  for j in 1:l
    cycle=nCycles[j]
    numElements=length(cycle)
    for i in 1:numElements
      if i==numElements
        x1=real(cycle[i])
        y1=imag(cycle[i])
        x2=real(cycle[1])
        y2=imag(cycle[1])
      else
        x1=real(cycle[i])
        y1=imag(cycle[i])
        x2=real(cycle[i+1])
        y2=imag(cycle[i+1])
      end
      x_values=[x1,x2]
      y_values=[y1,y2]
      PyPlot.plot(x_values,y_values,color=plotcolor[(numElements%7)+1],linestyle="-")
    end
  end
end


#Example
#=
tradPol=complex([-1.0,0.0,0.0,1.0])
pol2cycle=complex([2.0,-2.0,0.0,1.0])
pol3cycle=complex([-1.49175,1.0,0.0,-1.0])
anotherPol2cycle=complex([-5.29,0.0,-4.29,0.0,1.0])
plot_manyCycles_Lyapunov(pol2cycle,(-1.5,1.5),(-1.5,1.5),25,2,3,8)

gcf()

# clf(); can be useful to clear the current plot.
=#


# Plots the Newton fractal (i.e. the basins of attraction of every final point
# of the discrete semi-flux induced by a given complex polynomial and the 
# complex Newton-Raphson method) of a given polynomial, and plots all the 
# attracting n-cycles (with 2<=n<=maxcycle) detected in the 
# given rectangle of the complex plane, up to some given tolerance.
# Each cycle is represented in the rectangle as a closed polygonal curve
# connecting every element of the cycle. Also, for the shake of simplicity,
# every attracting n-cycle is plotted in cyan, a color that is not commongly
# used in the color map generated for the basins of attraction of fixed points,
# when the degree of the given polynomial is reasonably small. 
function plotWithCycles_BasinsOfAttraction_Lyapunov(coefficientlistnum::Array{Complex{Float64},1},coefficientlistden::Array{Complex{Float64},1},
  xinterval::Tuple{Float64,Float64},yinterval::Tuple{Float64,Float64},maxiterations::Int64=25,maxcycle::Int64=5,tolerance::Int64=8;inverted::Bool=false)
  # Pre:
  # Post: Plots the basin of attraction of the roots of the given polynomial when
  #       the complex Newton-Raphson method is applied, along with the detected
  #       attracting n-cycles (n>1) represented as polygonals, on the given region
  #       and up to some tolerance. Also, returns the number of colors used for
  #       representing the basins of attraction. 
  if inverted==true
    rect=invertedRectangle(xinterval,yinterval)
  else
    rect=rectangle(xinterval,yinterval)
  end
  a=size(rect)[1]
  b=size(rect)[2]
  tol=1.0/10^tolerance
  d=max(length(coefficientlistnum),length(coefficientlistden))-1
  fixedPoints=[]
  nCycles=[]
  if abs(evaluate_hpoly(coefficientlistnum,(complex(1.0),complex(0.0)),d))>tol && abs(evaluate_hpoly(coefficientlistden,(complex(1.0),complex(0.0)),d))<tol
    infFixedPoint=true # Checks whether the infinity point is a fixed point of the rational function given by Newton's method.
    push!(fixedPoints,Inf)
  else
    infFixedPoint=false
  end
  positer=Array{Vector{Int64}}(undef,a,b)
  @inbounds Threads.@threads for j in 1:b
    @inbounds for i in 1:a
      cycle=cycleDetector_Lyapunov(coefficientlistnum,coefficientlistden,rect[i,j],maxiterations,maxcycle,tolerance)
      if infFixedPoint==true && length(cycle[1])==1 && abs(cycle[1][1][2])<tol # In case the orbit of the point converges to infinity.
        positer[i,j]=[1.0,cycle[2]]
      elseif cycle[3]<tol && length(cycle[1])==1 && abs(cycle[1][1][2])>tol # In case the orbit of the point converges to another super-attracting fixed point.
        fPoint=cycle[1][1][1]/cycle[1][1][2]
        if real(fPoint)<tol
          fPoint=imag(fPoint)*im
        elseif imag(fPoint)<tol
          fPoint=real(fPoint)
        end
        pos=contains(fixedPoints,complex(fPoint),tolerance)
        if pos[1]==false
          push!(fixedPoints,complex(fPoint))
          positer[i,j]=[length(fixedPoints),cycle[2]]
        else
          positer[i,j]=[pos[2],cycle[2]]
        end
      else # In case the orbit of the point does not converge to any point (1-cycle).
        positer[i,j]=[0.0,cycle[2]]
        if length(cycle[1])>1 && length(cycle[1])<=maxcycle && cycle[3]<1.0 && IsInfIn(cycle[1])==false # If an attracting n-cycle (not including the infinity point) is detected.
          complexCycle=P1toComplex_cycle(cycle[1])
          if !degeneratedCycle(complexCycle,tolerance) && !contains_list(nCycles,complexCycle,6) # Checks that the detected cycle isn't degenerated and hasn't been already computed.
            push!(nCycles,complexCycle)
          end
        end
      end
    end
  end 
  # At this point, we have the positer matrix, a list with every fixed point (fixedPoints)
  # and a list with all detected attracting n-cycles with 2<=n<=maxcycle (nCycles).
  numFixedPoints=length(fixedPoints)
  numColors=numFixedPoints+1
  @inbounds colorMatrix=[positer[i,j][1]+0.25 for i=1:a,j=1:b]
  colorMatrix[1,1]=0.0
  colorMatrix[a,b]=numFixedPoints    
  if  infFixedPoint==false
    seed1=[Colors.RGB(0.0,0.0,0.0),Colors.RGB(0.5,0.5,0.5)]
  else
    seed1=[Colors.RGB(0.0,0.0,0.0),Colors.RGB(0.5,0.5,0.5),Colors.RGB(1.0,1.0,0.0)]
  end
  seed=union(seed1,[Colors.RGB(1.0,0.0,0.0),Colors.RGB(0.0,1.0,0.0),Colors.RGB(0.0,0.0,1.0)])
  gseed=Colors.distinguishable_colors(numColors,seed)
  fpcm=PyPlot.ColorMap(gseed)
  img=PyPlot.imshow(colorMatrix,cmap=fpcm,extent=[xinterval[1],xinterval[2],yinterval[1],yinterval[2]])
  numCycles=length(nCycles)
  for j in 1:numCycles
    cycle=nCycles[j]
    numElements=length(cycle)
    for i in 1:numElements
      if i==numElements
        x1=real(cycle[i])
        y1=imag(cycle[i])
        x2=real(cycle[1])
        y2=imag(cycle[1])
      else
        x1=real(cycle[i])
        y1=imag(cycle[i])
        x2=real(cycle[i+1])
        y2=imag(cycle[i+1])
      end
      x_values=[x1,x2]
      y_values=[y1,y2]
      PyPlot.plot(x_values,y_values,color="blue",linestyle="-")
    end
  end
  return img, numColors
end


#Example
#=
pol2cycle=complex([2.0,-2.0,0.0,1.0])
img=plotWithCycles_NewtonFractal_Lyapunov(normTradNum,normTradDen,(-1.5,1.5),(-1.5,1.5),25,3,8)[1]

gcf()

# clf(); can be useful to clear the current plot.
=#


# The philosophy of this method is exactly the same than positionIteration_superAttractingFixedPoints_givenFixedPoints;
# it works exactly in the same way than the previous one, is considerably more time-consuming, and its purpose is to
# preserve the color consistency between different graphics.
function plotWithCycles_BasinsOfAttraction_Lyapunov_givenFixedPoints(coefficientlistnum::Array{Complex{Float64},1},coefficientlistden::Array{Complex{Float64},1},
  xinterval::Tuple{Float64,Float64},yinterval::Tuple{Float64,Float64},fixedPoints::Array{Any,1},maxiterations::Int64=25,maxcycle::Int64=5,tolerance::Int64=8;inverted::Bool=false)
  # Pre:
  # Post: Plots the basin of attraction of the roots of the given polynomial when
  #       the complex Newton-Raphson method is applied, along with the detected
  #       attracting n-cycles (n>1) represented as polygonals, on the given region
  #       and up to some tolerance. Also, returns the number of colors used for
  #       representing the basins of attraction. 
  if inverted==true
    rect=invertedRectangle(xinterval,yinterval)
  else
    rect=rectangle(xinterval,yinterval)
  end
  a=size(rect)[1]
  b=size(rect)[2]
  tol=1.0/10^tolerance
  d=max(length(coefficientlistnum),length(coefficientlistden))-1
  nCycles=[]
  if abs(evaluate_hpoly(coefficientlistnum,(complex(1.0),complex(0.0)),d))>tol && abs(evaluate_hpoly(coefficientlistden,(complex(1.0),complex(0.0)),d))<tol
    infFixedPoint=true # Checks whether the infinity point is a fixed point of the rational function given by Newton's method.
  else
    infFixedPoint=false
  end
  positer=Array{Vector{Int64}}(undef,a,b)
  @inbounds Threads.@threads for j in 1:b
    @inbounds for i in 1:a
      cycle=cycleDetector_Lyapunov(coefficientlistnum,coefficientlistden,rect[i,j],maxiterations,maxcycle,tolerance)
      if infFixedPoint==true && length(cycle[1])==1 && abs(cycle[1][1][2])<tol # In case the orbit of the point converges to infinity.
        positer[i,j]=[1.0,cycle[2]]
      elseif cycle[3]<tol && length(cycle[1])==1 && abs(cycle[1][1][2])>tol # In case the orbit of the point converges to another super-attracting fixed point.
        fPoint=cycle[1][1][1]/cycle[1][1][2]
        if real(fPoint)<tol
          fPoint=imag(fPoint)*im
        elseif imag(fPoint)<tol
          fPoint=real(fPoint)
        end
        pos=contains(fixedPoints,complex(fPoint),tolerance)
        positer[i,j]=[pos[2],cycle[2]]
      else # In case the orbit of the point does not converge to any point (1-cycle).
        positer[i,j]=[0.0,cycle[2]]
        if length(cycle[1])>1 && length(cycle[1])<=maxcycle && cycle[3]<1.0 && IsInfIn(cycle[1])==false # If an attracting n-cycle (not including the infinity point) is detected.
          complexCycle=P1toComplex_cycle(cycle[1])
          if !degeneratedCycle(complexCycle,tolerance) && !contains_list(nCycles,complexCycle,6) # Checks that the detected cycle isn't degenerated and hasn't been already computed.
            push!(nCycles,complexCycle)
          end
        end
      end
    end
  end 
  # At this point, we have the positer matrix and a list with all detected attracting n-cycles with 2<=n<=maxcycle (nCycles).
  numFixedPoints=length(fixedPoints)
  numColors=numFixedPoints+1
  @inbounds colorMatrix=[positer[i,j][1]+0.25 for i=1:a,j=1:b]
  colorMatrix[1,1]=0.0
  colorMatrix[a,b]=numFixedPoints    
  if  infFixedPoint==false
    seed1=[Colors.RGB(0.0,0.0,0.0),Colors.RGB(0.5,0.5,0.5)]
  else
    seed1=[Colors.RGB(0.0,0.0,0.0),Colors.RGB(0.5,0.5,0.5),Colors.RGB(1.0,1.0,0.0)]
  end
  seed=union(seed1,[Colors.RGB(1.0,0.0,0.0),Colors.RGB(0.0,1.0,0.0),Colors.RGB(0.0,0.0,1.0)])
  gseed=Colors.distinguishable_colors(numColors,seed)
  fpcm=PyPlot.ColorMap(gseed)
  img=PyPlot.imshow(colorMatrix,cmap=fpcm,extent=[xinterval[1],xinterval[2],yinterval[1],yinterval[2]])
  numCycles=length(nCycles)
  for j in 1:numCycles
    cycle=nCycles[j]
    numElements=length(cycle)
    for i in 1:numElements
      if i==numElements
        x1=real(cycle[i])
        y1=imag(cycle[i])
        x2=real(cycle[1])
        y2=imag(cycle[1])
      else
        x1=real(cycle[i])
        y1=imag(cycle[i])
        x2=real(cycle[i+1])
        y2=imag(cycle[i+1])
      end
      x_values=[x1,x2]
      y_values=[y1,y2]
      PyPlot.plot(x_values,y_values,color="blue",linestyle="-")
    end
  end
  return img, numColors
end


# This method is specifically written to study the behaviour of cubic polynomials
# when the complex Newton-Raphson method is applied.
# In particular, we consider polynomials with two of its roots fixed (im and -im, in
# this case), letting its third root lambda take values in a complex parameter space.
# This way, we can study the behaviour of all general cubic polynomials (modulus a Moebius transformation),
# because we can always fix two of its roots to im and -im through a Moebius transformation.
# So, given a rectangle in the complex plane, and for every point in this rectangle,
# this method defines the polynomial whose third root is such point and checks
# whether if the rational function obtained by applying the N-R method to that polynomial
# has some sort of "cathastrophic behaviour". That is, if the rational function has any
# attracting n-cycles for which the N-R method never converges to a root, and "fails" in this sense.
# Then, for every point in the given rectangle, a color is assigned based on if the 
# associated polynomial presents this kind of cyclic behaviour. This process is followed
# in a reasonable execution time thanks to Fatou's Theorem. For more specific and
# detailed explanations of the implementation or the mathematical framework that supports it,
# we can only invite the reader to consult the source paper.
function detectCathastrophicBehaviour_cubicPolynomials_Newton(xinterval::Tuple{Float64,Float64},yinterval::Tuple{Float64,Float64},
  maxiterations::Int64=25,maxcycle::Int64=5,tolerance::Int64=8)
  # Pre:
  # Post: Plots the image, and returns the number of plotted colors.
  rect=rectangle(xinterval,yinterval)
  a=size(rect)[1]
  b=size(rect)[2]
  tol=1.0/10^tolerance
  posMatrix=Array{Float64}(undef,a,b)
  @inbounds Threads.@threads for j in 1:b
    @inbounds for i in 1:a
      lambda=rect[i,j][1]
      coefficientlistnum=complex([lambda,0.0,-1.0*lambda,2.0]) # Constructs the rational function given by the N-R method applied to the polynomial associated with this particular point. 
      coefficientlistden=complex([1.0,-2.0*lambda,3.0])
      cycle=cycleDetector_Lyapunov(coefficientlistnum,coefficientlistden,(lambda,complex(3.0)),maxiterations,maxcycle,tolerance)
      if cycle[3]<tol && length(cycle[1])==1 && abs(cycle[1][1][2])>tol # In case the orbit of the point converges to a super-attracting fixed point (a root of the polynomial).
        if abs(cycle[1][1][1]/cycle[1][1][2]-1.0*im)<tol # Checks if said point is im.
          posMatrix[i,j]=1.25
        elseif abs(cycle[1][1][1]/cycle[1][1][2]+1.0*im)<tol # Checks if said point is -im.
          posMatrix[i,j]=2.25
        else
          posMatrix[i,j]=4.25
        end
      else # In case the orbit of the point does not converge to any root of the polynomial (1-cycle).
        posMatrix[i,j]=0.25
      end
    end
  end
  posMatrix[1,1]=0
  posMatrix[a,b]=4
  seed=[Colors.RGB(1.0,0.0,0.0),Colors.RGB(0.0,1.0,0.0),Colors.RGB(0.0,0.0,1.0)]
  numColors=4
  gseed=Colors.distinguishable_colors(numColors,seed)
  fpcm=PyPlot.ColorMap(gseed)
  img=PyPlot.imshow(posMatrix,cmap=fpcm,extent=[xinterval[1],xinterval[2],yinterval[1],yinterval[2]])
  return img, numColors
end


#Examples
#detectCathastrophicBehaviour_cubicPolynomials_Newton((-2.0,2.0),(-2.0,2.0),40,5,8)[1]

#gcf()

# Note: In order to improve the resolution of certain regions of the plot
#       when zooming, it can be helpful to increase significantly the 
#       maximum number of iterations allowed up to 50 (for example).
#       It is also important to note that the more points in the 
#       selected region whose associated polynomial presents cathastrophic 
#       behaviour under the iteration of Newton-Raphson method, the higher 
#       is the computational cost of this method and the more time it'll 
#       take to execute it.


# Plots the same graphic as described in the previous method, but distinguishing
# the cyclic behaviour by the length of the cycle; that is, if the polynomial 
# associated with the point z in the parameter space presents 2-cyclic behaviour,
# this method chooses a color for z, distinct from any other colors, and does the
# same for every case in which n-cyclic behaviour appears, with 1<n<=maxcycle.
# As in the previous method, it is important to note that the computational cost
# of this method can escalate quickly if maxiterations or maxcycle are high, or 
# if most of the chosen region of the complex plane correspond to polynomials 
# that present catastrophic behaviour.
function detectCathastrophicBehaviour_distinguishCycles_cubicPolynomials_Newton(xinterval::Tuple{Float64,Float64},yinterval::Tuple{Float64,Float64},
  maxiterations::Int64=25,maxcycle::Int64=5,tolerance::Int64=8)
  # Pre:
  # Post: Plots the image, and returns the number of plotted colors.
  rect=rectangle(xinterval,yinterval)
  a=size(rect)[1]
  b=size(rect)[2]
  tol=1.0/10^tolerance
  posMatrix=Array{Float64}(undef,a,b)
  @inbounds Threads.@threads for j in 1:b
    @inbounds for i in 1:a
      lambda=rect[i,j][1]
      coefficientlistnum=complex([lambda,0.0,-1.0*lambda,2.0]) # Constructs the rational function given by the N-R method applied to the polynomial associated with this particular point. 
      coefficientlistden=complex([1.0,-2.0*lambda,3.0])
      cycle=cycleDetector_Lyapunov(coefficientlistnum,coefficientlistden,(lambda,complex(3.0)),maxiterations,maxcycle,tolerance)
      if cycle[3]<tol && length(cycle[1])==1 && abs(cycle[1][1][2])>tol # In case the orbit of the point converges to a super-attracting fixed point (a root of the polynomial).
        if abs(cycle[1][1][1]/cycle[1][1][2]-1.0*im)<tol # Checks if said point is im.
          posMatrix[i,j]=1.25
        elseif abs(cycle[1][1][1]/cycle[1][1][2]+1.0*im)<tol # Checks if said point is -im.
          posMatrix[i,j]=2.25
        else
          posMatrix[i,j]=3.25
        end
      else # In case the orbit of the point converges to an attracting n-cycle (n>1).
        convergesToCycle=false
        for n in 2:maxcycle
          if length(cycle[1])==n && cycle[3]<1.0 && IsInfIn(cycle[1])==false
            convergesToCycle=true
            posMatrix[i,j]=2.25+n*1.0
          end
        end
        if !convergesToCycle # In case the orbit of the point does not converge to any attracting n-cycle.
          posMatrix[i,j]=0.25
        end
      end
    end
  end
  posMatrix[1,1]=0
  posMatrix[a,b]=3+maxcycle
  seed=[Colors.RGB(1.0,0.0,0.0),Colors.RGB(0.0,1.0,0.0),Colors.RGB(0.0,0.0,1.0)]
  numColors=3+maxcycle
  gseed=Colors.distinguishable_colors(numColors,seed)
  fpcm=PyPlot.ColorMap(gseed)
  img=PyPlot.imshow(posMatrix,cmap=fpcm,extent=[xinterval[1],xinterval[2],yinterval[1],yinterval[2]])
  return img, numColors
end


#Examples
#detectCathastrophicBehaviour_cubicPolynomials_Newton((-2.0,2.0),(-2.0,2.0),40,5,8)

#gcf()

# Note: In order to improve the resolution of certain regions of the plot
#       when zooming, it can be helpful to increase significantly the 
#       maximum number of iterations allowed up to 50 (for example).
#       It is also important to note that the more points in the 
#       selected region whose associated polynomial presents cathastrophic 
#       behaviour under the iteration of Newton-Raphson method, the higher 
#       is the computational cost of this method and the more time it'll 
#       take to execute it.


# Computes the Lyapunov constants of the orbits of alpha/3 under the iteration
# of the rational map given by the application of Newton's method to the polynomial
# (z-alpha)(z^2+1) associated with each point alpha in the parameter space
# previously described, and plots it. A color is assigned to constants equal to
# 0.0, another to constants between 0.0 and 1.0 (both not included), another to
# constants equal to 1.0, and another to constants greater than 1.0. 
function plotLyapunovConstants_cubicPolynomials_Newton(xinterval::Tuple{Float64,Float64},yinterval::Tuple{Float64,Float64},
  maxiterations::Int64=25,maxcycle::Int64=5,tolerance::Int64=8)
  # Pre:
  # Post: Plots the image, and returns the number of plotted colors.
  rect=rectangle(xinterval,yinterval)
  a=size(rect)[1]
  b=size(rect)[2]
  tol=1.0/10^tolerance
  posMatrix=Array{Float64}(undef,a,b)
  @inbounds Threads.@threads for j in 1:b
    @inbounds for i in 1:a
      lambda=rect[i,j][1]
      coefficientlistnum=complex([lambda,0.0,-1.0*lambda,2.0]) # Constructs the rational function given by the N-R method applied to the polynomial associated with this particular point. 
      coefficientlistden=complex([1.0,-2.0*lambda,3.0])
      cycle=cycleDetector_Lyapunov(coefficientlistnum,coefficientlistden,(lambda,complex(3.0)),maxiterations,maxcycle,tolerance)
      if cycle[3]<tol # In case the orbit of the point converges to a super-attracting fixed point (a root of the polynomial) or a super-attracting n-cycle.
        posMatrix[i,j]=0.0
      elseif cycle[3]<1.0
        posMatrix[i,j]=1+cycle[3]
      elseif abs(cycle[3]-1.0)<tol
        posMatrix[i,j]=3.0
      else
        posMatrix[i,j]=4.0+cycle[3]/(cycle[3]+1)
      end
    end
  end
  posMatrix[1,1]=0
  posMatrix[a,b]=4
  seed=[Colors.RGB(1.0,0.0,0.0),Colors.RGB(0.0,1.0,0.0),Colors.RGB(0.0,0.0,1.0)]
  numColors=5
  gseed=Colors.distinguishable_colors(numColors,seed)
  fpcm=PyPlot.ColorMap(gseed)
  img=PyPlot.imshow(posMatrix,cmap=fpcm,extent=[xinterval[1],xinterval[2],yinterval[1],yinterval[2]])
  return img, numColors
end



end #module