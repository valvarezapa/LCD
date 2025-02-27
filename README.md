﻿# Lyapunov Cycle Detector

Welcome to the *LCD* module! We present a collection of algorithms written in Julia Language which you might find very useful when approaching to certain Numerical Analysis problems involving rational maps. This kind of problems arise very often in many different areas of Mathematics and Computing like Complex Dynamics, Fractal Geometry or Algorithmic Methods. For example, when you apply Newton's method (the most currently used method to approximate the solutions of a non-linear equation) to a polynomial, the resulting function you have to iterate is a rational map.
The basic functionality of this software is to detect attracting *n*-cycles when a complex rational map is iterated, and to plot its basins of attraction on the plane of complex numbers. In order to do this, a new method that relies in great mathematical theorems (such as Birkhoff's Ergodic Theorem and many other results given by Sullivan, recently awarded with the Abel prize) is implemented and described extensively in the article *Algorithms for computing basins of attraction associated with a rational self-map of the Hopf fibration based on Lyapunov exponents* (https://doi.org/10.1016/j.topol.2023.108578).
This collection of algorithms avoid many computational problems that often appears in Numerical Analysis, such as overflows caused by denominators close to zero or indeterminations that appear when the numerator and denominator are both equal to zero.
This new method (which we will often refer to as *Lyapunov's method*) is focused on detecting attracting *n*-cycles. This approach is more general than some traditional algorithms that compute the basin of attraction of the fixed points of a given rational map. Also, another advantage of the algorithms we present is that they do not need any previous knowledge about the fixed points (nor to compute them beforehand) of the rational map in order to compute the basins of attraction.
Some of these improvements come with a prize; this collection of algorithms might not be the fastest option to compute only the basins of attraction of fixed points of a rational map, since we want to detect not only the attracting *1*-cycles (fixed points) but also any attracting *n*-cycle. However, we benefit from some built-in Julia macros for Parallel Programming and Multithreading in order to make this algorithms work in a reasonable amount of time.

If you are unfamiliar with the mathematical concepts and technicallities behind these algorithms, you can check the article *Algorithms for computing basins of attraction associated with a rational self-map of the Hopf fibration based on Lyapunov exponents* (https://doi.org/10.1016/j.topol.2023.108578) and its references. If you aren't but still want to plot some cool-looking fractals or to study the basins of attractions of a rational map of your choice, you can follow the *[Quick User Guide](https://github.com/valvarezapa/LCD/blob/main/LCD%20Quick%20User%20Guide.ipynb)* in order to look for specific details of the usage of this package.

Let's see an illustrative example of what kind of results this software can get, and the mathematical context behind them.

![firstExample](https://github.com/valvarezapa/LCD/blob/main/Examples/Newton's%20method/z%5E3-1/z%5E3-1%20(origin).PNG "Basins of z^3-1")

A reader familiar with the theory of iteration of rational maps will soon notice that the previous image corresponds to the basins of attraction of the attracting *n*-cycles (in this case, just *1*-cycles) of the rational map given by the iteration of complex Newton-Raphson's method over the polynomial *z^3-1*. Attending to the colors on the colorbar, note that colors 2, 3 and 4 correspond to the basin of attraction of a root of z^3-1 each, and colors 0 and 1 correspond respectively to the basin of attraction of *n*-cycles (n>1) and the basin of attraction of infinity (repulsive fixed point).
This gray and black colors do not appear on the graphic. Of course, much more mathematically sophisticated examples can be studied using the algorithms presented here, as is exposed in the already mentioned reference. 

# How to install

How do I import this collection of algorithms in Julia? It's simple, just execute the following code in the Pkg> command line to install the *LCD* package.

~~~
add https://github.com/valvarezapa/LCD
~~~

You need to install the module only once. It is *REALLY IMPORTANT* to note that the presented module has 3 dependencies: *Polynomials.jl*, *Colors.jl* and *PyPlot*. Also, this package was developed with *Julia 1.6.4*, so with later versions some methods could not work. In order to use the algorithms that the module contains, just add the following code to load it

~~~
using LCD.jl
~~~

whenever you want to. Since the algorithms are collected in a module, to invoke a certain *LCD* method, it must always be preceded by the *LCD.* prefix, as we shall see next (and as it appears in the *Quick User Guide*).

Also, you will need to have the Python [Matplotlib](http://matplotlib.org/) library installed in order to use PyPlot, which is used to plot the graphics. You can either do inline plotting with [IJulia](https://github.com/JuliaLang/IJulia.jl), which doesn't require a GUI backend, or use the Qt, wx, or GTK+ backends of Matplotlib (as is shown in the *Quick User Guide*).

# How to use it

We strongly recommend to read the *LCD Quick User Guide* in order to get started with the functionalities of the *LCD* module. However, here we present a brief summary:

The functionality of this module is divided essenially in 3 different sections. The code of the first one (which is the central one) is used in the other 2 sections, but in terms of functionality all of them serve different purposes and can be used independently.

The first section is kind of a fixed-points-oriented one, and it consists of the basic methods that are used to detect said *n*-cyclic behaviour when a rational map is iterated. The main method of this section is the following:

~~~
LCD.plotBasinsOfAttraction_Lyapunov(coefficientlistnum,coefficientlistden,(-1.5,1.5),(-1.5,1.5),200,1,8;coloringStrat="position")
~~~

As the name of the method suggests, it is used to compute and plot the basins of attraction of the attracting *n*-cycles of a rational map (given the coefficients of its numerator and denominator) on the given rectangle of the complex plane. Of course, a maximum number of iterations of Lyapunov's method, the maximum period of the detected cycles and a tolerance have to be selected.
Every basin of attraction of each fixed point has a different color in the graphic, while the points whose orbit converges to an *n*-cycle (*n>1*) or diverges appear in another different color. There are a few different coloring strategies implemented in the code one might follow in order to generate this plots, also taking into account the number of iterations of Lyapunov's method it took for each point to converge, for example.
Let's see an example of this method's functionality. In the following graphic we can see the basins of attraction of the complex polynomial *z^5-1*, both in a neighborhood of the origin and in a neighborhood of infinity. The black areas are the points whose orbit does not converge to any fixed point.

![z^5-1_(origin)](https://github.com/valvarezapa/LCD/blob/main/Examples/Newton's%20method/z%5E5-1/z%5E5-1_(origin).png "Basins of z^5-1")

![z^5-1_(infinity)](https://github.com/valvarezapa/LCD/blob/main/Examples/Newton's%20method/z%5E5-1/z%5E5-1_(infinity).png "Basins of z^5-1 in a neighborhood of infinity")

Again, in this file only a brief review of the methods is given; for a more specific and in-depth explanation of the functionality of each method one can consult the article *Algorithms for computing basins of attraction associated with a rational self-map of the Hopf fibration based on Lyapunov exponents* (https://doi.org/10.1016/j.topol.2023.108578).

The second one of the sections in which this module is divided is devoted to plotting the attracting *n*-cycles as such. The main method of this section is the following:

~~~
LCD.plotWithCycles_BasinsOfAttraction_Lyapunov(coefficientlistnum,coefficientlistden,(-1.5,1.5),(-1.5,1.5),200,4,8)
~~~

It is used to plot the basins of attraction of each attracting *n*-cycle, with a different color for each fixed point (attracting *1*-cycle) and a different color for each attracting *n*-cycle (*n>1*).  This section also has method to detect and plot the *n*-cycles as such, representing each *n*-cycle with a polygonal on the complex plane whose vertices are the elements of the cycle. The plotting colors are selected in a way that the detected *n*-cycle is clearly visible, as we shall see in the following figure.

![origin](https://github.com/valvarezapa/LCD/blob/main/Examples/Newton's%20method/(z%5E2%2B1)(z-(0.845%2B0.81i))%20%5B2%20-%20CYCLE%5D/2cycle_colorbar.PNG "Basins in a neighborhood of the origin")

![infinity](https://github.com/valvarezapa/LCD/blob/main/Examples/Newton's%20method/(z%5E2%2B1)(z-(0.845%2B0.81i))%20%5B2%20-%20CYCLE%5D/neighborhood_infinity_colorbar.PNG "Basins in a neighborhood of infinity")

In this graphic we can see the basins of attraction of the complex polynomial *(z^2+1)(z-(0.845+0.81i))*, both in a neighborhood of the origin and in a neighborhood of infinity. The black areas are the points whose orbit does not converge to any fixed point. In fact, the black areas we see in these two graphics are the basins of an attracting *2*-cycle. As one can see in the graphics, the mentioned *2*-cycle is clearly visible in blue.

In the last section there are a few specific methods implemented to study the *n*-cyclic behaviour that appears in cubic polynomials when some iterative method is applied. This generates pretty interesting images, in which one can clearly see the fractal structure of the regions associated with those polynomials that behave "badly" (that is, that presents some attracting *n*-cycles when a given iterative method is applied). For a more technical explanation and mathematical justifications (which involve results given by Fatou and the Scaling Theorem) for the complex parameter space that is considered in the algorithms of this section in order to study the behaviour of cubic polynomials, we again refer the reader to the article *Algorithms for computing basins of attraction associated with a rational self-map of the Hopf fibration based on Lyapunov exponents* (https://doi.org/10.1016/j.topol.2023.108578). Also, similar studies of complex parameter spaces associated with cubic polynomials can be found in the literature for many different iterative methods. A few of them can be found in the references of the previously mentioned article.
Several iterative methods are supported in this software, such as Newton-Raphson's method, Halley's, Schroeder's, Chebyshev's,... so one can compare the behaviour of cubic polynomials under iteration of these methods.
An example of this is the following:

![cubicPolynomials_Newton](https://github.com/valvarezapa/LCD/blob/main/Examples/Newton's%20method/Parameter%20plane%20of%20cubic%20polynomials%20%5BDynamical%20Study%5D/parameter_plane_cubicPolynomials.PNG)

This method is specific for Newton-Raphson's method. As in the other methods, a rectangle, and the maximum number of iterations, maximum period of the detected *n*-cycles and a tolerance must be given. An example of how to obtain this kind of graphics is the following:

~~~
LCD.detectCathastrophicBehaviour_cubicPolynomials_Newton((-2.0,2.0),(-2.0,2.0),200,5,8)
~~~

As we can see in the following images, the areas whose associated polynomials behave "badly" (that is, that present *n*-cyclic behaviour when Newton's method is applied, or that the orbit of some point in the complex plane diverges), despite being relatively small and scattered, are in fact Mandelbrot-like sets.
Also, each one of this methods has an alternative version in which the areas corresponding to "bad" polynomials are colored by which *n*-cyclic behaviour they present; a color for attracting *2*-cycles, another for *3*-cycles, and so on.

![cubicPolynomials_cycles_Newton](https://github.com/valvarezapa/LCD/blob/main/Examples/Newton's%20method/Parameter%20plane%20of%20cubic%20polynomials%20%5BDynamical%20Study%5D/topRight_1bulb_2cycles.png)

The results are intriguing, as one can see in the figure above.

# Additional comments

It is important to note that the methods contained in this repository address the problem of computing the basins of attraction induced by a discrete dynamical system induced by a rational map. Other Julia repositories and packages of similar nature (like (https://juliadynamics.github.io/DynamicalSystems.jl/previews/PR156/chaos/basins/) also address the same problem, often in a more direct and more computationally efficient way. However, the methods presented in this repository provide an alternative way of computing basins of attraction based on topological/geometrical techniques that allow our approach to avoid numerical problems like overflows, underflows and mathematical indeterminations, and to obtain some other information regarding the dynamical system, as for example computing the induced Lyapunov exponents (see (https://doi.org/10.1007/s00208-002-0404-7)) or computing the basins in a neighborhood of the infinity point, all of this without any previous knowledge of the system (neither attractors, fixed points, etc,...). For a more detailed explanation on the differences between this method and the more classical approach often chosen by other repositories and packages devoted to computing basins of attraction, we refer the reader to both the *[Quick User Guide](https://github.com/valvarezapa/LCD/blob/main/LCD%20Quick%20User%20Guide.ipynb)* and the article *Algorithms for computing basins of attraction associated with a rational self-map of the Hopf fibration based on Lyapunov exponents* (https://doi.org/10.1016/j.topol.2023.108578).

# Cite this work

This collection of algorithms have been developed by Víctor Álvarez Aparicio, Luis Javier Hernández Paricio, María Teresa Rivas Rodríguez and José Manuel García Calcines. If you would like to cite our results in your work, we urge you to cite our paper *Algorithms for computing basins of attraction associated with a rational self-map of the Hopf fibration based on Lyapunov exponents* (https://doi.org/10.1016/j.topol.2023.108578). You can also download the citation to this repository, which can be found in the *About* section on the right, in the homepage.
