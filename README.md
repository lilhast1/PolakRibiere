Implementation of Polak-Ribiere algorithm in julia.
The algorithm finds local minima using conjugate gradients[1] with Polak-Ribiere formula.

For the purpose of single variable minimization is used a variaton of Brents algorithm[2] which makes use of the first derivative.

Usage:
PolakRibiere(x->-exp(-(x[1]*x[1]+x[2]*x[2])), [2.,1.]) # the function and the initial guess
PolakRibiere(x->x[1] + x[2]^2/(4*x[1]) + x[3]^2/x[2] + 2/ x[3], [.4, 1.1, 1.1]) 
dbrent(0.0, 3.14, x->sin(x), x->cos(x)) # the interval [a,b] and then the function, and its derivative


[1] Ž. Jurić, Numerički algoritmi, 1. izd. Sarajevo: Elektrotehnički fakultet, 2018, p. XI, 481.
[2] William H. Press, Saul A. Teukolsky, William T. Vetterling, and Brian P. Flannery. 2007. Numerical Recipes 3rd Edition: The Art of Scientific Computing (3rd. ed.). Cambridge University Press, USA.
