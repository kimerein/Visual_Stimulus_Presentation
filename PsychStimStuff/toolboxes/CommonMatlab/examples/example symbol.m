syms A B f
rho =int(1/(A + B*f^2),f)

syms m x s Pi
solve(-(m - x)^3/(3*s^2) - (m - x)*(-((m - x)^2/s^2) + log(exp((m - x)^2/s^2)/((2*Pi)^(1/2)* s)))-.98,x)
