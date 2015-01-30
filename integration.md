### Integration

#### The standard integration routine

The standard function for numerical integration in Julia Base is `quadgr()`. It implements an adaptive Gauss-Kronrod procedure of order (7,15) by default and is fast and extremely accurate, especially for smooth functions over finite intervals.

As an example, integrate the function $f(x) = e^{-x}\,\cos(x)$ from $0$ to $\pi$:

```julia
    f(x) = exp(-x) * cos(x);

    q = quadgk(f, 0, pi)
    (0.5216069591318859,4.937117381587086e-11)
```

`quadgk()` returns a vector, the first entry the computed value of the integral, the second the estimated absolute error. As we know the exact value of the integral is $\frac{1}{2} (1 + e^{-\pi})$, the true absolut error is:

    v = 0.5*(1+exp(-pi));

    abs(q[1] - v)
    2.220446049250313e-16

The following step function is a well-known test example on the interval `[0, 3]`:

    f24(x) = floor(exp(x));

    q24 = quadgk(f24, 0.0, 3.0)
    (17.690917276136695,2.4104866703467365e-7)

Because of its step function property the exact value can be easily calculated.

    q24_exact = sum(diff([log(1:20), 3]) .* [1:20])
    17.664383539246515

So the true absolute error is much larger than the estimated one. `quadgk()` can be forced to split the integration interval into smaller pieces, adding up the intermediate integral to return a more exact answer. The points at which `f24()` has discontinuities are `[log(1:20), 3];` and the call to `quadgk()` is now

    ds = [log(1:20), 3];

    quadgk(f24, ds...)
    (17.664383,7.771561172376096e-16)

which reduces the true absolute error to about 0.5e-06 -- still not optimal and certainly not as accurate as the estimated absolute error suggests.


#### Other integration routines

In package **Calculus** (and similar in **Calculus2**) there are implementations of an adaptive variant of _Simpson's rule_ and of an approach to _Monte Carlo integration_.

    using Calculus

    q = integrate(f, 0.0, pi)
    0.5216069591315619

with an absolute error smaller than 0.5e-12. Applying Monte Carlo integration will not give as good results,

   q = integrate(f, 0.0, pi, :monte_carlo)

results that will vary between 0.50 and 0.53 because of its stochastic nature. Monte-Carlo normally is applied only for higher dimensional integrals.

Compare this with a very simple approach of using the _trapezoidal rule_ on discrete points representing the function. An implementation is available in package **NumericalMath**.

    using NumericalMath

    xs = linspace(0, pi, 10000); ys = map(f, xs);

    q = trapz(xs, ys)
    0.5216069677136992                  # abs.err < 1e-08

We can improve this result by applying an "edge correction" term $-\frac{h^{12}}{12} (f'(b) - f'(a))$ with step length $h$ and discretized derivatives at end points $a$ and $b$.

    h = pi/(10000-1);

    q_corr = q - h * ((ys[end]-ys[end-1]) - (ys[2]-ys[1])) / 12;

    abs(q_corr - v)
    6.439293542825908e-15

With an absolute error of about 5e-15 a corrected trapezoid integration with 10000 nodes is more accurate than an adaptive Simpson integration.

**NumericalMath** also implements a version of *Romberg integration*. This method is based Richardson extrapolation by halving the step width and applying the trapezoidal rule.

    romberg(x -> exp(-x)*cos(x), 0, pi)
    (0.5216069591318865,1.0e-15)

Romberg integration works best with smooth and non-oscillatory functions. It needs the least number of function calls among all integration routines with comparable accuracies.


#### Singularities at interval boundaries

Gaussian integration rules never evaluate the integrand exactly at the interval end points. Therefore, `quadgk()` can handle singularities of integrand functions to a certain extent.

We will compute the integral of the _dilogarithm function_ $f(x) = log(1-x)/x$. This function has a singularity at 1.0, but not at 0.0 (where its value is -1.0).

    dilog(x) = log(1-x) / x;

    q, e = quadgk(dilog, 0, 1)
    (-1.644934063619709,1.8252236055171762e-8)

and the exact value of this integral is $\pi^2/6 = 1.6449340668482264\ldots$.

As another example, we will look at the following integral whose exact value we know:
$$
    \int_0^1 x^{a - 1} dx = \frac{1}{a}
$$
where $0 < a < 1$. These integrals are finite though there is a singulariy at $0$. We will compute it numerically for $a = 0.5$ and $a = 0.025$, the corresponding functions being $1/\sqrt(x)$ and $x^{-0.975}$.

    quadgk(x -> 1/sqrt(x), 0, 1)
    (1.9999999845983916,2.3762511591521646e-8)

    quadgk(x -> x^-0.975, 0, 1)
    (39.99998745356245,5.939527519982054e-7)    # abs.err = 1.25e-05

and the true absolute error is now larger than the estimated one. For $a \le 0.02$ we get a "DomainError".

    julia> quadgk(x -> x^-0.98, 0, 1)
    ERROR: DomainError
     in evalrule at quadgk.jl:75
     in do_quadgk at quadgk.jl:133
     in quadgk at quadgk.jl:174

Actually, no routine will reliably compute integrals for all kinds of singularities, the user always has to check for correctness or at least plausibility.


#### Integrals on infinite domains

As the help page of `quadgk()` says one or both of the endpoints of the real interval may be infinite (i.e., the interval may be infinite). In this case a "coordinate transformation is performed internally to map the infinite interval to a finite one". Let's test this with some examples, for example the well-known Gauss error integral
$$
    \frac{1}{\sqrt{2 \pi}} \int_{-\infty}^{\infty} e^{-\frac{1}{2} t^2} \mathrm{d}t
$$
whose value must be $1.0$.

    julia> fge1(t) =  exp(-t^2/2);

    julia> q = quadgk(fge1, -Inf, Inf);

    julia> q[1] / sqrt(2*pi)
    1.0000000000032585

But if we set the peak far outside, `quadgk()` has difficulties finding it there. The same will happen if the function is multimodal and the peaks are far away from each other.

    julia> mu = 1000; fge2(t) = exp(-(t-mu)^2/2);

    julia> quadgk(fge2, -Inf, Inf)
    (0.0,0.0)

If instead we integrate from 0 to 2000, with the paek at 1000 in the middle, then we get the correct result again.

    julia> quadgk(fge2, 0, 2000)
    (2.5066282746310016,3.417400287863345e-9)

Here are some more examples of infinite integrals whose exact values are known:

    julia> f(x) = sqrt(x) * exp(-x);

    julia> quadgk(f, 0, Inf)
    (0.8862269258634896,7.167163683171998e-9)   # sqrt(pi) / 2.0

    julia> f(x) = x * exp(-x^2)

    julia> quadgk(f, 0, Inf)
    (0.5,3.7291023073046645e-10)                # 0.5

    julia> f(x) = 1.0 / (1 + x^2)

    julia> quadgk(f, -Inf, Inf)
    (3.141592653589793,3.615143517876618e-9)    # pi

And here is a function that is difficult to integrate numerically on its infinite interval, mostly because it is obviously an oscillating function.

    julia> f(x) = sin(x) / x;

    julia> quadgk(f, 0.0, Inf)
    ERROR: DomainError
     in f at none:1
     in anonymous at quadgk.jl:107
     in evalrule at quadgk.jl:54
     in do_quadgk at quadgk.jl:134
     in do_quadgk at quadgk.jl:94

Nonetheless, the error message is disturbing as we would not expect a domain error. Even if the transformation maps `Inf` to 0 the integration routine should not evaluate the function at that point.

This integral is the sum of an alternating, decreasing sequence and is therefore convergent. An upper and lower bound can be obtained as follows (where we use `sinpi` instead of `sin` as it is sometimes said this function is more accurate for large input values):

    julia> N = 10^6; fpi(x) = sinpi(x) / x

    julia> quadgk(fpi, 0.0, N)[1], quadgk(fpi, 0.0, N+1)[1]
    (1.5707941572040225,1.570796613700503)

and the true value $\pi / 2$ fits perfectly into this interval.


#### Higher dimensional integrals


#### Complex line integrals

It is obvious how to integrate a complex-imaginary function on a real interval, i.e. a function $f: [a, b] \to C$. Split the function into real and imaginary parts, integrate each function as a real function on the interval, and combine the two results as real and imaginary part of the complex integral.

To integrate a complex function along a path in the complex plane is only slightly more difficult. Function `line_integral` in **NumericalMath** accepts a vector of points in the complex plane, computes the integrals on the straight lines between these points, and sums up a final result.

As an example, take the complex function $fz(z) = 1/z$ along two closed curves around the origin. The first is along the points $-1-i, 1-i, 1+i, -1+i, -1-i$, that is along a rectangle around $0+0i$.

    fz(z::Complex) = 1 ./ z;

    points = [-1.0-1.0im, 1.0-1.0im, 1.0+0im, -1.0+1.0im, -1.0-1.0im];

    line_integral(fz, points)
    -2.510715007782571e-15 + 6.283185307179583im

Because function $fz()$ (as a Laurent series) has a pole with residuum $1$ at the origin, every complex integral around $0+0i$ has the value $2 \pi i$. (Note that the path is going in the mathematically positive direction.)

    x = linspace(0.0, 2*pi);

    y = cos(x) + sin(x)*1.0im

    line_integral(fz, y)
    -2.1572909173345303e-17 + 6.283185307179586im

This is almost exactly the same result as above.
