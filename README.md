Download Link: https://assignmentchef.com/product/solved-math567-homework-2
<br>
. (Spectral differentiation). As discussed in class, the discrete Fourier series for a 2π periodic function u, sampled at the points tj = 2 Nπj, j = 0,1,…,N − 1, where N is odd, is given

,                                                    (1)

where the discrete Fourier coefficients are given as

.                             (2)

Here uj are the samples of u at tj. The formula for the coefficients arise from the trapezoidal rule approximation of the continuous Fourier coefficients. Recall that pN(t) interpolates uj at tj.

The forward and inverse discrete Fourier transforms (DFTs) are commonly defined as follows:

N−1

forward DFT: ˆuk = X uje−2πijk/N, k = 0,…,N − 1,                                   (3)

inverse DFT: .                      (4)

Letting uˆ denote the row-vector containing the N coefficients ˆuk, the discrete Fourier coefficients can be recovered from uˆ as follows:

uˆ  .

See the class notes for the definitions when N is even. The fast Fourier transform (FFT) and its inverse compute (3) and (4), respectively.

In this problem you will be approximating derivatives of the 2π periodic function

u(t) = ecos(5(t−0.1))

over [0,2π) using the FFT and finite differences. For all problems, you will be using the sample points  1, for an odd value of N.

(a)    Implement a code that uses the FFT and its inverse to approximate u0(t) at the sample points tj, j = 0,1,…,N − 1. Your code should take in the samples uj and compute the approximate derivatives u0j using the FFT library that is part of the language you are using:

•   In Matlab you can just call fft and ifft directly (these use FFTW)

•   In Juila you need to add the FFTW package

•   In Python you can use the FFT from NumPy, or better yet, use pyFFTW.

(b)   Using your code from part (a), approximate the derivative of u using N = 63 and N = 127 and compare this to correct derivative by plotting the difference between the correct derivative and the approximate derivative vs. the sample points tj, for each value of N. Include plots of the error for both values of N in your write-up.

(c)    Compute the approximate derivative of u at tj using the centered second order accurate finite difference formula

,

where h = 2π/N. In this formula you can use the fact that u is periodic to get values for u−1 and uN:

u−1 = uN−1         and      uN = u0.

You may wish to construct a sparse differentiation matrix DN for this computation, which takes the form, when N = 7,

.

The general matrix DN can be easily constructed using the toeplitz functions in Matlab, Julia, or NumPy. Similar to part (b), compute these approximations using N = 63 and N = 127 and make plots of the error. Include these plots in your write-up.

(d)   Comment on what you notice from the results from (b) and (c) in terms of which technique is more accurate and what happens to the errors when N is increased from 63 to 127.

2. (Discrete convolution). The discrete convolution of two vectors a   and x =   is defined as a new vector b with entries

N−1

bm = X am−jxj,               m = 0,…,N − 1,                                         (5)

j=0

or in matrix-vector form as

The structure exhibited in A gets a special name in linear algebra: circulant. Equation (5) can be viewed as a discretization of the continuous circular convolution (up to scaling) of two T-periodic functions a(t) and x(t) on the interval [t0,t0 + T]:

(7)

where am−j ≈ a(tm − tj), xj ≈ x(tj), and t` = t0 + (T/N)`. Circular convolutions and their discrete versions play a fundamental role in signal processing (your cell phones would not work without them!).

One of the most important results of the discrete Fourier transform (DFT) is that it can be used to relate the discrete convolution of two vectors to the DFT of the individual vectors in a beautiful way. The result is known as the discrete (or circular) convolution theorem and says the following. Let aˆ and xˆ denote the DFTs of a and x in (6), respectively. Then the DFT of the discrete convolution b of a and x is given by

or in matrix-vector form:ˆbm = aˆmxˆm,m = 0,…,N − 1,(8)(9)

Aˆ := diag(aˆ)                     xˆ                 bˆ

This theorem means that we can recover the result of the discrete convolution of a and x by the following procedure:

i Compute the DFTs of a and x. ii Obtain bˆ according to (8).

iii Compute the inverse DFT of bˆ to obtain b in (5).

By using the FFT to compute the DFT, this procedure requires O(N logN) operations instead of O(N2) as formula (6) might suggest.

Problem: Verify numerically that the discrete convolution theorem is correct by computing b directly using (6) and computing it using the procedure (i)–(iii) above with FFT. Use the following

to construct the vectors a and x. Choose N = 64 for the length of the vectors and σ = 0.1. Report the max-norm of the difference between b computed from the two methods.

If one scales b by 2π/N, then this gives the approximation to (7). Scale b in this fashion and plot it together with x and use t as the independent variable. Note that the scaled b is a smoothed out version of x. This is one of the main applications of the discrete convolution theorem in signal processing.

3. Fast solutions of tridiagonal linear systems. Consider the tridiagonal linear system of equations

,                    (10)

u                    fand suppose u0 = uN = f0 = fN = 0, and b 6= 0. This system can be solved efficiently using the Crout (or Thomas) algorithm for linear systems. However, in this problem you will show how to solve the system using the Discrete Sine Transform (DST).

(a)    Let uˆ and ˆf denote the DST of the vectors u and f in (10), i.e.

,

for k = 1,…,N − 1. Now, we can express the entries of u and f as

,

Substitute these expressions into the linear system (10) to show that row j of this systems simplifies to

,

in the case of j = 1 and j = N − 1 use the fact that u0 = uN = 0.

(b)   Use the results from part (a) to show that the DST of the solution to (10) is given by

.

(c)    Put parts (a) and (b) together to explain how to obtain the solution to (10).

(d)   Download the Matlab codes dst.m and idst.m on the course Blackboard page for computing the DST and inverse DST (these codes use the FFT to do the computation and are thus ‘fast’—although not as fast as they could be since they use complex arithemtic). Use these codes, or translate them into another language, and the procedure outlined in (a)–(c) to solve the linear system (10) for N = 100, a = 1, b = −2, and fj = h2 tanh(4sin(πj/N)), for j = 1,…,N − 1, where h = π/N. Make a plot of the solution u.

Note: This gives a second-order finite difference solution to the boundary value problem u00(x) = tanh(4sin(x)), 0 ≤ x ≤ π, u(0) = u(π) = 0, at the points xj = πj/N, j = 1,…,N − 1.

4. Implicit FD methods. Finite difference (FD) schemes can be made more accurate by adding more nodes to the formulas (i.e. increasing the stencil widths). For example the schemes

and

are the respective second and fourth-order centered FD approximations to the 1-D Poisson equation u00(x) = f. The problem with increasing the stencil width is that special treatment must be given to points near the boundary. This has the potential for altering the nice banded structure of the resulting linear systems, which can lead to a considerable increase in the computational cost involved in solving these systems. An additional problem, especially for the Poisson equation, is that the discrete maximum principle of the schemes can be lost (as illustrated in the second FD scheme above).

An alternative idea for increasing the order of accuracy of the FD schemes without increasing the stencil width is based on an idea first developed in the 1950s by L. Collatz and originally termed Mehrstellenverfahren (however, now the term “implicit”, “compact”, or “Hermite” is typically used, as in Homework 1). The basic idea is to use information from the underlying differential equation to increase the order. The following is an example of such a scheme for the 3-point, centered, approximation to the 1-D Poisson equation u00(x) = f(x)

(11)

which matches the formula (5) from problem 3 of Homework 1. In that homework you also showed that this formula is fourth-order accurate.

Fornberg [1, p. 67] describes the following simple technique for deriving the approximation (11). First note that

.                             (12)

By differentiating both sides of the differential equation with respect to x twice, we have u0000(x) = f00(x). Discretizing this equation with second-order accurate FD formulas for u0000(x) and f00(x) gives

.                  (13)

Combining equations (12) and (13) with the linear combination

results in precisely the fourth-order accurate scheme (11).

(a)    Generalize the technique above to derive an implicit, three-point, fourth-order accurate FD formula for the following BVP:

u00(x) + k2u(x) = f(x), a ≤ x ≤ b

which is referred to as the Helmholtz equation and k is interpreted as the ‘wave-number’.

(b)   Use the implicit formula from part (a) to solve Helmholtz equation for 0 ≤ x ≤ 1, with f(x) = 1, k = 150, and boundary conditions u(0) = 1, u(1) = 0. Compare your solution to the exact answer,

,

by plotting relative two-norm of the error in the approximate solution as a function of the grid spacing, h. Use h = 2−`, ` = 7,8,…,16 for this plot. Verify that the convergence is O(h4). Also plot the approximate solution for the finest grid. For full credit, use the tridiagonal solver that you developed in problem 3 for this problem.