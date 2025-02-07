# Multivariate Closed Skew Normal (MVCSN) Distribution 
The MVCSN repository provides MATLAB code for calculating the probability density of a *d*-dimensional, Closed Skew Normal Distribution at realization vector $\boldsymbol{x}$ with centralization vector $\boldsymbol{\mu}$, vector $\boldsymbol{v}$, scale matrix $\Sigma$, matrix $\Delta$, and matrix $\Gamma$. The probability density function defined in [1] is:  <br>

$$
\begin{gather}
    p(\bm{x})&= \frac{\Phi(\Gamma(\bm{x}-\boldsymbol{\mu})\,|\,\bm{v}, \Delta)}{\Phi(\boldsymbol{0}\,|\,\bm{v}, \Delta + \Gamma\Sigma\Gamma^T)}\,\mathcal{N}(\bm{x}\,|\boldsymbol{\mu}, \Sigma),
    \\ 
    \text{where}\Phi(\bm{a}\,|\,\bm{b}, C) = \int_{-\infty}^{\bm{a}}\mathcal{N}\,(\bm{x}\,|\,\bm{b}, C)d\bm{x} \text{  and  } \mathcal{N}(\bm{x}\,|\boldsymbol{\mu}, \Sigma) = \frac{1}{\sqrt{(2\pi)^k|\Sigma|}}\exp\Bigg(-\frac{1}{2}(\bm{x}-\boldsymbol{\mu})^T\Sigma^{-1}(\bm{x}-\boldsymbol{\mu})\Bigg)
\end{gather}
$$

When $\Gamma = 0$, *mvcsn.m* becomes *mvnpdf.m*. Please direct any questions to blhanson@ucsd.edu. <br><br>

## References
[1] Genton, M.G. (2004). Skew-elliptical distributions and their applications: A journey beyond normality. In G. Gonzalez-Farias J.A. Domínguez-Molina A.K. Gupta, (Eds.), The closed skew-normal distribution. Chapman and Hall/ CRC. 
