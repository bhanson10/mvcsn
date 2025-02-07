# Multivariate Closed Skew Normal (MVCSN) Distribution 
The MVCSN repository provides MATLAB code for calculating the probability density of a *d*-dimensional, Closed Skew Normal Distribution at realization vector $\boldsymbol{x}$ with centralization vector $\boldsymbol{\mu}$, vector $\boldsymbol{v}$, scale matrix $\Sigma$, matrix $\Delta$, and matrix $\Gamma$. The probability density function defined in [1] is:  <br>

$$
\begin{gather}
    p(\boldsymbol{x})&= \frac{\Phi(\Gamma(\boldsymbol{x}-\boldsymbol{\mu})\,|\,\boldsymbol{v}, \Delta)}{\Phi(\boldsymbol{0}\,|\,\boldsymbol{v}, \Delta + \Gamma\Sigma\Gamma^T)}\,\mathcal{N}(\boldsymbol{x}\,|\boldsymbol{\mu}, \Sigma),
    \\ 
    \text{where}\Phi(\boldsymbol{a}\,|\,\boldsymbol{b}, C) = \int_{-\infty}^{\boldsymbol{a}}\mathcal{N}\,(\boldsymbol{x}\,|\,\boldsymbol{b}, C)d\boldsymbol{x} \text{  and  } \mathcal{N}(\boldsymbol{x}\,|\boldsymbol{\mu}, \Sigma) = \frac{1}{\sqrt{(2\pi)^k|\Sigma|}}\exp\Bigg(-\frac{1}{2}(\boldsymbol{x}-\boldsymbol{\mu})^T\Sigma^{-1}(\boldsymbol{x}-\boldsymbol{\mu})\Bigg)
\end{gather}
$$

When $\Gamma = 0$, *mvcsn.m* becomes *mvnpdf.m*. Please direct any questions to blhanson@ucsd.edu. <br><br>

## References
[1] Genton, M.G. (2004). Skew-elliptical distributions and their applications: A journey beyond normality. In G. Gonzalez-Farias J.A. Domínguez-Molina A.K. Gupta, (Eds.), The closed skew-normal distribution. Chapman and Hall/ CRC. 
