# Multivariate Closed Skew Normal (MVCSN) Distribution 
The MVCSN repository provides MATLAB code for calculating the probability density of a *d*-dimensional, Closed Skew Normal Distribution at realization vector $\mathbf{x}$ with centralization vector $\mathbf{\mu}$, vector $\mathbf{v}$, scale matrix $\Sigma$, matrix $\Delta$, and matrix $\Gamma$. The probability density function defined in [1] is:  <br>

$$
\begin{gather}
    p(\mathbf{x})= \frac{\Phi(\Gamma(\mathbf{x}-\mathbf{\mu})|\mathbf{v}, \Delta)}{\Phi(\mathbf{0}|\mathbf{v}, \Delta + \Gamma\Sigma\Gamma^T)}\mathcal{N}(\mathbf{x}|\mathbf{\mu}, \Sigma),
    \\ 
    \text{where} \quad \Phi(\mathbf{a}|\mathbf{b}, C) = \int_{-\infty}^{\mathbf{a}}\mathcal{N}(\mathbf{x}|\mathbf{b}, C)d\mathbf{x} \quad \text{and} \quad \mathcal{N}(\mathbf{x}|\mathbf{\mu}, \Sigma) = \frac{1}{\sqrt{(2\pi)^k|\Sigma|}}\exp\Bigg(-\frac{1}{2}(\mathbf{x}-\mathbf{\mu})^T\Sigma^{-1}(\mathbf{x}-\mathbf{\mu})\Bigg)
\end{gather}
$$

When $\Gamma = 0$, *mvcsn.m* becomes *mvnpdf.m*. Please direct any questions to blhanson@ucsd.edu. <br><br>

## References
[1] Genton, M.G. (2004). Skew-elliptical distributions and their applications: A journey beyond normality. In G. Gonzalez-Farias J.A. Domínguez-Molina A.K. Gupta, (Eds.), The closed skew-normal distribution. Chapman and Hall/ CRC. 
