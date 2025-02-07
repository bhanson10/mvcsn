function y = mvcsn(X, Mu, V, Sigma, Delta, Gamma)
%MVCSN Multivariate Closed Skew Normal (CSN) Distribution 
%   Y = MVCSN(X) returns the probability density of the multivariate
%   closed skew normal distribution with zero mean, zero v, identity 
%   covariance matrix, zero Delta, and zero Gamma, evaluated at each row of 
%   X.  Rows of the N-by-D matrix X correspond to observations or points, 
%   and columns correspond to variables or coordinates.  Y is an N-by-1 vector.
%
%   Y = MVCSN(X,MU) returns the density of the multivariate CSN with 
%   centralization parameter MU, zero v, identity scale matrix, zero Delta,
%   and zero Gamma, evaluated at each row of X.  MU is a 1-by-D vector, a 
%   scalar value, which MVCSN replicates to match the size of X.
%
%   Y = MVCSN(X,MU,V) returns the density of the multivariate CSN with 
%   centralization parameter MU, V, identiy scale matrix, zero Delta, and
%   zero Gamma, evaluated at each row of X. V is a 1-by-D vector, or a 
%   scalar value, which MVCSN replicates to match the size of X.
%
%   Y = MVCSN(X,MU,V,SIGMA) returns the density of the multivariate CSN with 
%   centralization parameter MU, V, scale matrix SIGMA, zero Delta, and 
%   zero Gamma, evaluated at each row of X.  SIGMA is a D-by-D matrix.
%
%   Y = MVCSN(X,MU,V,SIGMA,DELTA) returns the density of the multivariate CSN 
%   with centralization parameter MU, V, scale matrix SIGMA, DELTA,
%   and zero Gamma, evaluated at each row of X.  DELTA is a D-by-D matrix.
%
%   Y = MVCSN(X,MU,V,SIGMA,DELTA,GAMMA) returns the density of the 
%   multivariate CSN with centralization parameter MU, V, scale matrix 
%   SIGMA, DELTA, and GAMMA, evaluated at each row of X. GAMMA is a D-by-D 
%   matrix.
%
%   If GAMMA=0, MVCSN becomes MVNPDF.
%
%   Example:
%
%     Mu = [0 0]; V = [0 0]; Sigma = [2 0; 0 2]; Delta = [3 0; 0 3]; Gamma = [-5 0; 0 -5]; 
%     [X1,X2] = meshgrid(linspace(-5,1,100)', linspace(-5,1,100)');
%     X = [X1(:) X2(:)];
%     y = mvcsn(X, Mu, V, Sigma, Delta, Gamma);
%     surf(X1,X2,reshape(y,100,100));
%   
%   References:
%      Genton, M.G. (2004). Skew-elliptical distributions and their 
%      applications: A journey beyond normality. In G. Gonzalez-Farias J.A. 
%      Domínguez-Molina A.K. Gupta, (Eds.), The closed skew-normal 
%      distribution. Chapman and Hall/ CRC. 
%
% Copyright 2025 by Benjamin L. Hanson, published under BSD 2-Clause License.

if nargin<1
    error("TooFewInputs");
elseif ~ismatrix(X)
    error("InvalidData");
end

% Get size of data.  Column vectors provisionally interpreted as multiple scalar data.
[~,d] = size(X);
if d<1
    error("TooFewDimensions");
end

% Assume zero mean, data are already centered
if nargin < 2 || isempty(Mu)
    X0 = X;  
    V = zeros(1, d); 
    Sigma = eye(d); 
    Delta = zeros(d); 
    Gamma = zeros(d); 
% Get scalar mean, and use it to center data
elseif nargin < 3
    X0 = X - Mu;
    V = zeros(1, d); 
    Sigma = eye(d); 
    Delta = zeros(d); 
    Gamma = zeros(d); 
elseif nargin < 4
    X0 = X - Mu;
    [~, d_V] = size(V);
    if d~=d_V
        error("BadVDimension");
    end
    Sigma = eye(d); 
    Delta = zeros(d); 
    Gamma = zeros(d); 
elseif nargin < 5
    X0 = X - Mu;
    [~, d_V] = size(V);
    if d~=d_V
        error("BadVDimension");
    end
    [d1_S, d2_S] = size(Sigma);
    if d1_S~=d2_S
        error("NotSquareSigma");
    elseif d1_S~=d
        error("BadSigmaDimension");
    end
    Delta = zeros(d); 
    Gamma = zeros(d); 
elseif nargin < 6
    X0 = X - Mu;
    [~, d_V] = size(V);
    if d~=d_V
        error("BadVDimension");
    end
    [d1_S, d2_S] = size(Sigma);
    if d1_S~=d2_S
        error("NotSquareSigma");
    elseif d1_S~=d
        error("BadSigmaDimension");
    end
    [d1_D, d2_D] = size(Delta);
    if d1_D~=d2_D
        error("NotSquareDelta");
    elseif d1_D~=d
        error("BadDeltaDimension");
    end
    Gamma = zeros(d); 
elseif nargin < 7
    X0 = X - Mu;
    [~, d_V] = size(V);
    if d~=d_V
        error("BadVDimension");
    end
    [d1_S, d2_S] = size(Sigma);
    if d1_S~=d2_S
        error("NotSquareSigma");
    elseif d1_S~=d
        error("BadSigmaDimension");
    end
    [d1_D, d2_D] = size(Delta);
    if d1_D~=d2_D
        error("NotSquareDelta");
    elseif d1_D~=d
        error("BadDeltaDimension");
    end
    [d1_G, d2_G] = size(Gamma);
    if d1_G~=d2_G
        error("NotSquareGamma");
    elseif d1_G~=d
        error("BadGammaDimension");
    end
else
    error("BadInputs");
end

y = (mvncdf(zeros(size(X)), V, Delta + Gamma * Sigma * Gamma')).^-1 ...
    .* mvncdf(pagemtimes(Gamma,'none', X0, 'transpose')', V, Delta)...
    .* mvnpdf(X0, zeros(1,d), Sigma); 
end