function y = mvcsn(X, Mu, Sigma, Gamma, Nu, Delta)
%MVCSN Multivariate Closed Skew Normal (CSN) Distribution 
%   Y = MVCSN(X) returns the probability density of the multivariate
%   closed skew normal distribution with zero mean, identity covariance 
%   matrix, zero Gamma, zero v, and identity Delta, evaluated at each row of 
%   X. Rows of the N-by-D matrix X correspond to observations or points, 
%   and columns correspond to variables or coordinates. Y is an N-by-1 vector.
%
%   Y = MVCSN(X,MU) returns the density of the multivariate CSN with 
%   centralization parameter MU, identity scale matrix, zero Gamma, 
%   zero Nu, and identity Delta, evaluated at each row of X. MU is a D-by-1 
%   vector or a scalar value, which MVCSN replicates to match the size of X.
%
%   Y = MVCSN(X,MU,SIGMA) returns the density of the multivariate CSN with 
%   centralization parameter MU, scale matrix SIGMA, zero Gamma, zero Nu, and
%   identity Delta, evaluated at each row of X. SIGMA is a D-by-D matrix or
%   a scalar value, which MVCSN multiplies by the identity matrix.
%
%   Y = MVCSN(X,MU,SIGMA,GAMMA) returns the density of the multivariate CSN 
%   with centralization parameter MU, scale matrix SIGMA, skew matrix GAMMA, 
%   zero Nu, and zero Delta, evaluated at each row of X. GAMMA is a D-by-D 
%   matrix or a scalar value, which MVCSN multiplies by the identity matrix,
%   that controls the skew of the distribution.
%
%   Y = MVCSN(X,MU,SIGMA,GAMMA,NU) returns the density of the multivariate CSN 
%   with centralization parameter MU, scale matrix SIGMA, skew matrix GAMMA,
%   NU, and zero Gamma, evaluated at each row of X. NU is a D-by-1 vector or 
%   a scalar value, which MVCSN replicates to match the size of X.
%
%   Y = MVCSN(X,MU,SIGMA,GAMMA,NU,DELTA) returns the density of the 
%   multivariate CSN with centralization parameter MU, scale matrix SIGMA, 
%   skew matrix GAMMA, NU, and DELTA, evaluated at each row of X. DELTA is 
%   a D-by-D matrix or a scalar value, which MVCSN multiplies by the identity 
%   matrix
%
%   If GAMMA=0, MVCSN becomes MVNPDF.
%
%   Example:
%
%     Mu = [0; 0]; Sigma = [2 0; 0 2]; Gamma = [-5 0; 0 -5]; Nu = [0; 0]; Delta = [3 0; 0 3];  
%     [X1,X2] = meshgrid(linspace(-5,1,100)', linspace(-5,1,100)');
%     X = [X1(:) X2(:)];
%     y = mvcsn(X, Mu, Sigma, Gamma, Nu, Delta);
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

if nargin < 2 || isempty(Mu)
    X0 = X;
    Sigma = eye(d); 
    Gamma = zeros(d); 
    Nu = zeros(1, d);
    Delta = eye(d); 
elseif nargin < 3
    if isscalar(Mu)
        Mu = Mu .* ones(1, d);
    else
        [d_Mu, ~] = size(Mu);
        if d~=d_Mu
            error("BadMuDimension");
        else
            Mu = Mu';
        end
    end
    X0 = X - Mu;
    Sigma = eye(d); 
    Gamma = zeros(d); 
    Nu = zeros(1, d); 
    Delta = eye(d);  
elseif nargin < 4
    if isscalar(Mu)
        Mu = Mu .* ones(1, d);
    else
        [d_Mu, ~] = size(Mu);
        if d~=d_Mu
            error("BadMuDimension");
        else
            Mu = Mu';
        end
    end
    X0 = X - Mu;
    if isscalar(Sigma)
        Sigma = Sigma .* eye(d); 
    else
        [d1_Sigma, d2_Sigma] = size(Sigma);
        if (d~=d1_Sigma)||(d~=d2_Sigma)
            error("BadSigmaDimension");
        end
    end
    Gamma = zeros(d); 
    Nu = zeros(1, d); 
    Delta = eye(d);  
elseif nargin < 5
    if isscalar(Mu)
        Mu = Mu .* ones(1, d);
    else
        [d_Mu, ~] = size(Mu);
        if d~=d_Mu
            error("BadMuDimension");
        else
            Mu = Mu';
        end
    end
    X0 = X - Mu;
    if isscalar(Sigma)
        Sigma = Sigma .* eye(d); 
    else
        [d1_Sigma, d2_Sigma] = size(Sigma);
        if (d~=d1_Sigma)||(d~=d2_Sigma)
            error("BadSigmaDimension");
        end
    end
    if isscalar(Gamma)
        Gamma = Gamma .* eye(d); 
    else
        [d1_Gamma, d2_Gamma] = size(Gamma);
        if (d~=d1_Gamma)||(d~=d2_Gamma)
            error("BadGammaDimension");
        end
    end
    Nu = zeros(1, d); 
    Delta = eye(d);  
elseif nargin < 6
    if isscalar(Mu)
        Mu = Mu .* ones(1, d);
    else
        [d_Mu, ~] = size(Mu);
        if d~=d_Mu
            error("BadMuDimension");
        else
            Mu = Mu';
        end
    end
    X0 = X - Mu;
    if isscalar(Sigma)
        Sigma = Sigma .* eye(d); 
    else
        [d1_Sigma, d2_Sigma] = size(Sigma);
        if (d~=d1_Sigma)||(d~=d2_Sigma)
            error("BadSigmaDimension");
        end
    end
    if isscalar(Gamma)
        Gamma = Gamma .* eye(d); 
    else
        [d1_Gamma, d2_Gamma] = size(Gamma);
        if (d~=d1_Gamma)||(d~=d2_Gamma)
            error("BadGammaDimension");
        end
    end
    if isscalar(Nu)
        Nu = Nu .* ones(1, d);
    else
        [d_Nu, ~] = size(Nu);
        if d~=d_Nu
            error("BadNuDimension");
        else
            Nu = Nu';
        end
    end
    Delta = eye(d);  
elseif nargin < 7
    if isscalar(Mu)
        Mu = Mu .* ones(1, d);
    else
        [d_Mu, ~] = size(Mu);
        if d~=d_Mu
            error("BadMuDimension");
        else
            Mu = Mu';
        end
    end
    X0 = X - Mu;
    if isscalar(Sigma)
        Sigma = Sigma .* eye(d); 
    else
        [d1_Sigma, d2_Sigma] = size(Sigma);
        if (d~=d1_Sigma)||(d~=d2_Sigma)
            error("BadSigmaDimension");
        end
    end
    if isscalar(Gamma)
        Gamma = Gamma .* eye(d); 
    else
        [d1_Gamma, d2_Gamma] = size(Gamma);
        if (d~=d1_Gamma)||(d~=d2_Gamma)
            error("BadGammaDimension");
        end
    end
    if isscalar(Nu)
        Nu = Nu .* ones(1, d);
    else
        [d_Nu, ~] = size(Nu);
        if d~=d_Nu
            error("BadNuDimension");
        else
            Nu = Nu';
        end
    end
    if isscalar(Delta)
        Delta = Delta .* eye(d); 
    else
        [d1_Delta, d2_Delta] = size(Delta);
        if (d~=d1_Delta)||(d~=d2_Delta)
            error("BadDeltaDimension");
        end
    end
else
    error("BadNumberOfInputs");
end

y = (mvncdf(zeros(size(X0)), Nu, Delta + Gamma * Sigma * Gamma')).^-1 ...
    .* mvncdf(pagemtimes(Gamma,'none', X0, 'transpose')', Nu, Delta)...
    .* mvnpdf(X0, zeros(1,d), Sigma); 
end