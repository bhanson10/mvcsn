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
%      Mu = [1 -1]; Sigma = [.9 .4; .4 .3]; Beta = 2.2; 
%      [X1,X2] = meshgrid(linspace(-2,5,25)', linspace(-3,1,25)');
%      X = [X1(:) X2(:)];
%      y = mvggd(X, Mu, Sigma, Beta);
%      surf(X1,X2,reshape(y,25,25));
%   
%   References:
%      Genton, M.G. (2004). Skew-elliptical distributions and their 
%      applications: A journey beyond normality. In G. Gonzalez-Farias J.A. 
%      Domínguez-Molina A.K. Gupta, (Eds.), The closed skew-normal 
%      distribution. Chapman and Hall/ CRC. 
%
% Copyright 2025 by Benjamin L. Hanson, published under BSD 2-Clause License.

if nargin<1
    error(message('stats:mvcsn:TooFewInputs'));
elseif ~ismatrix(X)
    error(message('stats:mvcsn:InvalidData'));
end

% Get size of data.  Column vectors provisionally interpreted as multiple scalar data.
[~,d] = size(X);
if d<1
    error(message('stats:mvcsn:TooFewDimensions'));
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
    [~, d_V] = size(V);
    if d~=d_V
        error(message('stats:mvcsn:BadVDimension'));
    end
    Sigma = eye(d); 
    Delta = zeros(d); 
    Gamma = zeros(d); 
elseif nargin < 4
    X0 = X - Mu;
    [~, d_V] = size(V);
    if d~=d_V
        error(message('stats:mvcsn:BadVDimension'));
    end
    [d1_S, d2_S] = size(Sigma);
    if d1_S~=d2_S
        error(message('stats:mvcsn:NotSquareSigma'));
    elseif d1_S~=d
        error(message('stats:mvcsn:BadSigmaDimension'));
    end
    Delta = zeros(d); 
    Gamma = zeros(d); 
elseif nargin < 5
    X0 = X - Mu;
    if d~=d_V
        error(message('stats:mvcsn:BadVDimension'));
    end
    [d1_S, d2_S] = size(Sigma);
    if d1_S~=d2_S
        error(message('stats:mvcsn:NotSquareSigma'));
    elseif d1_S~=d
        error(message('stats:mvcsn:BadSigmaDimension'));
    end
    [d1_D, d2_D] = size(Delta);
    if d1_D~=d2_D
        error(message('stats:mvcsn:NotSquareDelta'));
    elseif d1_D~=d
        error(message('stats:mvcsn:BadDeltaDimension'));
    end
    Gamma = zeros(d); 
elseif nargin < 6
    X0 = X - Mu;
    if d~=d_V
        error(message('stats:mvcsn:BadVDimension'));
    end
    [d1_S, d2_S] = size(Sigma);
    if d1_S~=d2_S
        error(message('stats:mvcsn:NotSquareSigma'));
    elseif d1_S~=d
        error(message('stats:mvcsn:BadSigmaDimension'));
    end
    [d1_D, d2_D] = size(Delta);
    if d1_D~=d2_D
        error(message('stats:mvcsn:NotSquareDelta'));
    elseif d1_D~=d
        error(message('stats:mvcsn:BadDeltaDimension'));
    end
    [d1_G, d2_G] = size(Gamma);
    if d1_G~=d2_G
        error(message('stats:mvcsn:NotSquareGamma'));
    elseif d1_G~=d
        error(message('stats:mvcsn:BadGammaDimension'));
    end
else
    error(message('stats:mvcsn:BadInputs'));
end

y = (mvncdf(zeros(size(X)), V, Delta + Gamma * Sigma * Gamma')).^-1 ...
    .* (mvncdf(Gamma * X0, V, Delta))...
    .* mvnpdf(X0, zeros(d,1), Sigma); 
end