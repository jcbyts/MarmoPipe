function Wrot = whiteningFromCovariance(CC, eps)
% takes as input the matrix CC of channel pairwise correlations
% outputs a symmetric rotation matrix (also Nchan by Nchan) that rotates
% the data onto uncorrelated, unit-norm axes

[E, D] 	= svd(CC); % covariance eigendecomposition (same as svd for positive-definite matrix)
D       = diag(D); % take the non-zero values from the diagonal
if nargin < 2
    eps 	= 1e-6;
end

Wrot 	= E * diag(1./(D + eps).^.5) * E'; % this is the symmetric whitening matrix (ZCA transform)
