function U = spin1_transform(c)
% SPIN1_TRANSFORM Generate a unitary matrix for spin-1 system transformation.
%   U = spin1_transform(c) returns a 3x3 unitary matrix U such that
%   U * [1; 0; 0] = c, where c is the target state vector (length 3).
%
%   Input:
%       c: Target state vector (3-element column or row vector, complex or real).
%          It should be normalized, but the function normalizes it if necessary.
%   Output:
%       U: 3x3 unitary matrix (complex if c is complex).
%
%   Example:
%       c = [1; 0; 0]; % Target state
%       U = spin1_transform(c); % Returns identity matrix
%
%       c = [1i; 0; 1]/sqrt(2); % Complex target
%       U = spin1_transform(c);
%       norm(U'*U - eye(3)) % Should be near zero (unitary check)
%       U * [1; 0; 0] % Should equal c

% Ensure c is a column vector
c = c(:);

% Normalize c to handle potential non-normalized input
n = norm(c);
if abs(n) < eps
    error('Input vector is zero or near zero.');
end
c = c / n;

% Compute the null space of c' (conjugate transpose of c)
% null(c') returns a 3x2 matrix whose columns form an orthonormal basis for the null space
N = null(c'); % N is 3x2, columns are orthonormal to c

% Construct the unitary matrix: first column is c, next columns are from N
U = [c, N];

% Optional: Verify unitarity (for debugging, not necessary in production)
% if max(abs(U'*U - eye(3)), [], 'all') > 1e-10
%     warning('Unitarity check failed. Numerical inaccuracies may be present.');
% end
end