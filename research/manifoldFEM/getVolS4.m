function volForm = getVolS4(dim, pt)
%  GETVOLS4 Compute volume form for S^d in stereographic coordinates
%
%   volForm = getVolS4(dim, pt) computes the volume form that arises
%   from the pullback of the standard volume form on S^d to the stereographic
%   projection chart.
%
%   The volume form in stereographic coordinates is:
%       volForm = 2^d * (1 + |x|^2)^(-d)
%
%   This accounts for the conformal distortion of the stereographic projection.
%
% Inputs:
%   dim - Embedding dimension
%   pt  - Element centers in stereographic coordinates [nElem x dim]
%
% Outputs:
%   volForm - Volume form at element centers [nElem x 1]

    x2 = dot(pt, pt, 2);
    volForm = 2^dim * (1 + x2).^(-dim);
end