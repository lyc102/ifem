function secondOrderCoeff = getMetricS4(dim, pt)
%  GETMETRIC Compute second-order metric coefficient for Laplace-Beltrami
%
%   metricTensor = getMetricS4(dim, coordinate) computes the metric-dependent
%   scaling factor that arises from the pullback of the Laplace-Beltrami operator
%   on S^d to the stereographic projection chart.
%
%   The conformal metric in stereographic coordinates is:
%       g_ij = (4/(1+|x|^2)^2) * delta_{ij}
%
%   The second-order coefficient is:
%       metricTensor = 2^(d-2) * (1 + |x|^2)^(2-d)
%
%   This transforms the flat Laplacian on the cube to the Laplace-Beltrami
%   operator on the sphere.
%
% Inputs:
%   dim- Embedding dimension
%   pt - Element centers in stereographic coordinates [nElem x dim]
%
% Outputs:
%   metricTensor - Metric coefficient at element centers [nElem x 1]

    ns = dot(pt, pt, 2);
    secondOrderCoeff = 2^(dim - 2) * (1 + ns).^(2 - dim);
end