function pde = S4LaplacianBeltramidata
%% S4LAPLACIANBELTRAMIDATA Data for Laplace-Beltrami equation on S^4
%
% This provides test data for the Laplace-Beltrami equation on the 4-sphere
% using stereographic projection coordinates.
%
% The equation is: -Delta_M u + b*u = f on S^4
%
% See also: sphereS4Poisson
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

dim = 4;
b = 1;

pde = struct('f', @f, 'exactu', @exactu, 'stereographicProj', @stereographicProj, 'dim', dim, 'b', b);
    %% Right-hand side function
    function val = f(center, b, chart)
        % Center: element centers in stereographic coordinates
        % b: L2 coefficient
        % chart: 'lower' or 'upper'
        switch chart
            case 'lower'
                u = lowerStgProjection(center);
            case 'upper'
                u = upperStgProjection(center);
            otherwise
                error('Chart must be ''lower'' or ''upper''');
        end
        val = (dim + b) * u;
    end

    %% Exact solution at nodes
    function u = exactu(node, r, h, chart)
        % node: node coordinate matrix
        % h: mesh size
        % chart: 'lower' or 'upper'
        pt = -r + node * h;
        
        switch chart
            case 'upper'
                u = upperStgProjection(pt);
            case 'lower'
                u = lowerStgProjection(pt);
            otherwise
                error('Chart must be lower or upper');
        end
    end

    %% Analytical solution wrapper
    function u = stereographicProj(pt, chart)
        % pt: points in stereographic coordinates
        % chart: 'lower' or 'upper'
        switch chart
            case 'upper'
                u = upperStgProjection(pt);
            case 'lower'
                u = lowerStgProjection(pt);
            otherwise
                error('Chart must be lower or upper');
        end
    end

    %% Stereographic projection functions
    function u=upperStgProjection(pt)
        ns = dot(pt,pt,2);
        u=(1-ns)./(1+ns);
    end

    function u=lowerStgProjection(pt)
        ns = dot(pt,pt,2);
        u=(-1+ns)./(1+ns);
    end


end