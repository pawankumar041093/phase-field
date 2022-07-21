function [p,dp_dxi,dp_deta] = orthopoly2D_deriv_xieta(x,n)
%
% [p,dp_dxi,dp_deta] = orthopoly2D_deriv_rst(x,n)
% Computes the ortogonal base of 2D polynomials of degree less 
% or equal to n at the point x=(r,s) in [-1,1]^2
%

% p = orthopoly2D(x,n)
% Computes the ortogonal base of 2D polynomials of degree less 
% or equal to n at the point x=(xi,eta) in the reference triangle
%

xi = x(1); eta = x(2); 

if eta==1 
    r = -1; s=1;
else
    r = 2*(1+xi)/(1-eta)-1;
    s = eta;
end

[p,dp_dxi,dp_deta] = orthopoly2D_deriv_rst([r,s],n);