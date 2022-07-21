function [P,dPdxi]=orthogonalPolynomialsAndDerivatives1D(degree,Xi)

x = Xi;
U = ones(size(x));
O = zeros(size(x));

switch degree
    case 1
        P = [U x];
        dPdxi = [O U];
    case 2
        P = [U, x, -U/2+(3/2)*x.^2];
        dPdxi = [O, U, 3*x];
    case 3
        P = [U, x, -U/2+(3/2)*x.^2, (5/2)*x.^3-x*3/2];
        dPdxi = [O, U, 3*x, (15/2)*x.^2-U*3/2];
    case 4
        P = [U, x, -U/2+(3/2)*x.^2, (5/2)*x.^3-x*3/2, U*3/8+(35/8)*x.^4-(15/4)*x.^2];
        dPdxi = [O, U, 3*x, (15/2)*x.^2-U*3/2, (35/2)*x.^3-(15/2)*x];
    case 5
        P = [U, x, -U/2+(3/2)*x.^2, (5/2)*x.^3-x*3/2, U*3/8+(35/8)*x.^4-(15/4)*x.^2, (63/8)*x.^5-(35/4)*x.^3+(15/8)*x];
        dPdxi = [O, U, 3*x, (15/2)*x.^2-U*3/2, (35/2)*x.^3-(15/2)*x, (315/8)*x.^4-(105/4)*x.^2+15/8];
    otherwise
        error('Degree not implemented in orthogonalPolynomialsAndDerivatives1D')
end

