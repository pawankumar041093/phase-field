function plotContinuosSolution(X,T,u,referenceElement,nDegRef)

% Check input
if nargin == 4
    nDegRef = 40;
end

% Plotting element (equal spaced points)
nodes = [];
h = 1/nDegRef;
for j = 0:nDegRef
    i = (0:nDegRef-j)';
    aux = j*ones(size(i));
    nodes = [nodes; [i aux]*h];
end
nodes = 2*nodes - 1;

% Delaunay triangulation of plotting element
elemTriRef = delaunayn(nodes);
coordRef=referenceElement.NodesCoord;
[shapeFunctions,~,~]=evaluateNodalBasisTri(nodes,coordRef,referenceElement.degree);

nOfElements = size(T,1);

Xplot = []; uplot = []; Tplot = [];
% Loop in elements
for ielem = 1:nOfElements
    Te = T(ielem,:);
    Tplot = [Tplot; elemTriRef + size(Xplot,1)];
    Xplot = [Xplot; shapeFunctions*X(Te,:)];
    uplot = [uplot; shapeFunctions*u(Te)];
end
trisurf(Tplot,Xplot(:,1),Xplot(:,2),uplot)
axis equal

view(0,90)
shading interp