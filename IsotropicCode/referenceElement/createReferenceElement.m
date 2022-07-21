function theReferenceElement=createReferenceElement(elementType,nOfElementNodes,handle,varargin)

% theReferenceElement=createReferenceElement(elementType,nOfElementNodes)
% Input:
%  elementType: 0 for quadrilateral, 1 for triangle
%  nOfElementNodes: number of nodes in the reference element
%  nOfGaussPoints (optional): n� of gauss points of the 1D quadrature. The
%                             default value is nDeg + 2 where nDeg is the
%                             degree of interpolation
% Output:
%  theReferenceElement: struct containing
%     .IPcoordinates: coordinates of the integration points for 2D elemens
%     .IPweights: weights of the integration points for 2D elements
%     .N: shape functions at the IP
%     .Nxi,.Neta: derivatives of the shape functions at the IP
%     .IPcoordinates1d: coordinates of the integration points for 1D boundary elemens
%     .IPweights1d: weights of the integration points for 1D boundary elements
%     .N1d: 1D shape functions at the IP
%     .N1dxi: derivatives of the 1D shape functions at the IP
%     .faceNodes: matrix [nOfFaces nOfNodesPerFace] with the edge nodes numbering
%     .innerNodes: vector [1 nOfInnerNodes] with the inner nodes numbering
%     .faceNodes1d: vector [1 nOfNodesPerElement] with the 1D nodes numbering
%     .NodesCoord: spatial coordinates of the element nodes
%     .NodesCoord1d: spatial coordinates of the 1D element nodes

if elementType == 1 %Triangles
    
    switch nOfElementNodes
        case 3 %P1
            nDeg = 1;
            faceNodes = [1 2; 2 3; 3 1];
            innerNodes = [];
            faceNodes1d = 1:2;
            coord2d = [-1 -1; 1 -1; -1 1];
            coord1d = [-1; 1];
        case 6 %P2
            nDeg = 2;
            faceNodes = [1 4 2; 2 5 3; 3 6 1];
            innerNodes = [];
            faceNodes1d = 1:3;
            coord2d = [-1 -1; 1 -1; -1 1; 0 -1; 0 0; -1 0];
            coord1d = [-1; 0; 1];
        case 10 %P3
            nDeg = 3;
            faceNodes = [1 4 5 2; 2 6 7 3; 3 8 9 1];
            innerNodes = 10;
            faceNodes1d = 1:4;
        case 15 %P4
            nDeg = 4;
            faceNodes = [1 4 5 6 2; 2 7 8 9 3; 3 10 11 12 1];
            innerNodes = 13:15;
            faceNodes1d = 1:5;
        case 21 %P5
            nDeg = 5;
            faceNodes = [1 4:7 2; 2 8:11 3; 3 12:15 1];
            innerNodes = 16:21;
            faceNodes1d = 1:6;
        case 28 %P6
            nDeg = 6;
            faceNodes = [1 4:8 2; 2 9:13 3; 3 14:18 1];
            innerNodes = 19:28;
            faceNodes1d = 1:7;
        case 36 %P7
            nDeg = 7;
            faceNodes = [1 4:9 2; 2 10:15 3; 3 16:21 1];
            innerNodes = 22:36;
            faceNodes1d = 1:8;
        case 45 %P8
            nDeg = 8;
            faceNodes = [1 4:10 2; 2 11:17 3; 3 18:24 1];
            innerNodes = 25:45;
            faceNodes1d = 1:9;
        case 55 %P9
            nDeg = 9;
            faceNodes = [1 4:11 2; 2 12:19 3; 3 20:27 1];
            innerNodes = 28:55;
            faceNodes1d = 1:10;
        case 66 %P10
            nDeg = 10;
            faceNodes = [1 4:12 2; 2 13:21 3; 3 22:30 1];
            innerNodes = 31:66;
            faceNodes1d = 1:11;
        case 78 %P11
            nDeg = 11;
            faceNodes = [1 4:13 2; 2 14:23 3; 3 24:33 1];
            innerNodes = 34:78;
            faceNodes1d = 1:12;
        otherwise
            setOutput({'??? Reference element not implemented yet'},handle)
            error('Error in Berkhoff GUI')
    end
    
    if nDeg >= 3 %EZ4U rules imply fekete nodes
        %coord2d = feketeNodesTri2D(nDeg,faceNodes,innerNodes);
        feketeFile = load('positionFeketeNodesTri2D_EZ4U.mat'); %EZ4U reference element
        coord2d = feketeFile.feketeNodesPosition.(['P' num2str(nDeg)]);
        coord1d = feketeNodes1D(nDeg,faceNodes1d);
%                         [coord2d coord1d] = equispacedNodesReferenceElement(nDeg);  % for equispaced nodes.. (20/12/2010 G.G.)
    end
    
elseif elementType == 0 %Quadrilaterals
    
    switch nOfElementNodes
        case 4 %Q1
            nDeg = 1;
            faceNodes = [1 2; 2 3; 3 4; 4 1];
            innerNodes = [];
            faceNodes1d = 1:2;
            coord2d = [-1 -1; 1 -1; 1 1; -1 1];
            coord1d = [-1; 1];
        case 9 %Q2
            nDeg = 2;
            faceNodes = [1 5 2; 2 6 3; 3 7 4; 4 8 1];
            innerNodes = 9;
            faceNodes1d = 1:3;
            coord2d = [-1 -1; 1 -1; 1 1; -1 1; 0 -1; 1 0; 0 1; -1 0; 0 0];
            coord1d = [-1; 0; 1];
        otherwise
            setOutput({'??? Reference element not implemented yet'},handle)
            error('Error in Berkhoff GUI')
    end
    
else
    error('Element not allowed')
end

%Compute shape functions and quadrature
if isempty(varargin)
    %     nOfGaussPoints = nDeg + 2;
    nOfGaussPoints = 2*nDeg + 1;
else
    nOfGaussPoints = varargin{:};
end
[shapeFun2d,gw2d,gp2d] = ...
    computeShapeFunctionsReferenceElement(nDeg,coord2d,nOfGaussPoints,elementType);
[shapeFun1d,gw1d,gp1d] = ...
    computeShapeFunctionsReferenceElement(nDeg,coord1d,nOfGaussPoints);
N = shapeFun2d(:,:,1)';
Nxi = shapeFun2d(:,:,2)';
Neta = shapeFun2d(:,:,3)';
N1 = shapeFun1d(:,:,1)';
Nxi1 = shapeFun1d(:,:,2)';

%Creating reference element structure
theReferenceElement = struct('IPcoordinates',gp2d,...
    'IPweights',gw2d,'N',N,'Nxi',Nxi,'Neta',Neta,...
    'IPcoordinates1d',gp1d,'IPweights1d',gw1d,...
    'N1d',N1,'N1dxi',Nxi1,'faceNodes',faceNodes,...
    'innerNodes',innerNodes,'faceNodes1d',faceNodes1d,...
    'NodesCoord',coord2d,'NodesCoord1d',coord1d,'degree',nDeg);

