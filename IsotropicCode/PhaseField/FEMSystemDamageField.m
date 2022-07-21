function [ind_i,ind_j,coef_K,f]=FEMSystemDamageField(X,T,referenceElement,Gc,l,H_g)

GaussWeights=referenceElement.IPweights;
N=referenceElement.N;
Nxi=referenceElement.Nxi;
Neta=referenceElement.Neta;

nOfNodes = size(X,1);
nOfElements = size(T,1);
nOfElementNodes = size(T,2);

f = zeros(nOfNodes,1);
%[mKe,nKe] = size(Ke)
mKe = nOfElementNodes;
nKe = nOfElementNodes;
ind_i = zeros(mKe*nKe,nOfElements);
ind_j = zeros(mKe*nKe,nOfElements);
coef_K = zeros(mKe*nKe,nOfElements);

%Loop in elements
for i = 1:nOfElements
    Te = T(i,:); %nodes in the element
    Xe = X(Te,:); %coordinates of the element nodes  
    He_g = H_g(:,i);
    [Ke,fe]=computeElementalMatrices(Xe,GaussWeights,N,Nxi,Neta,Gc,l,He_g);    
    
    coef_K(:,i) = Ke(:);
    [mi,mj] = meshgrid(Te,Te);
    ind_i(:,i) = mi(:); ind_j(:,i) = mj(:);
    f(Te) = f(Te) + fe;
end


%_______________________________________
%Computation of elemental matrix and vector
function [Ke,fe]=computeElementalMatrices(Xe,GaussWeights,N,Nxi,Neta,Gc,l,Hg)

J11 = Nxi*Xe(:,1); J12 = Nxi*Xe(:,2); 
J21 = Neta*Xe(:,1); J22 = Neta*Xe(:,2); 
detJ = J11.*J22-J12.*J21;
dvolu = GaussWeights.*detJ;
dvolu_diag = diag(dvolu);
invJ11 = diag(J22./detJ);
invJ12 = diag(-J12./detJ);
invJ21 = diag(-J21./detJ);
invJ22 = diag(J11./detJ);
% xy-derivatives
Nx = invJ11*Nxi + invJ12*Neta;
Ny = invJ21*Nxi + invJ22*Neta;
%Elemental matrices
fe = 2*N'*(dvolu_diag*Hg);
Kxxe = Nx'*(dvolu_diag*Nx);
Kyye = Ny'*(dvolu_diag*Ny);
Ke=Gc*l*(Kxxe+Kyye) + N'*(diag(dvolu.*(Gc/l+2*Hg))*N);
 