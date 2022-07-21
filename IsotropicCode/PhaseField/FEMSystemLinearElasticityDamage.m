function [ind_i,ind_j,coef_K,f] = FEMSystemLinearElasticityDamage(X,T,d,referenceElement,E,nu)

GaussWeights=referenceElement.IPweights;
N=referenceElement.N;
Nxi=referenceElement.Nxi;
Neta=referenceElement.Neta;

nOfNodes = size(X,1);
nOfElements = size(T,1);
nOfElementNodes = size(T,2);

f = zeros(2*nOfNodes,1);

%[mKe,nKe] = size(Ke)
mKe = 2*nOfElementNodes;
nKe = mKe;
ind_i  = zeros(mKe*nKe,nOfElements); ind_j  = zeros(mKe*nKe,nOfElements); coef_K = zeros(mKe*nKe,nOfElements); 
%Loop in elements
for i=1:nOfElements
    Te=T(i,:); %nodes in the element
    Xe=X(Te,:); %coordinates of the element nodes
    de = d(Te);
    Ke=computeElementalMatrices(Xe,de,GaussWeights,N,Nxi,Neta,E,nu);    
    
    
    ind = [Te,Te+size(X,1)];
    coef_K(:,i) = Ke(:);
    [mi,mj] = meshgrid(ind,ind);
    ind_i(:,i) = mi(:);
    ind_j(:,i) = mj(:);
end

end


%_______________________________________
%Computation of elemental matrix and vector
function Ke = computeElementalMatrices(Xe,de,GaussWeights,N,Nxi,Neta,E,nu)
J11 = Nxi*Xe(:,1); J12 = Nxi*Xe(:,2); 
J21 = Neta*Xe(:,1); J22 = Neta*Xe(:,2); 
detJ = J11.*J22-J12.*J21;
dvolu = GaussWeights.*detJ;
invJ11 = diag(J22./detJ);
invJ12 = diag(-J12./detJ);
invJ21 = diag(-J21./detJ);
invJ22 = diag(J11./detJ);
% xy-derivatives
Nx = invJ11*Nxi + invJ12*Neta;
Ny = invJ21*Nxi + invJ22*Neta;
dg = N*de;
dvolu_d = diag(dvolu.*((1-dg).^2));
%Elemental matrices
Kxxe = Nx'*(dvolu_d*Nx);
Kyye = Ny'*(dvolu_d*Ny);
Kxye = Nx'*(dvolu_d*Ny);
Kyxe = Ny'*(dvolu_d*Nx);
Ke=E/((1+nu)*(1-2*nu))*[(1-nu)*Kxxe+(1-2*nu)/2*Kyye,nu*Kxye+(1-2*nu)/2*Kyxe;
                        nu*Kyxe+(1-2*nu)/2*Kxye,(1-nu)*Kyye+(1-2*nu)/2*Kxxe];

end  