function H = computeH(ux,uy,H_previous,referenceElement,X,T,lambda,mu)

% H(:,iElem) = H evaluated in IP of element iElem

nOfElements = size(T,1);
nIP = length(referenceElement.IPweights);

Nxi=referenceElement.Nxi;
Neta=referenceElement.Neta;

H = zeros(nIP,nOfElements);

for iElem = 1:nOfElements
    Te = T(iElem,:);
    Xe=X(Te,:);
    ux_e = ux(Te);
    uy_e = uy(Te);

    [J1,J2,J3,J4] = computeJ_IP(Xe,ux_e,uy_e,Nxi,Neta);
    strain = [J1,J4,1/2*(J2+J3)]; %strain(i,:) = [eps_x, eps_y, eps_xy]
    trStrain = strain(:,1)+strain(:,2);
    trStrainSquared = strain(:,1).^2+2*(strain(:,3).^2)+strain(:,2).^2;
    elasticEnergyDensity = 1/2*lambda*trStrain.^2 + mu*trStrainSquared;
    H(:,iElem) = max(H_previous(:,iElem),elasticEnergyDensity);
end

end

function [J1,J2,J3,J4] = computeJ_IP(Xe,ux,uy,Nxi,Neta)
    J11 = Nxi*Xe(:,1); J12 = Nxi*Xe(:,2); 
    J21 = Neta*Xe(:,1); J22 = Neta*Xe(:,2); 
    detJ = J11.*J22-J12.*J21;
    invJ11 = diag(J22./detJ);
    invJ12 = diag(-J12./detJ);
    invJ21 = diag(-J21./detJ);
    invJ22 = diag(J11./detJ);
    Nx = invJ11*Nxi + invJ12*Neta;
    Ny = invJ21*Nxi + invJ22*Neta;  
    J1 = Nx*ux; J2 = Ny*ux; 
    J3 = Nx*uy; J4 = Ny*uy;
end