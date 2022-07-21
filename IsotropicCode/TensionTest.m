%% TENSION TEST
%% Isotropic phase-field model for fracture

clear all, close all
setpath

%% Parameters of the problem
E = 210; nu = 0.3; 
lambda = E*nu/((1+nu)*(1-2*nu)); mu = E/(2*(1+nu));
tol = 1.e-2; 
Gc = 2.7e-3;
l  = 0.015;
increments = 1e-4:1e-4:6.5e-3 ;
%% Mesh
degree = 1;
load('MeshTension.mat');
figure(1),clf,kk = plotMesh(X,T);
nOfElements = size(T,1);
x = X(:,1); y = X(:,2); 
nOfNodes = size(X,1);

%Linear elasticity boundary nodes
nodesTop = find(abs(y-0.5)<tol); 
nodesBottom = find(abs(y+0.5)<tol);
nodesCCD = [nodesBottom;nodesTop];

loads = zeros(1,length(increments));

referenceElement = createReferenceElementTri(1);
nGaussPoints = length(referenceElement.IPweights); 
d = zeros(nOfNodes,1);
H = zeros(nGaussPoints,nOfElements); 
Hprevious = H ;
%% Loop in load steps
for ind=1:length(increments)
    k = increments(ind) ;
    disp(sprintf('Load step %d -- Prescribed displacement %0.5g',ind,k))
    stoppingCriterion = 0; s = 0;

    %% Staggered iterations
    while (stoppingCriterion == 0)
        s = s+1;
        disp(sprintf('  Iteration %d',s))
        %% Linear elasticity
        [ind_i,ind_j,coef_K,f] = FEMSystemLinearElasticityDamage(X,T,d,referenceElement,E,nu);
        K = sparse(ind_i,ind_j,coef_K);
        
        % Dirichlet BC 
        uCCD = [zeros(length(nodesCCD)+length(nodesBottom),1); k+zeros(length(nodesTop),1)];
        dofCCD = [nodesCCD; nodesCCD+size(X,1)];

        %% BC by system reduction       
        notCCD= setdiff(1:(2*size(X,1)),dofCCD); %actual degrees of freedom (not boundary nodes)
        f = f(notCCD)-K(notCCD,dofCCD)*uCCD;
        Kdn = K(dofCCD,notCCD); Kdd = K(dofCCD,dofCCD);
        K = K(notCCD,notCCD);
        sol = K\f;
        %   Nodal values
        u = zeros(2*size(X,1),1);
        u(notCCD) = sol; u(dofCCD) = uCCD;

        %% History field H
        ux = u(1:size(X,1)); uy = u(size(X,1)+1:end);
        H = computeH(ux,uy,Hprevious,referenceElement,X,T,lambda,mu);
        %% Damage field equation           
        % Computation
        % Loop in elements
        [ind_i,ind_j,coef_K,f_damage] = FEMSystemDamageField(X,T,referenceElement,Gc,l,H);
        K_damage = sparse(ind_i,ind_j,coef_K);
        d = K_damage\f_damage;
        
        %% Stopping criterion for staggered iterations
        if (s == 1)
            dprevious = d;
            uprevious = u;
        else
            errorDisplacement = computeEuclideanNormRelative(u,uprevious);
            errorDamage = computeEuclideanNormRelative(d,dprevious);           
            if(errorDisplacement<tol && errorDamage<tol)
                stoppingCriterion=1;
            else
                stoppingCriterion=0;
            end      
            dprevious = d;
            uprevious = u;
        end  
    end
    
    %% Load
    lagrange_multipliers = Kdn*sol+Kdd*uCCD; %load == langrangemultipliers
    lm_top = lagrange_multipliers(length(nodesCCD)+length(nodesBottom)+1:end);
    load = sum(lm_top);
    loads(ind) = load;
    %% History variable
    Hprevious = H ;
end


%% Plot damage field
figure(2),clf,plotContinuosSolution(X,T,d,referenceElement,4),colorbar
set(gca,'FontSize',17);
title('FEM Damage field')

%% Plot force-displacement curve
figure(3)
plot(increments,loads,'ro-','LineWidth',2.5)
xlabel('Displacement','FontSize',24);
ylabel('Force','FontSize',24,'VerticalAlignment','bottom');
axis square ; 
set(gca,'Fontsize',24)


