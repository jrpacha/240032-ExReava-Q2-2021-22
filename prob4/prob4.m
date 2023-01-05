clearvars
close all
clc

fileName = 'prob4.txt';

%Material properties
%th=1.0;               %thickness (in mm)
E=4.0e+8;              %Young's modulus (in N/mm^2)
nu=0.18;               %Poisson's ratio (adimensional)
Ft = [0;-1026];        %N/m
elemVM = 198;

%Define the plane elasticity problem: 
%modelProblem=1; %plane stress
%modelProblem=2; %plane strain (2)

modelProblem=2; %Plane strain
          
eval('channel1'); %load the mesh k

[numNod,ndim]=size(nodes);
numElem=size(elem,1);

numbering=0; %=1 shows nodes and element numbering
plotElementsOld(nodes,elem,numbering);
hold on

%Find Boundary points
[indNodBd,~,~,~] = boundaryNodes(nodes, elem);
indLeft=find(nodes(:,1) < 0.01);
indRight=find(nodes(:,1) > 4.99);
indBottom = find(nodes(:,2) < 0.01);
indTop = setdiff(indNodBd,unique([indLeft;indRight;indBottom]));

plot(nodes(indBottom,1),nodes(indBottom,2),'o','color','black',...
    'markerFaceColor','magenta','markerSize',7)
plot(nodes(indLeft,1),nodes(indLeft,2),'o','color','black',...
    'markerFaceColor','green','markerSize',7)
plot(nodes(indRight,1),nodes(indRight,2),'o','color','black',...
    'markerFaceColor','blue','markerSize',7)
plot(nodes(indTop,1),nodes(indTop,2),'o','color','black',...
    'markerFaceColor','yellow','markerSize',7)
hold off

%% modelProblem
switch modelProblem
    case 1
        c11=E/(1-nu^2);
        c22=c11;
        c12=nu*c11;
        c21=c12;
        c33=E/(2*(1+nu));
        fprintf('Plane stress problem\n')
    case 2
        th=1.0;
        c11=E*(1-nu)/((1+nu)*(1-2*nu));
        c22=c11;
        c12=c11*nu/(1-nu);
        c21=c12;
        c33=E/(2*(1+nu));
        fprintf('Plane strain problem\n')
    otherwise
        error('modelProblem should be 1 (stress) or 2 (strain)');
end
C=[c11, c12, 0; c21, c22, 0; 0, 0, c33];

%Computation of the stiffness matrix
K=zeros(ndim*numNod);
F=zeros(ndim*numNod,1);
Q=zeros(ndim*numNod,1);

for e=1:numElem
    Ke=planeElastTriangStiffMatrix(nodes,elem,e,C,th);
    %
    % Assemble the stiffness matrices
    %
    row=[2*elem(e,1)-1; 2*elem(e,1); ...
         2*elem(e,2)-1; 2*elem(e,2); ...
         2*elem(e,3)-1; 2*elem(e,3)];
    col=row;
    K(row,col)=K(row,col)+Ke;
end
%Boundary conditions
%Natural BC: constant traction on the right boundary
nodLoads=indBottom'; %nodes where the load is applied (transposed)
Q=applyLoadsTriang(nodes,elem,nodLoads,Q,Ft);
%Essential BC: left boundary fixed
u=zeros(ndim*numNod,1);
fixedNod=[ndim*indLeft-1;ndim*indLeft;ndim*indRight];
u(fixedNod)=0.0;

%Reduced System
freeNod=setdiff(1:ndim*numNod,fixedNod);
%Fm=F(freeNod)-K(freeNod,fixedNod)*u(fixedNod); %Note: here u(fixedNod)=0,
                                                %so this is not necessary
                                                %in this case
%Fm=Fm+Q(freeNod);
Qm=Q(freeNod);
Km=K(freeNod,freeNod);

%Solve the reduced system
um=Km\Qm;
u(freeNod)=um;

%Post Process
%Compute the Stress & Strain
stress=zeros(numElem,3);
strain=zeros(numElem,3);
vonMisses=zeros(numElem,1);
for e=1:numElem
    v1=nodes(elem(e,1),:);
    v2=nodes(elem(e,2),:);
    v3=nodes(elem(e,3),:);
    beta=[v2(2)-v3(2),v3(2)-v1(2),v1(2)-v2(2)];
    gamma=-[v2(1)-v3(1),v3(1)-v1(1),v1(1)-v2(1)];
    Area=0.5*det([v1 1; v2 1; v3 1]);
    B=[beta(1), 0, beta(2), 0, beta(3), 0; 
       0, gamma(1), 0 gamma(2), 0 gamma(3);
       gamma(1), beta(1), gamma(2), beta(2), gamma(3), beta(3)]/(2*Area);   
    row=[2*elem(e,1)-1; 2*elem(e,1); ...
         2*elem(e,2)-1; 2*elem(e,2); ...
         2*elem(e,3)-1; 2*elem(e,3)];
    ue=u(row,:);
    strain(e,:)=B*ue;   %This is constant for each element
    stress(e,:)=C*B*ue; %This as constant for each element
    sxx=stress(e,1);
    syy=stress(e,2);
    sxy=stress(e,3);
    vonMisses(e)=sqrt(sxx^2+syy^2-sxx*syy+3*sxy^2);
end

%Output (fancy output: don't try this at the exams!)
%Displacements only the last 10 nodes:
displacements=[u(1:2:end),u(2:2:end)];
tableDispl=[(1:numNod)',nodes(:,1),nodes(:,2),displacements];
fprintf('\n%51s\n\n','Displacements (only for the last 10 nodes)')
fprintf('%7s%8s%12s%12s%11s\n',...
    'Num.Nod.','X','Y','U','V')
fprintf('%4d%16.4e%12.4e%12.4e%12.4e\n',tableDispl(end-9:end,:)');    
%Stress:
tableStress=[(1:numElem)',stress,vonMisses]; 
fprintf('\n%48s\n\n','Stress (only for the last 10 elements)')
fprintf('%5s%12s%12s%12s%15s\n','Elem.','SXX','SYY','SXY','vonMisses')
fprintf('%4d%16.4e%12.4e%12.4e%12.4e\n',tableStress(end-9:end,:)')

%Graphical output
esc=2500;
plotPlaneNodElemDespl(nodes, elem, u, esc);
valueToShow=u(1:2:end);
title='Desp.X';
colorScale='jet';
plotContourSolution(nodes,elem,valueToShow,title,colorScale);
plotStressVM(nodes,elem,vonMisses);

%Solutions
UX = u(1:2:end);
UY = u(2:2:end);
[~,indMaxNorm]=max(sqrt(UX.^2+UY.^2));

fOut = fopen(fileName,'w');
fprintf(fOut,'Prob 4\n');
fprintf(fOut,'(a) Mean of the y-component of the nodes in the top boundary\n');
fprintf(fOut,'    of the domain (without the two corners): %.4e\n',...
    sum(nodes(indTop,2))/length(indTop));
fprintf(fOut,'    Hint. The number of nodes (without the corners) in the\n');
fprintf(fOut,'          top boundary of the domain is %d\n',length(indTop));
fprintf(fOut,'(b) The maximum of the absolute value of the vertical\n');
fprintf(fOut,'    displacement of the nodes %.4e\n',max(abs(UY)));
fprintf(fOut,'    Hint. The maximum of the absolute value of the horizontal\n');
fprintf(fOut,'          displacement fo the nodes %.4e\n',max(abs(UX)));
fprintf(fOut,'(c) Node in which it reaches the maximum of the norm in the\n');
fprintf(fOut,'    displacement %d\n',indMaxNorm);
fprintf(fOut,'(d) The Von Misses stress of element %d is sigmaVM(%d) =\n %.4e\n',...
    elemVM,elemVM,vonMisses(elemVM)); 
type(fileName)
