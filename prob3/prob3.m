clearvars
close all

fileSols = 'prob3.txt';

kc = 3.0;
ff = 27.67;

tempBottom = 16.47;
tempTop = 29.08;

qn = 10.0; %on the circle centered at (6,1). Circle 2

beta1 = 7.0; TInf1 = 25.0; %on the circles centered at (2,1) and (6,3). Circle 1 amd 4
beta2 = 6.0; TInf2 = 24.0; %on the circles centered at (2,3). Circle 3

elemB = 50;
elemD = 37;

eval('meshPlaca4foratsTriang');

numNodes= size(nodes,1);
numElem= size(elem,1);

numbering= 0; %= 1 shows nodes and element numbering
plotElementsOld(nodes,elem,numbering);

% Select Boundary points

indT = find(nodes(:,2) > 3.99); %indices of the nodes at the top boundary
indB = find(nodes(:,2) < 0.01); %indices of the nodes at the bottom boundary
indCirc1 = find(sqrt((nodes(:,1)-2).^2 + (nodes(:,2)-1).^2) < 0.501); 
indCirc2 = find(sqrt((nodes(:,1)-6).^2 + (nodes(:,2)-1).^2) < 0.501);
indCirc3 = find(sqrt((nodes(:,1)-2).^2 + (nodes(:,2)-3).^2) < 0.501); 
indCirc4 = find(sqrt((nodes(:,1)-6).^2 + (nodes(:,2)-3).^2) < 0.501);

hold on
plot(nodes(indT,1),nodes(indT,2),'ok','lineWidth',1,'markerFaceColor',...
    'red','markerSize',5)
plot(nodes(indB,1),nodes(indB,2),'ok','lineWidth',1,'markerFaceColor',...
    'blue','markerSize',5)
plot(nodes(indCirc1,1),nodes(indCirc1,2),'ok','lineWidth',1,'markerFaceColor',...
    'green','markerSize',5)
plot(nodes(indCirc2,1),nodes(indCirc2,2),'ok','lineWidth',1,'markerFaceColor',...
    'magenta','markerSize',5)
plot(nodes(indCirc3,1),nodes(indCirc3,2),'ok','lineWidth',1,'markerFaceColor',...
    'yellow','markerSize',5)
plot(nodes(indCirc4,1),nodes(indCirc4,2),'ok','lineWidth',1,'markerFaceColor',...
    'green','markerSize',5)
hold off

%% Define the coefficients vector of the model equation
a11=kc;
a12=0;
a21=a12;
a22=a11;
a00=0;
f=ff;
coeff=[a11,a12,a21,a22,a00,f];

%Compute the global stiff matrix
K=zeros(numNodes);    %global stiff matrix
F=zeros(numNodes,1);  %global internal forces vector
Q=zeros(numNodes,1);  %global secondary variables vector

for e = 1:numElem
    [Ke, Fe] = linearTriangElement(coeff,nodes,elem,e);
    rows= [elem(e,1); elem(e,2); elem(e,3)];
    cols= rows;
    K(rows,cols)= K(rows,cols)+Ke;
    if (coeff(6) ~= 0)
        F(rows)= F(rows) + Fe;
    end
end %end for elements
Kini=K; %We save a copy of the initial K and F arrays
Fini=F; %for the post-process step 

%Booundary Conditions
fixedNodes= [indB',indT'];                 %fixed Nodes (global numbering)
freeNodes= setdiff(1:numNodes,fixedNodes); %free Nodes (global numbering)

% Natural BC
%------------- Constant natural BC
Q=applyConstantNaturalBC(nodes,elem,indCirc2',qn,Q);

%------------- Convetion BC
%[K,Q]=applyConvTriangJR(indCV,beta,Tinf,K,Q,nodes,elem); %<--DO NOT USE IT!
indCV=[indCirc1',indCirc4'];
[K,Q]=applyConvTriang(indCV,beta1,TInf1,K,Q,nodes,elem);
indCV=indCirc3';
[K,Q]=applyConvTriang(indCV,beta2,TInf2,K,Q,nodes,elem);

% Essential B.C.
u=zeros(numNodes,1);
u(indB)=tempBottom;  %Temperature at the bottom boundary
u(indT)=tempTop;     %Temperature at the bottom boundary
Fm = F(freeNodes) - K(freeNodes,fixedNodes)*u(fixedNodes);

%Reduced system
Km = K(freeNodes,freeNodes);
Fm = Fm + Q(freeNodes);

%Compute the solution
um = Km\Fm;
u(freeNodes)= um;

%PostProcess: Compute secondary variables and plot results
Q = Kini*u - Fini;
titol='Temperature Distribution';
colorScale='jet';
plotContourSolution(nodes,elem,u,titol,colorScale);

%Fancy output
%tableSol=[(1:numNodes)',nodes,u,Q];
%fprintf('%8s%9s%15s%15s%14s\n','Num.Nod','X','Y','T','Q')
%fprintf('%5d%18.7e%15.7e%15.7e%15.7e\n',tableSol')

%Solutions
[Ke, Fe] = linearTriangElement(coeff,nodes,elem,elemB);
interpTempP = [0.25,0.25,0.5]*u(elem(elemD,:));

fOut=fopen(fileSols,'w');
clc
fprintf(fOut,'Prob. 3\n');
fprintf(fOut,'(a) Number of nodes in the boundary of the hole centered at (6,1): %d\n',...
    length(indCirc2));
fprintf(fOut,'(b) The value of K^{50}(2,3) is K^{50}(2,3) = %.4e\n',...
    Ke(2,3));
fprintf(fOut,'    Hint. K(1,1) = %.4e\n',Kini(1,1));
fprintf(fOut,'(c) The mean value of the temperature in the nodes, <u> = %.4e\n',...
    sum(u)/numNodes);
fprintf(fOut,'    Hint. The maximum value of the temperature at the nodes is\n');
fprintf(fOut,'          max u = %.4e\n',max(u));
fprintf(fOut,'(d) The temperature for the point p of the element %d with barycentric\n',...
    elemD);
fprintf(fOut,'    coordinates (0.25,0.25,0.5) is T = %.4e\n',interpTempP);


%Compute the temperature for the point p=[4,2].
q= [4,2];

for e=1:numElem
    vertexs= nodes(elem(e,:),:);
    [alphas,isInside] = baryCoord(vertexs,q);
    if (isInside > 0)
        pElem = e;
        numNodElem= elem(e,:);
        break;
    end
end

interpTempQ = alphas*u(numNodElem);
fprintf(fOut,'    Hint. T(4,2) = %.4e\n',interpTempQ);
fclose(fOut);

type(fileSols);

