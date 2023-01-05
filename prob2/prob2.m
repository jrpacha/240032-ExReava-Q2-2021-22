clearvars
close all

fileSols = 'prob2.txt';
fOut = fileSols;
fOut=fopen(fileSols,'w');

Area = 3500; %mm^2
Y = 150;     %150GPa = 150 kN/mm^2

%Goemetry

nodes = [0, 0;
    1, 1;
    2, 1;
    3, 0;
    1, 2;
    2, 2];

nodes = 1000*nodes;
       
elem=[
    1, 2;
    2, 3;
    3, 4;
    2, 5;
    3, 6;
    1, 5;
    5, 6;
    6, 4
    ];

      
numNod=size(nodes,1);
numElem=size(elem,1);
ndim=size(nodes,2);

numbering=1;
%plotElements(nodes,elem,numbering);
plotElementsOld(nodes,elem,numbering);       
grid on;

%% Real constants
A=Area*ones(1,numElem);
E=Y*ones(1,numElem);

%Assembly
u=zeros(ndim*numNod,1);
Q=zeros(ndim*numNod,1);
K=zeros(ndim*numNod);

for e=1:numElem
    Ke=planeLinkStiffMatrix(nodes,elem,e,E,A);
    rows=[ndim*elem(e,1)-1,ndim*elem(e,1),...
          ndim*elem(e,2)-1,ndim*elem(e,2)];
    cols=rows;
    K(rows,cols)=K(rows,cols)+Ke; %Assembly
end

%Loads
%node 6:
nod=6;
Q(ndim*nod-1)=200.0; %Fx = 200kN
Q(ndim*nod)=0.0;     %Fy = 0kN

%node 3:
nod=3;
Q(ndim*nod-1)=0.0;   %Fx = 0kN
Q(ndim*nod)=200.0;   %Fy = 200kN

%Boundary Conditions
fixedNods=[];
%node 1:
nod=1;
fixedNods=[fixedNods,ndim*nod-1];   %(u1_x=0);
fixedNods=[fixedNods,ndim*nod];     %(u1_y=0);
u(ndim*nod-1)=0.0;
u(ndim*nod)=0.0;

%node 5:
nod=5;
fixedNods=[fixedNods,ndim*nod-1];     %(u5_x=0);
fixedNods=[fixedNods,ndim*nod];       %(u5_y=0);
u(ndim*nod-1)=0.0;
u(ndim*nod)=0.0;

%node 4:
nod=4;
fixedNods=[fixedNods,ndim*nod];      %(u4_y=0);
u(ndim*nod)=0.0;

%Reduced system
freeNods=setdiff(1:ndim*numNod,fixedNods);
Qm=Q(freeNods,1)-K(freeNods,fixedNods)*u(fixedNods);
Km=K(freeNods,freeNods);

%Solve the reduced system
um=Km\Qm;
u(freeNods)=um;

UX=u(1:2:end); UY=u(2:2:end);
X = nodes(:,1); Y=nodes(:,2); XF = X + UX; YF = Y + UY;

%Reaction forces
Fr=K*u-Q;

L = sqrt((X(elem(:,2))-X(elem(:,1))).^2+(Y(elem(:,2))-Y(elem(:,1))).^2);
LF = sqrt((XF(elem(:,2))-XF(elem(:,1))).^2+(YF(elem(:,2))-YF(elem(:,1))).^2);

%Solutions
clc
fprintf(fOut,'Prob 2\n');
fprintf(fOut,'(a) Km(5,6) = %.4e\n',Km(5,6));
fprintf(fOut,'    Hint. UX(6) = %.6e\n',UX(6));
fprintf(fOut,'(b) The X-component of the right extreme of the 3rd bar\n');
fprintf(fOut,'    X%d = %.4e\n',elem(3,2),XF(elem(3,2)));
fprintf(fOut,'    Hint. The Y-component of the right extreme of the\n');
fprintf(fOut,'          7th bar is on Y%d = %.4e\n',elem(7,2),YF(elem(7,2)));
fprintf(fOut,'(c) Maximum of the length deformation of the bars\n');
fprintf(fOut,'    max |LF-L| = %.4e\n',max(abs(L-LF)));
fprintf(fOut,'    Hint. The final lenght of the 2nd bar is LF(2) = %.4e\n',...
    LF(2));

fclose(fOut);
type(fileSols);

%Post-process: plot deformed structure
esc=10;
plotDeformedTruss(nodes, elem, u, esc);