function Ke=planeLinkStiffMatrix(nod,elem,e,E,A)
%
% Matriu de rigidesa per un element Barra en 2-dim
% 
% (c) Numerical factory
%
x1=nod(elem(e,1),1);
y1=nod(elem(e,1),2);
x2=nod(elem(e,2),1);
y2=nod(elem(e,2),2);

x21=x2-x1;
y21=y2-y1;
Le=sqrt(x21*x21+y21*y21);
coef=((E(e)*A(e))/Le^3);
    Ke=[x21*x21, x21*y21, -x21*x21, -x21*y21;
        x21*y21, y21*y21, -x21*y21, -y21*y21;
        -x21*x21, -x21*y21, x21*x21, x21*y21;
        -x21*y21, -y21*y21, x21*y21, y21*y21];
Ke=coef*Ke;