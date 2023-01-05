clearvars
close all

%f(x) = f1*delta(x-pi/2)+f2*delta(x-3*pi/2)
f1 = 8.0;  %We shall pass this as a parameter
f2 = 8.0;

syms x;

a1 = 3.0;
a0 = 30.0;
u0 = 3.0;
du1 = -3.0;

p = 7/4;

Psi11 = @(x) 1 - x;
Psi12 = @(x) x;

Psi21 = @(x) 2*(x-3/2).*(x-2);
Psi22 = @(x) -4*(x-1).*(x-2);
Psi23 = @(x) 2*(x-1).*(x-3/2);

Psi1 = @(x) [Psi11(x);Psi12(x)];
Psi2 = @(x) [Psi21(x);Psi22(x);Psi23(x)];

K = a1*int(diff(Psi1(x),x)*diff(Psi1(x),x)',x,0,1);
K11 = double(K);
K = a0*int(Psi1(x)*Psi1(x)',x,0,1);
K10 = double(K);

K = a1*int(diff(Psi2(x),x)*diff(Psi2(x),x)',x,1,2);
K21 = double(K);
K = a0*int(Psi2(x)*Psi2(x)',x,1,2);
K20 = double(K);

K1 = K11 + K10;
K2 = K21 + K20;

F = f1*Psi1(1/2);
F1 = double(F);

F = f2*Psi2(3/2);
F2 = double(F);

clear x K F;

K = zeros(4);
F = zeros(4,1);
Q = zeros(4,1);
u = zeros(4,1);

K(1:2,1:2) = K1;
K(2:4,2:4) = K(2:4,2:4) + K2;

F(1:2) = F1;
F(2:4) = F(2:4) + F2;

%Bondary conditions
fixedNods = 1;
freeNods = setdiff(1:4,fixedNods);

%Natural B.C.
Q(4) = a1*du1;

%Esential B.C.
u(1) = u0;

%Reduced system
Qm = F(freeNods) + Q(freeNods) - K(freeNods,fixedNods)*u(fixedNods);
Km = K(freeNods,freeNods);

um = Km\Qm;

u(freeNods) = um;

Q = K*u-F;

%Interpolated value of u at x = p
interpU = u(2:4)'*Psi2(p);

%Solutions
fprintf('Prob. 1\n')
fprintf('(a)       F(2) = %.4e\n',F(2))
fprintf('  Hint.   F(3) = %.4e\n',F(3))
fprintf('(b)     K(2,2) = %.4e\n',K(2,2))
fprintf('  Hint. K(3,2) = %.4e\n',K(3,2))
fprintf('(c)       Q(1) = %.4e\n',Q(1))
fprintf('(d)       U(%f) %c %.4e\n',p,char(8776),interpU)