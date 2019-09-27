\\ Load the package modules
\\read("install.gp");
read("Hyper2RR.gp");
read("galrep.gp");

f = x^6+ 8*x^5+ 22*x^4+ 22*x^3+ 5*x^2+ 6*x+ 1; \\ Equation for the curve (must be smooth)
print("C : y² = ",f);
d = poldegree(f);
p = 5; \\ We choose to get it 5-adically
e = 2; \\ Target p-adic accuracy is O(5^2)
a = 3; \\ Degree of the unramified extension of Qp over which the relevant points are defined
Lp = hyperellcharpoly(Mod(f,p)); \\ Local L factor of the curve at p, needed to know the number of points on the Jacobian mod p
P1=[0,1];P2=[-3,1]; \\ We need TODO
[f,g,d0,L,LL,L1,L2]=Hyper2RR(f,P1,P2); \\ Precomputation of a generic-enough equation of the curve and of some Ri
J=PicInit(f,g,d0,L,LL,1,p,a,e);
print("Hyper2RR done");
U=PicEvalInit(J,[L1,L2]);
J1 = PicRed(J,1); \\ Reduction mod p

T = JgetT(J);
Z = JgetZ(J);
pe = Jgetpe(J);

b = [0,1];
c = [-3,1];
P = Z[1];
Dc = PicSub(J,HyperPicInbedInf(J,P,1),HyperPicInbedInf(J,P,-1)); \\ = inf_+ - inf-

Dc41 = PicMul(J,Dc,41,2); \\ 41*Dc. For some reason, the flag has to be 2

\\print(HyperInbedJQzero(J,Dc41,[b,c]));
