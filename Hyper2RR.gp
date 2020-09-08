read("install.gp")

HyperRR(n,g,u,v)=
{
	concat(vector(n+1-poldegree(u,'x),i,'x^(i-1)*u),(y-v)*vector(n-g,i,x^(i-1)));
}

Hyper2RR(f0,P1,P2)= /* y^2=f0(x). P1,P2 rat pts, not conjugate by hyper invol. */
{
	my([x1,y1]=P1,[x2,y2]=P2,f,g,d0,L,LL,L1,L2);
	f = subst(f0,variable(f0),'x);
	d0 = poldegree(f);
	if(poldegree(f)%2,
		while(polcoef(f,0)==0 || x1==0 || x2==0,
			f=subst(f,x,x+1);
			x1-=1;
			x2-=1
		);
		f='x*polrecip(f);
		d0 += 1;
		x1 = 1/x1;
		x2 = 1/x2;
		y1 /= x1^(d0/2);
		y2 /= x2^(d0/2);
	);
	print(f,[x1,y1],[x2,y2]);
	g=(d0-2)/2;
	L=HyperRR(g+1,g,1,0);
	LL=HyperRR(2*g+2,g,1,0);
	if(g%2,
		L1=HyperRR(3*(g+1)/2,g,'x-x1,y1);
		L2=HyperRR(3*(g+1)/2,g,'x-x2,y2);
	,
		L1=HyperRR(3*g/2+2,g,('x-x1)*('x-x2),(y2-y1)/(x2-x1)*'x+(y1*x2-y2*x1)/(x2-x1));
		L2=HyperRR(3*g/2+1,g,1,0)
	);
	[y^2-f,g,d0,L,LL,L1,L2];
}

Hyper2RR_2(f0)= /* y^2=f0(x). */
{
	my(f,g,d0,L,LL);
	f = subst(f0,variable(f0),'x);
	d0 = poldegree(f);
	if(poldegree(f)%2,
		while(polcoef(f,0)==0,
			f=subst(f,x,x+1);
		);
		f='x*polrecip(f);
		d0 += 1;
	);
	print(f);
	g=(d0-2)/2;
	L=HyperRR(g+1,g,1,0);
	LL=HyperRR(2*g+2,g,1,0);
	[y^2-f,g,d0,L,LL];
}

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ From now on all code is by Pim Spelier                                     \\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

\\ Given a point P = [x,y] \in Z^2 with f(x,y) = 0 mod p where f(x,y) = y^2-f(x), find all p lifts of P
\\ to (Z/p^2 Z)^2. The first point is P itself if f(x,y) = 0 mod p^2
\\ (Assumes that f defines a smooth curve over p, and p > 2, and y^2 = f(x) in general)
HyperDeformPoint(P,p,f)={
	\\ Want solutions to (y + ap)^2 = f(x + bp)
	\\ Or equivalently y^2 - f(x) + 2yap = bf'(x)p
	
	my(b,y2fxp,dfx,dfxi,pts,y2i);
	
	y2fxp = subst(subst(f,x,P[1]),y,P[2])/p;
	dfx = -subst(deriv(f,x),x,P[1]);
	
	\\ If y = 0 mod p, then the lifts are exactly given by b = f'(x)^{-1}*(y^2-f(x))/p
	if(P[2]%p == 0,
		\\ dfx is invertible as the curve is smooth
		dfxi = lift(Mod(dfx,p)^(-1));
		eps = (dfxi*y2fxp)%p;
		pts = [[(P[1] + eps*p)%p^2,(P[2]+a*p)%p^2] | a <- [0..(p-1)]];
		return(pts);
	);
	
	y2i = lift(Mod(2*P[2],p)^(-1)); \\ The inverse of 2y modulo p, as an element of Z
	
	pts = [[(P[1] + b*p)%(p^2),(P[2]+ (b*dfx-y2fxp)*y2i*p)%(p^2)] | b <- [0..(p-1)]];
	return(pts);
}

\\ If flag = -1, compute P - Inf_- inside the jacobian
\\ If flag = 1, compute P - Inf_+
HyperPicInbedInf(J,P,flag)=
{
	/* Need to find E >= 0 sucht that P - Inf_- = E - D_0, and compute L(2D_0 - E)
	We know for F = (X^(g+1) + Y)/Z^(g+1), we have div F = -2D - infty_+ + 2g+1 random points
	So we can start by calculating L(D_0 + Inf_- P), and then multiply by F to get
	L(2D_0 - E) for E = P + 2g+1 random points.
	
	For this, we first calculate L(D_0 + Inf_-). We evaluate at all points of Z 
	and also at P. Note that D_0 = (g+1)*(Inf_- + Inf_+)*/
	
	my(T,p,e,pe,g);
	
	T = JgetT(J);
	p = Jgetp(J);
	e = Jgete(J);
	pe = Jgetpe(J);
	g = Jgetg(J);
	Z = JgetZ(J);
	Zp = concat(Z,[P]);
	z = length(Z);
	
	LD0 = HyperPicEval(Zp,(g+1),poldegree(Jgetf(J)),T,pe);
	
	/* We know that L(D_0 + Inf_-) has dimension one higher than L(D_0), so we
	just need to add a single element with poles D_0 + Inf_-; for example,
	X*(X^(g+1)-Y)/Z^(g+2) works.
	*/
	vars = Vecsmall([0,1]); \\ for x, y resp.
	
	M = PicFnsEvalAt([x*(x^(g+1)+flag*y)], Zp, vars, T, p, e, pe);
	
	LD0I = concat(LD0,M);
	LD0IP = FqM_mul(LD0I,matkerpadic(Mat(LD0I[z+1,]),T,p,e),T,pe);
	\\ Drop the last row, the one corresponding to P, to represent the functions
	LD0IP = LD0IP[1..z,];
	\\ in LD0IP normally
	
	F = PicFnsEvalAt([x^(g+1)-flag*y],Z,vars,T,p,e,pe)[,1];
	return(DivMul(F,LD0IP,T,pe));
}

Cinb(x,t) = Pol([x],t); 
CinbP(x,t) = [Cinb(x[1],t),Cinb(x[2],t)]; 


\\ This is previous code; it is replaced by HyperInbedJLocalChart(J,D)
\\ Given D a divisor on J (J define modulo p^2) that is zero modulo p
\\ and divisors D_i 
\\ return a vector mu_1,...  with D = sum_i mu_i * D_i
\\ or -1 if no such vector exists
\\ Works by brute force
HyperInbedJQzero(J,D,Ds)={
	my(T,Z,p,pe,f,t,n,P,Pmu,Pmt,a,Dso);
	T = JgetT(J);
	Z = JgetZ(J);
	p = Jgetp(J);
	pe = Jgetpe(J);
	f = Jgetf(J);

	t = variable(Z[1][1]);
	\\ To go from a point in Z^2 to a point on the curve, we need to view 
	\\ the coordinates as variables in t.
	
	n = length(Ds);
	
	\\ Base case
	if(n == 0,
		if(PicIsZero(J,D),return ([]),return (-1))
	);
	
	Dso = Ds[1..(n-1)];
	\\ Counterintuitively, PicSub(J,A,B) = B-A
  	for (i = 1, p,
  		a = HyperInbedJQzero(J,D,Dso);
  		if(a != -1, return(concat(a,[i-1])));
  		D = PicSub(J,Ds[n],D)
  	);
  	return(-1);
}

\\ There is an F_p linear map: J(Z/p^2 Z)_0 \to \F_p^k for some k >= g; evaluate D according to this map
\\ This works in polynomial time! (See the PicChart function for more information)
HyperInbedJLocalChart(J,D) = {
	return (((PicChart(J,D,1) - PicChart(J,JgetW0(J),1))~)/Jgetp(J))
}

\\ Given a point P in C(Z/p Z), return a divisor generating {Q-P} where
\\ Q ranges over the deformations of P
HyperDeformGenDiv(J,P) = {
	my(f,p,dfs,t);
	p = Jgetp(J); \\ the prime p
	f = Jgetf(J); \\ the polynomial defining the curve C
	t = variable(Z[1][1]);
	
	dfs = HyperDeformPoint(P,p,f);
	return(PicSub(J,HyperPicInbedInf(J,CinbP(dfs[1],t),-1),HyperPicInbedInf(J,CinbP(dfs[2],t),-1)));
} 

\\ Given bds = [D_1,...,D_r] the image of a basis of J(Q)_0 in the F_p vector space J(Z/p^2Z)_0
\\ And some points p_1,...,p_g in C(Z/p^2Z) with J(Z/p^Z)_0 = { sum_{i=1}^g p_i,mu - p_i } with p_i,mu
\\ deformations of p_i
\\ And finally a list of points candidates in C(Z/pZ), determine for each P in candidates whether there's a unique Fp-linear combination
\\ of D_1,...,D_r equal to P_mu - P_nu for some deformations P_mu,P_nu of P
HyperDeformLiftQ(J,bds,candidates) = {
	\\ Note that first of all 0 is equal to P - P, so there's always at least one solution
	\\ We have maps from F_p^r to F_p^k induced by the D_i's and a map from F_p to F_p^k given 
	\\ by mu |-> P_mu - P
	\\ There is exactly one solution iff the kernel of kappa: F_p^r + F_p to F_p^k is zero

	my(p,f,T,kappa);

	p = Jgetp(J); \\ the prime p
	f = Jgetf(J); \\ the polynomial defining the curve C
	T = JgetT(J);
	
	\\ We first construct the matrix corresponding to the map F_p^{r} to F_p^l; its columns are given
	\\ by HyperInbedJLocalChart on the elements of bds
	kappa = matconcat([HyperInbedJLocalChart(J,D)~ | D <- bds]);
	
	\\ Then the matrix corresponding to the map F_p^r + F_p -> F_p^g is given by
	\\ concatenating the column corresponding to HyperDeformGenDiv(J,P)
	print([matconcat([kappa,HyperInbedJLocalChart(J,HyperDeformGenDiv(J,P))~])| P <- candidates]);
	return([matsize(matkerpadic(matconcat([kappa,HyperInbedJLocalChart(J,HyperDeformGenDiv(J,P))~]),T,p,1))[2] == 0 | P <- candidates]);
	\\bds = concat(bds,[HyperDeformGenDiv(J,P)]);
	\\return([kappa,matsize(matkermod(kappa,p))[2] == 0])
}
