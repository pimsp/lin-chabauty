#include "linalg.h"
#include "exp.h"
#include "pic.h"

GEN Jgetf(GEN J) {return gel(J,1);}
long Jgetg(GEN J) {return itos(gel(J,2));}
long Jgetd0(GEN J) {return itos(gel(J,3));}
GEN JgetT(GEN J) {return gel(J,4);}
GEN Jgetp(GEN J) {return gel(J,5);}
long Jgete(GEN J) {return itos(gel(J,6));}
GEN Jgetpe(GEN J) {return gel(J,7);}
GEN JgetFrobMat(GEN J) {return gel(J,8);}
GEN JgetV(GEN J) {return gel(J,9);}
GEN JgetKV(GEN J) {return gel(J,10);}
GEN JgetW0(GEN J) {return gel(J,11);}
GEN JgetZ(GEN J) {return gel(J,12);}
GEN JgetFrobCyc(GEN J) {return gel(J,13);}
GEN JgetV3(GEN J) {return gel(J,14);}
GEN JgetKV3(GEN J) {return gel(J,15);}

void JgetTpe(GEN J, GEN* T, GEN* pe, GEN* p, long* e)
{
	*T = gel(J,4);
	*p = gel(J,5);
	*e = itos(gel(J,6));
	*pe = gel(J,7);
}


GEN PicRed(GEN J, ulong e)
{
	GEN Je,p,pe;
	if(Jgete(J)<e) pari_err(e_MISC,"Cannot perform this reduction");
	Je = cgetg(lgJ+1,t_VEC);
	gel(Je,1) = gcopy(Jgetf(J));
	gel(Je,2) = stoi(Jgetg(J));
	gel(Je,3) = stoi(Jgetd0(J));
	gel(Je,4) = gcopy(JgetT(J));
	gel(Je,5) = p = gcopy(Jgetp(J));
	gel(Je,6) = utoi(e);
	gel(Je,7) = pe = powiu(p,e);
	gel(Je,8) = FpM_red(JgetFrobMat(J),pe);
	gel(Je,9) = FpXM_red(JgetV(J),pe);
	gel(Je,10) = FpXM_red(JgetKV(J),pe);
	gel(Je,11) = FpXM_red(JgetW0(J),pe);
	gel(Je,12) = FpXT_red(JgetZ(J),pe);
	gel(Je,13) = gcopy(JgetFrobCyc(J));
	gel(Je,14) = FpXM_red(JgetV3(J),pe);
	gel(Je,15) = FpXM_red(JgetKV3(J),pe);
	return Je;
}

GEN DivMul(GEN f, GEN W, GEN T, GEN pe)
{
	ulong nW,nZ,i,j;
	GEN fW,col;
	nW = lg(W);
	nZ = lg(f);
	fW = cgetg(nW,t_MAT);
	for(j=1;j<nW;j++)
	{
		col = cgetg(nZ,t_COL);
		for(i=1;i<nZ;i++)
		{
			gel(col,i) = Fq_mul(gel(f,i),gcoeff(W,i,j),T,pe);
		}
		gel(fW,j) = col;
	}
	return fW;
}

GEN DivAdd0(GEN WA, GEN WB, ulong d, GEN T, GEN p, long e, GEN pe, ulong excess)
{
  pari_sp av,av1,av0=avma;
  unsigned long nZ,nA,nB,j,P,r;
  GEN WAB,s,t,st;
  nZ = lg(gel(WA,1));
	nA = lg(WA)-1;
	nB = lg(WB)-1;
  WAB = cgetg(d+excess+1,t_MAT);
  while(1)
  {
    av1 = avma;
    for(j=1;j<=d+excess;j++)
    {
      av = avma;
			s = gel(WA,1+random_Fl(nA)); /* random basis fn in WA */
			t = gel(WB,1+random_Fl(nB)); /* random basis fn in WB */
      st = cgetg(nZ,t_COL); /* Product */
      for(P=1;P<nZ;P++)
      {
        gel(st,P) = Fq_mul(gel(s,P),gel(t,P),T,pe);
      }
      gel(WAB,j) = gerepileupto(av,st);
    } 
    r = FqM_rank(WAB,T,p); /* TODO faut-il reduire WAB d'abord? */
    if(r==d)
    {
      if(excess)
      {
        WAB = gerepileupto(av0,matimagepadic(WAB,T,p,e));
      }
      return WAB;
    }
    printf("add0(%lu/%lu)",r,d);
    avma = av1;
  }
}

GEN DivAdd1(GEN WA, GEN WB, ulong d, GEN T, GEN p, long e, GEN pe, ulong excess)
{
  pari_sp av,av1,av0=avma;
  unsigned long nZ,j,P,r;
  GEN WAB,s,t,st;
  nZ = lg(gel(WA,1));
  WAB = cgetg(d+excess+1,t_MAT);
  while(1)
  {
    av1 = avma;
    for(j=1;j<=d+excess;j++)
    {
      av = avma;
      s = RandVec_1(WA,pe); /* random fn in WA */
      t = RandVec_1(WB,pe); /* random fn in WB */
      st = cgetg(nZ,t_COL); /* Product */
      for(P=1;P<nZ;P++)
      {
        gel(st,P) = Fq_mul(gel(s,P),gel(t,P),T,pe);
      }
      gel(WAB,j) = gerepileupto(av,st);
    }
    r = FqM_rank(WAB,T,p); /* TODO faut-il reduire WAB d'abord? */
    if(r==d)
    {
      if(excess)
      {
        WAB = gerepileupto(av0,matimagepadic(WAB,T,p,e));
      }
      return WAB;
    }
    printf("add0(%lu/%lu)",r,d);
    avma = av1;
  }
}

GEN DivAdd(GEN WA, GEN WB, ulong d, GEN T, GEN p, long e, GEN pe, ulong excess)
{
	pari_sp av,av1,av0=avma;
	unsigned long nZ,j,P,r;
	GEN WAB,s,t,st;
	nZ = lg(gel(WA,1));
	WAB = cgetg(d+excess+1,t_MAT);
	while(1)
	{
		av1 = avma;
		for(j=1;j<=d+excess;j++)
		{ 
			av = avma;
			s = RandVec_padic(WA,T,p,pe); /* random fn in WA */
			t = RandVec_padic(WB,T,p,pe); /* random fn in WB */
			st = cgetg(nZ,t_COL); /* Product */
			for(P=1;P<nZ;P++)
			{
				gel(st,P) = Fq_mul(gel(s,P),gel(t,P),T,pe);
			}
			gel(WAB,j) = gerepileupto(av,st);
		}
		r = FqM_rank(WAB,T,p); /* TODO faut-il reduire WAB d'abord? */
		if(r==d)
		{
			if(excess)
			{
				WAB = gerepileupto(av0,matimagepadic(WAB,T,p,e));
			}
			return WAB;
		}
		printf("add(%lu/%lu)",r,d);
		avma = av1;
	}
}

GEN DivSub(GEN WA, GEN WB, GEN KV, ulong d, GEN T, GEN p, long e, GEN pe, ulong nIGS)
{
	pari_sp av1,av = avma;
	unsigned long nZ,P,nE,E,nV,nB,n,r;
	GEN KB,K,col,s,res;
	nZ = lg(KV);
	nV = lg(gel(KV,1))-1;
	KB = mateqnpadic(WB,T,p,e);
	nB = lg(gel(KB,1))-1;
	/* Prepare a mat K of size a v stack of KV + nIGS copies of KB */
	/* and copy KV at the top */
	nE = nV + nIGS*nB;
	K = cgetg(nZ,t_MAT);
	for(P=1;P<nZ;P++)
	{
		col = cgetg(nE+1,t_COL);
		for(E=1;E<=nV;E++)
		{
			gel(col,E) = gcoeff(KV,E,P);
		}
		gel(K,P) = col;
	}
	av1 = avma;
	while(1)
	{
		/* nIGS times, take rand s in WA, and stack s.KB down K */
		for(n=1;n<=nIGS;n++)
		{
			s = RandVec_padic(WA,T,p,pe);
			for(E=1;E<=nB;E++)
			{
				for(P=1;P<nZ;P++)
				{
					gcoeff(K,nV+(n-1)*nB+E,P) = Fq_mul(gel(s,P),gcoeff(KB,E,P),T,pe);
				}
			}
		}
		res = FqM_ker(K,T,p); /* TODO faut-il reduire K d'abord? */
		/* TODO take rand subset of eqns */
		/* TODO write fn for that */
		r = lg(res)-1;
		if(r==d)
		{
			return gerepileupto(av,e==1?res:matkerpadic(K,T,p,e));
		}
		printf("sub(%lu/%lu)",r,d);
		avma = av1;
	}
}

GEN PicChord(GEN J, GEN WA, GEN WB, long flag)
{
	pari_sp av = avma;
	GEN WAWB,WAB,s,sV,WC,res;
	GEN V,KV,KV3,W0,T,p,pe;
	long g,d0,e;

	V = JgetV(J);
	KV = JgetKV(J);
	KV3 = JgetKV3(J);
	W0 = JgetW0(J);
	JgetTpe(J,&T,&pe,&p,&e);
	g = Jgetg(J);
	d0 = Jgetd0(J);

	WAWB = DivAdd(WA,WB,2*d0+1-g,T,p,e,pe,0);
	WAB = DivSub(W0,WAWB,KV3,d0+1-g,T,p,e,pe,2);
	/* TODO can free some memory here */
	if(flag & 1) s = RandVec_padic(WAB,T,p,pe);
	else s = gel(WAB,1);
	sV = DivMul(s,V,T,pe);
	WC = DivSub(WAB,sV,KV,d0+1-g,T,p,e,pe,2);

	if(flag & 2)
	{
		res = cgetg(3,t_VEC);
		gel(res,1) = gcopy(WC);
		gel(res,2) = gcopy(s);
		return gerepileupto(av,res);
	}
	else
	{
		return gerepileupto(av,WC);
	}
}

GEN PicAdd(GEN J, GEN WA, GEN WB)
{
	pari_sp av = avma;
	GEN W;
	W = PicChord(J,WA,WB,0);
	W = PicChord(J,W,JgetW0(J),0);
	return gerepileupto(av,W);
}

GEN PicSub(GEN J, GEN WA, GEN WB)
{
	pari_sp av = avma;
	GEN W;
	W = PicChord(J,WB,JgetW0(J),0);
	//printf("First chord done\n");
	W = PicChord(J,W,WA,0);
	return gerepileupto(av,W);
}

GEN PicNeg(GEN J, GEN W) { return PicChord(J,W,JgetW0(J),0); }

GEN PicMul(GEN J, GEN W, GEN n, long flag)
{
	pari_sp av = avma;
	GEN W0,C,Wlist,WA,WB;
	ulong nC,i;
	long a,b;

	W0 = JgetW0(J);
	if(gequal0(n)) return W0;
	if(gequal(n,gen_1)) return gcopy(W);
	C = AddChain(n,flag&2);
	nC = lg(C);
	if(flag&2)
		pari_printf("   PicMul : Mul by %Ps in %lu steps\n",n,nC-2);
	else
		pari_printf("   PicMul : Mul by ±%Ps in %lu steps\n",n,nC-2);
	Wlist = cgetg(nC,t_VEC);
	gel(Wlist,1) = W;
	for(i=2;i<nC;i++)
	{
		a = gmael(C,i,2)[1];
		WA = a?gel(Wlist,a):W0;
		b = gmael(C,i,2)[2];
		WB = b?gel(Wlist,b):W0;
		gel(Wlist,i) = PicChord(J,WA,WB,(i==nC-1)&&(flag&1));
	}
	return gerepileupto(av,gel(Wlist,nC-1));
}

GEN ZpXQ_FrobMat(GEN T, GEN p, long e, GEN pe)
{
	pari_sp av = avma;
	GEN F,M,col,Fj;
	long v = gvar(T),d = degpol(T),i,j;
	F = ZpX_Frobenius(T,p,e);
	M = cgetg(d+1,t_MAT);
	col = cgetg(d+1,t_COL);
	gel(col,1) = gen_1;
	for(i=2;i<=d;i++) gel(col,i) = gen_0;
	gel(M,1) = col;
	if(d==1) return gerepileupto(av,M);
	col = cgetg(d+1,t_COL);
	for(i=0;i<d;i++) gel(col,i+1) = polcoef(F,i,v);
	gel(M,2) = col;
	Fj = F;
	for(j=2;j<d;j++)
	{
		Fj = Fq_mul(Fj,F,T,pe);
		col = cgetg(d+1,t_COL);
		for(i=0;i<d;i++) gel(col,i+1) = polcoef(Fj,i,v);
		gel(M,j+1) = col;
	}
	return gerepilecopy(av,M);
}

GEN Frob(GEN x, GEN FrobMat, GEN T, GEN pe)
{
	pari_sp av = avma;
	GEN cx,cy,y;
	long var = gvar(T), d = degpol(T), i;
	cx = cgetg(t_COL,d+1);
	for(i=0;i<d;i++) gel(cx,i+1) = polcoef(x,i,var);
	cy = FpM_FpC_mul(FrobMat,cx,pe);
	y = cgetg(d+2,t_POL);
	setvarn(y,var);
	for(i=1;i<=d;i++) gel(y,i+1) = gel(cy,i);
	y = normalizepol(y);
	return gerepilecopy(av,y);
}

GEN PicFrob(GEN J, GEN W)
{
	GEN W2,T,pe,FrobMat,FrobCyc;
	ulong o,i,j,k,c,nW,nZ,nCyc;

	T = JgetT(J);
	pe = Jgetpe(J);
	FrobMat = JgetFrobMat(J);
	FrobCyc = JgetFrobCyc(J);
	nW = lg(W);
	nZ = lg(JgetZ(J));
	nCyc = lg(FrobCyc);

	W2 = cgetg(nW,t_MAT);
	for(j=1;j<nW;j++)
	{
		gel(W2,j) = cgetg(nZ,t_COL);
	}

	i = 0;
	for(o=1;o<nCyc;o++)
	{
		c = FrobCyc[o];
		for(k=1;k<c;k++)
		{
			for(j=1;j<nW;j++)
			{
				gcoeff(W2,i+k+1,j) = Frob(gcoeff(W,i+k,j),FrobMat,T,pe);
			}
		}
		for(j=1;j<nW;j++)
		{
			gcoeff(W2,i+1,j) = Frob(gcoeff(W,i+c,j),FrobMat,T,pe);
		}
		i += c;
	}
	return W2;
}

GEN PicFrobPoly(GEN J, GEN W, GEN F)
{
	pari_sp av = avma;
	ulong d,i;
	GEN n,FW,res;

	d = degree(F);
	FW = W;
	n = truecoeff(F,0);
	if(d&1L) n = negi(n);
	res = PicMul(J,W,n,2);
	for(i=1;i<=d;i++)
	{
		FW = PicFrob(J,FW);
		n = truecoeff(F,i);
		if((d+1-i)&1L) n = negi(n);
		res = PicChord(J,res,PicMul(J,FW,n,2),0);
	}
	return gerepileupto(av,res);
}

long PicEq(GEN J, GEN WA, GEN WB)
{
	pari_sp av = avma;
	long e,r;
	GEN s,sWB,KsWB,K,KV,col,T,p,pe;
	ulong P,i,j,nZ,nW,nKV,nKsB,nK;

	JgetTpe(J,&T,&pe,&p,&e);
	KV = JgetKV(J);

	s = gel(WA,1);
	nZ = lg(s)-1;
	nW = lg(WA)-1;
	nKV = lg(gel(KV,1))-1;
	nKsB = nZ-nW;
	nK = nKV+nW*nKsB;

	sWB = cgetg(nW+1,t_MAT);
	for(j=1;j<=nW;j++)
	{
		col = cgetg(nZ+1,t_COL);
		for(i=1;i<=nZ;i++)
		{
			gel(col,i) = Fq_mul(gel(s,i),gcoeff(WB,i,j),T,pe);
		}
		gel(sWB,j) = col;
	}

	KsWB = mateqnpadic(sWB,T,p,e);

	K = cgetg(nZ+1,t_MAT);
	for(j=1;j<=nZ;j++)
	{
		gel(K,j) = cgetg(nK+1,t_COL);
	}


	for(j=1;j<=nW;j++)
	{
		for(i=1;i<=nKsB;i++)
		{
			for(P=1;P<=nZ;P++)
			{
				gcoeff(K,(j-1)*nKsB+i,P) = Fq_mul(gcoeff(WA,P,j),gcoeff(KsWB,i,P),T,pe);
			}
		}
	}
	for(i=1;i<=nKV;i++)
	{
		for(P=1;P<=nZ;P++)
		{
			gcoeff(K,nW*nKsB+i,P) = gcoeff(KV,i,P);
		}
	}

	r = lg(matkerpadic(K,T,p,e))-1;

	avma = av;
	return r;
}

long PicIsZero(GEN J, GEN W)
{
	return PicEq(J,W,JgetW0(J));
}

GEN PicChart(GEN J, GEN W, ulong P0) /* /!\ Not Galois-equivariant ! */
{
	pari_sp av = avma;
	ulong d0,g,n1,n2,nV,nZ,nW;
	ulong j,P;
	long e;
	GEN V,KV,T,p,pe;
	GEN K,col,s,sV,U,res;

	g = Jgetg(J);
	d0 = Jgetd0(J);
	n1 = d0-g;
	n2 = d0-g;
	V = JgetV(J);
	KV = JgetKV(J);
	nV = lg(V)-1;
	nZ = lg(gel(V,1))-1;
	nW = lg(W)-1;
	JgetTpe(J,&T,&pe,&p,&e);

	K = cgetg(nW+1,t_MAT);
	for(j=1;j<=nW;j++)
	{
		col = cgetg(n1+1,t_COL);
		for(P=1;P<=n1;P++) gel(col,P) = gcoeff(W,P+P0,j);
		gel(K,j) = col;
	}
	K = matkerpadic(K,T,p,e);
	if(lg(K)!=2)
	{
		pari_printf("Genericity 1 failed in PicChart\n");
		avma = av;
		return NULL;
	}
	s = FqM_FqC_mul(W,gel(K,1),T,pe);

	sV = cgetg(nV+1,t_MAT);
	for(j=1;j<=nV;j++)
	{
		col = cgetg(nZ+1,t_COL);
		for(P=1;P<=nZ;P++) gel(col,P) = Fq_mul(gel(s,P),gcoeff(V,P,j),T,pe);
		gel(sV,j) = col;
	}
	U = DivSub(W,sV,KV,d0+1-g,T,p,e,pe,2);
	K = cgetg(d0+2-g,t_MAT);
  for(j=1;j<=d0+1-g;j++)
  { 
    col = cgetg(n2+1,t_COL);
    for(P=1;P<=n2;P++) gel(col,P) = gcoeff(U,P0+n1+P,j);
    gel(K,j) = col;
  }
	K = matkerpadic(K,T,p,e);
	if(lg(K)!=2) 
  {
    pari_printf("Genericity 2 failed in PicChart\n");
    avma = av;
    return NULL;
  }
	s = FqM_FqC_mul(U,gel(K,1),T,pe);
	res = cgetg(nZ-n1-n2,t_COL);
	for(j=1;j<=P0;j++) gel(res,j) = gel(s,j);
	for(j=P0+n1+n2+1;j<=nZ;j++) gel(res,j-n1-n2) = gel(s,j);
	return gerepilecopy(av,res);
}

GEN rand_subset(ulong n, ulong r)
{
	pari_sp av;
	GEN X,S;
	ulong m,i;
	S = cgetg(r+1,t_VECSMALL);
	av = avma;
	X = cgetg(n+1,t_VECSMALL);
	for(i=1;i<=n;i++) X[i] = 1;
	m = 0;
	while(m<r)
	{
		i = random_Fl(n)+1;
		if(X[i])
		{
			X[i] = 0;
			m++;
			S[m] = i;
		}
	}
	avma = av;
	return S;
}

GEN PicRand0(GEN J)
{
	pari_sp av = avma;
	ulong d0,nZ,nV;
	ulong i,j;
	long e;
	GEN T,p,pe,V;
	GEN S,col,K;

	d0 = Jgetd0(J);
	JgetTpe(J,&T,&pe,&p,&e);
	V = JgetV(J);
	nV = lg(V);
	nZ = lg(gel(V,1));

	K = cgetg(nV,t_MAT);
	S = rand_subset(nZ-1,d0);
	for(j=1;j<nV;j++)
	{
		col = cgetg(d0+1,t_COL);
		for(i=1;i<=d0;i++)
		{
			gel(col,i) = gcoeff(V,S[i],j);
		}
		gel(K,j) = col;
	}
	K = matkerpadic_hint(K,T,p,e,pe,nV-1-d0);
	K = FqM_mul(V,K,T,pe);
	return gerepileupto(av,K);
}
