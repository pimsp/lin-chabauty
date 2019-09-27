GEN HyperRandPt(GEN f, GEN T, GEN p, ulong e, GEN pe); // Get a random point on the curve modulo p^e
GEN RReval(GEN Ps, ulong n, ulong d, GEN T, GEN pe); // Get a basis for L(n*(infty_+ + infty_-)), where d = deg f
GEN HyperInit(GEN f, GEN p, ulong a, long e); // Construct the Jacobian
GEN HyperPicRand(GEN J,GEN f); // Get a random point on the Jacobian (by randomly constructing E)
GEN ordJ(GEN f, GEN p, ulong a); // Find the cardinality of the Jacobian corresponding to y^2 = f(x) in F_{p^a}
GEN HyperPicRandTors(GEN J, GEN f, GEN l, GEN C); // Get a random l torsion point? Not sure
GEN HyperPicEvalData(GEN J); // Unclear
GEN HyperPicEval(GEN J, GEN W, GEN U); // Unclear

GEN HyperPicInbedInfMin(GEN J, GEN P); // Represent a point P (for now assume to not be at infinity) in the Jacobian J as P - Inf_-
