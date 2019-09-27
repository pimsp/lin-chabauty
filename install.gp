WeiRed(f,h)=
{
 my(F=f+(h/2)^2);
 2^poldegree(F)*subst(F,x,x/2);
}

install("matkerpadic","GGGL",,"./liblinalg.so");
install("mateqnpadic","GGGL",,"./liblinalg.so");
install("VecSmallCompl","GU","VecSmallCompl","./liblinalg.so");
install("matF","GGGU","matF","./liblinalg.so");
install("FqM_mul","GGGG","FqM_mul","./liblinalg.so");

install("NAF","G","NAF","./libexp.so");
install("AddChain","GL","AddChain","./libexp.so");

install("PlaneZeta","GU","PlaneZeta","./libplanereg.so");

install("RRInit2","GUUGGGGUL","PicInit","./librr.so");
install("RREvalInit","GG","PicEvalInit","./librr.so");
install("RREval","GGG","PicEval","./librr.so");
install("OnePol","GGGG","OnePol","./librr.so");
install("AllPols","GGGLG","AllPols","./librr.so");
install("FnsEvalAt","GGGGGLG",PicFnsEvalAt,"./librr.so");
install("FnEvalAt","GGGGGLG",PicFnEvalAt,"./librr.so");

install("RReval","GUUGG","HyperPicEval","./libhyper.so");

install("DivMul","GGGG","DivMul","./libpic.so");
install("PicRed","GU","PicRed","./libpic.so");
install("PicChord","GGGL","PicChord","./libpic.so");
install("PicAdd","GGG","PicAdd","./libpic.so");
install("PicSub","GGG","PicSub","./libpic.so");
install("PicNeg","GG","PicNeg","./libpic.so");
install("PicMul","GGGL","PicMul","./libpic.so");
install("PicFrob","GG","PicFrob","./libpic.so");
install("PicFrobPoly","GGG","PicFrobPoly","./libpic.so");
install("PicEq","lGGG","PicEq","./libpic.so");
install("PicIsZero","lGG","PicIsZero","./libpic.so");
install("PicChart","GGU","PicChart","./libpic.so");
install("PicRand0","G","PicRand","./libpic.so");

install("JgetW0","G","JgetW0","./libpic.so");
install("Jgetg","lG","Jgetg","./libpic.so");
install("JgetT","G","JgetT","./libpic.so");
install("JgetZ","G","JgetZ","./libpic.so");
install("Jgetpe","G","Jgetpe","./libpic.so");
install("Jgetp","G","Jgetp","./libpic.so");
install("Jgetf","G","Jgetf","./libpic.so");
install("Jgete","lG","Jgete","./libpic.so");
install("JgetFrobMat","G","JgetFrobMat","./libpic.so");
install("Frob","GGGG","Frob","./libpic.so");

install("PicLift_worker","UUUGGGGGGG","PicLift_worker","./liblift.so");
install("PicLiftTors_worker","GGGGGGUUUUGGLGLGGGU","PicLiftTors_worker","./liblift.so");
install("PicLiftTors","GGLG","PicLiftTors","./liblift.so");

install("PicNorm","GGG","PicNorm","./libfreyruck.so");
install("PicFreyRuckMulti","GGGGGG","PicFreyRuckMulti","./libfreyruck.so");
install("PicTorsRels","GGGU","PicTorsRels","./libfreyruck.so");
install("Fq_zeta_l","GGG","Fq_zeta_l","./libfreyruck.so");
install("Fq_mu_l_log","GGGGG","Fq_mu_l_log","./libfreyruck.so");

