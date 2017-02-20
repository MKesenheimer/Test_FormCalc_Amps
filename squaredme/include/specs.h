#if 0
* specs.h
* process specifications
* generated by FormCalc 9.4 (7 Jun 2016) on 15-Dec-2016 16:37
#endif

#undef FCVERSION
#define FCVERSION "FormCalc 9.4 (7 Jun 2016)"

#undef PROCNAME
#define PROCNAME "Fb31iF31iV1V5i"

#undef SQUAREDME_FUNC
#define SQUAREDME_FUNC SquaredME

#undef KIN
#define KIN "2to2.F"

#undef IDENTICALFACTOR
#define IDENTICALFACTOR 1

#undef Compose
#define Compose(f,c,A,B,C,D) c(c(c(f(1,A,1),f(2,B,1)),f(3,C,2)),f(4,D,2))

#undef Generic
#define Generic(f,c) Compose(f,c,FERMION,FERMION,PHOTON,PHOTON)

#undef Anti
#define Anti(f,c) Compose(f,c,-1,1,1,1)

#undef Mass
#define Mass(f,c) Compose(f,c,MU,MU,0,0)

#undef Charge
#define Charge(f,c) Compose(f,c,-2/3.D0,2/3.D0,0,0)

#undef ColorCharge
#define ColorCharge(f,c) Compose(f,c,-2/Sqrt(3.D0),2/Sqrt(3.D0),0,Sqrt(3.D0))

