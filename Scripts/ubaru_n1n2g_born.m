(* ::Package:: *)

(*
generates the Fortran code for
p p -> weakino weakino jet in the MSSM
last modified Mar 2017
*)


Clear["Global`*"]
SetDirectory[NotebookDirectory[]];
<< FeynArts`
<< FeynArtsAdd`
<< FormCalc`
<< FormCalcAdd`
ClearProcess[]
<<"!rm *.frm"

time1 = SessionTime[]


(*You can now load the script with the command $ MathKernel -script ubaru_n1n2g_born.m "ubar" "u" "n1" "n2" "g"*)
Print[$CommandLine]
If[$CommandLine[[2]] === "-script",
	(p[1] = ToString[$CommandLine[[4]]];
	 p[2] = ToString[$CommandLine[[5]]];
	 p[3] = ToString[$CommandLine[[6]]];
	 p[4] = ToString[$CommandLine[[7]]];
	 p[5] = ToString[$CommandLine[[8]]];),
	(*Else*)
	(p[1] = "ubar";
	 p[2] = "u";
	 p[3] = "n1";
	 p[4] = "n2";
	 p[5] = "g";)
]

CalcProcess = p[1]<>p[2]<>"_"<>p[3]<>p[4]<>p[5];
name = CalcProcess;
Print[CalcProcess]

GluonLegs = {};
For[i=1, i<6, i++,
If[p[i] === "qu", P[i] = F[3],
If[p[i] === "qubar", P[i] = -F[3],
If[p[i] === "qd", P[i] = F[4],
If[p[i] === "qdbar", P[i] = -F[4],
If[p[i] === "nI", P[i] = F[11],
If[p[i] === "nJ", P[i] = F[11],
If[p[i] === "xI-", P[i] = F[12],
If[p[i] === "xI+", P[i] = -F[12],
If[p[i] === "xJ-", P[i] = F[12],
If[p[i] === "xJ+", P[i] = -F[12],
If[p[i] === "g", (GluonLegs = Join[GluonLegs, {i}]; P[i] = V[5]),

If[p[i] === "u", P[i] = F[3,{1}],
If[p[i] === "ubar", P[i] = -F[3,{1}],
If[p[i] === "c", P[i] = F[3,{2}],
If[p[i] === "cbar", P[i] = -F[3,{2}],
If[p[i] === "t", P[i] = F[3,{3}],
If[p[i] === "tbar", P[i] = -F[3,{3}],

If[p[i] === "d", P[i] = F[4,{1}],
If[p[i] === "dbar", P[i] = -F[4,{1}],
If[p[i] === "s", P[i] = F[4,{2}],
If[p[i] === "sbar", P[i] = -F[4,{2}],
If[p[i] === "b", P[i] = F[4,{3}],
If[p[i] === "bbar", P[i] = -F[4,{3}],

If[p[i] === "n1", P[i] = F[11,{1}],
If[p[i] === "n2", P[i] = F[11,{2}],
If[p[i] === "n3", P[i] = F[11,{3}],
If[p[i] === "n4", P[i] = F[11,{4}],

If[p[i] === "x1-", P[i] = F[12,{1}],
If[p[i] === "x1+", P[i] = -F[12,{1}],
If[p[i] === "x2-", P[i] = F[12,{2}],
If[p[i] === "x2+", P[i] = -F[12,{2}]
]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
]

process = {P[1], P[2]} -> {P[3], P[4], P[5]};
Print[process]


(*Neglect Masses (URL)*)
Neglect[MU] = Neglect[MU2] = 0;
Neglect[MC] = Neglect[MC2] = 0;
Neglect[MD] = Neglect[MD2] = 0;
Neglect[MS] = Neglect[MS2] = 0;
Neglect[MUC] = Neglect[MU2C] = 0;
Neglect[MCC] = Neglect[MC2C] = 0;
Neglect[MDC] = Neglect[MD2C] = 0;
Neglect[MSC] = Neglect[MS2C] = 0;
Neglect[_Mf] = Neglect[_Mf2] = 0;
Neglect[_MfC] = Neglect[_Mf2C] = 0;
(*Neglect[MB] = Neglect[MB2] = 0;
Neglect[MT] = Neglect[MT2] = 0;*)

(*DEBUG*)
(*Neglect[_MNeu] = Neglect[_MNeu2] = 0;*)
(*USf[i_,j_,k_,l_] := IndexDelta[i,j];
USfC[i_,j_,k_,l_] := IndexDelta[i,j];
ZNeu = IndexDelta;
ZNeuC = IndexDelta;
CKM = IndexDelta;
CKMC = IndexDelta;*)

(*particle widths*)
widths = {MGl2 -> MGl2 - I MGl WGl, MSf2[a] :> MSf2[a] - I MSf[a] WSf[a], MZ2 -> MZ2 - I MZ WZ, MW2 -> WZ2 - I WZ WW};
(*real widths*)
Scan[ (RealQ[#] = True)&, {WGl, _WSf, WW, WZ}];


(*Options*)
SetOptions[InsertFields, Model -> "MSSMCT",
           (*No Fermion-Higgs coupling*)
           Restrictions -> {NoLightFHCoupling},
           (*Exclude Top, Higgs, Neutrinos, massive Leptons, Sneutrinos, Sleptons*)
		   ExcludeParticles -> {S[1|2|3|4|5|6|11|12], F[1|2]},
		   (*no internal Weakinos*)
		   LastSelections -> {!F[11],!F[12]}];

SetOptions[Paint, PaintLevel -> {Classes}, ColumnsXRows -> {4, 5}, AutoEdit -> False];

(*Reduce tensor to scalar integrals and choose regularization scheme*)
(*D = dimensional regularization (default),*)
(*4 = constrained differential renormalization,*)
(*0 = keeps the whole amplitude D-dimensional*)
SetOptions[CalcFeynAmp,Dimension->D];

(*Save the Diagrams*)
$PaintSE = MkDir["Diagrams"];
DoPaint[diags_, type_, opt___] := Paint[diags, opt,
  DisplayFunction -> (Export[ToFileName[$PaintSE, name <> "_" <> type <> ".pdf"], #]&)]


Print["Born"]

tops = CreateTopologies[0, 2 -> 3];
ins = InsertFields[tops, process];
DoPaint[ins, "born"];

born = CalcFeynAmp[CreateFeynAmp[ins]];

(*insert the partice widths*)
born = born/.{Den[x_,y_]:>Den[x,y/.widths]}


(*Write files*)
amps = {born};
{born} = Abbreviate[amps, 6, Preprocess -> OnSize[100, Simplify, 500, DenCollect]];

col = ColourME[All, born];

abbr = OptimizeAbbr[Abbr[]]
subexpr = OptimizeAbbr[Subexpr[]]

dir = SetupCodeDir[name<>"_born", Drivers -> name <> "_drivers"];
WriteSquaredME[born, {}, col, abbr, subexpr, dir];


Print["time used: ", SessionTime[] - time1]
Exit[];
