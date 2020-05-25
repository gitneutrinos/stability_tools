(* ::Package:: *)

(* ::Subsection:: *)
(*Import data from a SINGLE file and a SINGLE radial bin - will be used inside of lists later, so everything is specific to one radius*)


(*Dialog prompt to fix in/out paths.  Can be replaced with direct path input later*)
id=ChoiceDialog[
"It's dangerous to go alone! Who are you and your computer companion?",{"Sam & Cuchulainn"-> "Sam", "Sherwood & Ganon"->"Sherwood"}];
(*Sets the paths based on the machine pick.  The paths are global*)
 Which[id=="Sam",
 inpath="G:\\My Drive\\Physics\\Neutrino Oscillation Research\\Fast Conversions\\lotsadata.tar\\lotsadata\\lotsadata\\";
 outpath="G:\\My Drive\\Physics\\Neutrino Oscillation Research\\Fast Conversions\\stability_data\\";
 ParallelEvaluate[$ProcessorCount];
 ,
 id=="Sherwood",
 inpath="/mnt/data/SamFlynn/lotsadata/";
 outpath="/mnt/data/SamFlynn/stability_data/";
 ];
 (*Want to implements toggle grid here to pick between the mass models and time slices.  It's easy, but definitely not a priority*)
filename = "112Msun_100ms_DO";

infile = inpath<>filename<>".h5";
outfolder = outpath<>filename;
(*Note: this used to contain a variable called "out_path". Be careful with the underscores; in Mathematica they are either function inputs or patterns.  If the font turn green, that's why.*)

ri = 1; (*Radial bin test case*)


ImportData[infile_,ri_]:=Association[
"lotsodo"->Import[infile,{"Datasets","distribution(erg|ccm,lab)"}][[ri]] (*distribution functions*),
"matters"->Import[infile,{"Datasets","rho(g|ccm,com)"}][[ri]], (*densities*)
"Yes"->Import[infile,{"Datasets","Ye"}][[ri]], (*electron fractions *)
"freqs"->Import[infile,{"Datasets", "distribution_frequency_grid(Hz,lab)"}], (*freq grid in hz*)
"freqmid"->Import[infile,{"Datasets", "distribution_frequency_mid(Hz,lab)"}], (*freq mid points*)
"muss"->Import[infile,{"Datasets", "distribution_costheta_grid(lab)"}], (*Cos\[Theta] grid*)
"mids"->Import[infile,{"Datasets", "distribution_costheta_mid(lab)"}], (*Cos\[Theta] bin midpoints*)
"radius"-> Import[infile,{"Datasets","r(cm)"}][[ri]]
]


(* ::Text:: *)
(*Set the following just so as not to immediately break Sam's code on merge*)
(*Thanks!*)


(*
ImportedData = ImportData[infile,ri];
lotsodo=ImportedData[["lotsodo"]];
matters=ImportedData[["matters"]];
Yes=ImportedData[["Yes"]];
freqs=ImportedData[["freqs"]];
freqmid=ImportedData[["freqmid"]];
muss=ImportedData[["muss"]];
mids=ImportedData[["mids"]];
*)


(*Constants*)
c=2.99792458 10^10; (* cm/s*)
h=6.6260755 10^-27; (*erg s*)
hbar = h/2 Pi; (*erg s*)
Gf=1.1663787 10^-5; (*GeV^-2*)
everg=1.60218 10^-12; (* convert eV to ergs*)
Geverg = everg*10^9; (* convert GeV to ergs *)
ergev=1.0/everg; (*convert ergs to eV*)
ergmev=ergev/10^6; (*convert erg to MeV*)
mp=1.6726219 10^-24; (*Proton mass in g*)
munits=Sqrt[2] Gf/Geverg^2 (hbar c)^3; (*Sqrt[2] Gf in erg cm^3*)
\[CapitalDelta]m12sq=7.59 10^-5;
\[Omega]EGev[En_]:=\[Omega]Eev[En]=(\[CapitalDelta]m12sq)/(2 En) everg/10^9;




(*Stability Matrix Construction*)

Com[A_,B_]:=Module[{a=A,b=B},
Return[A.B-B.A]
];


buildHamiltonians[infile_,ri_,testE_,hi_]:=Module[{n,\[Theta],name11,name12,name21,name22,\[Rho],\[Rho]b,A,Ab,Hm,Hvac,\[Mu],\[Mu]b,Hsi,H,Hb,\[Delta]H,\[Delta]Hb,data,nudensity,nubardensity,Ve,\[Omega],HsiRad},(

data=ImportData[infile,ri];

Ve=munits/mp *data["Yes"]  *data["matters"];
\[Omega]=\[Omega]EGev[testE];

name11="ee";
name12="ex";
name21="xe";
name22="xx";

n=Length[data["mids"]];
\[Theta]=ArcCos[data["mids"]];

Do[
\[Rho][i]={{ToExpression[StringJoin["\[Rho]",name11,ToString[i]]],ToExpression[StringJoin["\[Rho]",name12,ToString[i]]]},{ToExpression[StringJoin["\[Rho]",name21,ToString[i]]],ToExpression[StringJoin["\[Rho]",name22,ToString[i]]]}};
\[Rho]b[i]={{ToExpression[StringJoin["\[Rho]b",name11,ToString[i]]],ToExpression[StringJoin["\[Rho]b",name12,ToString[i]]]},{ToExpression[StringJoin["\[Rho]b",name21,ToString[i]]],ToExpression[StringJoin["\[Rho]b",name22,ToString[i]]]}};
A[i]={{0,ToExpression[StringJoin["A",name12,ToString[i]]]},{ToExpression[StringJoin["A",name21,ToString[i]]],0}};
Ab[i]={{0,ToExpression[StringJoin["Ab",name12,ToString[i]]]},{ToExpression[StringJoin["Ab",name21,ToString[i]]],0}};
,{i,1,n}];


nudensity[dt_]:= Sum[Sum[data["lotsodo"][[1,f,dt,dp]]/ (h (data["freqmid"][[f]]) (Abs[data["muss"][[dt+1]]-data["muss"][[dt]]])),{f,1,Length[data["freqs"]]-1}],{dp,1,2}];
nubardensity[dt_]:= Sum[Sum[data["lotsodo"][[2,f,dt,dp]]/ (h (data["freqmid"][[f]]) (Abs[data["muss"][[dt+1]]-data["muss"][[dt]]])),{f,1,Length[data["freqs"]]-1}],{dp,1,2}];

Hm={{Ve,0.},{0.,0.}};
Hvac=hi{{-\[Omega]/2,0.},{0.,\[Omega]/2}};
\[Mu]=munits Table[nudensity[i](data["muss"][[i+1]]-data["muss"][[i]]),{i,1,n},{j,1,n}];
\[Mu]b=munits Table[nubardensity[i](data["muss"][[i+1]]-data["muss"][[i]]),{i,1,n},{j,1,n}];


Do[
Hsi[j]=Sum[\[Mu][[k,j]]\[Rho][k](2Pi)(1-Cos[\[Theta][[j]]]Cos[\[Theta][[k]]]),{k,1,n}]+Sum[-\[Mu]b[[k,j]]\[Rho]b[k](2Pi)(1-Cos[\[Theta][[j]]]Cos[\[Theta][[k]]]),{k,1,n}];
,{j,1,n}];
Do[
H[i]=Hvac+Hm+Hsi[i];
Hb[i]=Hvac-Hm-Hsi[i];
\[Delta]H[i]=Sum[((D[H[i][[1,2]],\[Rho][j][[1,2]]]) A[j])+((D[H[i][[1,2]],\[Rho]b[j][[1,2]]])Ab[j]),{j,1,n}];
\[Delta]Hb[i]=Sum[((D[Hb[i][[1,2]],\[Rho][j][[1,2]]])A[j])+((D[Hb[i][[1,2]],\[Rho]b[j][[1,2]]])Ab[j]),{j,1,n}];
,{i,1,n}];

HsiRad=Hsi[n][[1,1]];

Return[{H,Hb,\[Rho],\[Rho]b,A,Ab,\[Delta]H,\[Delta]Hb,HsiRad}]
)
];


getEquations[infile_,ri_,testE_,hi_,k_]:=Module[{n,\[Theta],eqn,eqnb,hs,data},( 
hs=buildHamiltonians[infile,ri,testE,hi];
data=ImportData[infile,ri];
n=Length[data["mids"]];
\[Theta]=ArcCos[data["mids"]];
(*This could be replaced with a mapthread or with associations, but one step at time*)
With[{H=hs[[1]],Hb=hs[[2]],\[Rho]=hs[[3]],\[Rho]b=hs[[4]],A=hs[[5]],Ab=hs[[6]],\[Delta]H=hs[[7]],\[Delta]Hb=hs[[8]]}, 

Do[
eqn[j]=Com[H[j],A[j]][[1,2]]+ Com[\[Delta]H[j],\[Rho][j]][[1,2]]+(k Cos[\[Theta][[j]]] A[j][[1,2]]);
eqnb[j]=-Com[Hb[j],Ab[j]][[1,2]]- Com[\[Delta]Hb[j],\[Rho]b[j]][[1,2]]+(k Cos[\[Theta][[j]]] Ab[j][[1,2]]);
,{j,1,n}];

Return[{eqn,eqnb,A,Ab,hs[[9]]}]
](*Close with*)
)
];




rules[n_]:=Module[{r1,r2,r3,r4,rb1,rb2,rb3,rb4,rrules},( 
r1=Table[ToExpression["\[Rho]ee"<>ToString[i]]-> 1.,{i,1,n}];
r2=Table[ToExpression["\[Rho]ex"<>ToString[i]] -> 0.,{i,1,n}];
r3=Table[ToExpression["\[Rho]xe"<>ToString[i]]-> 0.,{i,1,n}];
r4=Table[ToExpression["\[Rho]xx"<>ToString[i]]-> 0.,{i,1,n}];

rb1=Table[ToExpression["\[Rho]bee"<>ToString[i]]-> 1.,{i,1,n}];
rb2=Table[ToExpression["\[Rho]bex"<>ToString[i]]-> 0.,{i,1,n}];
rb3=Table[ToExpression["\[Rho]bxe"<>ToString[i]]-> 0.,{i,1,n}];
rb4=Table[ToExpression["\[Rho]bxx"<>ToString[i]]-> 0.,{i,1,n}];

rrules=Flatten[{r1,r2,r3,r4,rb1,rb2,rb3,rb4}];
Return[rrules];
);
];

stabilityMatrix[infile_,ri_,testE_,hi_,k_]:=Module[{S1,S2,S3,S4,S,hs,ea,n,data,HsiRad}, 
ea=getEquations[infile,ri,testE,hi,k];
data=ImportData[infile,ri];
n=Length[data["mids"]];

With[{eqn=ea[[1]],eqnb=ea[[2]],A=ea[[3]],Ab=ea[[4]]},

S1=ParallelTable[Coefficient[eqn[l],A[m][[1,2]]],{l,1,n},{m,1,n}];
S2=ParallelTable[Coefficient[eqn[l],Ab[m][[1,2]]],{l,1,n},{m,1,n}];
S3=ParallelTable[Coefficient[eqnb[l],A[m][[1,2]]],{l,1,n},{m,1,n}];
S4=ParallelTable[Coefficient[eqnb[l],Ab[m][[1,2]]],{l,1,n},{m,1,n}];
S=ArrayFlatten[{{S1,S2},{S3,S4}}]/.rules[n];
HsiRad=ea[[5]]/.rules[n];
Return[{S,HsiRad}];
]
];


(*
build2bMatrix[wt_,kx_]:=Block[{S2ba,S2b},
S2ba=stabilityMatrix[2,{0.,Pi},wt,kx,0.,-1.];
S2b={{S2ba[[1,1]],S2ba[[1,4]]},{S2ba[[4,1]],S2ba[[4,4]]}};
Return[S2b]
];
*)




evscale[infile_,ri_,testE_,hi_,ktest_]:=Module[{\[Epsilon],A,As,kx0s,as,kx,kxs},
\[Epsilon]=$MachineEpsilon/2;
A=stabilityMatrix[infile,ri,testE,hi,kx][[1]];
As=Expand[(A/\[Epsilon])/.kx->\[Epsilon] kxs];
kx0s=ktest/\[Epsilon];
as=\[Epsilon] Eigenvalues[N[As]/.kxs->kx0s];
(*bs=\[Epsilon] Eigenvalues[As]/.kxs\[Rule]kx0s;*)
Return[as]
];



evscale[infile,ri,20,-1,10^-17]


SCalcScale[infile_,ri_,testE_,ktest_,hi_]:=SCalclotsScale[infile,ri,testE,ktest,hi]=Module[{evalsl,pot,mat,kt,ktarget},(
(*stabilityMatrix[infile_,ri_,testE,hi_,k_]*)
(*Create list of eigenvalues of S. Instability freqs*)
evalsl=Sort[Im[evscale[infile,ri,testE,hi,ktest]],Greater]; 
pot=stabilityMatrix[infile,ri,testE,hi,ktest][[2]];
mat=stabilityMatrix[infile,ri,testE,hi,ktest][[1]];
ktarget=kvar/.Solve[stabilityMatrix[infile,ri,testE,hi,kvar][[1]][[1,1]]==0,kvar][[1]];
Return[{evalsl,pot,mat,ktarget}]; 
);
];


SCalcScale[infile,ri,20,-1,10][[4]]



buildkGrid[infile_,ri_,testE_,hi_,nstep_]:=Module[{ktarget,kgrid,kvar,fspace},
fSpace[min_,max_,steps_,f_: Log]:=InverseFunction[ConditionalExpression[f[#],min<#<max]&]/@Range[f@min,f@max,(f@max-f@min)/(steps-1)];
ktarget=SCalcScale[infile,ri,testE,hi,0.][[4]];
kgrid=fSpace[ktarget*10^-2,ktarget*10^2,nstep];
Return[kgrid];
];


buildkGrid[infile,ri,20,-1,50]


(*Exports *)
kAdapt[infile_,rstr_,rend_,testE_,hi_,nstep_]:= kAdapt[infile,rstr,rend,testE,hi]=Module[{kl,evs1r,evout},
evout=Reap[Do[
kl=buildkGrid[infile,rx,testE,hi,nstep];
Sow[
evs1r=Reap[
Do[
Sow[{ImportData[infile,rx]["radius"],kl[[kx]],SCalcScale[infile,rx,testE,kl[[kx]],hi][[1]][[1]],SCalcScale[infile,rx,testE,kl[[kx]],hi][[2]]}]; 
,{kx,1,Length[kl]}
] (*close do over ktargets*)
][[2,1]] (*close reap over k*)
](*Close sow over r *)
,{rx,rstr,rend}
] (*close do over r*)
][[2,1]] (*Close reap over r*)
]; (*close module*)


(*Runs ok up to this point, stillr requires cross checks.  Seems slower.*)


kAdapt[infile,80,81,20,-1,50]
