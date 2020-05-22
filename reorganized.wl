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
"mids"->Import[infile,{"Datasets", "distribution_costheta_mid(lab)"}] (*Cos\[Theta] bin midpoints*)
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




(*Stability Matrix Construction*)

Com[A_,B_]:=Module[{a=A,b=B},
Return[A.B-B.A]
];


buildHamiltonians[infile_,ri_,\[Omega]_,Ve_,hi_]:=Module[{n,\[Theta],name11,name12,name21,name22,\[Rho],\[Rho]b,A,Ab,Hm,Hvac,\[Mu],\[Mu]b,Hsi,H,Hb,\[Delta]H,\[Delta]Hb,data,nudensity,nubardensity},(

data=ImportData[infile,ri];

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

Return[{H,Hb,\[Rho],\[Rho]b,A,Ab,\[Delta]H,\[Delta]Hb}]
)
];


getEquations[infile_,ri_,\[Omega]_,Ve_,hi_,k_]:=Module[{n,\[Theta],eqn,eqnb,hs,data},( 
hs=buildHamiltonians[infile,ri,\[Omega],Ve,hi];
data=ImportData[infile,ri];
n=Length[data["mids"]];
\[Theta]=ArcCos[data["mids"]];
(*This could be replaced with a mapthread or with associations, but one step at time*)
With[{H=hs[[1]],Hb=hs[[2]],\[Rho]=hs[[3]],\[Rho]b=hs[[4]],A=hs[[5]],Ab=hs[[6]],\[Delta]H=hs[[7]],\[Delta]Hb=hs[[8]]}, 

Do[
eqn[j]=Com[H[j],A[j]][[1,2]]+ Com[\[Delta]H[j],\[Rho][j]][[1,2]]+(k Cos[\[Theta][[j]]] A[j][[1,2]]);
eqnb[j]=-Com[Hb[j],Ab[j]][[1,2]]- Com[\[Delta]Hb[j],\[Rho]b[j]][[1,2]]+(k Cos[\[Theta][[j]]] Ab[j][[1,2]]);
,{j,1,n}];

Return[{eqn,eqnb,A,Ab}]
]; (*Close with*)
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

stabilityMatrix[infile_,ri_,\[Omega]_,Ve_,hi_,k_]:=Module[{S1,S2,S3,S4,S,hs,ea,n,data}, 
ea=getEquations[infile,ri,\[Omega],Ve,hi,k];
data=ImportData[infile,ri];
n=Length[data["mids"]];

With[{eqn=ea[[1]],eqnb=ea[[2]],A=ea[[3]],Ab=ea[[4]]},

S1=ParallelTable[Coefficient[eqn[l],A[m][[1,2]]],{l,1,n},{m,1,n}];
S2=ParallelTable[Coefficient[eqn[l],Ab[m][[1,2]]],{l,1,n},{m,1,n}];
S3=ParallelTable[Coefficient[eqnb[l],A[m][[1,2]]],{l,1,n},{m,1,n}];
S4=ParallelTable[Coefficient[eqnb[l],Ab[m][[1,2]]],{l,1,n},{m,1,n}];
S=ArrayFlatten[{{S1,S2},{S3,S4}}]/.rules[n];
Return[S];
]
];


(*
build2bMatrix[wt_,kx_]:=Block[{S2ba,S2b},
S2ba=stabilityMatrix[2,{0.,Pi},wt,kx,0.,-1.];
S2b={{S2ba[[1,1]],S2ba[[1,4]]},{S2ba[[4,1]],S2ba[[4,4]]}};
Return[S2b]
];
*)

\[Omega]Eev[En_]:=\[Omega]Eev[En]=(\[CapitalDelta]m12sq)/(2 En) everg;


stabilityMatrix[infile,ri,0.1,0,-1,0]//MatrixForm (*yay*)


evscale[A_,kx0_]:=Block[{\[Epsilon],As,kx0s},
\[Epsilon]=$MachineEpsilon;
As=Expand[(A/\[Epsilon])/.kx->\[Epsilon] kxs];
kx0s=kx0/\[Epsilon];
as=\[Epsilon] Eigenvalues[As/.kxs->kx0s];
(*bs=\[Epsilon] Eigenvalues[As]/.kxs\[Rule]kx0s;*)
Return[as]
];
evscale2[A_,kx0_]:=Block[{\[Epsilon],As,kx0s},
\[Epsilon]=$MachineEpsilon/2;
As=Expand[(A/\[Epsilon])/.kx->\[Epsilon] kxs];
kx0s=kx0/\[Epsilon];
as=\[Epsilon] Eigenvalues[As/.kxs->kx0s];
(*bs=\[Epsilon] Eigenvalues[As]/.kxs\[Rule]kx0s;*)
Return[as]
];

SCalclotsScale[rsrt_,rend_,testE_,ind_,ktest_,hi_]:=(
(*Saves output for particular set of inputs to avoid repeated computation*)
(*Piece of mind that absuing global variables isn't causing hold over values.*)
Clear[k,bangles,angles,nv,nvb,evalsl];
(*Create list of eigenvalues of S. Instability freqs*)
evalsl=Reap[ (* open reap *)
Do[ (* open do over radial bins*)
Do[ (* Open do over angluar bins*)
(* Populate matrix of neutrino numbers in each angular bin scaled by munits.  Used in SI part of Hamiltonian*)
nv[i]= Abs[dolots[ind,r,i](muss[[ind,i+1]]-muss[[ind,i]])];
nvb[i]= Abs[adolots[ind,r,i](muss[[ind,i+1]]-muss[[ind,i]])];
,{i,1,Length[mids[[ind]]]}]; (* Close angular bin do *)
(*Calls stabiltyMatrix for the current inputs.  Assigns a global value for the S matrix*)
Sow[Sort[Im[evscale[N[stabilityMatrix[Length[mids[[ind]]],ArcCos[mids[[ind]]],\[Omega]Eev[testE],kx,munits/mp *Yes[[ind,r]] *matters[[ind,r]],hi]],ktest]],Greater],eigen]; (* Calculates the numerical eignevalues of S, takes the imaginary parts, ranks them according to absolute value, then sows them with the tag "eigen"*)
Sow[Hsi[Length[mids[[ind]]]][[1,1]]/.rrules,pots] ;
(*Caculates the self interaction along the radial direction for the current nv/nvb and sows with tag pots*)
,{r,rsrt,rend}](*Close r do loop*)
,{eigen,pots}];(*Close reap and reap the tags eigen and pots*)
Return[evalsl]; 
);

SCalclotsScaleMat[rsrt_,rend_,testE_,ind_,ktest_,hi_]:=(
(*Saves output for particular set of inputs to avoid repeated computation*)
(*Piece of mind that absuing global variables isn't causing hold over values.*)
Clear[k,bangles,angles,nv,nvb,evalsl];
(*Create list of eigenvalues of S. Instability freqs*)
evalsl=Reap[ (* open reap *)
Do[ (* open do over radial bins*)
Do[ (* Open do over angluar bins*)
(* Populate matrix of neutrino numbers in each angular bin scaled by munits.  Used in SI part of Hamiltonian*)
nv[i]= Abs[dolots[ind,r,i](muss[[ind,i+1]]-muss[[ind,i]])];
nvb[i]= Abs[adolots[ind,r,i](muss[[ind,i+1]]-muss[[ind,i]])];
,{i,1,Length[mids[[ind]]]}]; (* Close angular bin do *)
(*Calls stabiltyMatrix for the current inputs.  Assigns a global value for the S matrix*)
Sow[Sort[Im[evscale2[N[stabilityMatrix[Length[mids[[ind]]],ArcCos[mids[[ind]]],\[Omega]Eev[testE],kx,munits/mp *Yes[[ind,r]] *matters[[ind,r]],hi]],ktest]],Greater],eigen]; (* Calculates the numerical eignevalues of S, takes the imaginary parts, ranks them according to absolute value, then sows them with the tag "eigen"*)
Sow[stabilityMatrix[Length[mids[[ind]]],ArcCos[mids[[ind]]],\[Omega]Eev[testE],ktest,munits/mp *Yes[[ind,r]] *matters[[ind,r]],hi],pots] ;
(*Caculates the self interaction along the radial direction for the current nv/nvb and sows with tag pots*)
,{r,rsrt,rend}](*Close r do loop*)
,{eigen,pots}];(*Close reap and reap the tags eigen and pots*)
Return[evalsl]; 
);
