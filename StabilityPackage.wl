(* ::Package:: *)

BeginPackage["StabililtyPackage`"]

ClearAll["StabililtyPackage`", "StabililtyPackage`"]

ImportData::usage =
	"ImportData[file] reads in the data from file."
buildHamiltonians::usage=
	"[data,test Energy, Hierarchy=1/-1]. Return elements of the hamiltonians"
getEquations::usage=
	"[data,test energy, hierarchy=-1/1, k]. Return the equations of motions"
stabilityMatrix::usage=
	"[data,Equations of motion (getEquations output)]. Returns the stability matrix"
buildkGrid::usage=
"[data, number of k steps, (optional) manual target k, 
(optional) k lowwer bound multiple; Default \!\(\*SuperscriptBox[\(10\), \(-3\)]\),(optional) k upper bound Default: 10]. Creates the k grid"
kAdapt::usage=
	"[file, radius start, radius end, test energy, hierarchy, number of k steps, (optional) manual target k, 
(optional) k lowwer bound multiple; Default \!\(\*SuperscriptBox[\(10\), \(-3\)]\),(optional) k upper bound Default: 10]. Runs the full routine"
GDValue::usage=
	"returns the maximum possible value of the frequency"
GDdata::usage=
	"Finds the gershgorin limits over a range of radii"
SelectSingleRadius::usage=
	"get single radius data"
boxFit::usage=
"[infile,species=1,2]. Box fits a file."
evscale::usage=
"[ktest,S,kx] Get eigenvalues of matrix with tiny scale."
Bfactor::usage=
"[data] finds the asymmetry factor for the x flavor"
rules::usage=
"replacement rules for density matrix values"
siPotential::usage=
"get SI potential for k targets"
evecscale::usage = 
"finds eigenvectors while scaling the matrix to reasonable numbers"
ndensities::usage = 
"finds the neutrino number densities for a species, and returns the matrix of densities in an angular bin"
esysscale::usage=
"returns the real and imaginary parts of eigenvalues and associated eigenvectors"
getInitialGuess::usage=
"[file,species]generates Initial Guesses for the first radius from expansion expresssion"
eBoxFitSingleRadius::usage=
"fit a single radius using the box fit parameters"
Com::usage=
"test"
ellipseMoments::usage=
"calculate the moments from an ellipse (box fit) with parameters a, \[Beta], \[Chi]"
eBoxFitToMoments::usage=
"given 3 input moments and set of guess parameters, fits ellipse parameters to moments"
getMoments::usage=
"takes moments out of moment data for a given file, radius, species. Energy integrated"




ClearAll["StabililtyPackage`", "StabililtyPackage`"]


c=2.99792458 10^10; (* cm/s*)
h=6.6260755 10^-27; (*erg s*)
hbar =h/(2 Pi); (*erg s*)
Gf=1.1663787 10^-5; (*GeV^-2*)
everg=1.60218 10^-12; (* convert eV to ergs*)
Geverg = everg*10^9; (* convert GeV to ergs *)
Meverg = everg*10^6; (*convert Mev to erg*) 
ergev=1.0/everg; (*convert ergs to eV*)
ergmev=ergev/10^6; (*convert erg to MeV*)
mp=1.6726219 10^-24; (*Proton mass in g*) 
munits=Sqrt[2] (Gf/Geverg^2 )(hbar c)^3; (*Sqrt[2] Gf in erg cm^3*)
\[CapitalDelta]m12sq=(7.59 10^-5) everg^2;
\[Omega]EMev[En_]:=(\[CapitalDelta]m12sq)/(2 (En Meverg));

Begin["`Private`"]


Com[A_,B_]:=Module[{a=A,b=B},Return[A.B-B.A]];





	
(*Have xflavor option  set the x do to 0 here*)

ImportData[infile_]:=
Association[
"Endensity"->Import[infile,{"Data","distribution(erg|ccm,lab)"}] (*distribution functions*),
"matters"->Import[infile,{"Data","rho(g|ccm,com)"}], (*densities*)
"Yes"->Import[infile,{"Data","Ye"}], (*electron fractions *)
"freqs"->Import[infile,{"Data", "distribution_frequency_grid(Hz,lab)"}], (*freq grid in hz*)
"freqmid"->Import[infile,{"Data", "distribution_frequency_mid(Hz,lab)"}], (*freq mid points*)
"muss"->Import[infile,{"Data", "distribution_costheta_grid(lab)"}], (*Cos\[Theta] grid*)
"mids"->Import[infile,{"Data", "distribution_costheta_mid(lab)"}], (*Cos\[Theta] bin midpoints*)
"radius"-> Import[infile,{"Data","r(cm)"}],
"phis"-> Import[infile,{"Data","/distribution_phi_grid(radians,lab)"}] (*"phi bin edges"*)
];



SelectSingleRadius[data_,ri_]:=
Association[
"Endensity"->data["Endensity"][[ri]] (*distribution functions*),
"matters"->data["matters"][[ri]], (*densities*)
"Yes"->data["Yes"][[ri]], (*electron fractions *)
"freqs"->data["freqs"], (*freq grid in hz*)
"freqmid"->data["freqmid"], (*freq mid points*)
"muss"->data["muss"], (*Cos\[Theta] grid*)
"mids"->data["mids"],
"phis"-> data["phis"] (*"phi bin edges"*) (*Cos\[Theta] bin midpoints*)
];


Options[ndensities]={"xflavor"-> True};
ndensities[data_,OptionsPattern[]]:=Module[{n,nudensity,nubardensity,nuxdensity,nd,ndb,ndx},

n=Length[data["mids"]];

nudensity[dt_]:= Sum[Sum[data["Endensity"][[1,f,dt,dp]]/ (h (data["freqmid"][[f]])) ,{f,1,Length[data["freqs"]]-1}],{dp,1,Length[data["phis"]]-1}];
nubardensity[dt_]:= Sum[Sum[data["Endensity"][[2,f,dt,dp]]/ (h (data["freqmid"][[f]])),{f,1,Length[data["freqs"]]-1}],{dp,1,Length[data["phis"]]-1}];

If[OptionValue["xflavor"],
nuxdensity[dt_]:= Sum[Sum[0.25 data["Endensity"][[3,f,dt,dp]]/ (h (data["freqmid"][[f]]) ),{f,1,Length[data["freqs"]]-1}],{dp,1,Length[data["phis"]]-1}]
,
nuxdensity[dt_]:=0.
];

nd= Table[nudensity[i],{i,1,n},{j,1,n}];
ndb= Table[nubardensity[i],{i,1,n},{j,1,n}];
ndx=Table[nuxdensity[i],{i,1,n},{j,1,n}];


Return[{nd,ndb,ndx}]
];


Options[siPotential]={"xflavor"-> True};
siPotential[data_,OptionsPattern[]]:=Module[{tot,m},
m=munits ndensities[data,"xflavor"-> OptionValue["xflavor"]];
tot=(Tr[m[[1]]]+Tr[m[[2]]]+2 Tr[m[[3]]]);

Return[tot]
]

(*This function finds the asymetry factor between the electron (anti)neutrinos and the x (anti)neutrinos *)

Options[Bfactor]={"xflavor"-> True};
Bfactor[data_,OptionsPattern[]]:=Module[{B,Bb,m},
m=munits ndensities[data,"xflavor"-> OptionValue["xflavor"]];
B=Tr[m[[3]]]/Tr[m[[1]]];
Bb=Tr[m[[3]]]/Tr[m[[2]]];
Return[{B,Bb}];
];



(*This function creates the hamiltonians based on the data set, and returns them in a "general" way that is readable.  
Returns 9 arguments with index,
1,3,5,7 = H,\[Rho],A,\[Delta]H
2,4,6,8 = Hb,\[Rho]b,Ab,\[Delta]Hb
9=HsiRadial
*)

Options[buildHamiltonians]={"xflavor"-> True};
buildHamiltonians[data_,testE_,hi_,OptionsPattern[]]:=Module[{n,\[Theta],name11,name12,name21,name22,\[Rho],\[Rho]b,A,Ab,Hm,Hvac,\[Mu],\[Mu]b,\[Mu]x,m,Hsi,H,Hb,\[Delta]H,\[Delta]Hb,Ve,\[Omega]},(

Ve=munits/mp *data["Yes"]  *data["matters"];
\[Omega]=\[Omega]EMev[testE];

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


Hm={{Ve,0.},{0.,0.}};
Hvac=hi{{-\[Omega]/2,0.},{0.,\[Omega]/2}};
m=munits ndensities[data,"xflavor"-> OptionValue["xflavor"]];
\[Mu]=m[[1]];
\[Mu]b=m[[2]];
\[Mu]x=m[[3]];


Do[
Hsi[j]=Sum[\[Mu][[k,j]]\[Rho][k](1-Cos[\[Theta][[j]]]Cos[\[Theta][[k]]]),{k,1,n}]+Sum[-\[Mu]b[[k,j]]\[Rho]b[k](1-Cos[\[Theta][[j]]]Cos[\[Theta][[k]]]),{k,1,n}];
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




(*Calculates the equations of motion by computing the relvant commutators. 
 Returns 5 arguments with index,
 1,3= neutrino equations of motion, A
 2,4 = Antineutrino equations of motion, Ab
 5=HsiRadial
 *)

Options[getEquations]={"xflavor"-> True,"inverse"-> False};

getEquations[data_,testE_,hi_,k_,OptionsPattern[]]:=Module[{n,\[Theta],eqn,eqnb,hs},
hs=buildHamiltonians[data,testE,hi,"xflavor"-> OptionValue["xflavor"]];
n=Length[data["mids"]];
\[Theta]=ArcCos[data["mids"]];

(*This could be replaced with a mapthread or with associations, but one step at time*)
With[{H=hs[[1]],Hb=hs[[2]],\[Rho]=hs[[3]],\[Rho]b=hs[[4]],A=hs[[5]],Ab=hs[[6]],\[Delta]H=hs[[7]],\[Delta]Hb=hs[[8]]}, 

If[OptionValue["inverse"],
Do[
eqn[j]=-(Com[H[j],A[j]][[1,2]]+ Com[\[Delta]H[j],\[Rho][j]][[1,2]]+(-k Cos[\[Theta][[j]]] A[j][[1,2]]));
eqnb[j]=-(-Com[Hb[j],Ab[j]][[1,2]]- Com[\[Delta]Hb[j],\[Rho]b[j]][[1,2]]+(-k Cos[\[Theta][[j]]] Ab[j][[1,2]]));
,{j,1,n}];
,
Do[
eqn[j]=Com[H[j],A[j]][[1,2]]+ Com[\[Delta]H[j],\[Rho][j]][[1,2]]+(k Cos[\[Theta][[j]]] A[j][[1,2]]);
eqnb[j]=-Com[Hb[j],Ab[j]][[1,2]]- Com[\[Delta]Hb[j],\[Rho]b[j]][[1,2]]+(k Cos[\[Theta][[j]]] Ab[j][[1,2]]);
,{j,1,n}];
];

Return[{eqn,eqnb,A,Ab}]
](*Close with*)

];


(*Substitution rules to change "named" density matrix components with initial flavor state*)

Options[rules]={"xflavor"-> True}; (*this is not needed as b will come back 0. if no x, and no option is needed noir is the data*)
rules[data_,OptionsPattern[]]:=Module[{r1,r2,r3,r4,rb1,rb2,rb3,rb4,rrules,n,B,Bb},
n=Length[data["mids"]];
B=Bfactor[data,"xflavor"-> OptionValue["xflavor"]][[1]];
Bb=Bfactor[data,"xflavor"-> OptionValue["xflavor"]][[2]];


r1=Table[ToExpression["\[Rho]ee"<>ToString[i]]-> 1./(1.+B),{i,1,n}];
r2=Table[ToExpression["\[Rho]ex"<>ToString[i]] -> 0.,{i,1,n}];
r3=Table[ToExpression["\[Rho]xe"<>ToString[i]]-> 0.,{i,1,n}];
r4=Table[ToExpression["\[Rho]xx"<>ToString[i]]-> B/(1.+B),{i,1,n}];

rb1=Table[ToExpression["\[Rho]bee"<>ToString[i]]-> 1./(1.+Bb),{i,1,n}];
rb2=Table[ToExpression["\[Rho]bex"<>ToString[i]]-> 0.,{i,1,n}];
rb3=Table[ToExpression["\[Rho]bxe"<>ToString[i]]-> 0.,{i,1,n}];
rb4=Table[ToExpression["\[Rho]bxx"<>ToString[i]]->Bb/(1.+Bb) ,{i,1,n}];

rrules=Flatten[{r1,r2,r3,r4,rb1,rb2,rb3,rb4}];
Return[rrules];
];

(*Generates the stability matrix for file infile, radial bin ri, vacuum energy scale testE, mass ordering hi=1 (NO) or -1 (IO), and wavenumber k.
Returns 2 arguments,
1=Stability Matrix
2=HsiRadial
*)
Options[stabilityMatrix]={"xflavor"-> True};
stabilityMatrix[data_,ea_,OptionsPattern[]]:=Module[{S1,S2,S3,S4,S,n}, 

n=Length[data["mids"]];

With[{eqn=ea[[1]],eqnb=ea[[2]],A=ea[[3]],Ab=ea[[4]]},

S1=Table[Coefficient[eqn[l],A[m][[1,2]]],{l,1,n},{m,1,n}];
S2=Table[Coefficient[eqn[l],Ab[m][[1,2]]],{l,1,n},{m,1,n}];
S3=Table[Coefficient[eqnb[l],A[m][[1,2]]],{l,1,n},{m,1,n}];
S4=Table[Coefficient[eqnb[l],Ab[m][[1,2]]],{l,1,n},{m,1,n}];
S=ArrayFlatten[{{S1,S2},{S3,S4}}]/.rules[data,"xflavor"-> OptionValue["xflavor"]];
Return[S];
]
];




(*This scales the stability matrix up to a more managable scale based on the machine prescision, Solves for the eigenvalues, and then scales backs.  Returns a list of eigenvalues*)
Options[evscale]={"output"-> "Eigenvalues"};
evscale[ktest_,S_,kx_,OptionsPattern[]]:=Module[{\[Epsilon],As,kx0s,as,kxs,evals,evecs},
\[Epsilon]=$MachineEpsilon/2;
As=Expand[(S/\[Epsilon])/.kx->\[Epsilon] kxs];
kx0s=ktest/\[Epsilon];

If[OptionValue["output"]=="Eigenvalues",
as=\[Epsilon] Eigenvalues[N[As]/.kxs->kx0s];
];

If[OptionValue["output"]=="Eigenvectors",
as=\[Epsilon] Eigenvectors[N[As]/.kxs->kx0s];
];
If[OptionValue["outout"]=="Eigensystem",
as={evals,evecs}=\[Epsilon] Eigensystem[N[As]/.kxs-> kx0s];
];

Return[as]
];




(*Constructs a nstep sized log spaced k grid based on the target k associated with the infile at radial bin r.  Currently the limits are 2 orders of magnitude above and below the target value, ignoring negatives for the moment *)

Options[buildkGrid]={"ktarget"-> 0.,"krange"-> {10.^-3,10.},"xflavor"-> True};
buildkGrid[data_,nstep_,OptionsPattern[]]:=Module[{kgrid,fSpace,ktarget,kblow,kbhigh},
If[OptionValue["ktarget"]== 0.,
ktarget=siPotential[data,"xflavor"-> OptionValue["xflavor"]]
,
ktarget=OptionValue["ktarget"]
];
fSpace[min_,max_,steps_,f_: Log]:=InverseFunction[ConditionalExpression[f[#],min<#<max]&]/@Range[f@min,f@max,(f@max-f@min)/(steps-1)];
kblow=OptionValue["krange"][[1]];
kbhigh=OptionValue["krange"][[2]];
kgrid=Join[fSpace[ktarget*kblow,ktarget*kbhigh,nstep],-fSpace[ktarget*kblow,ktarget*kbhigh,nstep/2]];
Return[kgrid];
];



(*Run buildkGrid and SCalcScale for several radial bins.*)

Options[kAdapt]={"xflavor"-> True,"ktarget"-> 0.,"krange"-> {10.^-3,10.},"koutput"-> Eigenvalues,"inverse"-> False};
kAdapt[infile_,rstr_,rend_,testE_,hi_,nstep_,OptionsPattern[]]:= Module[{kl,evout,data,singleRadiusData,ea,kvar,eout,pot,S},

data=ImportData[infile];
evout=
Reap[
	Do[
		singleRadiusData = SelectSingleRadius[data,rx];
		ea=getEquations[singleRadiusData,testE,hi,kvar,"xflavor"-> OptionValue["xflavor"],"inverse"-> OptionValue["inverse"]];

		S=stabilityMatrix[singleRadiusData,ea,"xflavor"-> OptionValue["xflavor"]];

		kl=buildkGrid[singleRadiusData,nstep,"ktarget"-> OptionValue["ktarget"],"krange"-> OptionValue["krange"],"xflavor"-> OptionValue["xflavor"]];
		Do[
		
			eout=evscale[kl[[kx]],S,kvar,"output"->OptionValue["koutput"]];
			
			pot=siPotential[singleRadiusData,"xflavor"-> OptionValue["xflavor"]];
			Sow[{data["radius"][[rx]],kl[[kx]],eout,pot}]; 
		,{kx,1,Length[kl]}] (*close do over ktargets*)
	,{rx,rstr,rend}] (*close do over r*)
][[2,1]];

Return[evout] (*Close reap over r*)
]; (*close module*)





GDValue[tm_]:=Module[{cs,rrow,rcol,drow,dcol,mv,reg,\[Epsilon],tms},
\[Epsilon]=$MachineEpsilon;
tms=Expand[(tm/\[Epsilon])];
cs=Diagonal[tms];
rrow=Total[Abs[tms],{2}]-Abs[cs];
rcol=Total[Abs[tms],{1}]-Abs[cs];
drow=MapThread[Disk,{ReIm[cs],rrow}];
dcol=MapThread[Disk,{ReIm[cs],rcol}];
reg=RegionIntersection[Apply[RegionUnion,drow],Apply[RegionUnion,dcol]];
mv=\[Epsilon] Max[MeshCoordinates[BoundaryDiscretizeRegion[reg,{{-1,1},{0,1}}]][[All,2]]];
Return[mv]];


GDdata[infile_,rstr_,rend_,testE_,hi_,nstep_]:= Module[{kl,evout,data,singleRadiusData,ea,kvar,evals,pot,ktarget,S,gdout,out,ea2,S2},
data=ImportData[infile];
out=
Reap[
	Do[
		singleRadiusData = SelectSingleRadius[data,rx];
		ea=getEquations[singleRadiusData,testE,hi,kvar];
		S=stabilityMatrix[data,ea];
		ktarget=(S/.kvar->0)[[1,1]];
		kl=buildkGrid[ktarget,nstep];
		Do[
		ea2=getEquations[singleRadiusData,testE,hi,kx];
		S2=stabilityMatrix[data,ea];
			gdout=GDValue[S2]; (*this isn't updating k correctly*)
			Sow[{data["radius"][[rx]],kl[[kx]],gdout}]; 
		,{kx,1,Length[kl]}] (*close do over ktargets*)
	,{rx,rstr,rend}] (*close do over r*)
][[2,1]];
Return[out] (*Close reap over r*)
		
];




getMoments[file_,r_,species_]:= Module[{data,datasr,moments},
data=ImportData[file]//Quiet;
datasr=SelectSingleRadius[data,r];
moments={Sum[datasr["Endensity"][[species,E,1]],{E,1,Length[data["freqs"]]-1}],Sum[datasr["Endensity"][[species,E,2]],{E,1,Length[data["freqs"]]-1}],Sum[datasr["Endensity"][[species,E,3]],{E,1,Length[data["freqs"]]-1}]};
Return[moments];
];


getInitialGuess[m0_,m1_,m2_]:=Module[{foc1234,ag,\[Beta]g,\[Chi]g,cg,bg,arg\[Beta],arg\[Chi]},
foc1234[x_,y_,z_]:=(1/(4 Pi ) )(m0 + 3 z m1+(5/2 (3 (m0 -m2 )/2 x^2+3 (m0 -m2 )/2 y^2+3 m2 z^2-m0)));

ag=0.5 (foc1234[Sin[ArcCos[-1.]],0.,-1.]+foc1234[Sin[ArcCos[1.]],0.,1.]); (*semi-major axis guess*)

cg=Abs[foc1234[Sin[ArcCos[1.]],0.,1.]-ag]; (*horizontal shift from the center *)

bg=Sqrt[foc1234[Sin[ArcCos[0.]],0.,0.]^2/(1-(cg/ag)^2)]; (*semi-minor axis*)

Assert[bg<= ag]; (*Assert the semi-major axis is larger than the semi-minor*)
If[bg> ag, Assert[bg/ag<= 1.001]]; (* If bg>ag, assert that the difference is small*)
If[bg> ag && bg/ag<= 1.001, ag=ag+2(bg-ag)]; (* If bg>ag, and the difference is small, switches them in the transform*)

arg\[Beta]=((2 bg/ag)-1);
arg\[Chi]=((2 cg/ag)-1);

\[Beta]g=ArcTanh[arg\[Beta]]; (*box transform b= a/2( tanh[\[Beta]]+1 *)
\[Chi]g=ArcTanh[arg\[Chi]]; (*box transform c=a/2( tanh[\[Chi]]+1 *)

Assert[Between[arg\[Beta],{-1.01,1.01}]];
Assert[Between[arg\[Chi],{-1.01,1.01}]];

Print[{ag,bg,cg,arg\[Beta],arg\[Chi]}];
Return[{ag,\[Beta]g,\[Chi]g}]
];



ellipseMoments[af_,\[Beta]f_,\[Chi]f_]:=Module[{ebox,esbox},

ebox[a_,\[Beta]_,\[Chi]_,m_]:=(a (1+Tanh[\[Beta]]) (1/4 a^2 m (1+Tanh[\[Beta]]) (1+Tanh[\[Chi]])+a Sqrt[-a^2 (-1+m^2)+1/4 a^2 m^2 (1+Tanh[\[Beta]])^2
+1/4 a^2 (-1+m^2) (1+Tanh[\[Chi]])^2]))/(2 (a^2+m^2 (-a^2+1/4 a^2 (1+Tanh[\[Beta]])^2)));

esbox[mom_]:=2 Pi NIntegrate[m^mom ebox[af,\[Beta]f,\[Chi]f,m],{m,-1.,1.},MinRecursion-> 8,MaxRecursion->16];

Return[{esbox[0],esbox[1],esbox[2]}]

];


(*Given 3 moments, fit parameters a, \[Beta], and \[Chi] so ellipseMoments match*)
eBoxFitToMoments[m0_,m1_,m2_,guesses_]:=Module[{emoments,br,g0=guesses,af,\[Beta]f,\[Chi]f},


br=FindRoot[{ellipseMoments[af,\[Beta]f,\[Chi]f][[1]]-m0,ellipseMoments[af,\[Beta]f,\[Chi]f][[2]]-m1,ellipseMoments[af,\[Beta]f,\[Chi]f][[3]]-m2},{{af,g0[[1]]},{\[Beta]f,g0[[2]]},{\[Chi]f,g0[[3]]}},Evaluated->False,MaxIterations-> 700];

Return[{af/.br,\[Beta]f/.br,\[Chi]f/.br}]


];


End[]

EndPackage[]








