(* ::Package:: *)

BeginPackage["StabililtyPackage`"]
ClearAll["StabililtyPackage`", "StabililtyPackage`"]


(* ::Subsection::Closed:: *)
(*Package Functions*)


ImportData::usage =
	"ImportData[file] reads in the data from file."
buildHamiltonians::usage=
	"[data,test Energy, Hierarchy=1/-1,ndens]. Return elements of the hamiltonians"
getEquations::usage=
	"[data,test energy, hierarchy=-1/1, k, ndens]. Return the equations of motions"
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
ndensities::usage = 
"finds the neutrino number densities for a species, and returns the matrix of densities in an angular bin"
getInitialGuess::usage=
"[file,species]generates Initial Guesses for the first radius from expansion expresssion"
eBoxFitSingleRadius::usage=
"fit a single radius using the box fit parameters"
Com::usage=
"test"
ellipseBoxMoments::usage=
"calculate the moments from an ellipse (box fit) with parameters a, \[Beta], \[Chi]"
ellipseSimpMoments::usage=
"calculate the moments from an ellipse (simple fit) with parameters a, b, cx"
eBoxFitToMoments::usage=
"given 3 input moments and set of guess parameters, fits ellipse parameters to moments"
eSimpFitToMoments::usage=
"given 3 input moments and set of guess parameters, fits ellipse parameters to moments"
getMoments::usage=
"takes moments out of moment data for a given file, radius, species. Energy integrated"
exportkadapt::usage=
"[outevs_,name_]"
ImportCalcGridData::usage=
"Imports the calcultion grid elements."
ImportCalcUniqueData::usage=
"Imports only unique values for each element, duplicates deleted."
ImportCalcOptions::usage=
"Imports 3 option settings, xflavor,inverse,krange"
ImportCalcInputs::usage=
"Imports the calc file, rsrt, rend,testE,hi,nstep used in the calculation"
ellipsefitfile::usage=
"gets ellipse parameters for a moments file"
ellipseparaerrors::usage=
"gives the errors in the ellisope fit approx for some given parameters"
boxToSimp::usage=
"converts from box transform to simple ellipses parameters"
simpToBox::usage=
"converts from simple ellipse parameters to box transform"
ellipseFitSingleSpecies::usage=
"fits a single species to an ellispe from rsrt to rend"
getInitialBoxGuess::usage=
"gets intiial guess for box transform coordinates"
getDOmoments::usage=
"calculate moments from a do file"
DOtoMoments::usage=
"convert a do file to a moment file"
getelipdata::usage=
"ellipse fits a moment file from rsrt to rend and returns association"
exportelipdata::usage=
"exports ellipse fit data in .h5 format"
exportDOasMoments::usage=
"exports DO data as moment style data"
MorinagaPlotter::usage=
"Make morinaga style plots"
makeThetaGrid::usage=
"gets a 10 bin theta grid refined (doubled) nref times"
compareDistributions::usage=
"Makes a plot given [srdodata_,srelipdata_,species_]"
ellipseEqnMoments::usage=
"analysi ellipse moment equations"
eEqnFitToMoments::usage=
"fit to moments with ellipse analytic equations"
ellipseEqnparaerrors::usage=
"ellise erros via equation"
esubfit::usage=
"fits just a and b when c has approached a"
eEqnMinFitToMoments::usage=
"fits to moments by least sqaures minimization of residuals"
ellipseBoxEqnMoments::usage=
"returns the moment equations in box transform coordinates"
eBoxEqnFitToMoments::usage=
"fits in the box transformed coordinates"
ellipseFitSingleSpeciesDO::usage=
"fits singel species given moments as an input"
kAdaptSlim::usage=
"for running datasets and exporting them in a particular txt file form"
averageEnergy::usage=
"compute the average energy for a file at a single radius"
buildS::usage=
"builds the Stability matrix given [infile,r,k,En,hier,ndens]"


(* ::Subsection::Closed:: *)
(*Units*)


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
\[CapitalDelta]m13sq=(2.4 10^-3 )everg^2;
\[Omega]EMev[En_]:=(\[CapitalDelta]m13sq)/(2 (En Meverg));
Begin["`Private`"]


(* ::Subsection::Closed:: *)
(*Import Functions*)


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
"Endensity"->data["Endensity"][[ri]] (*distribution functions [species,frequency,theta,phi]*),
"matters"->data["matters"][[ri]], (*densities*)
"Yes"->data["Yes"][[ri]], (*electron fractions *)
"freqs"->data["freqs"], (*freq grid in hz*)
"freqmid"->data["freqmid"], (*freq mid points*)
"muss"->data["muss"], (*Cos\[Theta] grid*)
"mids"->data["mids"],
"phis"-> data["phis"] (*"phi bin edges"*) (*Cos\[Theta] bin midpoints*)
];


ImportCalcGridData[infile_]:=Import[infile,"/grid_elements"];
ImportCalcOptions[infile_]:=Association[
"xflavor"-> ToExpression[Import[infile,{"Data","/settings/options"}]][[1]],
"inverse"->ToExpression[Import[infile,{"Data","/settings/options"}]][[2]],
"krange"->ToExpression[Import[infile,{"Data","/settings/options"}]][[3]]
];
ImportCalcInputs[infile_]:=Association[
"file"->ToString[ReadList[StringToStream[StringDelete[Import[infile,{"Data","settings/inputs"}],{"{","}"}]],Word,WordSeparators->{","}][[1]]],
"rsrt"-> ToExpression[ToString[ReadList[StringToStream[StringDelete[Import[infile,{"Data","settings/inputs"}],{"{","}"}]],Word,WordSeparators->{","}][[2]]]],
"rend"-> ToExpression[ToString[ReadList[StringToStream[StringDelete[Import[infile,{"Data","settings/inputs"}],{"{","}"}]],Word,WordSeparators->{","}][[3]]]],
"testE"->ToExpression[ToString[ReadList[StringToStream[StringDelete[Import[infile,{"Data","settings/inputs"}],{"{","}"}]],Word,WordSeparators->{","}][[4]]]],
"hi"-> ToExpression[ToString[ReadList[StringToStream[StringDelete[Import[infile,{"Data","settings/inputs"}],{"{","}"}]],Word,WordSeparators->{","}][[5]]]],
"nstep"-> ToExpression[ToString[ReadList[StringToStream[StringDelete[Import[infile,{"Data","settings/inputs"}],{"{","}"}]],Word,WordSeparators->{","}][[6]]]]
];


(* ::Subsection::Closed:: *)
(*Neutrino Densities and Potentials*)


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


siPotential[ndens_]:=Module[{tot,m},
m=munits ndens;
tot=(Tr[m[[1]]]-Tr[m[[3]]])-(Tr[m[[2]]]- Tr[m[[3]]]);
Return[tot]
]
(*This function finds the asymetry factor between the electron (anti)neutrinos and the x (anti)neutrinos *)
(*This is likely not needed and will need to be removed.*)
Options[Bfactor]={"xflavor"-> True};
Bfactor[data_,OptionsPattern[]]:=Module[{B,Bb,m},
(*m=munits ndensities[data,"xflavor"-> OptionValue["xflavor"]];*)
B=0.(*Tr[m[[3]]]/Tr[m[[1]]]*);
Bb=0.(*Tr[m[[3]]]/Tr[m[[2]]]*);
Return[{B,Bb}]
];


averageEnergy[dofile_,ri_]:=Module[{srdata,nuendensity,nubarendensity,Eav},
srdata=SelectSingleRadius[ImportData[dofile],ri];
nuendensity= Sum[Sum[Sum[srdata["Endensity"][[1,f,dt,dp]] ,{f,1,Length[srdata["freqs"]]-1}],{dp,1,Length[srdata["phis"]]-1}],{dt,1,Length[srdata["mids"]]}];
nubarendensity= Sum[Sum[Sum[srdata["Endensity"][[2,f,dt,dp]],{f,1,Length[srdata["freqs"]]-1}],{dp,1,Length[srdata["phis"]]-1}],{dt,1,Length[srdata["mids"]]}];
Eav=(nuendensity+nubarendensity)/((ndensities[srdata][[1]]//Tr)+(ndensities[srdata][[2]]//Tr));
Return[Eav]
]


(* ::Subsection::Closed:: *)
(*Stability Matrix Functions*)


(*This function creates the hamiltonians based on the data set, and returns them in a "general" way that is readable.  
Returns 9 arguments with index,
1,3,5,7 = H,\[Rho],A,\[Delta]H
2,4,6,8 = Hb,\[Rho]b,Ab,\[Delta]Hb
9=HsiRadial
*)
Options[buildHamiltonians]={"xflavor"-> True};
buildHamiltonians[data_,testE_,hi_,ndens_,OptionsPattern[]]:=Module[{n,\[Theta],name11,name12,name21,name22,\[Rho],\[Rho]b,A,Ab,Hm,Hvac,\[Mu],\[Mu]b,\[Mu]x,m,Hsi,H,Hb,\[Delta]H,\[Delta]Hb,Ve,\[Omega]},(
Ve=(*munits/mp *data["Yes"]  *data["matters"];*) 0.;
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
m=munits ndens;
\[Mu]=m[[1]]-m[[3]]; (*With x flavor, \[Mu]\[Rule] (Subscript[\[Mu], e]-Subscript[\[Mu], x])*)
\[Mu]b=m[[2]]-m[[3]]; (*Same for anti-neutrinos*)

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


Com[A_,B_]:=Module[{a=A,b=B},Return[A . B-B . A]];


(*Calculates the equations of motion by computing the relvant commutators. 
 Returns 5 arguments with index,
 1,3= neutrino equations of motion, A
 2,4 = Antineutrino equations of motion, Ab
 5=HsiRadial
 *)
Options[getEquations]={"xflavor"-> True,"inverse"-> False};
getEquations[data_,testE_,hi_,k_,ndens_,OptionsPattern[]]:=Module[{n,\[Theta],eqn,eqnb,hs},
hs=buildHamiltonians[data,testE,hi,ndens,"xflavor"-> OptionValue["xflavor"]];
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
]
(*Close with*)
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
S=ArrayFlatten[{{S1,S2},{S3,S4}}]/.Dispatch[rules[data,"xflavor"-> OptionValue["xflavor"]]];
Return[S];
]
];


(* ::Subsection::Closed:: *)
(*Stability Matrix Entry - by - Entry Method*)


(*This function builds the stability matrix by computing each component exactly, and then building the total table.*)
buildS[datasr_,k_,En_,hier_,ndens_]:=Module[{n,nb,Ve,\[Omega],\[Phi]0,\[Phi]1,\[Phi],cos\[Theta],nwhich,\[Omega]which,S,Smat,\[Mu]which},
cos\[Theta]=datasr["mids"];
\[Omega]=hier \[Omega]EMev[En];
Ve=0. ;(*munits/mp *datasr["Yes"]*datasr["matters"];*)
n=munits Diagonal[ndens[[1]]-ndens[[3]]];
nb=munits Diagonal[ndens[[2]]-ndens[[3]]];
\[Omega]which[i_]:=Which[i<= Length[cos\[Theta]], -1.,i>Length[cos\[Theta]],1.];
nwhich[j_]:=Which[j<= Length[cos\[Theta]], n[[j]],j>Length[cos\[Theta]], -nb[[(j-Length[cos\[Theta]])]] ];
\[Mu]which[i_]:= Which[i<= Length[cos\[Theta]],0.,i>Length[cos\[Theta]],Length[cos\[Theta]]];
\[Phi][rank_,mu_,mubar_,cos\[Theta]_]:=\[Phi][rank,mu,mubar,cos\[Theta]]=Sum[(mu[[i]]-mubar[[i]])cos\[Theta][[i]]^rank,{i,1,Length[cos\[Theta]]}];
S[i_,j_]:= (\[Omega]which[i] \[Omega] KroneckerDelta[i,j])+(Ve+\[Phi][0,n,nb,cos\[Theta]]+cos\[Theta][[i-\[Mu]which[i]]] (k-\[Phi][1,n,nb,cos\[Theta]])) KroneckerDelta[i,j]
-(1-(cos\[Theta][[i-\[Mu]which[i]]]cos\[Theta][[j-\[Mu]which[j]]]))nwhich[j];
Smat=Table[S[i,j],{i,1,2*Length[cos\[Theta]]},{j,1,2*Length[cos\[Theta]]}];
Return[Smat];
];


(* ::Subsection::Closed:: *)
(*Machine Scaled Eigensystem (evscale)*)


(*This scales the stability matrix up to a more managable scale based on the machine prescision, Solves for the eigenvalues, and then scales backs.  Returns a list of eigenvalues*)
Options[evscale]={"output"-> "RankedEigenvalues"};
evscale[ktest_,S_,kx_,OptionsPattern[]]:=Module[{\[Epsilon],As,kx0s,as,kxs,evals,evecs},
\[Epsilon]=$MachineEpsilon/2;
As=Expand[(S/\[Epsilon])/.kx->\[Epsilon] kxs];
kx0s=ktest/\[Epsilon];
If[OptionValue["output"]=="RankedEigenvalues",
as=Max[Im[\[Epsilon] Eigenvalues[N[As]/.kxs->kx0s]]];
];
If[OptionValue["output"]=="Eigenvalues",
as=\[Epsilon] Eigenvalues[N[As]/.kxs->kx0s];
];
If[OptionValue["output"]=="Eigenvectors",
as=\[Epsilon] Eigenvectors[N[As]/.kxs->kx0s];
];
If[OptionValue["output"]=="Eigensystem",
as=\[Epsilon] Eigensystem[N[As]/.kxs-> kx0s];
];
Return[as]
];


(* ::Subsection::Closed:: *)
(*k Grid and Adaptive k Solver*)


(*Constructs a nstep sized log spaced k grid based on the target k associated with the infile at radial bin r.  Currently the limits are 2 orders of magnitude above and below the target value, ignoring negatives for the moment *)
Options[buildkGrid]={"ktarget"-> "auto","krange"-> {10.^(-3),10.}};
buildkGrid[ndens_,nstep_,OptionsPattern[]]:=Module[{kgrid,fSpace,ktarget,kblow,kbhigh},
Which[OptionValue["ktarget"]== "auto",
ktarget=Abs[siPotential[ndens]],
OptionValue["ktarget"]!= "auto",
ktarget=OptionValue["ktarget"]
];
fSpace[min_,max_,steps_,f_: Log]:=InverseFunction[ConditionalExpression[f[#],min<#<max]&]/@Range[f@min,f@max,(f@max-f@min)/(steps-1)];
kblow=OptionValue["krange"][[1]];
kbhigh=OptionValue["krange"][[2]];

Which[OptionValue["ktarget"]!=  0.,
kgrid=Join[fSpace[ktarget*kblow,ktarget*kbhigh,nstep],-fSpace[ktarget*kblow,ktarget*kbhigh,nstep]];
,
OptionValue["ktarget"== 0.],
kgrid={0.};
];
Return[kgrid];
];


(*Run buildkGrid and SCalcScale for several radial bins.*)
Options[kAdapt]={"xflavor"-> True,"ktarget"-> "auto","krange"-> {10.^(-3),10.},"koutput"-> "Eigensystem","inverse"-> False,"timing"-> False};
kAdapt[infile_,rstr_,rend_,testE_,hi_,nstep_,OptionsPattern[]]:= Module[{kl,evout,data,singleRadiusData,ea,kvar,eout,pot,S,ndens,stimes,totaltime,dtime,dptime,pc,ttab,Sks,Sk,Skl},
stimes[1]=SessionTime[];
data=ImportData[infile];
stimes[2]=SessionTime[];

evout=
Reap[
	Do[
		stimes[3]=SessionTime[];
		singleRadiusData = SelectSingleRadius[data,rx];
		stimes[4]=SessionTime[];
		ndens = ndensities[singleRadiusData,"xflavor"-> OptionValue["xflavor"]];
		stimes[5]=SessionTime[];	
		S=buildS[singleRadiusData,kvar,testE,hi,ndens];
		stimes[6]=SessionTime[];
		kl=buildkGrid[ndens,nstep,"ktarget"-> OptionValue["ktarget"],"krange"-> OptionValue["krange"]];
		stimes[7]=SessionTime[];
		pot=siPotential[ndens];
		stimes[8]=SessionTime[];
		Skl=S/.Thread[{kvar}->#]&/@kl;
		Do[
	(*Sk=S/.kvar\[Rule] kl[[kx]];*)
	Sk=N[Skl[[kx]]];
	Sks=Sk/Min[Sk];
		Which[OptionValue["koutput"]== "Eigensystem",
	eout=Eigensystem[N[Sks]]*Min[Sk];
	,
	OptionValue["koutput"]== "RankedEigenvalues",
	eout=Max[Im[Eigenvalues[N[Sks]]*Min[Sk]]];
			];
			Sow[{rx,data["radius"][[rx]],kl[[kx]],eout,pot}]; 
		,{kx,1,Length[kl]}]; (*close do over ktargets*)
		stimes[9]=SessionTime[];
	,{rx,rstr,rend}] (*close do over r*)
][[2,1]];
stimes[10]=SessionTime[];
Which[OptionValue["timing"]== False,
Return[{evout,{OptionValue["xflavor"],OptionValue["inverse"],OptionValue["krange"]//N},{infile,rstr,rend,testE,hi,nstep}}]; (*Close reap over r*)
,
OptionValue["timing"]==True,
totaltime=stimes[10]-stimes[1];
dtime[i_,j_]:=stimes[j]-stimes[i];
dptime[i_,j_]:=100 (dtime[i,j]/totaltime);

pc=PieChart[{dtime[1,2],dtime[3,4],dtime[4,5],dtime[5,6],dtime[6,7],dtime[7,8],dtime[8,9]},ChartLegends-> {"Import Data","Select Single Radius","ndensities","buildS","buildkgrid","siPotential","evscale loop"}
,ChartLabels->Placed[ {dptime[1,2] ,dptime[3,4],dptime[4,5],dptime[5,6], dptime[6,7],dptime[7,8],(dptime[8,9]),dptime[9,10]},"RadialCallout"],ImageSize-> Scaled[0.25],ImagePadding-> All];

ttab=Grid[{{"Total Time (s)",totaltime,""},{"Function Call","Time Taken (s)","Percent of Total Time"},{"Import Data",dtime[1,2],dptime[1,2]},{"Select Single Radius",dtime[3,4],dptime[3,4]},{"ndensities",dtime[4,5],dptime[4,5]},
{"build S",dtime[5,6],dptime[5,6]},{"Build k Grid",dtime[6,7],dptime[6,7]},{"siPotential",dtime[7,8],dptime[7,8]},{"k loop",dtime[8,9],dptime[8,9]},{"Average evscale per k",dtime[9,10]/Length[kl],dptime[9,10]/Length[kl]}},Frame-> All];

Return[Row[{ttab,pc}]];
];

]; (*close module*)


Options[kAdaptSlim]={"rrange"-> {1,384},"xflavor"-> True,"ktarget"-> 0.,"krange"-> {10.^(-3),10.},"inverse"-> False};
kAdaptSlim[file_,nk_,OptionsPattern[]]:=Module[{max,kdat,etab,etab2,pts,evpts,rpts,tab,tab2},
kdat=kAdapt[file,OptionValue["rrange"][[1]],OptionValue["rrange"][[2]],Infinity,-1.,nk,"koutput"-> "RankedEigenvalues","krange"-> OptionValue["krange"],"ktarget"->OptionValue["ktarget"],"xflavor"-> OptionValue["xflavor"] ];
rpts=kdat[[1,All,2]]//DeleteDuplicates; (* Looks at all radial points, returns a list of all radii used*)
evpts=Map[Max,Partition[kdat[[1,All,4]],2*nk]]; (* For each r, picks out the eigenvalue for all k at that r*)
tab=Table[{rpts[[i]],evpts[[i]]},{i,1,Length[rpts]}];
tab2=PrependTo[tab,{{kdat[[2,All]]},{"radius (cm)","Im[\[CapitalOmega]]"}}];
Export[StringJoin[StringDrop[file,-3],"_kadapted_nk", ToString[nk],".txt"],tab2,"Table"];
]


exportkadapt[outevs_,name_]:=
Export[ToString[name]<>".h5",  {
"/grid_elements/ri"->{"Data"-> outevs[[1,All,1]]},
"/grid_elements/radius"-> {"Data"-> outevs[[1,All,2]],"Attributes"-> {"Units"-> "Centimeters"}},
"/grid_elements/k"-> {"Data"-> outevs[[1,All,3]],"Attributes"-> {"Units"-> "Ergs"}},
"/grid_elements/evs_Re"-> {"Data"-> Re[outevs[[1,All,4,1]]],"Attributes"-> {"Units"-> "Ergs"}},
"/grid_elements/evs_Im"-> {"Data"-> Im[outevs[[1,All,4,1]]],"Attributes"-> {"Units"-> "Ergs"}},
"/grid_elements/evecs_nu_Re"-> {"Data"-> Re[outevs[[1,All,4,2]][[All,All,1;;1/2 Length[outevs[[1,All,4,2]][[All,All]][[1]]] ]] ],"Attributes"-> {"Norm?"-> "Unnormalized"}},
"/grid_elements/evecs_nu_Im"-> {"Data"-> Im[outevs[[1,All,4,2]][[All,All,1;;1/2 Length[outevs[[1,All,4,2]][[All,All]][[1]]]]]],"Attributes"-> {"Norm?"-> "Unnormalized"}},
"/grid_elements/evecs_nubar_Re"-> {"Data"-> Re[outevs[[1,All,4,2]][[All,All,1+(1/2 Length[outevs[[1,All,4,2]][[All,All]][[1]]]);;Length[outevs[[1,All,4,2]][[All,All]][[1]]]]]],"Attributes"-> {"Norm?"-> "Unnormalized"}},
"/grid_elements/evecs_nubar_Im"-> {"Data"-> Im[outevs[[1,All,4,2]][[All,All,1+(1/2 Length[outevs[[1,All,4,2]][[All,All]][[1]]]);;Length[outevs[[1,All,4,2]][[All,All]][[1]]]]]],"Attributes"-> {"Norm?"-> "Unnormalized"}},
"/grid_elements/Vsi"-> {"Data"-> outevs[[1,All,5]],"Attributes"-> {"Units"-> "Ergs"}},
"/settings/options"-> {"Data"-> ToString[outevs[[2]]],"Attributes"-> {"Order"-> "xflavor value, inverse value,krange"}},
"/settings/inputs"-> {"Data"-> ToString[outevs[[3]]],"Attributes"-> {"Order"-> "file,rsrt,rend,testE,hi,nstep"}}
}
];
(*
---Indicies Guide for Editing---
kAdapt calculates data via 2 loops, one over radial position, and one over k.  For each r and k, kAdapt outputs the following, called "evout";
{Radial Index, Radial position in cm, k in ergs, eigensystem information, sipotential (\[Sqrt]2}Subscript[G, f](total neutrinos +total anti-neutrinos)}
---
The Eigensystem information is itself a list containing {eigenvalues, eigenvectors}, where the eigenvalues are a list of 2Subscript[n, \[Theta]] eigenvalues, an the eigenvectors contain 2Subscript[n, \[Theta]] vectors each of length 2Subscript[n, \[Theta].]
---
Thus, evout is a list of length number of Subscript[\[CapitalDelta]r, i] x nstep(1+1/2), where nstep is the number of k points in the positive direction (plus half as much in the negative k direction),
with each entry of evout being of length 5, with the length of the 4th entry (the eigensystem) having 2 parts, with the first part holding a 2Subscript[n, \[Theta]] list, and the second part holding a 2Subscript[n, \[Theta]] list with each entry being a 2Subscript[n, \[Theta]] list.
---
kAdapt itself exports a list of 3 entries; {evout, {List of Options},{list of inputs}.
The list of options is length 3: {xflavor, inverse, krange}. 
The length of the list of inputs is 6: {infile,rstrt,rend,testE,hi,nstep}.
---
To export grid elements, we export [[1,All,i]]; 1 for the evout entry, all for all r/k combinations, and i for a particular part of evout.  Example; For radial positions, the second entry in kadapt, [[1,All,2]].
To export eigenvalues, we jsut export all of the first entries of the 4th entry of evout. [[1,All,4,1]]
To export eigenvectors, we export the second part of the 4th entry, [[1,All,4,2]].  
To pick out just the first (or last) 10 entries of each eigenvector (the neutrino entries), we need to additionally pick out the first ten entries of each vector, for each eigenvalue So we need [[1,All,4,2]] then [[All,All,1;;10]].
---
Finally, to export outputs we jsut take the second list, and for the inputs the third.
*)
(*This currently outputs several datasets, each containing the unique elements from a run of kadapt, belonging to the group unique_elements.  Will add in groups that are, for instance, sorted by r. i.e. all k and omega combinations for a given radial index.*)




(* ::Subsection::Closed:: *)
(*Gershgorin Disks*)


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


(* ::Subsection::Closed:: *)
(*Ellipse Fitting Approximations*)


getMoments[file_,r_,species_]:= Module[{data,datasr,moments},
data=ImportData[file]//Quiet;
datasr=SelectSingleRadius[data,r];
moments={
Sum[(datasr["Endensity"][[species,E,1]])/(h data["freqmid"][[E]]),{E,1,Length[data["freqs"]]-1}],
Sum[(datasr["Endensity"][[species,E,2]])/(h data["freqmid"][[E]]),{E,1,Length[data["freqs"]]-1}],
Sum[(datasr["Endensity"][[species,E,3]])/(h data["freqmid"][[E]]),{E,1,Length[data["freqs"]]-1}]
};
If[moments[[2]]< 0., moments[[2]]=10^-8];
Return[moments//N];
];


getInitialGuess[m0_,m1_,m2_]:=Module[{foc1234,ag,\[Beta]g,\[Chi]g,cg,bg,arg\[Beta],arg\[Chi]},
ag=0.5 m0; 
bg=0.5 m0 3/2 (1- m2/m0);
cg=0.5 m1;
(*(*semi-major axis guess*)
cg=Abs[foc1234[Sin[ArcCos[1.]],0.,1.]-ag]; (*horizontal shift from the center *)
bg=Sqrt[foc1234[Sin[ArcCos[0.]],0.,0.]^2/(1-(cg/ag)^2)]; (*semi-minor axis*)
If[bg> ag && bg/ag<= 1.001, ag=ag+2(bg-ag)]; (* If bg>ag, and the difference is small, switches them in the transform*)
*)
Return[{ag,bg,cg}//N]
];


ellipseEqnMoments[af_,bf_,cxf_]:=Module[{em,fm,pm},
(*
em[a_,b_,cx_]:= (b (-2 Sqrt[b^2] cx ArcTanh[cx/a]-I a Sqrt[a^2-b^2-cx^2] (Log[Sqrt[b^2]-I Sqrt[a^2-b^2-cx^2]]-Log[Sqrt[b^2]+I Sqrt[a^2-b^2-cx^2]])))/(-a^2+b^2);
fm[a_,b_,cx_]:=-((2 b^2 cx (Sqrt[(a-b) (a+b)]-a ArcCoth[a/Sqrt[a^2-b^2]]))/((a-b) (a+b))^(3/2));
pm[a_,b_,cx_]:=(a b (2 Sqrt[b^2] (-a^2+b^2)+4 a Sqrt[b^2] cx ArcTanh[cx/a]+(I (a^4-b^2 cx^2-a^2 (b^2+cx^2)) (Log[Sqrt[b^2]-I Sqrt[a^2-b^2-cx^2]]-Log[Sqrt[b^2]+I Sqrt[a^2-b^2-cx^2]]))/Sqrt[a^2-b^2-cx^2]))/(2 (a^2-b^2)^2);
*)
em[a_,b_,cf_]:= (2 b (a Sqrt[a^2-b^2-cf^2] ArcCot[b/Sqrt[a^2-b^2-cf^2]]+b cf ArcTanh[cf/a]))/(a^2-b^2);
fm[a_,b_,cf_]:=-((2 b^2 cf (Sqrt[a^2-b^2]-a ArcCoth[a/Sqrt[a^2-b^2]]))/((a-b) (a+b))^(3/2));
pm[a_,b_,cf_]:=(a b (-a^2 b+b^3+((a^4-b^2 cf^2-a^2 (b^2+cf^2)) ArcCot[b/Sqrt[a^2-b^2-cf^2]])/Sqrt[a^2-b^2-cf^2]-a b cf Log[(a-cf)/(a+cf)]))/(a^2-b^2)^2;
Return[{em[af,bf,cxf],fm[af,bf,cxf],pm[af,bf,cxf]}];
];


ellipseBoxEqnMoments[a_,\[Beta]f_,\[Chi]f_]:=Module[{em,fm,pm},
em[af_,\[Beta]_,\[Chi]_]:=(af (1+Tanh[\[Beta]]) (1/4 af^2 ArcTanh[1/2 (1+Tanh[\[Chi]])] (1+Tanh[\[Beta]]) (1+Tanh[\[Chi]])+af ArcCot[(af (1+Tanh[\[Beta]]))/(2 Sqrt[af^2-1/4 af^2 (1+Tanh[\[Beta]])^2-1/4 af^2 (1+Tanh[\[Chi]])^2])] Sqrt[af^2-1/4 af^2 (1+Tanh[\[Beta]])^2-1/4 af^2 (1+Tanh[\[Chi]])^2]))/(af^2-1/4 af^2 (1+Tanh[\[Beta]])^2);
fm[af_,\[Beta]_,\[Chi]_]:=-((af^3 (1+Tanh[\[Beta]])^2 (-af ArcCoth[af/Sqrt[af^2-1/4 af^2 (1+Tanh[\[Beta]])^2]]+Sqrt[af^2-1/4 af^2 (1+Tanh[\[Beta]])^2]) (1+Tanh[\[Chi]]))/(4 ((af-1/2 af (1+Tanh[\[Beta]])) (af+1/2 af (1+Tanh[\[Beta]])))^(3/2)));
pm[af_,\[Beta]_,\[Chi]_]:=(af (1+Tanh[\[Beta]]) (-4 (1+Tanh[\[Beta]])+(1+Tanh[\[Beta]])^3-2 Log[-1+4/(3+Tanh[\[Chi]])] (1+Tanh[\[Beta]]) (1+Tanh[\[Chi]])+(af (4+2 E^(-2 \[Beta])+2 E^(-2 \[Chi])+E^(-2 (\[Beta]+\[Chi]))-2 E^(2 (\[Beta]+\[Chi]))) ArcCot[(af (1+Tanh[\[Beta]]))/Sqrt[-af^2 (-2+Tanh[\[Beta]] (2+Tanh[\[Beta]])+Tanh[\[Chi]] (2+Tanh[\[Chi]]))]] Sech[\[Beta]]^2 Sech[\[Chi]]^2)/Sqrt[-af^2 (-2+Tanh[\[Beta]] (2+Tanh[\[Beta]])+Tanh[\[Chi]] (2+Tanh[\[Chi]]))]))/((-1+Tanh[\[Beta]])^2 (3+Tanh[\[Beta]])^2);
Return[{em[a,\[Beta]f,\[Chi]f],fm[a,\[Beta]f,\[Chi]f],pm[a,\[Beta]f,\[Chi]f]}];
];


(*Converts the 3 box transform parameters to the simple ellipse parameters*)
boxToSimp[a_,\[Beta]_,\[Chi]_]:=Module[{b,cx},
b=a/2 (Tanh[\[Beta]]+1);
If[b==0.,
b=10^-10;
]; 
cx=a/2 (Tanh[\[Chi]]+1);
Return[{a,b,cx}]];


(*Converts the 3 simple elipse parameters into the box transform parameters*)
simpToBox[a_,b_,cx_]:=Module[{\[Beta],\[Chi]},
\[Beta]=ArcTanh[(2 b)/a-1]//Re;
\[Chi]=ArcTanh[(2 cx)/a-1]//Re;
Return[{a,\[Beta],\[Chi]}]];


eEqnFitToMoments[m0_,m1_,m2_,guesses_]:=Module[{emoments,br,g0=guesses,af,bf,cf,momeqns},
momeqns=ellipseEqnMoments[af,bf,cf]//N;
br=FindRoot[{(momeqns[[1]]-m0)/m0,(momeqns[[2]]-m1)/m0,(momeqns[[3]]-m2)/m0},{{af,g0[[1]]},{bf,g0[[2]]},{cf,g0[[3]]}},Evaluated->False,MaxIterations->20,AccuracyGoal-> 6];
Return[{Re[af/.br],Re[bf/.br],Re[cf/.br]}]
];


eEqnMinFitToMoments[m0_,m1_,m2_,guesses_]:=Module[{emoments,br,g0=guesses,af,bf,cf,momeqns},
momeqns=ellipseEqnMoments[af,bf,cf]//N;
br=FindMinimum[{1/m0 Total[momeqns^2],af>=cf},{{af,g0[[1]]},{bf,g0[[2]]},{cf,g0[[3]]}},AccuracyGoal->10,MaxIterations->1000];
Return[{Re[af/.br[[2]]],Re[bf/.br[[2]]],Re[cf/.br[[2]]]}]
];


eBoxEqnFitToMoments[m0_,m1_,m2_,guesses_]:=Module[{emoments,br,g0=guesses,g1,af,bf,cf,\[Beta]f,\[Chi]f,momeqns,ans},
momeqns=ellipseBoxEqnMoments[af,\[Beta]f,\[Chi]f]//N;
g1=Apply[simpToBox,guesses];
br=FindRoot[{(momeqns[[1]]-m0)/m0,(momeqns[[2]]-m1)/m0,(momeqns[[3]]-m2)/m0},{{af,g1[[1]]},{\[Beta]f,g1[[2]],-Infinity,Infinity},{\[Chi]f,g1[[3]],-Infinity,Infinity}},Evaluated->False,MaxIterations->20,AccuracyGoal-> 6];
ans=Apply[boxToSimp,{Re[af/.br],Re[\[Beta]f/.br],Re[\[Chi]f/.br]}];
Return[ans]
];


esubfit[m0_,m1_,ag_,bg_,cxg_]:=Module[{af,bf,momeqns,br},
momeqns=ellipseEqnMoments[af,bf,cxg]//N;
br=FindRoot[{(momeqns[[1]]-m0)/m0,(momeqns[[2]]-m1)/m0},{{af,ag},{bf,bg}},Evaluated->False,MaxIterations->20,AccuracyGoal-> 6];
Return[{Re[af/.br],Re[bf/.br],cxg//Re}]
]


ellipseEqnparaerrors[a_,b_,cx_,m0_,m1_,m2_]:=Module[{er0,er1,er2,fits},
(*again evaluate only once*)
er0=(ellipseEqnMoments[a,b,cx][[1]]-m0)/m0;
er1=(ellipseEqnMoments[a,b,cx][[2]]-m1)/m0;
er2=(ellipseEqnMoments[a,b,cx][[3]]-m2)/m0;
Return[{er0,er1,er2}//Re];
];


ellipseFitSingleSpecies[file_,species_,rsrt_,rend_]:=Module[{moms,igs,paras,errs,out,simpparas,boxparas,simperrs,boxerrs},
out=Reap[
	Do[
		moms=getMoments[file,ri,species];
		igs=Apply[getInitialGuess,moms];
		simpparas=Apply[eEqnFitToMoments,Join[moms,{igs}]]//Re;
		simperrs=Apply[ellipseEqnparaerrors,Join[simpparas,moms]];
		Sow[{simpparas,simperrs}];
	,{ri,rsrt,rend}](*Clsoe Do *)
][[2,1]];(*Close Reap*)
Return[out];
];


(*Special ellipse fit function to take DO moments as an input instead of a file*)
ellipseFitSingleSpeciesDO[moments_,species_,rsrt_,rend_]:=Module[{moms,igs,paras,errs,out,simpparas,boxparas,simperrs,boxerrs},
out=Reap[
	Do[
		igs=Apply[getInitialGuess,moments[[ri,All,species]]];
		simpparas=Apply[eEqnFitToMoments,Join[moments[[ri,All,species]],{igs}]]//Re;
		simperrs=Apply[ellipseEqnparaerrors,Join[simpparas,moments[[ri,All,species]]]];
		Sow[{simpparas,simperrs}];
	,{ri,rsrt,rend}](*Clsoe Do *)
][[2,1]];(*Close Reap*)
Return[out];
];


(* Caclulates the 3 angular moments from a DO file and returns a list of the 3 moments from radial index rsrt to radial index rend. *)
getDOmoments[dofile_,rsrt_,rend_]:=Module[{dat,srdat,out,emom,fmom,pmom},
dat=ImportData[dofile];
out=Reap[
	Do[
srdat=SelectSingleRadius[dat,ri];
emom={Sum[Sum[Sum[srdat["Endensity"][[1,f,dt,dp]] / (h (srdat["freqmid"][[f]])),{f,1,Length[dat["freqs"]]-1}],{dp,1,Length[dat["phis"]]-1}],{dt,1,Length[srdat["mids"]]}],
Sum[Sum[Sum[srdat["Endensity"][[2,f,dt,dp]] / (h (srdat["freqmid"][[f]])),{f,1,Length[dat["freqs"]]-1}],{dp,1,Length[dat["phis"]]-1}],{dt,1,Length[srdat["mids"]]}],
Sum[Sum[Sum[srdat["Endensity"][[3,f,dt,dp]] / (h (srdat["freqmid"][[f]])),{f,1,Length[dat["freqs"]]-1}],{dp,1,Length[dat["phis"]]-1}],{dt,1,Length[srdat["mids"]]}]
};
fmom={Sum[Sum[Sum[srdat["Endensity"][[1,f,dt,dp]] srdat["mids"][[dt]]/ (h (srdat["freqmid"][[f]])),{f,1,Length[dat["freqs"]]-1}],{dp,1,Length[dat["phis"]]-1}],{dt,1,Length[srdat["mids"]]}],
Sum[Sum[Sum[srdat["Endensity"][[2,f,dt,dp]] srdat["mids"][[dt]]/ (h (srdat["freqmid"][[f]])),{f,1,Length[dat["freqs"]]-1}],{dp,1,Length[dat["phis"]]-1}],{dt,1,Length[srdat["mids"]]}],
Sum[Sum[Sum[srdat["Endensity"][[3,f,dt,dp]] srdat["mids"][[dt]]/ (h (srdat["freqmid"][[f]])),{f,1,Length[dat["freqs"]]-1}],{dp,1,Length[dat["phis"]]-1}],{dt,1,Length[srdat["mids"]]}]
};
pmom={Sum[Sum[Sum[srdat["Endensity"][[1,f,dt,dp]] srdat["mids"][[dt]]^2/ (h (srdat["freqmid"][[f]])) ,{f,1,Length[dat["freqs"]]-1}],{dp,1,Length[dat["phis"]]-1}],{dt,1,Length[srdat["mids"]]}],
Sum[Sum[Sum[srdat["Endensity"][[2,f,dt,dp]] srdat["mids"][[dt]]^2/ (h (srdat["freqmid"][[f]])) ,{f,1,Length[dat["freqs"]]-1}],{dp,1,Length[dat["phis"]]-1}],{dt,1,Length[srdat["mids"]]}],
Sum[Sum[Sum[srdat["Endensity"][[3,f,dt,dp]] srdat["mids"][[dt]]^2/ (h (srdat["freqmid"][[f]])) ,{f,1,Length[dat["freqs"]]-1}],{dp,1,Length[dat["phis"]]-1}],{dt,1,Length[srdat["mids"]]}]
};
Sow[{emom,fmom,pmom}];
,{ri,rsrt,rend}]
][[2,1]];
Return[out];
];


(*Get a 10 bin theta grid refined nref times*)
makeThetaGrid[nref_]:=Module[{walls,middles},
walls=Import["walls.m"];(*Walls indicate where bin walls should be for a 10,20,40,80,160, or 320 bin grid. This is currently hardcoded in a .m file*)
middles=Import["middles.m"];(*Middles works the same as walls but for the bin centers.*)
Return[{walls[[nref+1]],middles[[nref+1]]}]
]


(*Take a discrete ordinates file, pull out it's moments, and returns data ready to be exported.  
nref is a integer \[GreaterEqual] 0 indicating the number of times the angular bins should be doubled. nref=0 means no grid refinement.
!!!This ASSUMES a nref=0 is a 10 angular bin grid.!!
*)
DOtoMoments[dofile_,rsrt_,rend_,nref_]:=Module[{data,nup,nbp,nxp,dodat,domom,walls,middles},
{walls,middles}=makeThetaGrid[nref];
dodat=ImportData[dofile];
domom=getDOmoments[dofile,rsrt,rend];
data= Association[
"muss"-> walls[[nref+1]],
"matters"-> dodat["matters"],
"Yes"-> dodat["Yes"],
"mids"-> middles[[nref+1]],
"freqs"->{0,2},
"Endensity"->Table[domom[[r,s,m]],{r,1,Length[domom]},{s,1,3},{f,1,1},{m,1,3}], (*Indicies; r is radial bin, s is species, f is freq bin, m is moments *)
 "freqmid"-> {1/h},
"phis"-> {0,2},
"radius"-> Table[dodat["radius"][[i]],{i,rsrt,rend}]
];
Return[data];
];


(*Exports a DOtoMoments dodata file in h5 format.*)
exportDOasMoments[dodata_,name_]:= 
Export[name<>".h5",{
"distribution(erg|ccm,lab)"-> {"Data"-> dodata["Endensity"]},
"rho(g|ccm,com)"-> {"Data"-> dodata["matters"]},
"Ye"-> {"Data"-> dodata["Yes"]},
"distribution_frequency_grid(Hz,lab)"-> {"Data"->dodata["freqs"]},
"distribution_frequency_mid(Hz,lab)"-> {"Data"-> dodata["freqmid"]},
"distribution_costheta_grid(lab)"-> {"Data"-> dodata["muss"]},
"distribution_costheta_mid(lab)"-> {"Data"-> dodata["mids"]},
"r(cm)"-> {"Data"-> dodata["radius"]},
"/distribution_phi_grid(radians,lab)"-> {"Data"-> dodata["phis"]}
}
]


(*Given a moment file, ellipse fits and returns an association from rsrt to rend*)
getelipdata[momentfile_,rsrt_,rend_,nref_]:=Module[{data,modat,efits,esimp,walls,middles,efitsord},
esimp[a_,b_,cx_,m_]:=(b (b m cx+a Sqrt[b^2 m^2-a^2 (-1+m^2)+(-1+m^2) cx^2]))/(a^2+(-a^2+b^2) m^2);
{walls,middles}=makeThetaGrid[nref];
modat=ImportData[momentfile];
efits={ellipseFitSingleSpecies[momentfile,1,rsrt,rend][[All,1]],ellipseFitSingleSpecies[momentfile,2,rsrt,rend][[All,1]],ellipseFitSingleSpecies[momentfile,3,rsrt,rend][[All,1]]};(*List of ellipse fits for the 3 species. 
part[[All,1]] takes the parameters for all radii as ellipseFitSingleSpecies has dimmensions {parameters, errors}*)
efitsord=Table[efits[[All,i]],{i,1,(rend-rsrt+1)}];(*Reorders so the endicies are in the same order as moment data*)
data= Association[
"muss"-> walls,
"matters"-> modat["matters"],
"Yes"-> modat["Yes"],
"mids"-> middles,
"freqs"->{0,2}, (*This is arbitrary*)
"Endensity"->
(*Indicies; s is species, r is radius, dt is theta bin, f is frew bin, dp is phi bin.*)
Abs[Table[esimp[efitsord[[r,s,1]],efitsord[[r,s,2]],efitsord[[r,s,3]],middles[[dt]]](walls[[dt+1]]-walls[[dt]]),{r,1,Length[efits[[1]]]},{s,1,3},{f,1,1},{dt,1,Length[middles]},{dp,1,1}]],
 "freqmid"-> {1/h},
"phis"-> {0,2}, (*This is arbitrary*)
"radius"-> Table[modat["radius"][[i]],{i,rsrt,rend}]
];
Return[data];
];


(*Exports elkipdata from getelipdata in the .h5 format*)
exportelipdata[name_,elipdata_]:= 
Export[name<>".h5",{
"distribution(erg|ccm,lab)"-> {"Data"-> elipdata["Endensity"]},
"rho(g|ccm,com)"-> {"Data"-> elipdata["matters"]},
"Ye"-> {"Data"-> elipdata["Yes"]},
"distribution_frequency_grid(Hz,lab)"-> {"Data"->elipdata["freqs"]},
"distribution_frequency_mid(Hz,lab)"-> {"Data"-> elipdata["freqmid"]},
"distribution_costheta_grid(lab)"-> {"Data"-> elipdata["muss"]},
"distribution_costheta_mid(lab)"-> {"Data"-> elipdata["mids"]},
"r(cm)"-> {"Data"-> elipdata["radius"]},
"/distribution_phi_grid(radians,lab)"-> {"Data"-> elipdata["phis"]}
}
]


(* ::Subsection::Closed:: *)
(*Plotting Tools*)


(*Calculates a plot in morinaga units for a dofile that and an output of kadapt (for an ellipsefit) on dofile called kadaptfile at radius ri.*)
MorinagaPlotter[dofile_,kadaptfile_,rido_,rielip_]:=Module[{datain,datainsr,cos\[Theta],datapos,evp,reimpairs,kspos,reevpos,imevpos,zeropos,nonzeropos,mu,mubar,\[Phi]0,\[Phi]1,Vmatter,kp,pts1,pts2,pts3,sp1,kdata},
datain=ImportData[dofile]; (*Import DO datafile*)
datainsr=SelectSingleRadius[datain,rido];
kdata=ImportCalcGridData[kadaptfile]; (*Import kadapt data file*)
cos\[Theta]=datain["mids"];
datapos=Position[kdata["ri"],rielip]; (*Finds the datapoints associated with radius ri; i.e. all ks, eigenvalues, eigenvectors etc*)
evp=Table[Transpose@{Extract[kdata["evs_Re"],datapos[[i]]],Extract[kdata["evs_Im"],datapos[[i]]]},{i,1,Length[datapos]}]; (*Extracts the position of the relavent eigenvalues at position datapos*)
reimpairs=Table[Sort[evp[[i]],#1[[2]]>#2[[2]]&][[1]],{i,1,Length[datapos]}]; (*Creates a list of Re and Im eigenvalue pairs, sorted by the largest imaginary part.*)
kspos=Table[Extract[kdata["k"],datapos[[i]]],{i,1,Length[datapos]}]; (*Extracts the positions of the relevant k values at datapos*)
reevpos=Transpose@{kspos,reimpairs[[All,1]]}; (*Ordered pairs of k and the re eigenvalues*)
imevpos=Transpose@{kspos,reimpairs[[All,2]]}; (*Ordered pairs of k and Im eigenvalues*)
zeropos=Position[imevpos[[All,2]],0.]//Flatten; (*Poisitions of eigenvalues with 0 Im part*)
nonzeropos=Delete[datapos//Flatten,Position[imevpos[[All,2]],0.]]; (*Positions of the eigenvalues with nonzero Im parts*)
Do[reevpos[[zeropos[[i]],2]]=0.,{i,1,Length[zeropos]}]; (*Sets the real parts of eigenvalues with no Im part to 0 for plotting purposes*)
(*Begin calculation of quantities for calculating k' and \[CapitalOmega]'*)
mu=munits Diagonal[ndensities[datainsr,"xflavor"->False][[1]]];
mubar=munits Diagonal[ndensities[datainsr,"xflavor"->False][[2]]];
\[Phi]0= Sum[(mu[[i]]-mubar[[i]])         ,{i,1,Length[cos\[Theta]]}];
\[Phi]1=Sum[(mu[[i]]-mubar[[i]])cos\[Theta][[i]],{i,1,Length[cos\[Theta]]}];
Vmatter =munits datainsr["Yes"] datainsr["matters"]/mp;
(*End calculating prime quantities*)
kp = Table[kspos[[i]]-\[Phi]1,{i,1,Length[kspos]}]; (*Table of k's*)
pts1=Table[{(1/(hbar c 10^-4))kp[[i]],(1/(hbar c 10^-4))imevpos[[i,2]]},{i,1,Length[datapos]}]; (*Ordered pairs of k' and Im eigenvalues in morinaga units*)
pts2=Table[{(1/(hbar c 10^-4))kp[[i]],(1/(hbar c 10^-4))(0.05)(reevpos[[i,2]]-\[Phi]0-Vmatter)},{i,1,Length[datapos]}]; (*Ordered pairs of k' and Re|\[CapitalOmega]'| in morinaga units*)
pts3=Table[{pts2[[i,1]],(pts2[[i+1,2]]-pts2[[i,2]])/(pts2[[i+1,1]]-pts2[[i,1]])},{i,1,Length[pts2]-1}]; (*Ordered pairs of k' and dRe|\[CapitalOmega]'|/dk'*)
sp1=ListPlot[{pts1,pts2,pts3},PlotRange->{{-30,85},{-1,1.5}},Filling-> {1-> Axis,2-> Axis},AxesLabel-> {Style["k' (\!\(\*SuperscriptBox[\(10\), \(-4\)]\) \!\(\*SuperscriptBox[\(cm\), \(-1\)]\))",FontSize-> 14,Bold],Style["Im[\[CapitalOmega]] (\!\(\*SuperscriptBox[\(10\), \(-4\)]\) \!\(\*SuperscriptBox[\(cm\), \(-1\)]\))",FontSize-> 14,Bold]},
ImageSize-> Scaled[0.65],Frame-> False,PlotStyle->{Blue,Red,Darker[Green,0.4],Purple},PlotLegends-> {"IM|\[CapitalOmega]|","0.05 Re|\[CapitalOmega]'|","\!\(\*FractionBox[\(d\), \(dk\)]\)Re|\[CapitalOmega]'|"}]; (*The plot*)
Return[sp1];
]


(*Takes in single radius selected DO data and single radius selected elipse fit data and makes a distribution plot for a species*)
compareDistributions[srdodata_,srelipdata_,species_]:=Module[{dopts,efpts,plot},
dopts=Transpose@{srdodata["mids"],ndensities[srdodata][[species]]//Diagonal};
efpts=Transpose@{srdodata["mids"],ndensities[srelipdata][[species]]//Diagonal};
plot=ListPlot[{dopts,efpts},Joined-> {False,True},PlotRange-> All,PlotLegends-> {"DO","Elipse"},AxesLabel-> {Style["cos\[Theta]",FontSize-> 14,Bold],Style["#/\!\(\*SuperscriptBox[\(cm\), \(3\)]\)",FontSize-> 14,Bold]},ImageSize-> Scaled[0.65]];
Return[plot]
]


End[]
EndPackage[]
