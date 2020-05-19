(* ::Package:: *)

(* ::Section::Closed:: *)
(*Preamble Loading*)


(* ::Input::Initialization:: *)
CloudDeploy[EvaluationNotebook[]] (*Saves to the wolfram cloud instead of locally. *)
ParallelEvaluate[$ProcessorCount] (*Open a kernel on all available cores.  Cuchulainn: 8   Kneller06: 4 *)
(*
file="/home/sdflynn2/Documents/Wolfram Mathematica/MC_moments_SR.h5";
filedo="/home/sdflynn2/Documents/Wolfram Mathematica/DO_SR_v2.h5";
fileang="/home/sdflynn2/Documents/Wolfram Mathematica/MC_angular_SR.h5";
*)

SetDirectory[] (*Sets working directory to local default*)



(*Import the radial grid, and ellipse fit parameter grid from the cloud.  Radial grid is the only thing used here. *)
rs = CloudImport["rgrid.m"];
paragrid=CloudImport["paragrid.m"];
paragridnoSR=CloudImport["paras_noSR.m"];


(*Import from lots 'o data*)

pathlots="G:\\My Drive\\Physics\\Neutrino Oscillation Research\\Fast Conversions\\lotsadata.tar\\lotsadata\\lotsadata";
SetDirectory[pathlots]; (*Sets working directory to the folder containing "lotsadata"*)
filenames=FileNames[]
lotsodo=ParallelTable[Import[StringJoin[pathlots,"\\",ToString[filenames[[i]]]],{"Datasets","distribution(erg|ccm,lab)"}],{i,1,Length[filenames]}]; (*distribution functions*)
matters=ParallelTable[Import[StringJoin[pathlots,"\\",ToString[filenames[[i]]]],{"Datasets","rho(g|ccm,com)"}],{i,1,Length[filenames]}]; (*densities*)
Yes=ParallelTable[Import[StringJoin[pathlots,"\\",ToString[filenames[[i]]]],{"Datasets","Ye"}],{i,1,Length[filenames]}]; (*electron fractions *)
freqs=ParallelTable[Import[StringJoin[pathlots,"\\",ToString[filenames[[i]]]],{"Datasets", "distribution_frequency_grid(Hz,lab)"}],{i,1,Length[filenames]}]; (*freq grid in hz*)
freqmid=ParallelTable[Import[StringJoin[pathlots,"\\",ToString[filenames[[i]]]],{"Datasets", "distribution_frequency_mid(Hz,lab)"}],{i,1,Length[filenames]}]; (*freq mid points*)
muss=ParallelTable[Import[StringJoin[pathlots,"\\",ToString[filenames[[i]]]],{"Datasets", "distribution_costheta_grid(lab)"}],{i,1,Length[filenames]}]; (*Cos\[Theta] grid*)
mids=ParallelTable[Import[StringJoin[pathlots,"\\",ToString[filenames[[i]]]],{"Datasets", "distribution_costheta_mid(lab)"}],{i,1,Length[filenames]}]; (*Cos\[Theta] bin midpoints*)
SetDirectory[] (*sets working directory back to local default*)




(* ::Section::Closed:: *)
(*Constants, Stability Modules, Scaling Modules, SCalc modules, and GD Modules*)


(* ::Input::Initialization:: *)
(*some constants*)
c=2.99792458 10^10; (* cm/s*)
h=6.6260755 10^-27; (*erg s*)
hbar = h/2 Pi; (*erg s*)
Gf=1.1663787 10^-5; (*GeV^-2*)
everg=1.60218 10^-12; (* convert eV to ergs*)
Geverg = everg*10^9 (* convert GeV to ergs *)
ergev=1.0/everg; (*convert ergs to eV*)
cmkm=10^-5; (*convert cm to km*)
ergmev=ergev/10^6; (*convert erg to MeV*)
mp=1.6726219 10^-24; (*Proton mass in g*)
munits=Sqrt[2] Gf/Geverg^2 (hbar c)^3; (*Sqrt[2] Gf in erg cm^3*)

buildInitials[n_]:=( 
name11="ee";
name12="ex";
name21="xe";
name22="xx";
Do[
\[Rho][i]={{ToExpression[StringJoin["\[Rho]",name11,ToString[i]]],ToExpression[StringJoin["\[Rho]",name12,ToString[i]]]},{ToExpression[StringJoin["\[Rho]",name21,ToString[i]]],ToExpression[StringJoin["\[Rho]",name22,ToString[i]]]}};
\[Rho]b[i]={{ToExpression[StringJoin["\[Rho]b",name11,ToString[i]]],ToExpression[StringJoin["\[Rho]b",name12,ToString[i]]]},{ToExpression[StringJoin["\[Rho]b",name21,ToString[i]]],ToExpression[StringJoin["\[Rho]b",name22,ToString[i]]]}};
A[i]={{0,ToExpression[StringJoin["A",name12,ToString[i]]]},{ToExpression[StringJoin["A",name21,ToString[i]]],0}};
Ab[i]={{0,ToExpression[StringJoin["Ab",name12,ToString[i]]]},{ToExpression[StringJoin["Ab",name21,ToString[i]]],0}};
,{i,1,n}];
);

defineFunctions[n_]:=( 
Com[A_,B_]:=Module[{a=A,b=B},
Return[A.B-B.A]
];
);

buildHamiltonians[n_,\[Theta]_]:=(
Hm={{Ve,0.},{0.,0.}};
Hvac=hi{{-\[Omega]/2,0.},{0.,\[Omega]/2}};
\[Mu]=munits Table[nv[i],{i,1,n},{j,1,n}];
\[Mu]b=munits Table[nvb[i],{i,1,n},{j,1,n}];
Do[
Hsi[j]=Sum[\[Mu][[k,j]]\[Rho][k](2Pi)(1-Cos[\[Theta][[j]]]Cos[\[Theta][[k]]]),{k,1,n}]+Sum[-\[Mu]b[[k,j]]\[Rho]b[k](2Pi)(1-Cos[\[Theta][[j]]]Cos[\[Theta][[k]]]),{k,1,n}];
,{j,1,n}];
Do[
H[i]=Hvac+Hm+Hsi[i];
Hb[i]=Hvac-Hm-Hsi[i];
\[Delta]H[i]=Sum[((D[H[i][[1,2]],\[Rho][j][[1,2]]]) A[j])+((D[H[i][[1,2]],\[Rho]b[j][[1,2]]])Ab[j]),{j,1,n}];
\[Delta]Hb[i]=Sum[((D[Hb[i][[1,2]],\[Rho][j][[1,2]]])A[j])+((D[Hb[i][[1,2]],\[Rho]b[j][[1,2]]])Ab[j]),{j,1,n}];
,{i,1,n}];
)
getEquations[n_,\[Theta]_]:=( 
Do[
eqn[j]=Com[H[j],A[j]][[1,2]]+ Com[\[Delta]H[j],\[Rho][j]][[1,2]]+(k Cos[\[Theta][[j]]] A[j][[1,2]]);
eqnb[j]=-Com[Hb[j],Ab[j]][[1,2]]- Com[\[Delta]Hb[j],\[Rho]b[j]][[1,2]]+(k Cos[\[Theta][[j]]] Ab[j][[1,2]]);
,{j,1,n}];
);

giveValues[freq_,ve_,n_,hir_,kx_]:=( 
\[Omega]=freq;
Ve=ve;
hi=hir;
k=kx;

r1=Table[ToExpression["\[Rho]ee"<>ToString[i]]-> 1.,{i,1,n}];
r2=Table[ToExpression["\[Rho]ex"<>ToString[i]] -> 0.,{i,1,n}];
r3=Table[ToExpression["\[Rho]xe"<>ToString[i]]-> 0.,{i,1,n}];
r4=Table[ToExpression["\[Rho]xx"<>ToString[i]]-> 0.,{i,1,n}];

rb1=Table[ToExpression["\[Rho]bee"<>ToString[i]]-> 1.,{i,1,n}];
rb2=Table[ToExpression["\[Rho]bex"<>ToString[i]]-> 0.,{i,1,n}];
rb3=Table[ToExpression["\[Rho]bxe"<>ToString[i]]-> 0.,{i,1,n}];
rb4=Table[ToExpression["\[Rho]bxx"<>ToString[i]]-> 0.,{i,1,n}];

rrules=Flatten[{r1,r2,r3,r4,rb1,rb2,rb3,rb4}];


);

buildMatrix[n_]:=Module[{S1,S2,S3,S4,S}, 
S1=ParallelTable[Coefficient[eqn[l],A[m][[1,2]]],{l,1,n},{m,1,n}];
S2=ParallelTable[Coefficient[eqn[l],Ab[m][[1,2]]],{l,1,n},{m,1,n}];
S3=ParallelTable[Coefficient[eqnb[l],A[m][[1,2]]],{l,1,n},{m,1,n}];
S4=ParallelTable[Coefficient[eqnb[l],Ab[m][[1,2]]],{l,1,n},{m,1,n}];
S=ArrayFlatten[{{S1,S2},{S3,S4}}]/.rrules;
Return[S];
];

stabilityMatrix[z_,t_,w_,k_,ve_,hi_]:=( 
buildInitials[z];
defineFunctions[z];
buildHamiltonians[z,t];
getEquations[z,t];
giveValues[w,ve,z,hi,k];
Return[buildMatrix[z]];
);

build2bMatrix[wt_,kx_]:=Block[{S2ba,S2b},
S2ba=stabilityMatrix[2,{0.,Pi},wt,kx,0.,-1.];
S2b={{S2ba[[1,1]],S2ba[[1,4]]},{S2ba[[4,1]],S2ba[[4,4]]}};
Return[S2b]
];

\[Omega]Eev[En_]:=\[Omega]Eev[En]=(7 10^-5)/(2 En) everg; (*Constant Vacuum Frequency*)
(*For file of number ind in filenames list, at radial bin rta, and in angular bin with walls dt to dt+1, returns the distribution function in that big, "integrated" over energies, with the average energy and anglular bin width divided out*)
(* dolots multiplied by the bin width returns the number of neutrinos in that angular bin across all energies*)

dolots[ind_,rta_,dt_]:=
Sum[Sum[lotsodo[[ind,rta,1,f,dt,dp]] /(h(freqmid[[ind,f]]) (Abs[muss[[ind,dt+1]]-muss[[ind,dt]]])),{f,1,Length[freqs[[ind]]]-1}],{dp,1,2}];
adolots[ind_,rta_,dt_]:=
Sum[Sum[lotsodo[[ind,rta,2,f,dt,dp]] /(h(freqmid[[ind,f]]) (Abs[muss[[ind,dt+1]]-muss[[ind,dt]]])),{f,1,Length[freqs[[ind]]]-1}],{dp,1,2}];


(*Performs stability analysis and returns ranked eigenvalues, and potentials, for various inputs.
From radial bin # (not value) rsrt, to rend, 
with vacuum freqency for a testE energy neutrino,
for file # ind in filenames list,
for wave number ktest,
for hierarchy hi =1 for standard, hi=-1 for inverted
*)


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



Clear[GDPlot,GDValues];
GDPlot[tm_]:=GDPlot[tm]=Block[{cs,rrow,rcol,drow,dcol,reg,tmv,cp,cpc,cpr,cpreg},cs=Diagonal[tm];
rrow=Total[Abs[tm],{2}]-Abs[cs];
rcol=Total[Abs[tm],{1}]-Abs[cs];
drow=MapThread[Disk,{ReIm[cs],rrow}];
dcol=MapThread[Disk,{ReIm[cs],rcol}];
reg=MeshPrimitives[DiscretizeRegion[RegionIntersection[Apply[RegionUnion,drow],Apply[RegionUnion,dcol]]],2];
tmv=Max[Join[rrow,rcol]];
cp=ComplexListPlot[cs,PlotStyle->{Directive[GrayLevel[0],PointSize[0.02]]},Sequence[Prolog->{Opacity[0.35],EdgeForm[Black],Gray,reg,Blue,drow,Red,dcol},PlotRange->N[tmv*2],ImageSize->Large,PlotLegends->{"Disk center"},AxesLabel->{"Real Part","Imaginary Party"}]];
cpr=ComplexListPlot[cs,PlotStyle->{Directive[GrayLevel[0],PointSize[0.02]]},Sequence[Prolog->{Opacity[0.35],EdgeForm[Black],Red,drow},PlotRange->N[tmv*2],ImageSize->Large,PlotLegends->{"Disk center"},AxesLabel->{"Real Part","Imaginary Party"}]];
cpc=ComplexListPlot[cs,PlotStyle->{Directive[GrayLevel[0],PointSize[0.02]]},Sequence[Prolog->{Opacity[0.35],EdgeForm[Black],Cyan,dcol},PlotRange->N[tmv*2],ImageSize->Large,PlotLegends->{"Disk center"},AxesLabel->{"Real Part","Imaginary Party"}]];
cpreg=ComplexListPlot[cs,PlotStyle->{Directive[GrayLevel[0],PointSize[0.02]]},Sequence[Prolog->{Opacity[0.35],EdgeForm[Black],Gray,reg},PlotRange->N[tmv*2],ImageSize->Large,PlotLegends->{"Disk center"},AxesLabel->{"Real Part","Imaginary Party"}]];
Return[{cp,cpr,cpc,cpreg}]];
GDValues[tm_]:=GDValues[tm]=Block[{cs,rrow,rcol,drow,dcol,tmv,mv,reg,evm,ev,As,\[Epsilon],kxs,tms},
\[Epsilon]=$MachineEpsilon;
tms=Expand[(tm/\[Epsilon])];
cs=Diagonal[tms];
rrow=Total[Abs[tms],{2}]-Abs[cs];
rcol=Total[Abs[tms],{1}]-Abs[cs];
drow=MapThread[Disk,{ReIm[cs],rrow}];
dcol=MapThread[Disk,{ReIm[cs],rcol}];
VerificationTest[Count[Join[rrow,rcol],0]==0,True];
reg=RegionIntersection[Apply[RegionUnion,drow],Apply[RegionUnion,dcol]];
tmv=\[Epsilon] Max[Join[rrow,rcol]];
mv=\[Epsilon] Max[MeshCoordinates[BoundaryDiscretizeRegion[reg]][[All,2]]];
(*evm=\[Epsilon] Max[Im[ Eigenvalues[tms]]]*);
Return[{tmv,mv}]];




ktabr=Table[x,{x,-2 10^-17,4 10^-17,10^-18}];
rtabr=rs[[180;;250]];


(* ::Section::Closed:: *)
(*Non - adaptive automated calculations*)


(* ::Input:: *)
(*(**)
(*(*Calculate Data Section*)*)
(*indss={12,13,29};*)
(*Do[*)
(*CreateDirectory["rscan"<> ToString[indss[[m]]]];*)
(*SetDirectory["C:\\Users\\Sam\\rscan"<>ToString[indss[[m]]]];*)
(*Do[*)
(*scantab=Table[SCalclotsScale[179+i,179+i,20 10^6,indss[[m]],ktabr[[j]],-1],{j,1,Length[ktabr]}];*)
(*scanmat=Table[SCalclotsScaleMat[179+i,179+i,20 10^6,indss[[m]],ktabr[[j]],-1],{j,1,Length[ktabr]}];*)
(*(*ScanGD=Table[GDValues[scanmat[[i]]],{i,1,Length[scanmat]}];*)*)
(*Export["rscan"<>ToString[indss[[m]]]<>"r"<> ToString[i+179]<>".m",scantab];*)
(*Export["scanmat"<>ToString[indss[[m]]]<>"r"<> ToString[i+179]<>".m",scanmat];*)
(*(*Export["scanGD"<>ToString[indss[[m]]]<>"r"<> ToString[i+179]<>".m",scanGD];*)*)
(*,{i,1,Length[rtabr]}*)
(*];*)
(*SetDirectory[];*)
(*Print[ToString[m]<>" has completed"]*)
(*,{m,1,Length[indss]}*)
(*];*)
(**)
(*indss=Table[i,{i,3,19}];*)
(*Do[*)
(*SetDirectory["C:\\Users\\Sam\\rscan"<>ToString[indss[[m]]]];*)
(*rfilenames=FileNames["rscan*.m"];*)
(*rfile=Reap[Do[*)
(*raws=Import[rfilenames[[i]]];*)
(*Sow[raws[[All,2,1,1,1,1]]];*)
(*,{i,1,Length[rfilenames]}*)
(*]*)
(*][[2,1]];*)
(**)
(*potfile=Reap[Do[*)
(*rawspoy=Import[rfilenames[[i]]];*)
(*Sow[rawspoy[[1,2,2,1,1]]];*)
(*,{i,1,Length[rfilenames]}*)
(*]*)
(*][[2,1]];*)
(*(*GDfile=Reap[Do[*)
(*gdfilenames=FileNames["scanGD*.m"];*)
(*rawgd=Import[gdfilenames[[i]]];*)
(*Sow[rawgd[[2]]];*)
(*,{i,1,Length[gdfilenames[[i]]]}*)
(*]*)
(*][[2,1]];*)
(**)*)
(**)
(*SetDirectory[];*)
(*Export["rscan"<>ToString[indss[[m]]]<>".m",rfile];*)
(*Export["rpot"<>ToString[indss[[m]]]<>".m",potfile];*)
(*(*Export["rgd"<>ToString[indss[[m]]]<>".m",GDfile];*)*)
(*,{m,1,Length[indss]}]*)
(**)
(**)
(**)*)


(* ::Input:: *)
(*(**)
(*indss={5,3};*)
(*Do[*)
(*SetDirectory["C:\\Users\\Sam\\rscan"<>ToString[indss[[m]]]];*)
(*matfilenames=FileNames["scanmat*"];*)
(*Do[*)
(*matraw=Import[matfilenames[[i]]];*)
(*scanGD=GDValues[matraw];*)
(*Export["scanGD"<>ToString[indss[[m]]]<>"r"<> ToString[i+179]<>".m",scanGD];*)
(*,{i,1,Length[rtabr]}*)
(*];*)
(*SetDirectory[];*)
(*Print[ToString[m]<>" has completed"]*)
(*,{m,1,Length[indss]}*)
(*];*)
(**)
(**)
(*Do[*)
(*SetDirectory["C:\\Users\\Sam\\rscan"<>ToString[indss[[m]]]];*)
(*GDfile=Reap[Do[*)
(*gdfilenames=FileNames["scanGD*.m"];*)
(*rawgd=Import[gdfilenames[[i]]];*)
(*Sow[rawgd[[2]]];*)
(*,{i,1,Length[gdfilenames[[i]]]}*)
(*]*)
(*][[2,1]];*)
(*SetDirectory[];*)
(*Export["rgd"<>ToString[indss[[m]]]<>".m",GDfile];*)
(*,{m,1,Length[indss]}]*)
(**)*)


(* ::Section::Closed:: *)
(*k adaptive automated calculations*)


(* ::Input:: *)
(*  (* This is the run all for adaptive k steps *)*)


(* ::Input:: *)
(*indss={2,31,36,37,4,5,6,7,8,3,9,23};*)
(*SetDirectory[];*)
(*fSpace[min_,max_,steps_,f_: Log]:=InverseFunction[ConditionalExpression[f[#],min<#<max]&]/@Range[f@min,f@max,(f@max-f@min)/(steps-1)]*)
(*kadapt[ind_]:=Table[Join[fSpace[ktarget[ind][[i]] 10^-1,ktarget[ind][[i]] 10^1,50],-fSpace[ktarget[ind][[i]] 10^-1,ktarget[ind][[i]] 10^1,20]],{i,1,71}];*)
(*Do[*)
(*SetDirectory["C:\\Users\\Sam\\rscan"<> ToString[indss[[m]]]];*)
(*matfilenames=FileNames["scanmat*"];*)
(*ktarget[indss[[m]]]=Reap[*)
(*Do[*)
(*testm=Import[matfilenames[[i]]];*)
(*Sow[kx/.Solve[testm[[1,2,2,1,1,1,1]]==0,kx][[1]]]*)
(*,{i,1,Length[matfilenames]}]*)
(*][[2,1]];*)
(*Export["ktarget"<>ToString[indss[[m]]]<>".m",ktarget[indss[[m]]]];*)
(*Export["kadapt"<>ToString[indss[[m]]]<>".m",kadapt[indss[[m]]]];*)
(*,{m,1,Length[indss]}*)
(*];*)
(**)
(*SetDirectory[];*)
(**)
(*Do[*)
(*CreateDirectory["rscan_targeted_adapt"<> ToString[indss[[m]]]];*)
(*SetDirectory["C:\\Users\\Sam\\rscan_targeted_adapt"<>ToString[indss[[m]]]];*)
(*Do[*)
(*scantab=Table[SCalclotsScale[179+i,179+i,20 10^6,indss[[m]],kadapt[indss[[m]]][[i,j]],-1],{j,1,Length[kadapt[indss[[m]]][[1]]]}];*)
(*Export["rscan_adapt"<>ToString[indss[[m]]]<>"r"<> ToString[i+179]<>".m",scantab];*)
(*,{i,1,Length[rtabr]}*)
(*];*)
(*Print["rscans_"<>ToString[indss[[m]]]<>" has completed at time: "<>DateString["Time"]];*)
(*rfilenames=FileNames["rscan_adapt"<>ToString[indss[[m]]]<>"r*.m"];*)
(*rfile=Reap[Do[*)
(*raws=Import[rfilenames[[i]]];*)
(*Sow[raws[[All,2,1,1,1,1]]];*)
(*,{i,1,Length[rfilenames]}*)
(*]*)
(*][[2,1]];*)
(*SetDirectory[];*)
(*Export["rscan_adapt"<>ToString[indss[[m]]]<>".m",rfile];*)
(*,{m,1,Length[indss]}*)
(*];*)
(**)
(**)


(* ::Section::Closed:: *)
(*File collection and collating routines*)


(* ::Input:: *)
(*scanfilenames=FileNames["rscan_adapt*.m"];*)
(*potfilenames=FileNames["rpot*.m"];*)
(*(*gdfilenames=FileNames["rgd*.m"];*)*)
(*rscanall=Table[Import[scanfilenames[[i]]]*)
(*,{i,1,Length[scanfilenames]}];*)
(*potscanall=Table[Import[potfilenames[[i]]]*)
(*,{i,1,Length[potfilenames]}];*)
(*(*gdscanall=Table[Import[gdfilenames[[i]]]*)
(*,{i,1,Length[gdfilenames]}];*)*)
(*key=Grid[*)
(*Table[*)
(*{i,scanfilenames[[i]],filenames[[FromDigits[ToExpression[StringCases[scanfilenames[[i]],DigitCharacter]],10]]]}*)
(*,{i,1,Length[scanfilenames]}]*)
(*,Frame-> All]*)


(* ::Input:: *)
(*Do[*)
(*Print[i];*)
(*rscanps[i]=Table[{rtabr[[x]],ktabr[[y]],rscanall[[i,x,y]]},{x,1,71},{y,1,61}];*)
(*,{i,1,Length[scanfilenames]}*)
(*];*)
(**)
(**)


(* ::Input:: *)
(*Do[*)
(*If[*)
(*Position[rscanps[i],Null]!= {},*)
(*Print[i,filenames[[i]]]*)
(*]*)
(*,{i,1,Length[filenames]}*)
(*]*)


(* ::Section::Closed:: *)
(*Crossing finder*)


(* ::Input:: *)
(**)


(* ::Input:: *)
(*(**)
(*crosslook[ind_,r_]:=crosslook[ind,r]=ParallelTable[Sign[dolots[ind,r,i]Abs[(muss[[ind,i+1]]-muss[[ind,i]])]-adolots[ind,r,i]Abs[(muss[[ind,i+1]]-muss[[ind,i]])]],{i,1,Length[mids[[ind]]]}];*)
(*ctab=Reap[*)
(*Do[*)
(*ind=12;*)
(*crosslook[ind,r];*)
(*Do[If[crosslook[ind,r][[i]]\[NotEqual] crosslook[ind,r][[i+1]],Sow[{r}]]*)
(*,{i,1,Length[mids[[ind]]]-1}*)
(*]*)
(*,{r,190,300}*)
(*]*)
(*][[2,1]]//Flatten//DeleteDuplicates;*)
(**)*)


(* ::Section::Closed:: *)
(*Auto Plotting*)


(* ::Input:: *)
(*(**)
(*Do[*)
(*SetDirectory[];*)
(*plotname=filenames[[FromDigits[ToExpression[StringCases[scanfilenames[[pind[[j]]]],DigitCharacter]],10]]];*)
(*CreateDirectory["Plots_"<>plotname];*)
(*SetDirectory["Plots_"<>plotname];*)
(**)
(*Export["potential"<>plotname<>".pdf",Rasterize@ListLogPlot[Transpose@{rs[[180;;250]],potscanall[[j]]},ImageSize\[Rule] Large]];*)
(*Export["density"<>plotname<>".pdf",Rasterize@ListLogPlot[Table[{rs[[r]],matters[[pind[[j]],r]]},{r,180,250}],ImageSize\[Rule] Large]];*)
(*(*Export["GDs"<>plotname<>".pdf",Rasterize@ListLogPlot[Transpose@{rs[[180;;250]],gdscanall[[j]]},ImageSize\[Rule] Large]];*)*)
(**)
(*allplot=ListPointPlot3D[Flatten[rscanps[pind[[j]]],1],AxesLabel\[Rule] {"radius (cm)","k (erg)","Im[\[CapitalOmega]] (erg)"},PlotRange\[Rule] All,ImageSize\[Rule]Scaled[0.8],TicksStyle\[Rule] Directive[16],LabelStyle\[Rule] Directive[18],Filling\[Rule] Bottom,ColorFunction\[Rule] "Rainbow",RegionFunction\[Rule]Function[ {x,y,z},z> 10^-23],PlotLegends\[Rule] Automatic,PlotLabel\[Rule] Style[plotname,Bold,Directive[20]],PlotTheme\[Rule] "Scientific"];*)
(*views={{Front,Left,Top}, {Back,Right,Top}, Back,Left,Right,Front};*)
(*Do[*)
(*Export["3Dv"<>ToString[j]<>"_"<>plotname<>".pdf",Rasterize@Show[allplot,PlotRange\[Rule] {{6 10^6,1.4 10^7},{-2 10^-17,4 10^-17},All},ViewPoint\[Rule]views[[k]]]]*)
(*,{k,1,6}];*)
(*SetDirectory[];*)
(*,{j,1,Length[pind]}]*)
(**) *)


(* ::Section::Closed:: *)
(*Debugging, testing, and scratch work*)


(* ::Input:: *)
(*SetDirectory["C:\\Users\\Sam\\rscan3"]*)
(*matfilenames=FileNames["scanmat*"];*)
(*ktarget3=Reap[*)
(*Do[*)
(*testm=Import[matfilenames[[i]]];*)
(*Sow[kx/.Solve[testm[[1,2,2,1,1,1,1]]==0,kx][[1]]]*)
(*,{i,1,Length[matfilenames]}]*)
(*][[2,1]];*)


(* ::Input:: *)
(*SetDirectory["C:\\Users\\Sam\\rscan3"]*)
(*matfilenames=FileNames["scanmat*"]*)
(*ktarget[3]=Reap[*)
(*Do[*)
(*testm=Import[matfilenames[[i]]];*)
(*Sow[kx/.Solve[testm[[1,2,2,1,1,1,1]]==0,kx][[1]]]*)
(*,{i,1,Length[matfilenames]}]*)
(*][[2,1]];*)
(*SetDirectory["C:\\Users\\Sam\\rscan22"]*)
(*matfilenames=FileNames["scanmat*"]*)
(*ktarget[22]=Reap[*)
(*Do[*)
(*testm=Import[matfilenames[[i]]];*)
(*Sow[kx/.Solve[testm[[1,2,2,1,1,1,1]]==0,kx][[1]]]*)
(*,{i,1,Length[matfilenames]}]*)
(*][[2,1]];*)


(* ::Input:: *)
(*ktarget[22]*)


(* ::Input:: *)
(*fSpace[min_,max_,steps_,f_: Log]:=InverseFunction[ConditionalExpression[f[#],min<#<max]&]/@Range[f@min,f@max,(f@max-f@min)/(steps-1)]*)
(*kadapt[ind_]:=Table[Join[fSpace[ktarget[ind][[i]] 10^-1,ktarget[ind][[i]] 10^1,50],-fSpace[ktarget[ind][[i]] 10^-1,ktarget[ind][[i]] 10^1,20]],{i,1,71}];*)
(**)
(**)


(* ::Input:: *)
(*kadapt[22]//Dimensions*)


(* ::Input:: *)
(*NumberLinePlot[{kadapt[22][[20]],ktabr}]*)


(* ::Input:: *)
(*SetDirectory[];*)
(*indss={22};*)
(*Do[*)
(*CreateDirectory["rscan_targeted_adapt"<> ToString[indss[[m]]]];*)
(*SetDirectory["C:\\Users\\Sam\\rscan_targeted_adapt"<>ToString[indss[[m]]]];*)
(*Do[*)
(*scantab=Table[SCalclotsScale[179+i,179+i,20 10^6,indss[[m]],kadapt[indss[[m]]][[i,j]],-1],{j,1,Length[kadapt[indss[[m]]][[1]]]}];*)
(*Export["rscan_adapt"<>ToString[indss[[m]]]<>"r"<> ToString[i+179]<>".m",scantab];*)
(*,{i,1,Length[rtabr]}*)
(*];*)
(*SetDirectory[];*)
(*Print[ToString[m]<>" has completed"]*)
(*,{m,1,Length[indss]}*)
(*];*)
(**)


(* ::Input:: *)
(*SetDirectory["C:\\Users\\Sam\\rscan_targeted_adapt22"]*)
(*rfilenames=FileNames["rscan_adapt22r*.m"];*)
(*rfile=Reap[Do[*)
(*raws=Import[rfilenames[[i]]];*)
(*Sow[raws[[All,2,1,1,1,1]]];*)
(*,{i,1,Length[rfilenames]}*)
(*]*)
(*][[2,1]];*)
(*Export["rscan_adapt22.m",rfile];*)


(* ::Input:: *)
(*rscan22adapt=Import["C:\\Users\\Sam\\rscan_targeted_adapt22\\rscan_adapt22.m"];*)


(* ::Input:: *)
(*rscan22adapt[[1,1]]*)


(* ::Input:: *)
(*radapt22=Table[{rtabr[[i]],kadapt[22][[i,j]],rscan22adapt[[i,j]]},{i,1,71},{j,1,60}];*)


(* ::Input:: *)
(*rfile[[1,1,1,1,71]]//Dimensions*)


(* ::Input:: *)
(*radapt22//Dimensions*)


(* ::Input:: *)
(*rscan22adapt[[1,1]]*)


(* ::Input:: *)
(*Flatten[radapt22,1][[1]]*)


(* ::Input:: *)
(*ListPointPlot3D[Flatten[radapt22,1],AxesLabel-> {"radius (cm)","k (erg)","Im[\[CapitalOmega]] (erg)"},PlotRange-> {All,{-2 10^-17,4 10^-17},All},ImageSize->Scaled[0.6],TicksStyle-> Directive[16],LabelStyle-> Directive[18],Filling-> Bottom,ColorFunction-> "Rainbow",RegionFunction->Function[ {x,y,z},z> 10^-23],PlotLegends-> Automatic,PlotLabel-> Style[plotname,Bold,Directive[20]],PlotTheme-> "Scientific"]*)


(* ::Input:: *)
(*ListPointPlot3D[Flatten[radapt22,1],AxesLabel-> {"radius (cm)","k (erg)","Im[\[CapitalOmega]] (erg)"},PlotRange-> {All,All,All},ImageSize->Scaled[0.6],TicksStyle-> Directive[16],LabelStyle-> Directive[18],Filling-> Bottom,ColorFunction-> "Rainbow",RegionFunction->Function[ {x,y,z},z> 10^-20],PlotLegends-> Automatic,PlotLabel-> Style[plotname,Bold,Directive[20]],PlotTheme-> "Scientific"]*)


(* ::Input:: *)
(*allplot=ListPointPlot3D[Flatten[rscanps[22],1],AxesLabel-> {"radius (cm)","k (erg)","Im[\[CapitalOmega]] (erg)"},PlotRange-> All,ImageSize->Scaled[0.6],TicksStyle-> Directive[16],LabelStyle-> Directive[18],Filling-> Bottom,ColorFunction-> "Rainbow",RegionFunction->Function[ {x,y,z},z> 10^-23],PlotLegends-> Automatic,PlotLabel-> Style[plotname,Bold,Directive[20]],PlotTheme-> "Scientific"]*)


(* ::Input:: *)
(*NumberLinePlot[kadapt[22]]*)
(*NumberLinePlot[ktarget22]*)


(* ::Input:: *)
(*ktarget[22]*)
