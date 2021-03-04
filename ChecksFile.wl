(* ::Package:: *)

(* ::Subsection::Closed:: *)
(*User Initialization*)


(* ::Input::Initialization:: *)
(*Dialog prompt to fix in/out paths.  Can be replaced with direct path input later*)
id=ChoiceDialog[
"It's dangerous to go alone! Who are you and your computer companion?",{"Sam & Cuchulain"-> "Sam", "Sherwood & Ganon"->"Sherwood","Someone new"->"Other"}];
(*Sets the paths based on the machine pick.  The paths are global*)
 Which[id=="Sam",
 inpath="G:\\My Drive\\Physics\\Neutrino Oscillation Research\\Fast Conversions\\lotsadata.tar\\lotsadata\\lotsadata\\";
 outpath="G:\\My Drive\\Physics\\Neutrino Oscillation Research\\Fast Conversions\\stability_data\\";
SetDirectory["C:\\Users\\Sam\\Documents\\GitHub\\stability_tools"];
Get["StabilityPackage`",Path->"."];
SetOptions[EvaluationNotebook[],
DockedCells-> 
Cell[BoxData[ToBoxes[
Grid[{{Item[Style["Stability Tools",FontFamily->"Helvetica",12,Bold],Alignment->Left],Item[ButtonBar[{Style["Refresh Package",10]:>Get["StabilityPackage`",Path->"."],Style["Abort Evaluation",10]:>FrontEndTokenExecute["EvaluatorAbort"],Style["Quit Kernel",10]:>FrontEndTokenExecute["EvaluatorQuit"]}],Alignment->Right]}},ItemSize->{{Scaled[0.3],Scaled[0.7]}}]
]],"DockedCell",Background->LightBlue]
];
 ,
 id=="Sherwood",
 inpath="/mnt/data/SamFlynn/lotsadata/";
 outpath="/mnt/data/SamFlynn/stability_data/";
Get["StabilityPackage`"];
 ,
 id=="Other",
 inpath = SystemDialogInput["Directory",WindowTitle-> "Choose the folder containing the data set"];
 outpath = inpath<>"out\\";
Get["StabilityPackage`"];
 ];
 (*Want to implements toggle grid here to pick between the mass models and time slices.  It's easy, but definitely not a priority*)
(*
filename = "15Msun_50ms_DO";
infile = inpath<>filename<>".h5";
outfolder = outpath<>filename;
*)
(*Note: this used to contain a variable called "out_path". Be careful with the underscores; in Mathematica they are either function inputs or patterns.  If the font turn green, that's why.*)
(*Constants*)






(* ::Subsection::Closed:: *)
(*Regression test for growth rate over a range of wavenumbers. The new results (blue) should match the old results (orange).*)


(*Test on kAdapt. Comparing to old data (below). Generates a plot which should look similar (different exact k values tested*)
OldData={{3.943351801881396`*^-20,0},{6.403122930332929`*^-20,0},{1.0397242072440586`*^-19,6.475185304707499`*^-23},{1.6882799828943322`*^-19,0},{2.7413897654617433`*^-19,1.2805250574405962`*^-19},
{4.451404934206825`*^-19,2.149878747929582`*^-19},{7.228087788874928`*^-19,2.877126337699488`*^-19},{1.173680082038912`*^-18,3.418162146974343`*^-19},{1.9057944164639052`*^-18,3.5031848436345418`*^-19},
{3.0945846431298252`*^-18,2.6089021800081186`*^-19},{5.02491456096483`*^-18,7.0760561598301275`*^-22},{8.159339380505439`*^-18,7.514263464205905`*^-22},{1.3248945493211215`*^-17,7.789552045160011`*^-22},
{2.1513329510655564`*^-17,7.961350592910255`*^-22},{3.493284404190471`*^-17,8.068054415566711`*^-22},{5.672313958895247`*^-17,8.134119405806582`*^-22},{9.210571463829571`*^-17,8.174940729310282`*^-22},
{1.4955911697601837`*^-16,8.200132465051396`*^-22},{2.4285061527926315`*^-16,8.215666277345043`*^-22},{3.9433518018813854`*^-16,8.225240355557669`*^-22},{-3.943351801881396`*^-20,1.2187126603736879`*^-22},
{-6.403122930332929`*^-20,7.691141781968944`*^-22},{-1.0397242072440586`*^-19,2.607457007692258`*^-20},{-1.6882799828943322`*^-19,5.587693993721637`*^-20},{-2.7413897654617433`*^-19,6.185616237315653`*^-20},
{-4.451404934206825`*^-19,5.810765742537203`*^-20},{-7.228087788874928`*^-19,1.6819680254319893`*^-20},{-1.173680082038912`*^-18,1.1963573929349437`*^-21},{-1.9057944164639052`*^-18,1.1670703933149374`*^-21},
{-3.0945846431298252`*^-18,1.0322217805945945`*^-21},{-5.02491456096483`*^-18,9.494282367992625`*^-22},{-8.159339380505439`*^-18,9.001232424984445`*^-22},{-1.3248945493211215`*^-17,8.704722289137772`*^-22},
{-2.1513329510655564`*^-17,8.524818134192752`*^-22},{-3.493284404190471`*^-17,8.415032243261679`*^-22},{-5.672313958895247`*^-17,8.347797496201637`*^-22},{-9.210571463829571`*^-17,8.306532421569712`*^-22},
{-1.4955911697601837`*^-16,8.281172359470442`*^-22},{-2.4285061527926315`*^-16,8.265574524679777`*^-22},{-3.9433518018813854`*^-16,8.255976606410335`*^-22}};
kdebug=Module[{data,ri=200,testE=20,hi=-1,kvar,nstep=20},
file=inpath<>"1D_withV_withPairBrems_DO.h5";
(*SCalcScale[ImportData[inpath<>file<>".h5"],ri,testE,hi,0][[3]]//MatrixForm*)
(*buildkGrid[ImportData[inpath<>file<>".h5"],ri,testE,hi,40]*)
kAdapt[file,ri,ri,testE,hi,nstep,"xflavor"-> False]
]
debugdata=Transpose@{kdebug[[1,All,3]],Table[kdebug[[1,All,4,1]][[i]]//Im//Max,{i,1,40}]};
plot1=ListLogPlot[{debugdata},ImageSize-> Scaled[0.25]]; (*Current data as calculated by kadapt*)
plot2=ListLogPlot[{OldData},ImageSize-> Scaled[0.25]]; (*Old reference data set*)
diffplot=ListPlot[debugdata[[All,2]]-OldData[[All,2]],Joined-> True,ImageSize->Scaled[0.25]]; (*Plot of the difference between OldData and current data*)
rowplot=GraphicsRow[{plot1,plot2,diffplot},Frame-> True]





(* ::Subsection::Closed:: *)
(*Two-beam test. Initialize data with neutrinos moving right and antineutrinos moving left. The real and imaginary parts of the eigenvalues should match the theoretical results from Chakraborty+2016 (Self-induced neutrino flavor conversion without flavor mixing)*)


(*Artifical data set which defines 2 beams only, with 2 angular bins, and have chakraborty asymmetry parameters "a" *)
get2bdata[]:=Module[{S2ba,S2b,data2b},
data2b=
 Association[
"muss"-> {-1,0,1},
"matters"-> 0.,
"Yes"-> 0.,
"mids"-> {-1,1},
"freqs"->{0,2},
"Endensity"-> {{{{0., 0.},{(1+a)/(munits/( h)),0.}}},{{{(1-a)/(munits/(h)),0.},{0.,0.}}}},
 "freqmid"-> {1},
 "munits"-> munits
]
];
get2bdatax[]:=Module[{S2ba,S2b,data2b},
data2b=
 Association[
"muss"-> {-1,0,1},
"matters"-> 0.,
"Yes"-> 0.,
"mids"-> {-1,1},
"freqs"->{0,2},
"Endensity"->{{{{0., 0.},{(1+a)/(munits/( h)),0.}}},{{{(1-a)/(munits/(h)),0.},{0.,0.}}}},
 "freqmid"-> {1},
 "munits"-> munits
]
];
(*builds a 4x4 matrix with the two beam data, then reduces the size of the matrix *)
build2bMatrix[En_,k_]:=Module[{fakeEn,S2ba,S2b},
S2ba=stabilityMatrix[get2bdata[],getEquations[get2bdata[],En,-1.,k,ndensities[get2bdata[],"xflavor"-> False],"xflavor"-> False],"xflavor"-> False];
S2b={{S2ba[[2,2]],S2ba[[2,3]]},{S2ba[[3,2]],S2ba[[3,3]]}};
Return[S2b]
];


(*Chakraborty's analytic expression for the 2x2 case, and a verification test that they're the same within one part in 1000*)
\[CapitalOmega]ch[k_,\[Mu]ch_,a_,\[Omega]_]:=Sort[{2 a \[Mu]ch+Sqrt[(2 a \[Mu]ch)^2+(\[Omega]+k)((\[Omega]+k)-4 \[Mu]ch)],2 a \[Mu]ch-Sqrt[(2 a \[Mu]ch)^2+(\[Omega]+k)((\[Omega]+k)-4 \[Mu]ch)]}];
cm[k_,\[Mu]ch_,w_]:=DiagonalMatrix[{w+k,-w-k,w-k,-w+k}]+2 \[Mu]ch{{l+lb,-lb,-l,0.},{-r,r+rb,0,-rb},{-r,0,r+rb,-rb},{0.,-lb,-l,l+lb}};
cma[k_,\[Mu]ch_,a_,w_]:=cm[k,\[Mu]ch,w]/.{rb-> 0.,l-> 0.,r-> (1+a),lb-> -(1-a)};


(* ::Subsection::Closed:: *)
(*Preliminaries for Dispersion Checks*)


SIpotential[ndens_]:=munits{Diagonal[ndens[[1]]],Diagonal[ndens[[2]]]}


Phi[rank_,mu_,mubar_,cos\[Theta]_]:=Sum[(mu[[i]]-mubar[[i]])cos\[Theta][[i]]^rank,{i,1,Length[cos\[Theta]]}];


IdisShifts[data_,\[Phi]0_,\[Phi]1_,\[CapitalOmega]_,k_]:=Module[{\[CapitalOmega]p,kp,\[Omega],Vmatter},

(* "Shifted" Eigenvalue and k*)
Vmatter = munits data["Yes"] data["matters"]/mp;
\[CapitalOmega]p = N[\[CapitalOmega]-Vmatter-\[Phi]0];
kp = k-\[Phi]1;

Return[{\[CapitalOmega]p,kp}];
]



IdisBottom[\[CapitalOmega]p_,kp_,cos\[Theta]_]:=\[CapitalOmega]p-kp cos\[Theta]


(*Calculates and Returns the nth I for the dispersion check.  Returns a single value of In*)
Idis[data_,cos\[Theta]_,mu_,mubar_,\[CapitalOmega]_,k_,En_,n_]:=Module[{\[Phi]0,\[Phi]1,\[CapitalOmega]p,kp,\[Omega],Vmatter,\[CapitalOmega]minuskpcos\[Theta],result},

(*Defined in Gail's Blue equation 30 and 31 *)
\[Phi]0 = Phi[0,mu,mubar,cos\[Theta]];
\[Phi]1 = Phi[1,mu,mubar,cos\[Theta]];

\[Omega]=\[Omega]EMev[En];
{\[CapitalOmega]p,kp}=IdisShifts[data,\[Phi]0,\[Phi]1,\[CapitalOmega],k];
\[CapitalOmega]minuskpcos\[Theta]=IdisBottom[\[CapitalOmega]p,kp,cos\[Theta]];

(* make sure the denominator is not tiny *)
(*On[Assert];
Do[Assert[(\[CapitalOmega]p-(kp cos\[Theta][[i]]))/(\[CapitalOmega]p+(kp cos\[Theta][[i]]))>= 1,"\[CapitalOmega]p' k'Cos[\[Theta]] percent difference less than 1"],{i,1,Length[cos\[Theta]]}]; 
*)
(*Definition of I from Gail's equation (41)*)
result= Sum[cos\[Theta][[i]]^n (((mu[[i]]-mubar[[i]])\[CapitalOmega]minuskpcos\[Theta][[i]]) + (mu[[i]]+mubar[[i]])\[Omega]) / (\[CapitalOmega]minuskpcos\[Theta][[i]]^2-\[Omega]^2) ,{i,1,Length[cos\[Theta]]}];
Return[result];
];


dispersionCheck[data_,cos\[Theta]_,\[CapitalOmega]_,k_,En_,xflavor_]:=Module[{I0,I1,I2,ndens,mu,mubar},
(* neutrino number densities disguised as SI potentials *)
ndens = ndensities[data,"xflavor"->xflavor];
{mu,mubar}=SIpotential[ndens];

(*The condition is that Equatrion (43), below, should be 0 if the vacuum is*)
I0 = Idis[data,cos\[Theta],mu,mubar,\[CapitalOmega],k,En,0];
I1 = Idis[data,cos\[Theta],mu,mubar,\[CapitalOmega],k,En,1];
I2 = Idis[data,cos\[Theta],mu,mubar,\[CapitalOmega],k,En,2];
Return[Abs[(I0+1.)(I2-1.)-I1^2]]
];


test2bdispersionCheck[k_,En_,atest_,xflavor_]:=Module[{cos\[Theta],I0,I1,I2,\[CapitalOmega],data},
data = get2bdata[]/.a-> atest;
cos\[Theta]=data["mids"];
\[CapitalOmega]=evscale[k,build2bMatrix[En,k]/.a-> 0.,kvar,"output"->"Eigenvalues"][[1]];
Return[dispersionCheck[data,cos\[Theta],\[CapitalOmega],k,En,xflavor]]
];


test4bdispersionCheck[k_,En_,atest_,xflavor_]:=Module[{cos\[Theta],I0,I1,I2,\[CapitalOmega],data,equations},
data = get2bdata[]/.a-> atest;
cos\[Theta]=data["mids"];
equations = getEquations[data,En,-1.,k,ndensities[data,"xflavor"-> xflavor],"xflavor"->xflavor];
\[CapitalOmega]=evscale[k,stabilityMatrix[data,equations,"xflavor"->xflavor],kx,"output"-> "Eigenvalues"][[1]];
Return[dispersionCheck[data,cos\[Theta],\[CapitalOmega],k,En,xflavor]]
];


(* ::Subsection::Closed:: *)
(*Real data dispersion check*)


(*Calculates the dispersion relation based on a hdf file of a stability calculation. Then, returns ordered pairs of the dispersion relation together with the percent difference of \[CapitalOmega]p-kpcos\[Theta]*)
realdatadispcalc[infile_,hdffile_,ri_]:=Module[{data,datasr,pos,testpos,griddata,testk,test\[CapitalOmega]s,evecchecks,dischecks,ops,ins,bottoms,cos\[Theta],sumbottoms,percentdiff,mindiff,\[CapitalOmega]p,kp,\[Phi]0,\[Phi]1,ndens,mu,mubar},
data=ImportData[infile];
datasr=SelectSingleRadius[data,ri];
cos\[Theta]=data["mids"];
griddata=ImportCalcGridData[hdffile];
ops=ImportCalcOptions[hdffile]; (*imports options*)
ins=ImportCalcInputs[hdffile]; (*imports the inputs*)
pos=Position[griddata["ri"],ri]; (*A list of positions in the grid which are at radius ri*)
testpos=pos[[1,1]];(*RandomChoice[pos]; Pulls out a random choice of these positions, i.e. the random k to use*)
testk=griddata["k"][[testpos]];
test\[CapitalOmega]s=griddata["evs_Re"][[testpos]]+I griddata["evs_Im"][[testpos]]; (*Reconstructs the evs*)
(*
Print[testpos];
Print[testk];
Print[test\[CapitalOmega]s];
Print[ins["testE"]];
*)
(*Returns a list of ordered pairs of the maximum component of each of the {neutrino, antineutrino} eigenvectors*)
ndens = ndensities[datasr,"xflavor"->ops["xflavor"]];
{mu,mubar} = SIpotential[ndens];
\[Phi]0 = Phi[0,mu,mubar,cos\[Theta]];
\[Phi]1 = Phi[1,mu,mubar,cos\[Theta]];
{\[CapitalOmega]p,kp}=Transpose@Table[IdisShifts[datasr,\[Phi]0,\[Phi]1,test\[CapitalOmega]s[[i]],testk],{i,1,Length[test\[CapitalOmega]s]}];
bottoms=Table[IdisBottom[\[CapitalOmega]p[[i]],kp[[i]],cos\[Theta]],{i,1,Length[test\[CapitalOmega]s]}]; (*Values of \[CapitalOmega]p-kpcos\[Theta] for each eigenvalue, and each angle.*)
sumbottoms=Table[2 test\[CapitalOmega]s[[i]]- bottoms[[i]],{i,1,Length[test\[CapitalOmega]s]}]; (*value of \[CapitalOmega]p+kpcos[\[Theta]] done via 2\[CapitalOmega]p-(\[CapitalOmega]ps-kpcos\[Theta])*)
percentdiff=Table[bottoms[[i]]/sumbottoms[[i]],{i,1,Length[test\[CapitalOmega]s]}]; (*percent difference; (\[CapitalOmega]p-kpcos\[Theta])/(\[CapitalOmega]p+kpcos\[Theta]) for each \[CapitalOmega]p and angle*)
mindiff=Table[Min[percentdiff[[i]]//Re],{i,1,Length[percentdiff[[1]]]}]; (*For each eigenvalue, takes the minimum percent difference across the 10 angular bins*)
(*Performs a dispersion checks for each of the 20 eigenvalues, which should be 0. Collects all of the dispchecks*)
dischecks=Table[
	Abs[dispersionCheck[datasr,cos\[Theta],test\[CapitalOmega]s[[i]],testk,ins["testE"],ops["xflavor"]]],
	{i,1,2(*Length[griddata["evs_Re"][[testpos]]]*)}];
Return[Transpose@{dischecks,mindiff[[1;;2]],test\[CapitalOmega]s[[1;;2]]}];
]; (*Returns a list of ordered pairs for each of the 20 eigenvalues tests where the first component is the result of the dispersion check, and the second contains a "" or error message*)


(*Calls realdatadispcalc and performs the logic for the check on whether the disp checks are passing*)
(*Returns a list of true/false for whether the test passed*)
realdatadispersioncheck[infile_,hdffile_,ri_]:=Module[{dispouts,checks,ans,dispcond},
dispouts=realdatadispcalc[infile,hdffile,ri];
(*Check each pair, returns true if the disp passes OR the percent difference is very small.*)
(*This collects whether or not the check should pass, conditionally or naturally*)
checks=Table[dispouts[[j,1]]<10^(-3) \[Or] dispouts[[j,2]]< 10^(-10), {j,1,Length[dispouts]}];
ans=Apply[And,checks];	(*Checks whether all elements of checks are true*)
(*Discond is using for determining the condition by which the check passes.*)
dispcond=Reap[
	Do[
		Sow[{dispouts[[j,3]],
				Which[dispouts[[j,1]]<10^-3, Style["Natural Pass",Darker[Green,0.4]],
				dispouts[[j,2]]< 10^-10 \[And] dispouts[[j,1]]>10^-3, Style["Conditional Pass",Darker[Yellow,0.4]],
				True,Style["Fail",Darker[Red,0.4]](*Close Which*)
			} (*Close ordered pair in sow*)
		] (*Close sow*)
	,{j,1,Length[dispouts]}
	] 
][[2,1]];
(*Print a grid with first row of eignevalues and second row pass conditions.*)
Grid[{{hdffile,SpanFromLeft},dispcond[[All,1]],dispcond[[All,2]]},Frame-> All]//Print;
Return[ans];
];


(* ::Subsection::Closed:: *)
(*Check that ellipse construction results in the correct moments given hand-chosen moments*)


(*Ellipse fits to artifical parameters*)
eCheckArtificial[]:= Module[{pans,ans,mom,igs},
mom={1., 1. 10^-8, N[1/3]};
igs=Apply[getInitialGuess,mom];
pans=Apply[eSimpFitToMoments,Join[mom,{igs}]];
ans=Apply[ellipseparaerrors,Join[pans,mom]];
Return[ans]
];

ellipsefitrealdatacheck[]:=Module[{file,moms,igs,paras,errs},
file=inpath<>"4timesHigh_1D_withV_withPairBrems_MC_moments.h5";
moms=getMoments[file,220,1];
igs=Apply[getInitialGuess,moms];
paras=Apply[eSimpFitToMoments,Join[moms,{igs}]];
errs=Apply[ellipseparaerrors,Join[paras,moms]];
Return[errs]
];
(*Imports real CSSN data and then calls ellipse fit errors for the tests file*)


(* ::Subsection::Closed:: *)
(*xflavor tests*)


ncheck[file_,ri_]:=Module[{ncheckdata},
ncheckdata=SelectSingleRadius[ImportData[file],ri];
Return[siPotential[ndensities[ncheckdata,"xflavor"->True]]>siPotential[ndensities[ncheckdata,"xflavor"-> False],"xflavor"->False]];
];

netleptoncheck[file_,ri_]:=Module[{ncheckdata},
ncheckdata=SelectSingleRadius[ImportData[file],ri];
Return[((Tr[ndensities[ncheckdata,"xflavor"-> True][[1]]]+Tr[ndensities[ncheckdata,"xflavor"-> True][[3]]])
-(Tr[ndensities[ncheckdata,"xflavor"-> True][[2]]]+Tr[ndensities[ncheckdata,"xflavor"-> True][[3]]]))==
((Tr[ndensities[ncheckdata,"xflavor"->False][[1]]]+Tr[ndensities[ncheckdata,"xflavor"-> False][[3]]])
-(Tr[ndensities[ncheckdata,"xflavor"-> False][[2]]]+Tr[ndensities[ncheckdata,"xflavor"-> False][[3]]]))]
];


(* ::Subsection::Closed:: *)
(*Morinaga Plot*)


Options[MorinagaTestPlot]={"xflavor"->False};
MorinagaTestPlot[ri_,\[Omega]_,res_,OptionsPattern[]]:=Module[{data,datain,cos\[Theta],datainsr,costheta,datapos,evspos,evsposre,kspos,mu,mubar,\[Phi]0,\[Phi]1,Vmatter,\[CapitalOmega]p,kp,pts,ptsre,p1,g1},
exportkadapt[kAdapt[inpath<>"112Msun_100ms_DO.h5",ri,ri,\[Omega],-1,res,"xflavor"-> OptionValue["xflavor"],"ktarget"->4.22 10^-19,"krange"-> {0.6,1.4}],"MorinagaTestPlot"];
data=ImportCalcGridData["MorinagaTestPlot.h5"];
datain=ImportData[inpath<>"112Msun_100ms_DO.h5"];
datainsr=SelectSingleRadius[datain,ri];
cos\[Theta]=datain["mids"];
datapos=Position[data["ri"],ri];
evspos=Table[Extract[data["evs_Im"],datapos[[i]]]//Max,{i,1,Length[datapos]}];
evsposre=Table[Extract[data["evs_Re"],datapos[[i]]]//Max,{i,1,Length[datapos]}];
kspos=Table[Extract[data["k"],datapos[[i]]],{i,1,Length[datapos]}];
{mu,mubar}=SIpotential[ndensities[datainsr,"xflavor"->OptionValue["xflavor"]]];
\[Phi]0 = Phi[0,mu,mubar,cos\[Theta]];
\[Phi]1 = Phi[1,mu,mubar,cos\[Theta]];
{\[CapitalOmega]p,kp}=Transpose@Table[IdisShifts[datainsr,\[Phi]0,\[Phi]1,evsposre[[i]],kspos[[i]]],{i,1,Length[evspos]}];

pts=Transpose@{1/(hbar c 10^-4) kp,1/(hbar c 10^-4) evspos};
ptsre=Transpose@{1/(hbar c 10^-4) kp,1/(hbar c 10^-4) \[CapitalOmega]p};
p1=ListPlot[pts,PlotRange-> {{-30,60},All},Filling-> Bottom,AxesLabel-> {Style["k' (\!\(\*SuperscriptBox[\(10\), \(-4\)]\) \!\(\*SuperscriptBox[\(cm\), \(-1\)]\))",FontSize-> 14,Bold],Style["Im[\[CapitalOmega]] (\!\(\*SuperscriptBox[\(10\), \(-4\)]\) \!\(\*SuperscriptBox[\(cm\), \(-1\)]\))",FontSize-> 14,Bold]},ImageSize-> Scaled[0.45],Frame-> False,PlotStyle-> Black];
Return[p1]
];


(* ::Subsection::Closed:: *)
(*Ellipse Fit a DO File*)


(*Convert a DO file to moments and export it *)
exportDOasMoments[DOtoMoments[inpath<>"112Msun_100ms_DO.h5",264,264,0],"testfiledomom"];
(*Ellipse fit that data and export it*)
exportelipdata["testfileelipdata",getelipdata["testfiledomom.h5",1,1,0]];
(*Search for unstable modes and export it*)
exportkadapt[kAdapt["morelipfit1.h5",1,1,Infinity,-1.,10,"xflavor"-> False],"testfilekadapt"];
(*Make a test plot of that dataset*)
tdata=ImportCalcGridData["testfilekadapt.h5"];
(*Check that it returned SOMETHING*)
{tdata["radius"],tdata["k"],tdata["evs_Im"]};


(* ::Subsection::Closed:: *)
(*Test Report*)


(* Generate the data for the dispersion check*)
dispersionCheckRi=250;
(*Note, I changed this to take 2 ksteps, as otherwise there's an error from the way the formula for the log spacing is calculated (1/0 error)*)
outevs = kAdapt[inpath <> "112Msun_100ms_DO.h5", dispersionCheckRi, dispersionCheckRi, Infinity, -1., 2, "xflavor" -> False]; 
exportkadapt[outevs,"112Msun_100ms_r200_r300_now_nox_test" ]
outevs = kAdapt[inpath <> "112Msun_100ms_DO.h5", dispersionCheckRi, dispersionCheckRi, 20., -1., 2,"xflavor" -> False];
exportkadapt[outevs,"112Msun_100ms_r200_r300_nox_test"]


Timing[tr=TestReport["testfiles.wlt"];]
Table[tr["TestResults"][i],{i,1,10}]//MatrixForm









