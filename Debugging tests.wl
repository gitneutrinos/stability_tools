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


SIpotential[ndens_]:=munits{Diagonal[ndens[[1]]-ndens[[3]]],Diagonal[ndens[[2]]-ndens[[3]]]}


Phi[rank_,mu_,mubar_,cos\[Theta]_]:=Sum[(mu[[i]]-mubar[[i]])cos\[Theta][[i]]^rank,{i,1,Length[cos\[Theta]]}];


IdisShifts[data_,\[Phi]0_,\[Phi]1_,\[CapitalOmega]_,k_]:=Module[{\[CapitalOmega]p,kp,\[Omega],Vmatter},

(* "Shifted" Eigenvalue and k*)
Vmatter = 0.;(* munits data["Yes"] data["matters"]/mp; *)
\[CapitalOmega]p = N[\[CapitalOmega]-Vmatter-\[Phi]0];
kp = k-\[Phi]1;

Return[{\[CapitalOmega]p,kp}];
]



IdisBottom[\[CapitalOmega]p_,kp_,cos\[Theta]_]:=\[CapitalOmega]p-kp cos\[Theta]


(*Calculates and Returns the nth I for the dispersion check.  Returns a single value of In*)
Options[Idis]={"hier"-> 1.};
Idis[data_,cos\[Theta]_,mu_,mubar_,\[CapitalOmega]_,k_,En_,n_,OptionsPattern[]]:=Module[{\[Phi]0,\[Phi]1,\[CapitalOmega]p,kp,\[Omega],Vmatter,\[CapitalOmega]minuskpcos\[Theta],result},

(*Defined in Gail's Blue equation 30 and 31 *)
\[Phi]0 = Phi[0,mu,mubar,cos\[Theta]];
\[Phi]1 = Phi[1,mu,mubar,cos\[Theta]];

\[Omega]=OptionValue["hier"]\[Omega]EMev[En];
{\[CapitalOmega]p,kp}=IdisShifts[data,\[Phi]0,\[Phi]1,\[CapitalOmega],k];
\[CapitalOmega]minuskpcos\[Theta]=IdisBottom[\[CapitalOmega]p,kp,cos\[Theta]];

(* make sure the denominator is not tiny *)
(*On[Assert];
Do[Assert[(\[CapitalOmega]p-(kp cos\[Theta][[i]]))/(\[CapitalOmega]p+(kp cos\[Theta][[i]]))>= 1,"\[CapitalOmega]p' k'Cos[\[Theta]] percent difference less than 1"],{i,1,Length[cos\[Theta]]}]; 
*)
(*Definition of I from Gail's equation (41)*)
result= Sum[cos\[Theta][[i]]^n (((mu[[i]]-mubar[[i]])\[CapitalOmega]minuskpcos\[Theta][[i]]) - (mu[[i]]+mubar[[i]])\[Omega]) / (\[CapitalOmega]minuskpcos\[Theta][[i]]^2-\[Omega]^2) ,{i,1,Length[cos\[Theta]]}];
Return[result];
];


Options[dispersionCheck]={"hier"-> -1.};
dispersionCheck[data_,cos\[Theta]_,\[CapitalOmega]_,k_,En_,xflavor_,OptionsPattern[]]:=Module[{I0,I1,I2,ndens,mu,mubar},
(* neutrino number densities disguised as SI potentials *)
ndens = ndensities[data,"xflavor"->xflavor];
{mu,mubar}=SIpotential[ndens];

(*The condition is that Equatrion (43), below, should be 0 if the vacuum is*)
I0 = Idis[data,cos\[Theta],mu,mubar,\[CapitalOmega],k,En,0,"hier"->OptionValue["hier"]];
I1 = Idis[data,cos\[Theta],mu,mubar,\[CapitalOmega],k,En,1,"hier"->OptionValue["hier"]];
I2 = Idis[data,cos\[Theta],mu,mubar,\[CapitalOmega],k,En,2,"hier"->OptionValue["hier"]];
Return[Abs[(I0+1.)(I2-1.)-I1^2]]
];


test2bdispersionCheck[k_,En_,atest_,xflavor_]:=Module[{cos\[Theta],I0,I1,I2,\[CapitalOmega],data},
data = get2bdata[]/.a-> atest;
cos\[Theta]=data["mids"];
\[CapitalOmega]=evscale[k,build2bMatrix[En,k]/.a-> atest,kvar,"output"->"Eigenvalues"][[1]];
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
realdatadispcalc[infile_,hdffile_,ri_]:=Module[{data,datasr,pos,\[Omega],testpos,griddata,testk,test\[CapitalOmega]s,evecchecks,dischecks,ops,ins,bottoms,cos\[Theta],sumbottoms,percentdiff,mindiff,\[CapitalOmega]p,kp,\[Phi]0,\[Phi]1,ndens,mu,mubar},
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
\[Omega]=ins[[5]]*\[Omega]EMev[ins[[4]]];
(*Returns a list of ordered pairs of the maximum component of each of the {neutrino, antineutrino} eigenvectors*)
ndens = ndensities[datasr,"xflavor"->ops["xflavor"]];
{mu,mubar} = SIpotential[ndens];
\[Phi]0 = Phi[0,mu,mubar,cos\[Theta]];
\[Phi]1 = Phi[1,mu,mubar,cos\[Theta]];
{\[CapitalOmega]p,kp}=Transpose@Table[IdisShifts[datasr,\[Phi]0,\[Phi]1,test\[CapitalOmega]s[[i]],testk],{i,1,Length[test\[CapitalOmega]s]}];
bottoms=Table[IdisBottom[\[CapitalOmega]p[[i]],kp[[i]],cos\[Theta]]^2-\[Omega]^2,{i,1,Length[test\[CapitalOmega]s]}]; (*Values of \[CapitalOmega]p-kpcos\[Theta] for each eigenvalue, and each angle.*)
sumbottoms=Table[(2 \[CapitalOmega]p[[i]]-IdisBottom[\[CapitalOmega]p[[i]],kp[[i]],cos\[Theta]])^2+\[Omega]^2,{i,1,Length[test\[CapitalOmega]s]}]; 
percentdiff=Table[bottoms[[i]]/sumbottoms[[i]],{i,1,Length[test\[CapitalOmega]s]}]; 
mindiff=Table[Min[percentdiff[[i]]//Abs],{i,1,Length[percentdiff]}]; (*For each eigenvalue, takes the minimum percent difference across the 10 angular bins*)
(*Performs a dispersion checks for each of the 20 eigenvalues, which should be 0. Collects all of the dispchecks*)
dischecks=Table[
	Abs[dispersionCheck[datasr,cos\[Theta],test\[CapitalOmega]s[[i]],testk,ins["testE"],ops["xflavor"],"hier"-> ins[[5]]]],
	{i,1,Length[griddata["evs_Re"][[testpos]]]}];
Return[Transpose@{dischecks,mindiff,test\[CapitalOmega]s}];
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
				True,Style["Fail",Darker[Red,0.4]]
				](*Close Which*)
			} (*Close ordered pair in sow*)
		] (*Close sow*)
	,{j,1,Length[dispouts]}
	] 
][[2,1]];
(*Print a grid with first row of eignevalues and second row pass conditions.*)
Grid[Prepend[Prepend[Table[{dispcond[[i,1]],dispcond[[i,2]],dispouts[[i,1]],dispouts[[i,2]]},{i,1,Length[dispcond]}],{"\[CapitalOmega]","Result","Disp. Relation","Denominator Percent Diff"}],{hdffile,SpanFromLeft,SpanFromLeft,SpanFromLeft}],Frame-> All]//Print;
(*Grid[{{hdffile,SpanFromLeft},dispcond[[All,1]],dispcond[[All,2]]},Frame-> All]//Print;*)
Return[ans];
];


Options[realdatadispcheck2]={"species"-> 1,"check"->1};
realdatadispcheck2[infile_,hdffile_,ri_,mi_,OptionsPattern[]]:=Module[
{data,datasr,pos,\[Omega],testpos,griddata,testk,test\[CapitalOmega]s,evecchecks,dischecks,ops,ins,bottoms,cos\[Theta],\[CapitalOmega]p,kp,\[Phi]0,\[Phi]1,ndens,mu,mubar,a0,a1,a0b,a1b,testAs,testAbs,lhs,rhs,Ve,testAsn,testAbsn,lvec,rvec,dp,tab,nA,nAb,vsi},
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
testAs=griddata["evecs_nu_Re"][[testpos]]+I griddata["evecs_nu_Im"][[testpos]];
testAbs=griddata["evecs_nubar_Re"][[testpos]]+I griddata["evecs_nubar_Im"][[testpos]];
\[Omega]=ins[[5]]*\[Omega]EMev[ins[[4]]];
(*Returns a list of ordered pairs of the maximum component of each of the {neutrino, antineutrino} eigenvectors*)
ndens = ndensities[datasr,"xflavor"->ops["xflavor"]];
{mu,mubar} = SIpotential[ndens];
\[Phi]0 = Phi[0,mu,mubar,cos\[Theta]];
\[Phi]1 = Phi[1,mu,mubar,cos\[Theta]];
{\[CapitalOmega]p,kp}=Transpose@Table[IdisShifts[datasr,\[Phi]0,\[Phi]1,test\[CapitalOmega]s[[i]],testk],{i,1,Length[test\[CapitalOmega]s]}];
a0=Sum[mu[[j]] testAs[[mi,j]],{j,1,Length[cos\[Theta]]}];
a1=Sum[mu[[j]] cos\[Theta][[j]] testAs[[mi,j]],{j,1,Length[cos\[Theta]]}];
a0b=Sum[mubar[[j]] testAbs[[mi,j]],{j,1,Length[cos\[Theta]]}];
a1b=Sum[mubar[[j]] cos\[Theta][[j]] testAbs[[mi,j]],{j,1,Length[cos\[Theta]]}];
Ve =0.;(* munits datasr["Yes"] datasr["matters"]/mp;*)
testAsn=Normalize[testAs[[mi]]];
testAbsn=Normalize[testAbs[[mi]]];
Which[OptionValue["check"]==1, (*This is an explcit check of equation 37.5 for a single eigenmode mi and one directional bin*)

Which[OptionValue["species"]==1,
lhs=test\[CapitalOmega]s[[mi]] testAs[[mi,1]];
rhs=(-\[Omega]+Ve+\[Phi]0+cos\[Theta][[1]](testk-\[Phi]1))testAs[[mi,1]]-Sum[(1-cos\[Theta][[1]]cos\[Theta][[j]])((mu[[j]]testAs[[mi,j]])-(mubar[[j]] testAbs[[mi,j]])),{j,1,Length[cos\[Theta]]}];
Return[{Between[Abs[lhs/rhs],{0.9999,1.0001}],lhs,rhs}];
,
OptionValue["species"]==2,
lhs=test\[CapitalOmega]s[[mi]] testAbs[[mi,1]];
rhs=(\[Omega]+Ve+\[Phi]0+cos\[Theta][[1]](testk-\[Phi]1))testAbs[[mi,1]]-Sum[(1-cos\[Theta][[1]]cos\[Theta][[j]])((mu[[j]]testAs[[mi,j]])-(mubar[[j]] testAbs[[mi,j]])),{j,1,Length[cos\[Theta]]}];
Return[{Between[Abs[lhs/rhs],{0.9999,1.0001}],lhs,rhs}];
];
,
OptionValue["check"]==2, (*This is an 'intermediary step' towards the dispersion relations. Equations 52 and 53*)
Which[OptionValue["species"]==1,
lhs=testAs[[mi,1]];
rhs=-(((a0-a0b)-cos\[Theta][[1]](a1-a1b))/(\[CapitalOmega]p[[mi]]-(cos\[Theta][[1]] kp[[1]])+\[Omega]));
Return[{Between[Abs[lhs/rhs],{0.9999,1.0001}],lhs,rhs}];
,
OptionValue["species"]==2 ,
lhs=testAbs[[mi,1]];
rhs=-(((a0-a0b)-cos\[Theta][[1]](a1-a1b))/(\[CapitalOmega]p[[mi]]-(cos\[Theta][[1]] kp[[1]])-\[Omega]));
Return[{Between[Abs[lhs/rhs],{0.9999,1.0001}],lhs,rhs}];
];
,
OptionValue["check"]==3, (*This is a check for subtractive cancellation*)
Grid[{{"\[CapitalOmega]",test\[CapitalOmega]s[[mi]],"\!\(\*SubscriptBox[\(A\), \(1\)]\)",testAs[[mi,1]],"\!\(\*SubscriptBox[\(Ab\), \(1\)]\)",testAbs[[mi,1]],"Eqn 37.5 Sum", Sum[(1-cos\[Theta][[1]]cos\[Theta][[j]])((mu[[j]]testAs[[mi,j]])-(mubar[[j]] testAbs[[mi,j]])),{j,1,Length[cos\[Theta]]}]},
{"a0",Sum[mu[[j]] testAs[[mi,j]],{j,1,Length[cos\[Theta]]}],"a0 terms",Table[mu[[j]] testAs[[mi,j]],{j,1,Length[cos\[Theta]]}]},
{"a1",Sum[mu[[j]]cos\[Theta][[1]] testAs[[mi,j]],{j,1,Length[cos\[Theta]]}],"a1 terms",Table[mu[[j]]cos\[Theta][[1]] testAs[[mi,j]],{j,1,Length[cos\[Theta]]}]},
{"a0b",Sum[mubar[[j]] testAbs[[mi,j]],{j,1,Length[cos\[Theta]]}],"a0b terms",Table[mubar[[j]] testAbs[[mi,j]],{j,1,Length[cos\[Theta]]}]},
{"a1b",Sum[mubar[[j]] cos\[Theta][[1]] testAbs[[mi,j]],{j,1,Length[cos\[Theta]]}],"a1b terms",Table[mubar[[j]] cos\[Theta][[1]] testAbs[[mi,j]],{j,1,Length[cos\[Theta]]}]}}
,Frame-> All]//Print
,
OptionValue["check"]==4, (*This is check of equation 37.5 for all values of cos\[Theta]*)

Grid[Prepend[Prepend[Table[{cos\[Theta][[l]],test\[CapitalOmega]s[[mi]] testAsn[[l]],
(-\[Omega]+Ve+\[Phi]0+cos\[Theta][[l]](testk-\[Phi]1))testAsn[[l]]-Sum[(1-cos\[Theta][[l]]cos\[Theta][[j]])((mu[[j]]testAsn[[j]])-(mubar[[j]] testAbsn[[j]])),{j,1,Length[cos\[Theta]]}],
testAsn[[l]],testAbsn[[l]]}
,{l,1,Length[cos\[Theta]]}],{"cos\[Theta]","lhs","rhs","\!\(\*SubscriptBox[\(A\), \(\[Theta]\)]\)","\!\(\*SubscriptBox[\(Ab\), \(\[Theta]\)]\)"}],{"\[CapitalOmega]=",test\[CapitalOmega]s[[mi]],SpanFromLeft,SpanFromLeft,SpanFromLeft}],Frame-> All]//Print
,
OptionValue["check"]==5, (*This is a better way of checking that check #4 is true for all angles by taking a vector dot product for a single mi*)
lvec=Table[test\[CapitalOmega]s[[mi]] testAsn[[l]],{l,1,Length[cos\[Theta]]}];
rvec=Table[(-\[Omega]+Ve+\[Phi]0+cos\[Theta][[l]](testk-\[Phi]1))testAsn[[l]]-Sum[(1-cos\[Theta][[l]]cos\[Theta][[j]])((mu[[j]]testAsn[[j]])-(mubar[[j]] testAbsn[[j]])),{j,1,Length[cos\[Theta]]}],{l,1,Length[cos\[Theta]]}];
Return[(lvec . rvec)/(lvec . lvec)];
,
OptionValue["check"]==6, (*This is a check on the dot product for all mi*)
nA[k_]:=Normalize[testAs[[k]]];
nAb[k_]:=Normalize[testAbs[[k]]];
lvec[k_]:=Table[test\[CapitalOmega]s[[k]] nA[k][[l]],{l,1,Length[cos\[Theta]]}];
rvec[k_]:=Table[(-\[Omega]+Ve+\[Phi]0+cos\[Theta][[l]](testk-\[Phi]1))nA[k][[l]]-Sum[(1-cos\[Theta][[l]]cos\[Theta][[j]])((mu[[j]]nA[k][[j]])-(mubar[[j]] nAb[k][[j]])),{j,1,Length[cos\[Theta]]}],{l,1,Length[cos\[Theta]]}];
dp[k_]:=(lvec[k] . rvec[k])/(lvec[k] . lvec[k]);
vsi=Sum[mu[[j]]-mubar[[j]],{j,1,Length[cos\[Theta]]}];
tab=Table[{\[Omega],vsi,test\[CapitalOmega]s[[i]],dp[i]},{i,1,Length[test\[CapitalOmega]s]}];
Grid[Prepend[tab,{"\[Omega]","\!\(\*SubscriptBox[\(V\), \(si\)]\)(\!\(\*SubscriptBox[\(a\), \(0\)]\)-\!\(\*SubscriptBox[\(ab\), \(0\)]\))","\[CapitalOmega]","dot product"}],Frame-> All]//Print
,
OptionValue["check"]==7, (*This prints some values for debugging*)
Return[{testk,\[Phi]0,\[Phi]1,testAs,testAbs,\[Omega],test\[CapitalOmega]s,mu,mubar}]
];
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


Timing[tr=TestReport["testfiles.wlt"];]
Table[tr["TestResults"][i],{i,1,10}]//MatrixForm


(* ::Section::Closed:: *)
(*Possible Code Changes to Watch Out For*)


(* ::Item::Closed:: *)
(*Ve : Is the Matter Potential on or off? *)


(* ::Subitem:: *)
(*  Check the StabilityTools package, IdisShifts, and realdatadispcheck2*)


(* ::Item::Closed:: *)
(*Is there an "extra" factor of 1/2 on \[Omega]?*)


(* ::Subitem:: *)
(*Check Idis, realdatadispcalc, and realdatadispcheck 2*)


(* ::Item::Closed:: *)
(*Is the hierarchy consistent?*)


(* ::Subitem:: *)
(*Check that the input hierarchy is accounted for. *)


(* ::Item:: *)
(*Check the neutrino energy used to generate the file*)
(* *)


(* ::Section::Closed:: *)
(*Generate Files for Checks on Results and/or Code-Paper Agreement*)



(* Generate the data for the dispersion check*)
dispersionCheckRi=250;
(*Note, I changed this to take 2 ksteps, as otherwise there's an error from the way the formula for the log spacing is calculated (1/0 error)*)
outevs = kAdapt[inpath <> "112Msun_100ms_DO.h5", dispersionCheckRi, dispersionCheckRi, Infinity, -1., 2, "xflavor" -> True]; 
exportkadapt[outevs,"112Msun_100ms_r200_r300_now_nox_test" ];
outevs = kAdapt[inpath <> "112Msun_100ms_DO.h5", dispersionCheckRi, dispersionCheckRi, 20., -1., 2,"xflavor" -> True];
exportkadapt[outevs,"112Msun_100ms_r200_r300_nox_test"];

(*Should one want to re-run the hdf files that are used to check 15msun datasets, these can be run to get "fresh" files, or commented out*)
outevs = kAdapt["G:\\.shortcut-targets-by-id\\1tib_pcM4sTpxj42Hb2g-kKJKo1pG_Qeq\\SamFlynn_fast_flavor\\40bin\\DO_data\\15Msun_400ms_DO_40_minerbo.h5", 255, 255, Infinity, -1., 4000, "xflavor" -> True]; 
exportkadapt[outevs,"15Msun_400ms_r255_minerbo_now_test" ]
outevs = kAdapt["G:\\.shortcut-targets-by-id\\1tib_pcM4sTpxj42Hb2g-kKJKo1pG_Qeq\\SamFlynn_fast_flavor\\40bin\\DO_data\\15Msun_400ms_DO_40_minerbo.h5", 255, 255, 20., -1., 4000, "xflavor" -> True]; 
exportkadapt[outevs,"15Msun_400ms_r255_minerbo_w/w_test" ]
outevs = kAdapt["G:\\.shortcut-targets-by-id\\1tib_pcM4sTpxj42Hb2g-kKJKo1pG_Qeq\\SamFlynn_fast_flavor\\40bin\\DO_data\\15Msun_400ms_DO_40.h5", 255, 255, Infinity, -1., 4000, "xflavor" -> True]; 
exportkadapt[outevs,"15Msun_400ms_r255_DO_now_test" ]
outevs = kAdapt["G:\\.shortcut-targets-by-id\\1tib_pcM4sTpxj42Hb2g-kKJKo1pG_Qeq\\SamFlynn_fast_flavor\\40bin\\DO_data\\15Msun_400ms_DO_40.h5", 255, 255, 20., -1., 4000, "xflavor" -> True]; 
exportkadapt[outevs,"15Msun_400ms_r255_DO_w/w_test" ]



(* ::Section::Closed:: *)
(*11.2 Msun_100ms*)


file11="G:\\My Drive\\Physics\\Neutrino Oscillation Research\\Fast Conversions\\lotsadata.tar\\lotsadata\\lotsadata\\112Msun_100ms_DO.h5";
hdf11now="C:\\Users\\Sam\\Documents\\GitHub\\stability_tools\\112Msun_100ms_r200_r300_now_nox_test.h5";
hdf11ww="C:\\Users\\Sam\\Documents\\GitHub\\stability_tools\\112Msun_100ms_r200_r300_nox_test.h5";
r11=250;


(* ::Subsection::Closed:: *)
(*\[Omega] = 0*)


(* ::Subsubsection::Closed:: *)
(*'Dispersion Relation' Check: (Subscript[I, 0]+1)(Subscript[I, 2]-1)-Subscript[I, 1]^2=0*)


realdatadispersioncheck[file11,hdf11now,r11]


(* ::Subsubsection::Closed:: *)
(*Check via Dot Products: Subscript[\[CapitalOmega]A, m]=(-\[Omega] + Ve +Subscript[\[Phi], 0]+Subscript[\[Mu], m](k-Subscript[\[Phi], 1]))Subscript[A, m]-Subscript[\[CapitalSigma], j](1-Subscript[\[Mu], m] Subscript[\[Mu], j])(Subscript[n, j] Subscript[A, j]-Subscript[nb, j] Subscript[Ab, j])*)


(* ::Text:: *)
(*The dot product checks for all direction components ' m' for a single eigenvalue eigenvector pair . Passes if the dot product is 1.*)


realdatadispcheck2[file11,hdf11now,r11,1,"check"-> 6]


(* ::Subsubsection::Closed:: *)
(*Subtractive cancellations check*)


(* ::Text:: *)
(*There is no issue of subtractive cancellation if the sum of equation 37.5 (final column) is less than a' s*)


Do[realdatadispcheck2[file11,hdf11now,r11,mi,"check"-> 3],{mi,1,20}]


(* ::Subsection::Closed:: *)
(*\[Omega] != 0*)


(* ::Subsubsection::Closed:: *)
(*'Dispersion Relation' Check: (Subscript[I, 0]+1)(Subscript[I, 2]-1)-Subscript[I, 1]^2=0*)


realdatadispersioncheck[file11,hdf11ww,r11]


(* ::Subsubsection::Closed:: *)
(*Check via Dot Products: Subscript[\[CapitalOmega]A, m]=(-\[Omega] + Ve +Subscript[\[Phi], 0]+Subscript[\[Mu], m](k-Subscript[\[Phi], 1]))Subscript[A, m]-Subscript[\[CapitalSigma], j](1-Subscript[\[Mu], m] Subscript[\[Mu], j])(Subscript[n, j] Subscript[A, j]-Subscript[nb, j] Subscript[Ab, j])*)


(* ::Text:: *)
(*The dot product checks for all direction components ' m' for a single eigenvalue eigenvector pair . Passes if the dot product is 1.*)


realdatadispcheck2[file11,hdf11ww,r11,1,"check"-> 6]


(* ::Subsubsection::Closed:: *)
(*Subtractive cancellations check*)


(* ::Text:: *)
(*There is no issue of subtractive cancellation if the sum of equation 37.5 (final column) is less than a' s*)


Do[realdatadispcheck2[file11,hdf11ww,r11,mi,"check"-> 3],{mi,1,20}]


(* ::Section::Closed:: *)
(*15 Msun_400ms (Minerbo)*)


file15="G:\\.shortcut-targets-by-id\\1tib_pcM4sTpxj42Hb2g-kKJKo1pG_Qeq\\SamFlynn_fast_flavor\\40bin\\DO_data\\15Msun_400ms_DO_40_minerbo.h5";
hdf15mnow="C:\\Users\\Sam\\Documents\\GitHub\\stability_tools\\15Msun_400ms_r255_minerbo_now_test.h5";
hdf15mww="C:\\Users\\Sam\\Documents\\GitHub\\stability_tools\\15Msun_400ms_r255_minerbo_w/w_test.h5";
r15=255;


(* ::Subsection::Closed:: *)
(*\[Omega] = 0*)


(* ::Subsubsection::Closed:: *)
(*'Dispersion Relation' Check: (Subscript[I, 0]+1)(Subscript[I, 2]-1)-Subscript[I, 1]^2=0*)


realdatadispersioncheck[file15,hdf15mnow,r15]


(* ::Subsubsection::Closed:: *)
(*Check via Dot Products: Subscript[\[CapitalOmega]A, m]=(-\[Omega] + Ve +Subscript[\[Phi], 0]+Subscript[\[Mu], m](k-Subscript[\[Phi], 1]))Subscript[A, m]-Subscript[\[CapitalSigma], j](1-Subscript[\[Mu], m] Subscript[\[Mu], j])(Subscript[n, j] Subscript[A, j]-Subscript[nb, j] Subscript[Ab, j])*)


(* ::Text:: *)
(*The dot product checks for all direction components ' m' for a single eigenvalue eigenvector pair . Passes if the dot product is 1.*)


realdatadispcheck2[file15,hdf15mnow,r15,1,"check"-> 6]


(* ::Subsubsection::Closed:: *)
(*Subtractive cancellations check*)


(* ::Text:: *)
(*There is no issue of subtractive cancellation if the sum of equation 37.5 (final column) is less than a' s*)


Do[realdatadispcheck2[file15,hdf15mnow,r15,mi,"check"-> 3],{mi,1,20}]


(* ::Subsection::Closed:: *)
(*\[Omega] != 0*)


(* ::Subsubsection::Closed:: *)
(*'Dispersion Relation' Check: (Subscript[I, 0]+1)(Subscript[I, 2]-1)-Subscript[I, 1]^2=0*)


realdatadispersioncheck[file15,hdf15mww,r15]


(* ::Subsubsection::Closed:: *)
(*Check via Dot Products: Subscript[\[CapitalOmega]A, m]=(-\[Omega] + Ve +Subscript[\[Phi], 0]+Subscript[\[Mu], m](k-Subscript[\[Phi], 1]))Subscript[A, m]-Subscript[\[CapitalSigma], j](1-Subscript[\[Mu], m] Subscript[\[Mu], j])(Subscript[n, j] Subscript[A, j]-Subscript[nb, j] Subscript[Ab, j])*)


(* ::Text:: *)
(*The dot product checks for all direction components ' m' for a single eigenvalue eigenvector pair . Passes if the dot product is 1.*)


realdatadispcheck2[file15,hdf15mww,r15,1,"check"-> 6]


(* ::Subsubsection::Closed:: *)
(*Subtractive cancellations check*)


(* ::Text:: *)
(*There is no issue of subtractive cancellation if the sum of equation 37.5 (final column) is less than a' s*)


Do[realdatadispcheck2[file15,hdf15mww,r15,mi,"check"-> 3],{mi,1,20}]


(* ::Section::Closed:: *)
(*15 Msun_400ms (DO)*)


file15do="G:\\.shortcut-targets-by-id\\1tib_pcM4sTpxj42Hb2g-kKJKo1pG_Qeq\\SamFlynn_fast_flavor\\40bin\\DO_data\\15Msun_400ms_DO_40.h5";
hdf15donow="C:\\Users\\Sam\\Documents\\GitHub\\stability_tools\\15Msun_400ms_r255_DO_now_test.h5";
hdf15doww="C:\\Users\\Sam\\Documents\\GitHub\\stability_tools\\15Msun_400ms_r255_DO_w/w_test.h5";


(* ::Subsection::Closed:: *)
(*\[Omega] = 0*)


(* ::Subsubsection::Closed:: *)
(*'Dispersion Relation' Check: (Subscript[I, 0]+1)(Subscript[I, 2]-1)-Subscript[I, 1]^2=0*)


realdatadispersioncheck[file15do,hdf15donow,r15]


(* ::Subsubsection::Closed:: *)
(*Check via Dot Products: Subscript[\[CapitalOmega]A, m]=(-\[Omega] + Ve +Subscript[\[Phi], 0]+Subscript[\[Mu], m](k-Subscript[\[Phi], 1]))Subscript[A, m]-Subscript[\[CapitalSigma], j](1-Subscript[\[Mu], m] Subscript[\[Mu], j])(Subscript[n, j] Subscript[A, j]-Subscript[nb, j] Subscript[Ab, j])*)


(* ::Text:: *)
(*The dot product checks for all direction components ' m' for a single eigenvalue eigenvector pair . Passes if the dot product is 1.*)


realdatadispcheck2[file15do,hdf15donow,r15,1,"check"-> 6]


(* ::Subsubsection::Closed:: *)
(*Subtractive cancellations check*)


(* ::Text:: *)
(*There is no issue of subtractive cancellation if the sum of equation 37.5 (final column) is less than a' s*)


Do[realdatadispcheck2[file15do,hdf15donow,r15,mi,"check"-> 3],{mi,1,20}]


(* ::Subsection::Closed:: *)
(*\[Omega] != 0*)


(* ::Subsubsection::Closed:: *)
(*'Dispersion Relation' Check: (Subscript[I, 0]+1)(Subscript[I, 2]-1)-Subscript[I, 1]^2=0*)


realdatadispersioncheck[file15do,hdf15doww,r15]


(* ::Subsubsection::Closed:: *)
(*Check via Dot Products: Subscript[\[CapitalOmega]A, m]=(-\[Omega] + Ve +Subscript[\[Phi], 0]+Subscript[\[Mu], m](k-Subscript[\[Phi], 1]))Subscript[A, m]-Subscript[\[CapitalSigma], j](1-Subscript[\[Mu], m] Subscript[\[Mu], j])(Subscript[n, j] Subscript[A, j]-Subscript[nb, j] Subscript[Ab, j])*)


(* ::Text:: *)
(*The dot product checks for all direction components ' m' for a single eigenvalue eigenvector pair . Passes if the dot product is 1.*)


realdatadispcheck2[file15do,hdf15doww,r15,1,"check"-> 6]


(* ::Subsubsection::Closed:: *)
(*Subtractive cancellations check*)


(* ::Text:: *)
(*There is no issue of subtractive cancellation if the sum of equation 37.5 (final column) is less than a' s*)


Do[realdatadispcheck2[file15do,hdf15doww,r15,mi,"check"-> 3],{mi,1,20}]


(* ::Section::Closed:: *)
(*Energy Dependence*)


(* ::Subsection::Closed:: *)
(*Generation of files*)


energyd[n_]:=N[20 10^-n];
Do[
exportkadapt[kAdapt[inpath <> "112Msun_100ms_DO.h5", 250,250, energyd[n], -1., 2, "xflavor" -> True],StringJoin["112Msun_100ms_r200_r300_E20n",ToString[n],"_nox_test"]]
,{n,0,6}]



(* ::Subsection::Closed:: *)
(*'Dispersion Relation' Check: (Subscript[I, 0]+1)(Subscript[I, 2]-1)-Subscript[I, 1]^2=0*)


Do[realdatadispersioncheck[file11,StringJoin["C:\\Users\\Sam\\Documents\\GitHub\\stability_tools\\112Msun_100ms_r200_r300_E20n",ToString[n],"_nox_test.h5"],r11]
,{n,0,6}]


(* ::Subsection::Closed:: *)
(*Check via Dot Products: Subscript[\[CapitalOmega]A, m]=(-\[Omega] + Ve +Subscript[\[Phi], 0]+Subscript[\[Mu], m](k-Subscript[\[Phi], 1]))Subscript[A, m]-Subscript[\[CapitalSigma], j](1-Subscript[\[Mu], m] Subscript[\[Mu], j])(Subscript[n, j] Subscript[A, j]-Subscript[nb, j] Subscript[Ab, j])*)


Do[realdatadispcheck2[file11,StringJoin["C:\\Users\\Sam\\Documents\\GitHub\\stability_tools\\112Msun_100ms_r200_r300_E20n",ToString[n],"_nox_test.h5"],r11,1,"check"-> 6]
,{n,0,6}]





(* ::Section::Closed:: *)
(*Stability Matrix Checks*)


(*This function builds the stability matrix by computing each component exactly, and then building the total table.*)
Options[buildS]={"xflavor"-> True};
buildS[file_,r_,k_,En_,hier_,OptionsPattern[]]:=Module[{datasr,n,nb,ndens,Ve,\[Omega],\[Phi]0,\[Phi]1,\[Phi],cos\[Theta],nwhich,\[Omega]which,S,Smat,\[Mu]which},
datasr=SelectSingleRadius[ImportData[file],r];
cos\[Theta]=datasr["mids"];
\[Omega]=hier \[Omega]EMev[En];
Ve=0. ;(*munits/mp *datasr["Yes"]*datasr["matters"];*)
ndens=ndensities[datasr,"xflavor"-> OptionValue["xflavor"]];
n=munits Diagonal[ndens[[1]]-ndens[[3]]];
nb=munits Diagonal[ndens[[2]]-ndens[[3]]];
\[Phi][rank_,mu_,mubar_,cos\[Theta]_]:=Sum[(mu[[i]]-mubar[[i]])cos\[Theta][[i]]^rank,{i,1,Length[cos\[Theta]]}];
\[Omega]which[i_]:=Which[i<= Length[cos\[Theta]], -1.,i>Length[cos\[Theta]],1.];
nwhich[j_]:=Which[j<= Length[cos\[Theta]], n[[j]],j>Length[cos\[Theta]], -nb[[(j-Length[cos\[Theta]])]] ];
\[Mu]which[i_]:= Which[i<= Length[cos\[Theta]],0.,i>Length[cos\[Theta]],Length[cos\[Theta]]];
S[i_,j_]:= (\[Omega]which[i] \[Omega] KroneckerDelta[i,j])+(Ve+\[Phi][0,n,nb,cos\[Theta]]+cos\[Theta][[i-\[Mu]which[i]]] (k-\[Phi][1,n,nb,cos\[Theta]])) KroneckerDelta[i,j]
-(1-(cos\[Theta][[i-\[Mu]which[i]]]cos\[Theta][[j-\[Mu]which[j]]]))nwhich[j];
Smat=Table[S[i,j],{i,1,2*Length[cos\[Theta]]},{j,1,2*Length[cos\[Theta]]}];
Return[Smat];

];


(* ::Subsection::Closed:: *)
(*Confirmation S from paper is S from package*)


ktest=6.0555782595072755`*^-21;


Clear[St,Stest]


(*S built entry-by entry*)
Stest=buildS[file11,250,ktest,20.,-1.];


(*S built by package*)
tdat=SelectSingleRadius[ImportData[file11],250];
ndens=ndensities[tdat];
ea=getEquations[tdat,20.,-1.,ktest,ndens];
St=stabilityMatrix[tdat,ea];


(*Test that Stest and St are not different by examining the matrix of relative differences*)
(Abs[Stest-St]/Abs[Stest+St])//Max
((Abs[Stest]-Abs[St])/Abs[Stest+St])//Max


(*Explicitly check that the diagonals of Stest and St are the same*)
Row[{Stest//Diagonal//MatrixForm,
St//Diagonal//MatrixForm}]


(*Confirm that Stest and St agree on eigenvalues, calculated 3 different ways*)

Row[
{Sort[Eigenvalues[Stest],Greater]//MatrixForm,
Sort[evscale[ktest,Stest,kt,"output"-> "Eigenvalues"],Greater]//MatrixForm,
Sort[Eigenvalues[Stest*Min[Abs[Stest]]^-1],Greater]*Min[Abs[Stest]]//MatrixForm,
Sort[Eigenvalues[St],Greater]//MatrixForm,
Sort[evscale[ktest,St,kt,"output"-> "Eigenvalues"],Greater]//MatrixForm,
Sort[Eigenvalues[St*Min[Abs[St]]^-1],Greater]*Min[Abs[St]]//MatrixForm,
}]
(*Output*)
(*| Eigenvalues[Stest] | evscale[Stest] | Eigenvalues[(Stest scaled by Stest//Min)] |Eigenvalues[St] | evscale[St] | Eigenvalues[(St scaled by St//Min)] |  *)


(*Calculate eigenvalues and eigenvectors of Stest*)
{test1val,test1vec}=Eigensystem[Stest];
{test2val,test2vec}=evscale[ktest,Stest,kt,"output"-> "Eigensystem"];
{test3val,test3vec}=Eigensystem[Stest*Min[Abs[Stest]]^-1]*Min[Abs[Stest]];


(*Check that the eigenvalues of Stest are the same when using scaled Stest's*)
Row[{test1val//MatrixForm,
test2val//MatrixForm,
test3val//MatrixForm}]
(*Test that the eigenvectors are the same*)
Row[{test1vec[[1]]//Normalize//MatrixForm,
test2vec[[1]]//Normalize//MatrixForm,
test3vec[[1]]//Normalize//MatrixForm}]





(* ::Subsection::Closed:: *)
(*Other consistency checks*)


(*Prints out k,\[Phi]0,\[Phi]1,As,Abs,\[Omega],\[CapitalOmega]s,mu,mubar used in debugging checks*)
realdatadispcheck2[file15,hdf15mww,r15,1,"check"-> 7]


(* ::Item:: *)
(*tested that the outputs of check-> 7 match those used in building S entry by entry*)


(* ::Item:: *)
(*Checked that the output of kadapt gets the same eigenvectors and values when called on S*)


(* ::Item:: *)
(*Checked that the hdf5 output file is correctly written by taking the values directly from kadapt and comparing to those taken from the hdf5 file itself.*)
