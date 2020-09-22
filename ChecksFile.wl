(* ::Package:: *)

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



(* ::Subsection:: *)
(*Regression test for growth rate over a range of wavenumbers. The new results (blue) should match the old results (orange).*)


(*Test on kAdapt. Comparing to old data (below). Generates a plot which should look similar (different exact k values tested*)
OldData={{3.943351801881396`*^-20,0},{6.403122930332929`*^-20,0},{1.0397242072440586`*^-19,6.475185304707499`*^-23},{1.6882799828943322`*^-19,0},
{2.7413897654617433`*^-19,1.2805250574405962`*^-19},{4.451404934206825`*^-19,2.149878747929582`*^-19},{7.228087788874928`*^-19,2.877126337699488`*^-19},
{1.173680082038912`*^-18,3.418162146974343`*^-19},{1.9057944164639052`*^-18,3.5031848436345418`*^-19},{3.0945846431298252`*^-18,2.6089021800081186`*^-19},
{5.02491456096483`*^-18,7.0760561598301275`*^-22},{8.159339380505439`*^-18,7.514263464205905`*^-22},{1.3248945493211215`*^-17,7.789552045160011`*^-22},
{2.1513329510655564`*^-17,7.961350592910255`*^-22},{3.493284404190471`*^-17,8.068054415566711`*^-22},{5.672313958895247`*^-17,8.134119405806582`*^-22},
{9.210571463829571`*^-17,8.174940729310282`*^-22},{1.4955911697601837`*^-16,8.200132465051396`*^-22},{2.4285061527926315`*^-16,8.215666277345043`*^-22},
{3.9433518018813854`*^-16,8.225240355557669`*^-22},{-3.943351801881396`*^-20,1.2187126603736879`*^-22},{-1.0972610632535474`*^-19,3.174638959109932`*^-20},
{-3.053194088231941`*^-19,5.832887809619426`*^-20},{-8.495693916972988`*^-19,1.851851392271892`*^-21},{-2.3639772986947037`*^-18,1.1000927210473065`*^-21},
{-6.577907259087188`*^-18,9.189580059662592`*^-22},{-1.830341769061953`*^-17,8.575191336733617`*^-22},{-5.0930346987557535`*^-17,8.360035843249597`*^-22},
{-1.4171671586789936`*^-16,8.283420407361107`*^-22},{-3.9433518018813854`*^-16,8.255976606410335`*^-22}};
kdebug=Module[{data,ri=200,testE=20,hi=-1,kvar},
file=inpath<>"1D_withV_withPairBrems_DO.h5";
(*SCalcScale[ImportData[inpath<>file<>".h5"],ri,testE,hi,0][[3]]//MatrixForm*)
(*buildkGrid[ImportData[inpath<>file<>".h5"],ri,testE,hi,40]*)
kAdapt[file,ri,ri,testE,hi,20,"xflavor"-> False]
];
plot1=ListLogPlot[{Transpose@{kdebug[[All,3]],kdebug[[All,4]]}},ImageSize-> Scaled[0.25]];
plot2=ListLogPlot[{OldData},ImageSize-> Scaled[0.25]];
VerificationTest[Rasterize[plot1]==Rasterize[plot2],TestID-> "kAdapt test; current calculation = old calculation"]





(* ::Subsection:: *)
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
(*builds a 4x4 matrix with the two beam data, then reduces the size of the matrix *)
build2bMatrix[En_,k_]:=Module[{fakeEn,S2ba,S2b},
S2ba=stabilityMatrix[get2bdata[],getEquations[get2bdata[],En,-1.,k,"xflavor"-> False],"xflavor"-> False];
S2b={{S2ba[[2,2]],S2ba[[2,3]]},{S2ba[[3,2]],S2ba[[3,3]]}};
Return[S2b]
];


(*Chakraborty's analytic expression for the 2x2 case, and a verification test that they're the same within one part in 1000*)
\[CapitalOmega]ch[k_,\[Mu]ch_,a_,\[Omega]_]:=Sort[{2 a \[Mu]ch+Sqrt[(2 a \[Mu]ch)^2+(\[Omega]+k)((\[Omega]+k)-4 \[Mu]ch)],2 a \[Mu]ch-Sqrt[(2 a \[Mu]ch)^2+(\[Omega]+k)((\[Omega]+k)-4 \[Mu]ch)]}];
cm[k_,\[Mu]ch_,w_]:=DiagonalMatrix[{w+k,-w-k,w-k,-w+k}]+2 \[Mu]ch{{l+lb,-lb,-l,0.},{-r,r+rb,0,-rb},{-r,0,r+rb,-rb},{0.,-lb,-l,l+lb}};
cma[k_,\[Mu]ch_,a_,w_]:=cm[k,\[Mu]ch,w]/.{rb-> 0.,l-> 0.,r-> (1+a),lb-> -(1-a)};



evtest=Sort[Eigenvalues[build2bMatrix[Infinity,2.]/.{a-> 0.}]]//Chop;
VerificationTest[
Re[\[CapitalOmega]ch[2.,1.,0.,0]]===Re[evtest] 
,TestID-> "2 Beam Growth Rate (Real part)"]
VerificationTest[
Im[\[CapitalOmega]ch[2.,1.,0.,0]]=== Im[evtest]
,TestID-> "2 Beam Growth Rate (Imaginary part)"]
Print[Re[\[CapitalOmega]ch[2.,1.,0.,0]],Re[evtest]];
Print[Im[\[CapitalOmega]ch[2.,1.,0.,0]],Im[evtest]];


(* ::Subsection:: *)
(*Preliminaries for dispersion checks*)


(*Calculates and Returns the nth I for the dispersion check.  Returns a single value of In*)
Idis[data_,\[CapitalOmega]_,k_,n_,xflavor_]:=Module[{cos\[Theta],\[Phi]0,\[Phi]1,\[CapitalOmega]p,kp,\[Omega],mu,mubar,Vmatter,\[CapitalOmega]minuskpcos\[Theta],result},

(* direction cosines *)
cos\[Theta]=data["mids"];

(* neutrino number densities disguised as SI potentials *)
mu[i_]:=   munits ndensities[data,"xflavor"->xflavor][[1,i,i]];
mubar[i_]:=munits ndensities[data,"xflavor"->xflavor][[2,i,i]];

(*Defined in Gail's Blue equation 30 and 31 *)
\[Phi]0 = Sum[(mu[i]-mubar[i])         ,{i,1,Length[cos\[Theta]]}];
\[Phi]1 = Sum[(mu[i]-mubar[i])cos\[Theta][[i]],{i,1,Length[cos\[Theta]]}];

(* "Shifted" Eigenvalue and k*)
Vmatter = munits data["Yes"] data["matters"]/mp;
\[CapitalOmega]p = N[\[CapitalOmega]-Vmatter-\[Phi]0];
kp = k-\[Phi]1;
\[Omega]=\[Omega]EMev[En];
\[CapitalOmega]minuskpcos\[Theta][i_]:=\[CapitalOmega]p-(kp cos\[Theta][[i]]);

(* make sure the denominator is not tiny *)
Do[Assert[(\[CapitalOmega]p-(kp cos\[Theta][[i]]))/(\[CapitalOmega]p+(kp cos\[Theta][[i]]))>= 1,"\[CapitalOmega]p' k'Cos[\[Theta]] percent difference less than 1"],{i,1,Length[\[Theta]]}]; 

(*Definition of I from Gail's equation (41)*)
result= Sum[cos\[Theta][[i]]^n (((mu[i]-mubar[i])\[CapitalOmega]minuskpcos\[Theta][i]) + (mu[i]+mubar[i])\[Omega]) / (\[CapitalOmega]minuskpcos\[Theta][i]^2-\[Omega]^2) ,{i,1,Length[cos\[Theta]]}];
Return[result]
];


dispersionCheck[data_,\[CapitalOmega]_,k_,En_,xflavor_]:=Module[{I0,I1,I2},
(*The condition is that Equatrion (43), below, should be 0 if the vacuum is*)
I0 = Idis[data,\[CapitalOmega],k,0,xflavor];
I1 = Idis[data,\[CapitalOmega],k,1,xflavor];
I2 = Idis[data,\[CapitalOmega],k,2,xflavor];
Return[(I0+1.)(I2-1.)-I1^2]
];


test2bdispersionCheck[k_,En_,atest_,xflavor_]:=Module[{I0,I1,I2,\[CapitalOmega],data},
data = get2bdata[]/.a-> atest;
\[CapitalOmega]=evscale[k,build2bMatrix[En,k]/.a-> 0.,kvar,"output"->"Eigenvalues"][[1]];
Return[dispersionCheck[data,\[CapitalOmega],k,En,xflavor]]
];

(*This exists so that I can run the testfile thing without having to worry about context paths for the package.*)
(*Will require similar modules in order to run in the testfile*)



test4bdispersionCheck[k_,En_,atest_,xflavor_]:=Module[{I0,I1,I2,\[CapitalOmega],data,equations},
data = get2bdata[]/.a-> atest;
equations = getEquations[data,En,-1.,k,"xflavor"->xflavor];
\[CapitalOmega]=evscale[k,stabilityMatrix[data,equations,"xflavor"->xflavor],kx,"output"-> "Eigenvalues"][[1]];
Return[dispersionCheck[data,\[CapitalOmega],k,En,xflavor]]
];




(* ::Subsection::Closed:: *)
(*Two-beam \[Omega]=0 dispersion check*)


(* Inputs *)
data = get2bdata[]/.a-> 0;
k=1;
En=Infinity; (* MeV *)

(* Calculations *)
\[CapitalOmega]=evscale[k,build2bMatrix[En,k]/.a-> 0.,kvar,"output"->"Eigenvalues"][[1]]
xflavor=False;
zero=dispersionCheck[data,\[CapitalOmega],k,En,xflavor]

(* Test *)
VerificationTest[Between[Abs[zero],{-0.01,0.01}],TestID-> "2 beam Dispersion Check"]


(* ::Subsection::Closed:: *)
(*2 Beam \[Omega]!=0 Dispersion Check*)


(* Inputs *)
data = get2bdata[]/.a-> 0.;
k=1.;
En=20.; (* MeV *)

(* Calculations *)
\[CapitalOmega] = evscale[k,build2bMatrix[En,k]/.a-> 0.,kvar,"output"->"Eigenvalues"][[1]]
zero = dispersionCheck[data,\[CapitalOmega],k,En,xflavor]

(* Test *)
VerificationTest[Between[Abs[zero],{-0.01,0.01}],TestID-> "2 Beam \[Omega]\[NotEqual]0 Dispersion Check"]


(* ::Subsection::Closed:: *)
(*Four-beam \[Omega]=0 dispersion check*)


(* Inputs *)
data = get2bdata[]/.a-> 0.;
k=1.;
En=Infinity;
hierarchy=-1;
xflavor=False;

(* Calculations *)
equations = getEquations[data,En,hierarchy,k,"xflavor"->xflavor];
\[CapitalOmega]=evscale[k,stabilityMatrix[data,equations,"xflavor"->xflavor],kx,"output"-> "Eigenvalues"][[1]]
zero = dispersionCheck[data,\[CapitalOmega],k,En,xflavor]

(* Test *)
VerificationTest[Between[Abs[zero],{-0.01,0.01}],TestID-> "4 Beam Dispersion Check"]


(* ::Subsection:: *)
(*Four - beam \[Omega] != 0 dispersion check*)


(* Inpus *)
data = get2bdata[]/.a-> 0.;
k=1.;
En=20;
hierarchy=-1;
xflavor=False;

(* Calculations *)
equations = getEquations[data,En,hierarchy,k,"xflavor"->xflavor];
\[CapitalOmega]=evscale[k,stabilityMatrix[data,equations,"xflavor"->xflavor],kx,"output"-> "Eigenvalues"][[1]]
zero = dispersionCheck[data,\[CapitalOmega],k,En,xflavor]

(* Test *)
VerificationTest[Between[Abs[zero],{-0.01,0.01}],TestID-> "4 Beam Dispersion Check"]


 tr=TestReport["C:\\Users\\Sam\\Documents\\GitHub\\stability_tools\\testfiles.wlt"]
 tr["TestResults"]


(* ::Subsection:: *)
(*Real data dispersion check*)


(* Inputs *)
data=ImportData[inpath<>"112Msun_100ms_DO.h5"];
datasr=SelectSingleRadius[data,250];
datasr["matters"]=0.; (*matter set to 0*)
k=0.;
En=Infinity; (* MeV *)
hierarchy=-1;
xflavor=False;

(* Calculations *)
equations = getEquations[datasr,En,hierarchy,k,"xflavor"-> False];
\[CapitalOmega]=evscale[k,stabilityMatrix[datasr,equations,"xflavor"-> False],kx,"output"-> "Eigenvalues"][[2]]
zero=dispersionCheck[datasr,\[CapitalOmega],k,En,xflavor]

(* Test *)
VerificationTest[Between[zero,{-0.01,0.01}],TestID-> "Real Data \[Omega]\[NotEqual]0 Dispersion Check"]


(* ::Subsection:: *)
(*Check that ellipse construction results in the correct moments given hand-chosen moments*)


ellipsefiterrors[m0_,m1_,m2_]:=Module[{er0,er1,er2,fits},
fits=eBoxFitToMoments[m0,m1,m2,getInitialGuess[m0,m1,m2]];
er0=(ellipseMoments[fits[[1]],fits[[2]],fits[[3]]][[1]]-m0)/m0;
er1=(ellipseMoments[fits[[1]],fits[[2]],fits[[3]]][[2]]-m1)/m1;
er2=(ellipseMoments[fits[[1]],fits[[2]],fits[[3]]][[3]]-m2)/m2;
(*Print["Initial Guess: ", getInitialGuess[m0,m1,m2]];*)
Return[{er0,er1,er2}];
]
(*Ellipse fits to 3 moments and returns the percent errors*)

ellipsefitrealdatacheck[]:=Module[{file,moms},
file=inpath<>"4timesHigh_1D_withV_withPairBrems_MC_moments.h5";
moms=Quiet[Quiet[getMoments[file,1,1],{Import::general}],{Import::noelem}]; (*quiets only the import complaint that there are no midpoints for moments*)
Return[ellipsefiterrors[moms[[1]],moms[[2]]//Abs,moms[[3]]]]
];
(*Imports real CSSN data and then calls ellipse fit errors for the tests file*)


m0=1.;
m1=10.^-8;
m2=1./3.;
fits=eBoxFitToMoments[m0,m1,m2,getInitialGuess[m0,m1,m2]];
er0=(ellipseMoments[fits[[1]],fits[[2]],fits[[3]]][[1]]-m0)/m0;
er1=(ellipseMoments[fits[[1]],fits[[2]],fits[[3]]][[2]]-m1)/m1;
er2=(ellipseMoments[fits[[1]],fits[[2]],fits[[3]]][[3]]-m2)/m2;
Print["Initial Guess: ", getInitialGuess[m0,m1,m2]];

VerificationTest[Abs[er0]<10^-5 && Abs[er1]< 10^-5 && Abs[er2]< 10^-5,TestID-> "Ellipse error check"]
{er0,er1,er2}


(* ::Subsection:: *)
(*Check that ellipse construction results in the correct moments given moments from CCSN data*)


file=inpath<>"4timesHigh_1D_withV_withPairBrems_MC_moments.h5";
moms=Quiet[Quiet[getMoments[file,1,1],{Import::general}],{Import::noelem}]; (*quiets only the import complaint that there are no midpoints for moments*)
m0=moms[[1]];
m1=moms[[2]]//Abs;
m2=moms[[3]];
fits=eBoxFitToMoments[m0,m1,m2,getInitialGuess[m0,m1,m2]];
er0=(ellipseMoments[fits[[1]],fits[[2]],fits[[3]]][[1]]-m0)/m0;
er1=(ellipseMoments[fits[[1]],fits[[2]],fits[[3]]][[2]]-m1)/m1;
er2=(ellipseMoments[fits[[1]],fits[[2]],fits[[3]]][[3]]-m2)/m2;

{VerificationTest[Abs[er0]<10^-4 && Abs[er1]< 10^-3 && Abs[er2]< 10^-4,TestID-> "Ellipse Data error check"],er0,er1,er2}












