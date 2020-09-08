(* ::Package:: *)

(************************************************************************)
(* This file was generated automatically by the Mathematica front end.  *)
(* It contains Initialization cells from a Notebook file, which         *)
(* typically will have the same name as this file except ending in      *)
(* ".nb" instead of ".m".                                               *)
(*                                                                      *)
(* This file is intended to be loaded into the Mathematica kernel using *)
(* the package loading commands Get or Needs.  Doing so is equivalent   *)
(* to using the Evaluate Initialization Cells menu command in the front *)
(* end.                                                                 *)
(*                                                                      *)
(* DO NOT EDIT THIS FILE.  This entire file is regenerated              *)
(* automatically each time the parent Notebook file is saved in the     *)
(* Mathematica front end.  Any changes you make to this file will be    *)
(* overwritten.                                                         *)
(************************************************************************)



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



On[Assert]
$MinPrecision=30;


(*Test on kAdapt. Comparing to old data (below). Generates a plot which should look similar (different exact k values tested*)
OldData={{2.61019*10^-17, 1.12835*10^-18}, {2.8674*10^-17, 
  7.11518*10^-19}, {3.14996*10^-17, 1.13532*10^-21}, {3.46036*10^-17, 
  1.12664*10^-21}, {3.80135*10^-17, 1.11669*10^-21}, {4.17594*10^-17, 
  1.10578*10^-21}, {4.58744*10^-17, 1.09417*10^-21}, {5.03949*10^-17, 
  1.08207*10^-21}, {5.53608*10^-17, 1.06968*10^-21}, {6.08162*10^-17, 
  1.05717*10^-21}, {6.68091*10^-17, 1.04469*10^-21}, {7.33925*10^-17, 
  1.03234*10^-21}, {8.06247*10^-17, 1.02023*10^-21}, {8.85695*10^-17, 
  1.00844*10^-21}, {9.72973*10^-17, 9.97028*10^-22}, {1.06885*10^-16, 
  9.86047*10^-22}, {1.17418*10^-16, 9.75532*10^-22}, {1.28988*10^-16, 
  9.78299*10^-22}, {1.41699*10^-16, 9.85668*10^-22}, {1.55662*10^-16, 
  9.92358*10^-22}, {1.71001*10^-16, 9.98434*10^-22}, {1.87852*10^-16, 
  1.00395*10^-21}, {2.06363*10^-16, 1.00897*10^-21}, {2.26698*10^-16, 
  1.01353*10^-21}, {2.49037*10^-16, 1.01767*10^-21}, {2.73577*10^-16, 
  1.02144*10^-21}, {3.00536*10^-16, 1.02487*10^-21}, {3.30151*10^-16, 
  1.02798*10^-21}, {3.62685*10^-16, 1.03082*10^-21}, {3.98424*10^-16, 
  1.03339*10^-21}, {4.37685*10^-16, 1.03574*10^-21}, {4.80815*10^-16, 
  1.03787*10^-21}, {5.28195*10^-16, 1.03981*10^-21}, {5.80244*10^-16, 
  1.04157*10^-21}, {6.37422*10^-16, 1.04318*10^-21}, {7.00235*10^-16, 
  1.04464*10^-21}, {7.69236*10^-16, 1.04597*10^-21}, {8.45038*10^-16, 
  1.04718*10^-21}, {9.28309*10^-16, 1.04828*10^-21}, {1.01979*10^-15, 
  1.04928*10^-21}, {1.12028*10^-15, 1.05019*10^-21}, {1.23067*10^-15, 
  1.05102*10^-21}, {1.35194*10^-15, 1.05178*10^-21}, {1.48516*10^-15, 
  1.05247*10^-21}, {1.63151*10^-15, 1.05309*10^-21}, {1.79228*10^-15, 
  1.05366*10^-21}, {1.9689*10^-15, 1.05418*10^-21}, {2.16292*10^-15, 
  1.05465*10^-21}, {2.37605*10^-15, 1.05508*10^-21}, {2.61019*10^-15, 
  1.05547*10^-21}, {-2.61019*10^-17, 1.44396*10^-21}, {-3.3261*10^-17,
   1.36214*10^-21}, {-4.23837*10^-17, 
  1.29781*10^-21}, {-5.40084*10^-17, 
  1.24715*10^-21}, {-6.88216*10^-17, 
  1.20722*10^-21}, {-8.76976*10^-17, 
  1.17574*10^-21}, {-1.11751*10^-16, 
  1.15093*10^-21}, {-1.42401*10^-16, 
  1.13139*10^-21}, {-1.81459*10^-16, 1.116*10^-21}, {-2.31228*10^-16, 
  1.10389*10^-21}, {-2.94648*10^-16, 
  1.09436*10^-21}, {-3.75463*10^-16, 
  1.08687*10^-21}, {-4.78443*10^-16, 
  1.08098*10^-21}, {-6.09668*10^-16, 
  1.07635*10^-21}, {-7.76884*10^-16, 
  1.07272*10^-21}, {-9.89964*10^-16, 
  1.06986*10^-21}, {-1.26149*10^-15, 
  1.06762*10^-21}, {-1.60748*10^-15, 
  1.06586*10^-21}, {-2.04837*10^-15, 
  1.06448*10^-21}, {-2.61019*10^-15, 1.06339*10^-21}};
kdebug=Module[{data,ri=200,testE=20,hi=-1,k=0,kvar},
file=inpath<>"1D_withV_withPairBrems_DO.h5";
(*SCalcScale[ImportData[inpath<>file<>".h5"],ri,testE,hi,0][[3]]//MatrixForm*)
(*buildkGrid[ImportData[inpath<>file<>".h5"],ri,testE,hi,40]*)
kAdapt[file,ri,ri,testE,hi,10,"xflavor"-> False]
];
ListLogPlot[{Transpose@{kdebug[[All,2]],kdebug[[All,3]]},OldData},ImageSize-> Scaled[0.25]]


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


(*Check the dispersion relation from Gail's paper*)
dispersionCheck[data_,\[CapitalOmega]_,k_]:=Module[{\[Phi]0,\[Phi]1,\[CapitalOmega]p,kp,check,Idis,\[Theta]},
\[Theta]=data["mids"];
(*Defined in Gail's Blue equation 30 and 31 *)
\[Phi]0=Sum[munits ndensities[data,"xflavor"-> False][[1,i,i]]-munits ndensities[data,"xflavor"-> False][[2,i,i]],{i,1,Length[\[Theta]]}];
\[Phi]1=Sum[(munits ndensities[data,"xflavor"-> False][[1,i,i]]-munits ndensities[data,"xflavor"-> False][[2,i,i]])\[Theta][[i]],{i,1,Length[\[Theta]]}];
(* "Shifted" Eigenvalue and k*)
\[CapitalOmega]p=N[\[CapitalOmega]-((munits/mp) *data["Yes"] *data["matters"])-\[Phi]0];
kp=k-\[Phi]1;
(*Definition of I from Gail's equation (41)*)
Idis[n_]:= Sum[(( munits(ndensities[data,"xflavor"-> False][[1]][[i,i]]- ndensities[data,"xflavor"-> False][[2]][[i,i]]))/(\[CapitalOmega]p-(kp \[Theta][[i]]))) \[Theta][[i]]^n,{i,1,Length[\[Theta]]}]//Chop;
(*The condition is that Equatrion (43), below, should be 0 if the vacuum is*)
check=((Idis[0]+1)(Idis[2]-1))-(Idis[1]^2)//Chop;
Return[check]
];


wdispersionCheck[data_,\[CapitalOmega]_,k_,En_]:=Module[{\[Phi]0,\[Phi]1,\[CapitalOmega]p,kp,check,Idis,\[Theta],\[Omega]},
\[Theta]=data["mids"];
(*Defined in Gail's Blue equation 30 and 31 *)
\[Phi]0=Sum[munits ndensities[data,"xflavor"-> False][[1,i,i]]-munits ndensities[data,"xflavor"-> False][[2,i,i]],{i,1,Length[\[Theta]]}];
\[Phi]1=Sum[(munits ndensities[data,"xflavor"-> False][[1,i,i]]-munits ndensities[data,"xflavor"-> False][[2,i,i]])\[Theta][[i]],{i,1,Length[\[Theta]]}];
(* "Shifted" Eigenvalue and k*)
\[CapitalOmega]p=N[\[CapitalOmega]-((munits/mp) *data["Yes"] *data["matters"])-\[Phi]0];
kp=k-\[Phi]1;
\[Omega]=\[Omega]EMev[En];
(*Definition of I from Gail's equation (41)*)
Idis[n_]:= Sum[(( munits(ndensities[data,"xflavor"-> False][[1]][[i,i]]- ndensities[data,"xflavor"-> False][[2]][[i,i]]))(\[CapitalOmega]p-(kp \[Theta][[i]])+( munits(ndensities[data,"xflavor"-> False][[1]][[i,i]]+ ndensities[data,"xflavor"-> False][[2]][[i,i]]))\[Omega]))/((\[CapitalOmega]p-(kp \[Theta][[i]]))^2-\[Omega]^2) \[Theta][[i]]^n,{i,1,Length[\[Theta]]}]//Chop;
(*The condition is that Equatrion (43), below, should be 0 if the vacuum is*)
check=((Idis[0]+1)(Idis[2]-1))-(Idis[1])^2//Chop;
Return[check]
];


allDispersions[]:=Module[{data,dc2,dc4,dcdata,datasr,t1,t2,t3,wdc2b,wdcdata,t4,t5,kx},
dc2=dispersionCheck[get2bdata[]/.a-> 0.,Eigenvalues[build2bMatrix[Infinity,2.]/.a-> 0.][[1]],2.];
dc4=dispersionCheck[get2bdata[]/.a-> 0.,Eigenvalues[stabilityMatrix[get2bdata[],getEquations[get2bdata[],Infinity,-1.,2.,"xflavor"-> False],"xflavor"-> False]/.a-> 0.][[1]],2.];
data=ImportData["G:\\My Drive\\Physics\\Neutrino Oscillation Research\\Fast Conversions\\lotsadata.tar\\lotsadata\\lotsadata\\112Msun_100ms_DO.h5"];
datasr=SelectSingleRadius[data,250];
dcdata=dispersionCheck[datasr,
evscale[0.,stabilityMatrix[datasr,getEquations[datasr,Infinity,-1.,0.,"xflavor"-> False],"xflavor"-> False],kx][[1]]
,0.];
wdc2b=wdispersionCheck[get2bdata[]/.a-> 0.,Eigenvalues[build2bMatrix[20.,2.]/.a-> 0.][[1]],2.,20.];
wdcdata=wdispersionCheck[datasr,Eigenvalues[stabilityMatrix[datasr,getEquations[datasr,20.,-1.,0.,"xflavor"-> False],"xflavor"-> False]][[1]],0.,20.];
t1=VerificationTest[
Between[dc2//Chop,{-0.01,0.01}],
TestID-> "2 beam Dispersion Check"];
t2=VerificationTest[
Between[dc4//Chop,{-0.01,0.01}],
TestID-> "4 Beam Dispersion Check"];
t3=VerificationTest[
Between[dcdata//Chop,{-0.01,0.01}],
TestID-> "Real Data Dispersion Check"];
t4=VerificationTest[
Between[wdc2b//Chop,{-0.01,0.01}],
TestID-> "Real Data Dispersion Check"];
t5=VerificationTest[
Between[wdcdata//Chop,{-0.01,0.01}],
TestID-> "Real Data Dispersion Check"];
Return[{{t1},{t2},{t3},{t4},{t5}}//MatrixForm]
];


allDispersions[]


(*Ellipse Check Section*)
ellipseCheck[]:=Module[{m0,m1,m2,er0,er1,er2,fits},
m0=1.;
m1=10.^-8;
m2=1./3.;
fits=eBoxFitToMoments[m0,m1,m2,getInitialGuess[m0,m1,m2]];
er0=(ellipseMoments[fits[[1]],fits[[2]],fits[[3]]][[1]]-m0)/m0;
er1=(ellipseMoments[fits[[1]],fits[[2]],fits[[3]]][[2]]-m1)/m1;
er2=(ellipseMoments[fits[[1]],fits[[2]],fits[[3]]][[3]]-m2)/m2;
Print["Initial Guess: ", getInitialGuess[m0,m1,m2]];
Return[{VerificationTest[Abs[er0]<10^-5 && Abs[er1]< 10^-5 && Abs[er2]< 10^-5,TestID-> "Ellipse error check"],er0,er1,er2}]
];


dataEllipseCheck[]:=Module[{m0,m1,m2,er0,er1,er2,fits,moms,file},
file="G:\\My Drive\\Physics\\Neutrino Oscillation Research\\Fast Conversions\\lotsadata.tar\\lotsadata\\lotsadata\\4timesHigh_1D_withV_withPairBrems_MC_moments.h5";
moms=Quiet[Quiet[getMoments[file,1,1],{Import::general}],{Import::noelem}]; (*quiets only the import complaint that there are no midpoints for moments*)
m0=moms[[1]];
m1=moms[[2]]//Abs;
m2=moms[[3]];
fits=eBoxFitToMoments[m0,m1,m2,getInitialGuess[m0,m1,m2]];
er0=(ellipseMoments[fits[[1]],fits[[2]],fits[[3]]][[1]]-m0)/m0;
er1=(ellipseMoments[fits[[1]],fits[[2]],fits[[3]]][[2]]-m1)/m1;
er2=(ellipseMoments[fits[[1]],fits[[2]],fits[[3]]][[3]]-m2)/m2;
Return[{VerificationTest[Abs[er0]<10^-4 && Abs[er1]< 10^-3 && Abs[er2]< 10^-4,TestID-> "Ellipse Data error check"],er0,er1,er2}]
];


ellipseCheck[]


dataEllipseCheck[]



