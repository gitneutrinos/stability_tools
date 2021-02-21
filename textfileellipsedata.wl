(* ::Package:: *)

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



(*Ellipse fits the entire file, calculates the parameters and errors, and then exports a text file*)
ellipsefitout[infile_,name_]:=Module[{moments,efitneutrinos,efitanti,efitx,tab,rs,rend,tab2},
rend=384; (*for debugging*)
(*1:Calculate the 3 moments from the distribution function of a DO file*)
moments=getDOmoments[infile,1,rend]; (*Returns a list with Subscript[n, r] entries, each entry is a list of the 3 moments*)
(*2: Ellipse fit each species for all r*)
{efitneutrinos,efitanti,efitx}={ellipseFitSingleSpeciesDO[moments,1,1,rend],ellipseFitSingleSpeciesDO[moments,2,1,rend],ellipseFitSingleSpeciesDO[moments,3,1,rend]}; 
	(*each efit contains a Subscript[n, r] entries with two parts in each entry. The first part is the list of 3 parameters, the second part is the list of 3 moment errors*)
(*3: Format into a table*)
rs=ImportData[infile]["radius"]; (*get Radii values*)
tab=Table[{rs[[r]],
efitneutrinos[[r,1,1]],efitneutrinos[[r,1,2]],efitneutrinos[[r,1,3]],efitneutrinos[[r,2,1]],efitneutrinos[[r,2,2]],efitneutrinos[[r,2,3]],
efitanti[[r,1,1]],efitanti[[r,1,2]],efitanti[[r,1,3]],efitanti[[r,2,1]],efitanti[[r,2,2]],efitanti[[r,2,3]],
efitx[[r,1,1]],efitx[[r,1,2]],efitx[[r,1,3]],efitx[[r,2,1]],efitx[[r,2,2]],efitx[[r,2,3]]},{r,1,rend}];
(*4: Add header row*)
tab2=Prepend[tab,{"radius (cm)","a(nu)","b(nu)","c(nu)","m0 err(nu)","m1 err (nu)","m2 err (nu)","a(nub)","b(nub)","c(nub)","m0 err(nub)","m1 err (nub)","m2 err (nub)","a(nux)","b(nux)","c(nux)","m0 err(nux)","m1 err (nux)","m2 err (nux)"}];
(*5: Export as text file*)
Export[name<>".txt",tab2,"Table"]
]


ellipsefitout["G:\\My Drive\\Physics\\Neutrino Oscillation Research\\Fast Conversions\\lotsadata.tar\\lotsadata\\lotsadata\\15Msun_1ms_DO.h5","test"]
