(* ::Package:: *)

(* ::Subsection:: *)
(*Import data from a SINGLE file and a SINGLE radial bin - will be used inside of lists later, so everything is specific to one radius*)


filename = "112Msun_100ms_DO";
infile = "/mnt/data/SamFlynn/lotsadata/"<>filename<>".h5"
out_path = "/mnt/data/SamFlynn/stability_data/"<>filename
ri = 1


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


ImportedData = ImportData[infile,ri];
lotsodo=ImportedData[["lotsodo"]];
matters=ImportedData[["matters"]];
Yes=ImportedData[["Yes"]];
freqs=ImportedData[["freqs"]];
freqmid=ImportedData[["freqmid"]];
muss=ImportedData[["muss"]];
mids=ImportedData[["mids"]];
