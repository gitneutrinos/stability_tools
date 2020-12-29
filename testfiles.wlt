BeginTestSection["testfiles_notebook"]

VerificationTest[(* 1 *)
	Equal[Select[Abs[Times[Plus[Cases[Part[OldData, All, 2], Except[0]], Times[-1, Cases[Part[debugdata, All, 2], Except[0]]]], Power[Plus[Cases[Part[OldData, All, 2], Except[0]], Cases[Part[debugdata, All, 2], Except[0]]], -1]]], GreaterThan[Power[10, -6]]], List[]]
	,
	True	
	,
	TestID->"kAdapt test; current calculation = old calculation"
]

VerificationTest[(* 2 *)
	Re[Chop[Sort[Eigenvalues[ReplaceAll[build2bMatrix[Infinity, 2.`], List[Rule[a, 0.`]]]]]]]
	,
	Re[\[CapitalOmega]ch[2.`, 1.`, 0.`, 0]]	
	,
	TestID->"2 Beam Growth Rate (Real part)"
]

VerificationTest[(* 3 *)
	Im[Chop[Sort[Eigenvalues[ReplaceAll[build2bMatrix[Infinity, 2.`], List[Rule[a, 0.`]]]]]]]
	,
	Im[\[CapitalOmega]ch[2.`, 1.`, 0.`, 0]]	
	,
	TestID->"2 Beam Growth Rate (Imaginary part)"
]

VerificationTest[(* 4 *)
	Between[Chop[Abs[test2bdispersionCheck[1, Infinity, 0, False]]], List[Times[-1, Power[10, -10]], Power[10, -10]]]
	,
	True	
	,
	TestID-> "2 beam \[Omega]=0 dispersion check"
]

VerificationTest[(* 5 *)
	Between[Chop[Abs[test2bdispersionCheck[1, 20, 0, False]]], List[Times[-1, Power[10, -10]], Power[10, -10]]]
	,
	True	
	,
	TestID-> "2 beam \[Omega]!=0 dispersion check"
]

VerificationTest[(* 6 *)
	Between[Chop[Abs[test4bdispersionCheck[1, Infinity, 0, False]]], List[Times[-1, Power[10, -10]], Power[10, -10]]]
	,
	True	
	,
	TestID->"4 beam \[Omega]=0 dispersion check"
]

VerificationTest[(* 7 *)
	Between[Chop[Abs[test4bdispersionCheck[1, 20, 0, False]]], List[Times[-1, Power[10, -10]], Power[10, -10]]]
	,
	True	
	,
	TestID->"4 beam \[Omega]!=0 dispersion check"
]

VerificationTest[(* 8 *)
	realdatadispersioncheck[StringJoin[inpath, "112Msun_100ms_DO.h5"], "112Msun_100ms_r200_r300_now_nox_test.h5", 250]
	,
	True	
	,
	TestID-> "Real data Dispersion Check; \[Omega]=0, no x flavor"
]

VerificationTest[(* 9 *)
	realdatadispersioncheck[StringJoin[inpath, "112Msun_100ms_DO.h5"], "112Msun_100ms_r200_r300_nox_test.h5", 250]
	,
	True	
	,
	TestID-> "Real data Dispersion Check; \[Omega]!=0, no x flavor"
]

VerificationTest[(* 10 *)
	And[LessEqual[Part[eCheckArtificial[], 1], Power[10, -4]], LessEqual[Part[eCheckArtificial[], 2], Power[10, -4]], LessEqual[Part[eCheckArtificial[], 3], Power[10, -4]]]
	,
	True	
	,
	TestID->"Ellipse Fitting: Check ellipse fit error in moments for artificial data"
]

VerificationTest[(* 11 *)
	And[Less[Part[ellipsefitrealdatacheck[], 1], Power[10, -3]], Less[Part[ellipsefitrealdatacheck[], 2], Power[10, -3]], Less[Part[ellipsefitrealdatacheck[], 3], Power[10, -3]]]
	,
	True	
	,
	TestID-> "Ellipse Fitting: Check ellipse fit error in moments for CSSN data: ri=220"
]

VerificationTest[(* 12 *)
	CompoundExpression[Set[ncheckdata, SelectSingleRadius[ImportData[StringJoin[inpath, "112Msun_100ms_DO.h5"]], 250]], Greater[siPotential[ncheckdata, Rule["xflavor", True]], siPotential[ncheckdata, Rule["xflavor", False]]]]
	,
	True	
	,
	TestID-> "Check that the total number of neutrinos+antineutrinos increases with xflavor on"
]

VerificationTest[(* 13 *)
	Equal[Plus[Plus[Tr[Part[ndensities[ncheckdata, Rule["xflavor", True]], 1]], Tr[Part[ndensities[ncheckdata, Rule["xflavor", True]], 3]]], Times[-1, Plus[Tr[Part[ndensities[ncheckdata, Rule["xflavor", True]], 2]], Tr[Part[ndensities[ncheckdata, Rule["xflavor", True]], 3]]]]], Plus[Plus[Tr[Part[ndensities[ncheckdata, Rule["xflavor", False]], 1]], Tr[Part[ndensities[ncheckdata, Rule["xflavor", False]], 3]]], Times[-1, Plus[Tr[Part[ndensities[ncheckdata, Rule["xflavor", False]], 2]], Tr[Part[ndensities[ncheckdata, Rule["xflavor", False]], 3]]]]]]
	,
	True	
	,
	TestID-> "Check that the net lepton current is not changed with xflavor on"
]

EndTestSection[]
