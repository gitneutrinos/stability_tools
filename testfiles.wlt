BeginTestSection["testfiles_notebook"]

VerificationTest[(* 1 *)
	ListLogPlot[List[Transpose[List[Part[kdebug, All, 3], Part[kdebug, All, 4]]]], Rule[ImageSize, Scaled[0.25`]]]
	,
	ListLogPlot[List[OldData], Rule[ImageSize, Scaled[0.25`]]]	
	,
	SameTest->(Rasterize[#1]==Rasterize[#2]&), TestID->"kAdapt test; current calculation = old calculation"
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
]

VerificationTest[(* 5 *)
	Between[Chop[Abs[test2bdispersionCheck[1, 20, 0, False]]], List[Times[-1, Power[10, -10]], Power[10, -10]]]
	,
	True	
]

EndTestSection[]
