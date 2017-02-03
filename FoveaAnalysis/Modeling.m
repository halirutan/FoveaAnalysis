(* Mathematica Package *)

BeginPackage["Modeling`"];


(* Functions *)
FoveaModel::usage = "FoveaModel[\[Mu], \[Sigma], \[Gamma], \[Alpha], x] returns the model function for the fovea.";

InterpolateFovea::usage = "InterpolateFovea[prop, opt] uses the information of a foveal fit property file to interpolate \
the retinal surface that is used in the modelfit. For this, the properties \"VolFile\", \"Center\" and \
\"CentralPixelHeight\" must be present in the property file.";

ProcessOCTFile::usage = "ProcessOCTFile  ";

IPCUFoveaCenterSegmentation::usage = "IPCUFoveaCenterSegmentation";
InterpolateParallel::usage = "InterpolateParallel is an option for InterpolateFovea.";
IPCUInterpolateLayer::usage;

(* Options *)

Begin[ "`Private`" ];

(* ::Subsection:: *)
(* The mathematical model used *)

model[\[Mu]_, \[Sigma]_, \[Gamma]_, \[Alpha]_, x_] := (x^\[Gamma] * \[Mu] * \[Sigma]^2) * Exp[- x^\[Gamma] * \[Mu]] + \[Alpha] (1 - Exp[- x^\[Gamma] * \[Mu]]);


FoveaModel[m_, s_, g_, a_, r_] := model[m, s, g, a, r];
FoveaModel[{m_, s_, g_, a_}, r_] := model[m, s, g, a, r];

DistributeDefinitions[model];

(* For fast processing I compile the model function *)
modelC := (modelC = Compile[{{m, _Real, 0}, {s, _Real, 0}, {g, _Real,
    0}, {a, _Real, 0}, {x, _Real, 0}}, # ,
    CompilationTarget -> "C",
    Parallelization -> True,
    RuntimeAttributes -> {Listable},
    RuntimeOptions -> "Speed"
]&[model[m, s, g, a, x]]);

(* some wrappers and the target function for fitting the model onto a given 
   InterpolatingFunction. See IPCUProcessOCTFile[] 
*)
modelCFunc[args__?NumericQ] := modelC[args];

(* ::Subsection:: *)
(*Creating the target function to use the model in a NMinmize fit*)

$DifferenceSamplingPoints = 255;

With[{expr = model[m, s, g, a, x]},
    ParallelEvaluate[
        SetSystemOptions["ParallelOptions" -> "ParallelThreadNumber" -> $ProcessorCount];
        sqrDiff := (sqrDiff = Compile[
            {{m, _Real, 0}, {s, _Real, 0}, {g, _Real, 0}, {a, _Real, 0}, {x, _Real, 0}, {fov, _Real, 0}},
            (expr - fov)^2,
            CompilationTarget -> "C",
            Parallelization -> True,
            RuntimeAttributes -> {Listable},
            RuntimeOptions -> "Speed"])
    ];
];

getTargetFuncPoints[phi_, fovea_, xend_, nSamplingPoints_] := Transpose[
    Table[{x, fovea[x * Sin[phi], x * Cos[phi]]}, {x, 0, xend, N[xend / (nSamplingPoints - 1)]}]
];

ParallelEvaluate[
    Module[{targetFuncMem, xs, ys, nSamples},
        initTargetFunc[pts : {xValues_, yValues_}] := (
            {xs, ys} = pts;
            nSamples = Length[xs];
            targetFuncMem[m_, s_, g_, a_] := targetFuncMem[m, s, g, a] = Total@sqrDiff[m, s, g, a, xs, ys] / N[nSamples];
        );
        callTargetFunc[m_?NumericQ, s_?NumericQ, g_?NumericQ, a_?NumericQ] := targetFuncMem[m, s, g, a];
        clearTargetFunc[] := ClearAll[targetFuncMem, xs, ys, nSamples];
    ]
];

DistributeDefinitions[modelCFunc, $DifferenceSamplingPoints];


(* ::Subsection:: *)
(* Interpolating a fovea *)

heyexRawFileQ[file_String] := TrueQ[FileExistsQ[file] && FileExtension[file] === "vol"];
heyexRawFileQ[___] := False;

With[{selector = If[Abs[#1] < 1000 && Abs[#2] > 1000, #1, #2] &},
    prepareSegmentationDataCompiled := prepareSegmentationDataCompiled =
            Compile[{{ilmm, _Real, 1}, {rpee, _Real, 1}, {centralValue, _Real, 0}, {scaleZ, _Real, 0}},
                Module[{rpe = rpee, ilm = ilmm, v = 0.0},
                    rpe = FoldList[selector, Reverse[FoldList[selector, Reverse[rpe]]]];
                    ilm = FoldList[selector, Reverse[FoldList[selector, Reverse[ilm]]]];
                    (If[Abs[#] > 496, 0, #] & /@ ((rpe - ilm) - centralValue)) * scaleZ
                ],
                CompilationTarget -> "C",
                RuntimeAttributes -> {Listable},
                RuntimeOptions -> "Speed",
                Parallelization -> True
            ];
];

GarwayHeathRescale[corneaAntR1_, ametropia_] := 1 / (17.21 / corneaAntR1 + 1.247 + ametropia / 17.455);


Options[InterpolateFovea] = {
    InterpolateParallel -> False
};


InterpolateFovea[prop_Association, opts : OptionsPattern[]] := With[{file = prop["VolFile"]},
    InterpolateFovea[file, prop, opts] /; Not[MissingQ[file]]
];
InterpolateFovea::missV = "Missing center or central height value in property file.";
InterpolateFovea[file_?heyexRawFileQ, prop_Association, OptionsPattern[]] := Module[
    {
        center = prop["Center"],
        rescaleQ = prop["RescaleOCTMagnification"],
        centralHeight = prop["CentralPixelHeight"],
        header = Import[file, {"Heyex", "FileHeader"}],
        data = {"ILM", "RPE"} /. Import[file, {"Heyex", "SegmentationData"}],
        optParallel, scaleZ, sx, sy, scanFocus
    },

    If[MissingQ[center] || MissingQ[centralHeight],
        Message[InterpolateFovea::missV];
        Abort[];
    ];
    scaleZ = "ScaleZ" /. header;
    scanFocus = "ScanFocus" /. header;
    data = (prepareSegmentationDataCompiled[##, centralHeight, scaleZ])& @@ Transpose[data];
    (* We need to be able to give the correct scaling of the OCT scan by ourselves *)
    (* possibly introducing a better altorithm for calculating the correct projection size on the retina *)
    If[TrueQ[rescaleQ] && NumericQ[prop["CorneaAntR1"]] && NumericQ[scanFocus],
        Module[{q = GarwayHeathRescale[prop["CorneaAntR1"], scanFocus]},
        (* Attention: there is no way to know how large the scanned OCT region in degree was, because this is not *)
        (* stored in the file. Need to hardcode 20 degree here TODO: I actually CAN calculate the degree *)
            {sx, sy} = 1.02302 * q * 20 / (Reverse[Dimensions[data]] - {0, 1});
        ],
        {sx, sy} = { "ScaleX" , "Distance"} /. header
    ];
    optParallel = OptionValue[InterpolateParallel];
    If[ optParallel =!= False && Head[optParallel] === Symbol,
        InterpolateFovea[data, center, {sx, sy}, optParallel],
        InterpolateFovea[data, center, {sx, sy}]
    ]
];

InterpolateFovea[data_, {cx_, cy_}, {sx_, sy_}] := Module[{nx, ny},
    {ny, nx} = Dimensions[data];
    ListInterpolation[data, {sy ({1, ny} - cy),
        sx ({1, nx} - cx)}, Method -> "Spline"]
];

InterpolateFovea[ddata_, {ccx_, ccy_}, {ssx_, ssy_}, sym_Symbol] := Module[{nnx, nny},
    {nny, nnx} = Dimensions[ddata];
    Function[{nx, ny, sx, sy, cx, cy, data},
        ParallelEvaluate[
            sym = ListInterpolation[data, {sy ({1, ny} - cy), sx ({1, nx} - cx)}, Method -> "Spline"]]
    ][nnx, nny, ssx, ssy, ccx, ccy, ddata];
    sym
];



(* ::Subsection:: *)
(*Processing the fitting of an a fovea OCT dataset*)

Options[ProcessOCTFile] = {
    "PropertiesPath" -> Automatic,
    "ParameterRanges" -> {{0.01, 12.0}, {0.01, 2.0}, {1.0, 10.0}, {-1, 1}},
    "AngleStepSize" -> Pi / 2,
    "MaxRadius" -> 2.0,
    "PreferPropertyFile" -> False,
    "InitialPoints" -> Automatic,
    "OCTPath" -> Automatic
};

ProcessOCTFile::noOCT = "Could not find OCT file";
ProcessOCTFile::moreOCT = "More than one matching OCT file was found. Please specify which one you want to use.";
ProcessOCTFile::nocentf = ",The properties are missing the center of the fovea.";
ProcessOCTFile::nocent = "`1` is neither a string denoting the path to the \
center file nor a tuple {cx,cy} of integers defining the center directly.";
ProcessOCTFile::wrongout = "`1` is not a valid output directory.";
ProcessOCTFile::pathNoString = "The path `1` should be a string. Using the default location.";
ProcessOCTFile::noCenter = "The properties file `1` does not contain a valid \"Center\" required for fitting a fovea.";
ProcessOCTFile::wrongParam = "Value for parameter `1` cannot be `2`";
ProcessOCTFile::noWrite = "Unable to write property file `1` to disk.";

ProcessOCTFile[prop_Association, opts: OptionsPattern[]] := Module[{file, octPath},
    octPath = OptionValue["OCTPath"];
    file = prop["VolFile"];
    If[FileExistsQ[file] && heyexRawFileQ[file],
        ProcessOCTFile[file, prop, opts],
        If[DirectoryQ[octPath],
            With[{filesearch = FileNames[FileBaseName[file] <> ".vol", {octPath}, Infinity]},
                Switch[Length[filesearch],
                    0, Message[processoctfile::nooct]; Abort[],
                    1, ProcessOCTFile[First[filesearch], prop, opts],
                    _, Message[ProcessOCTFile::moreOCT]; Abort[]
                ]
            ]
        ]
    ]
];

ProcessOCTFile[
    file_String /; FileExistsQ[file] && StringMatchQ[file, __ ~~ ".vol"],
    prop_Association, opts : OptionsPattern[]] := Module[
    {
        maxRadius,
        center,
        angleStepSize,
        parameterRanges,
        preferPropQ = OptionValue["PreferPropertyFile"],
        initPoints
    },

    initPoints = OptionValue["InitialPoints"];

    If[MissingQ[prop["Center"]],
        Message[ProcessOCTFile::noCenter];
        Return[$Failed]
    ];
    maxRadius = If[preferPropQ, prop["MaxRadius"], OptionValue["MaxRadius"]];
    If[Not[NumericQ[maxRadius]] || maxRadius <= 0,
        Message[ProcessOCTFile::wrongParam, "MaxRadius", maxRadius];
        Return[$Failed]
    ];
    center = prop["Center"];
    If[Not[MatchQ[center, {_Integer, _Integer}]],
        Message[ProcessOCTFile::wrongParam, "Center", center];
        Return[$Failed]
    ];
    angleStepSize = If[preferPropQ, prop["AngleStepSize"], OptionValue["AngleStepSize"]];
    If[Not[NumericQ[angleStepSize]] || angleStepSize <= 0 || angleStepSize >= 2Pi,
        Message[ProcessOCTFile::wrongParam, "AngleStepSize", angleStepSize];
        Return[$Failed]
    ];
    parameterRanges = If[preferPropQ, prop["ParameterRanges"], OptionValue["ParameterRanges"]];
    If[Not[MatrixQ[parameterRanges]] || Dimensions[parameterRanges] != {4, 2},
        Message[ProcessOCTFile::wrongParam, "ParameterRanges", parameterRanges];
        Return[$Failed]
    ];

    internalProcessOCTFile[file, center, angleStepSize, maxRadius, parameterRanges, prop, initPoints]

];


ProcessOCTFile::wrres = "Fit could not be calculated correctly for file `1`.";
ProcessOCTFile::winit = "Optionvalue for \"InitialPoints\" should be either Automatic or a list of points {{_,_,_,_}..} in the parameter space.";

internalProcessOCTFile[
    file_String /; FileExistsQ[file] && StringMatchQ[file, __ ~~ ".vol"],
    center : {_Integer, _Integer},
    dphi_?NumericQ,
    xe_ /; 0 < xe < 5,
    parameterRanges_,
    prop_Association,
    deInit_] :=
        Module[
            {
                foveaInterpol,
                minm, maxm, mins, maxs, ming, maxg, mina, maxa, (* min and max values for the parameters *)
                xend = xe,
                initialPoints = prop["Parameters"],
                dataPoints, data, phi, res, method
            },

            foveaInterpol = InterpolateFovea[file, prop];
            {{minm, maxm}, {mins, maxs}, {ming, maxg}, {mina, maxa}} = parameterRanges;

            (* Differential Evolution needs some good starting agents. If the prop file already contains an old fit, we will use 20 *)
            (* of them. If there are less then 20 or none in the prop-file, we will use random values which fulfill our constrains *)
            (* of the parameter ranges. *)
            If[ Not[MatchQ[deInit, {{_, _, _, _}..}]],
                With[{randPoints = Transpose[RandomReal[#, 20]& /@ parameterRanges]},
                    If[MatchQ[initialPoints, {{_, _, _, _}..}],
                        initialPoints = Join[initialPoints, Take[randPoints, Max[0, 20 - Length[initialPoints]]]],
                        initialPoints = randPoints
                    ]
                ], (* else *)
                initialPoints = deInit;
            ];

            method = {"DifferentialEvolution", "InitialPoints" -> initialPoints, "ScalingFactor" -> .8, "CrossProbability" -> .1};

            DistributeDefinitions[minm, mins, ming, mina, maxm, maxs, maxg, maxa, method];
            Block[{m, s, g, a, outputFile},
                dataPoints = Table[getTargetFuncPoints[phi, foveaInterpol, xend, $DifferenceSamplingPoints], {phi, 0, 2 Pi - dphi, N[dphi]}];
                res = ParallelTable[
                    Module[{optimum},
                        initTargetFunc[data];
                        optimum = NMinimize[{callTargetFunc[m, s, g, a],
                            And[
                                minm < m < maxm,
                                mins < s < maxs,
                                ming < g < maxg,
                                mina < a < maxa
                            ]
                        }, {m, s, g, a}, Method -> method ];
                        clearTargetFunc[];
                        optimum
                    ], {data, dataPoints}];
                (* Clean the subkernels from debris *)
                ParallelEvaluate[Remove[minm, mins, ming, mina, maxm, maxs, maxg, maxa]];
                prop[{
                    "Parameters" -> res[[All, 2, All, 2]],
                    "Errors" -> res[[All, 1]],
                    "Center" -> center,
                    "AngleStepSize" -> N[dphi],
                    "MaxRadius" -> N[xend],
                    "ParameterRanges" -> parameterRanges
                }][Overwrite];
                {prop["FileName"], res[[All, 1]]}
            ]
        ];



(* Takes one segmentation of the RPE and one of the ILM and calculates their difference *)
(* With this, we get rid of any global slope in the data *)
$foveaHeight := $foveaHeight = Compile[{{ilm, _Real, 1}, {rpe, _Real, 1}},
    With[{vector = rpe - ilm}, If[Abs[#] > 496, 0, #] & /@ vector ],
    CompilationTarget -> "C",
    RuntimeAttributes -> {Listable},
    RuntimeOptions -> "Speed",
    Parallelization -> True
];

$foveaHeighAsPixelValue := $foveaHeighAsPixelValue = Compile[
    {
        {ilm, _Real, 1},
        {rpe, _Real, 1}
    },
    With[{vector = rpe - ilm}, If[Abs[#] > 496, {0, 0, 0}, {Mod[10 # / 496, 1], 1, 1}] & /@ vector ],
    CompilationTarget -> "C",
    RuntimeAttributes -> {Listable},
    RuntimeOptions -> "Speed",
    Parallelization -> True
];

$foveaHeighAsPixelValueWithCenter := $foveaHeighAsPixelValueWithCenter = Compile[
    {
        {ilm, _Real, 1},
        {rpe, _Real, 1},
        {centValue, _Real, 0}
    },
    With[{vector = (rpe - ilm) - centValue}, If[Abs[#] > 496, {0, 0, 0}, {Mod[10 # / 496, 1], 1, 1}] & /@ vector ],
    CompilationTarget -> "C",
    RuntimeAttributes -> {Listable},
    RuntimeOptions -> "Speed",
    Parallelization -> True
];

End[];

EndPackage[];
