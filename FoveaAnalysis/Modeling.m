(* Mathematica Package *)

BeginPackage["FoveaAnalysis`Modeling`", {"HSF`"}];


(* Functions *)
foveaModel::usage = "foveaModel[\[Mu], \[Sigma], \[Gamma], \[Alpha], x] returns the model function for the fovea.";

interpolateFovea::usage = "interpolateFovea[prop, opt] uses the information of a foveal fit property file to interpolate \
the retinal surface that is used in the modelfit. For this, the properties \"VolFile\", \"Center\" and \
\"CentralPixelHeight\" must be present in the property file.";

processOCTFile::usage = "processOCTFile  ";

findCenter::usage = "findCenter[volFile] tries to calculate the center of the fovea.";
monotonicInterpolation::usage = "monotonicInterpolation[vector, opts] interpolates the vector using Steffen monotonic interpolation.";

(* Options *)

Begin[ "`Private`" ];

(* ::Subsection:: *)
(* The mathematical model used *)

model[\[Mu]_, \[Sigma]_, \[Gamma]_, \[Alpha]_, x_] := (x^\[Gamma] * \[Mu] * \[Sigma]^2) * Exp[- x^\[Gamma] * \[Mu]] + \[Alpha] (1 - Exp[- x^\[Gamma] * \[Mu]]);


foveaModel[m_, s_, g_, a_, r_] := model[m, s, g, a, r];
foveaModel[{m_, s_, g_, a_}, r_] := model[m, s, g, a, r];

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

$samplingPoints = 255;

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

DistributeDefinitions[modelCFunc, $samplingPoints];


(* ::Subsection:: *)
(* Interpolating a fovea *)

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

garwayHeathRescale[corneaAntR1_, ametropia_] := 1 / (17.21 / corneaAntR1 + 1.247 + ametropia / 17.455);


Options[interpolateFovea] = {
    "InterpolateParallel" -> False
};


interpolateFovea[prop_Association, opts : OptionsPattern[]] := With[{file = prop["VolFile"]},
    interpolateFovea[file, prop, opts] /; Not[MissingQ[file]]
];
interpolateFovea::missV = "Missing center or central height value in property file.";
interpolateFovea[file_?HSFFileQ, prop_Association, OptionsPattern[]] := Module[
    {
        center = prop["Center"],
        rescaleQ = prop["RescaleOCTMagnification"],
        centralHeight = prop["CentralPixelHeight"],
        header = HSFInfo[file],
        data = {"ILM", "RPE"} /. HSFLayerSegmentation[file],
        optParallel, scaleZ, sx, sy, scanFocus
    },

    If[MissingQ[center] || MissingQ[centralHeight],
        Message[interpolateFovea::missV];
        Abort[];
    ];
    {scaleZ, scanFocus} = {"ScaleZ", "ScanFocus"} /. header;
    data = (prepareSegmentationDataCompiled[##, centralHeight, scaleZ])& @@ Transpose[data];
    (* We need to be able to give the correct scaling of the OCT scan by ourselves *)
    (* possibly introducing a better altorithm for calculating the correct projection size on the retina *)
    If[TrueQ[rescaleQ] && NumericQ[prop["CorneaAntR1"]] && NumericQ[scanFocus],
        Module[{q = garwayHeathRescale[prop["CorneaAntR1"], scanFocus]},
        (* Attention: there is no way to know how large the scanned OCT region in degree was, because this is not *)
        (* stored in the file. Need to hardcode 20 degree here TODO: I actually CAN calculate the degree *)
            {sx, sy} = 1.02302 * q * 20 / (Reverse[Dimensions[data]] - {0, 1});
        ],
        {sx, sy} = { "ScaleX" , "Distance"} /. header
    ];
    optParallel = OptionValue["InterpolateParallel"];
    If[ optParallel =!= False && Head[optParallel] === Symbol,
        interpolateFovea[data, center, {sy, sx}, optParallel],
        interpolateFovea[data, center, {sy, sx}]
    ]
];

interpolateFovea[data_, {cy_, cx_}, {sy_, sx_}] := Module[
    {
        nx, ny
    },
    {ny, nx} = Dimensions[data];
    ListInterpolation[data, {sy ({1, ny} - cy), sx ({1, nx} - cx)}, Method -> "Spline"]
];

interpolateFovea[ddata_, {ccy_, ccx_}, {ssy_, ssx_}, sym_Symbol] := Module[
    {
        nnx, nny
    },
    {nny, nnx} = Dimensions[ddata];
    Function[{nx, ny, sx, sy, cx, cy, data},
        ParallelEvaluate[
            sym = ListInterpolation[data, {sy ({1, ny} - cy), sx ({1, nx} - cx)}, Method -> "Spline"]]
    ][nnx, nny, ssx, ssy, ccx, ccy, ddata];
    sym
];

(* ::Section:: *)
(* Finding the foveal center *)

retinalHeight := retinalHeight = Compile[
    {
        {rpe, _Real, 1}, {ilm, _Real, 1}
    },
    With[
        {
            vector = Transpose[{rpe, ilm, rpe - ilm}]
        },
        If[Abs[Compile`GetElement[#, 1]] > 10^3 ||
            Abs[Compile`GetElement[#, 2]] > 10^3 ||
            Abs[Compile`GetElement[#, 3]] > 496, -1,
            Compile`GetElement[#, 3]] & /@ vector],
    CompilationTarget -> "C",
    RuntimeAttributes -> {Listable},
    RuntimeOptions -> "Speed",
    Parallelization -> True
];

findCenter[vol_?HSFFileQ] := findCenter @@ Transpose[{"RPE", "ILM"} /.
    HSFLayerSegmentation[vol]
];

findCenter[rpe_?(MatrixQ[#, NumberQ] &), ilm_?(MatrixQ[#, NumberQ] &)] := Module[
    {
        x, y, min = Infinity, res = $Failed, xmin = 1, xmax, ymin = 1, ymax,
        data = GaussianFilter[retinalHeight[rpe, ilm], Dimensions[rpe] / 50.]
    },
    {ymax, xmax} = Dimensions[data];
    {xmin, xmax} = Round[{0.25 * (xmax - xmin), 0.75 * (xmax - xmin)}];
    Do[
        If[min > data[[y, x]],
            min = data[[y, x]];
            res = {y, x}
        ],
        {y, ymin, ymax, 1},
        {x, xmin, xmax, 1}
    ];
    If[Not[MatchQ[res, {_Integer, _Integer}]],
        Throw[findCenter],
        <|"Center" -> res, "CentralPixelHeight" ->  Extract[rpe, res] - Extract[ilm, res]|>
    ]
];


(* ::Subsection:: *)
(*Processing the fitting of an a fovea OCT dataset*)

Options[processOCTFile] = {
    "PropertiesPath" -> Automatic,
    "ParameterRanges" -> {{0.01, 12.0}, {0.01, 2.0}, {1.0, 10.0}, {-1, 1}},
    "AngleStepSize" -> Pi / 2,
    "MaxRadius" -> 2.0,
    "PreferPropertyFile" -> False,
    "InitialPoints" -> Automatic,
    "OCTPath" -> Automatic
};

processOCTFile::noOCT = "Could not find OCT file";
processOCTFile::moreOCT = "More than one matching OCT file was found. Please specify which one you want to use.";
processOCTFile::nocentf = ",The properties are missing the center of the fovea.";
processOCTFile::nocent = "`1` is neither a string denoting the path to the \
center file nor a tuple {cx,cy} of integers defining the center directly.";
processOCTFile::wrongout = "`1` is not a valid output directory.";
processOCTFile::pathNoString = "The path `1` should be a string. Using the default location.";
processOCTFile::noCenter = "The properties file `1` does not contain a valid \"Center\" required for fitting a fovea.";
processOCTFile::wrongParam = "Value for parameter `1` cannot be `2`";
processOCTFile::noWrite = "Unable to write property file `1` to disk.";

processOCTFile[prop_Association, opts : OptionsPattern[]] := Module[{file, octPath},
    octPath = OptionValue["OCTPath"];
    file = prop["VolFile"];
    If[FileExistsQ[file] && HSFFileQ[file],
        processOCTFile[file, prop, opts],
        If[DirectoryQ[octPath],
            With[{filesearch = FileNames[FileBaseName[file] <> ".vol", {octPath}, Infinity]},
                Switch[Length[filesearch],
                    0, Message[processoctfile::nooct]; Abort[],
                    1, processOCTFile[First[filesearch], prop, opts],
                    _, Message[processOCTFile::moreOCT]; Abort[]
                ]
            ]
        ]
    ]
];

processOCTFile[
    file_?HSFFileQ,
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
        Message[processOCTFile::noCenter];
        Return[$Failed]
    ];
    maxRadius = If[preferPropQ, prop["MaxRadius"], OptionValue["MaxRadius"]];
    If[Not[NumericQ[maxRadius]] || maxRadius <= 0,
        Message[processOCTFile::wrongParam, "MaxRadius", maxRadius];
        Return[$Failed]
    ];
    center = prop["Center"];
    If[Not[MatchQ[center, {_Integer, _Integer}]],
        Message[processOCTFile::wrongParam, "Center", center];
        Return[$Failed]
    ];
    angleStepSize = If[preferPropQ, prop["AngleStepSize"], OptionValue["AngleStepSize"]];
    If[Not[NumericQ[angleStepSize]] || angleStepSize <= 0 || angleStepSize >= 2Pi,
        Message[processOCTFile::wrongParam, "AngleStepSize", angleStepSize];
        Return[$Failed]
    ];
    parameterRanges = If[preferPropQ, prop["ParameterRanges"], OptionValue["ParameterRanges"]];
    If[Not[MatrixQ[parameterRanges]] || Dimensions[parameterRanges] != {4, 2},
        Message[processOCTFile::wrongParam, "ParameterRanges", parameterRanges];
        Return[$Failed]
    ];

    internalProcessOCTFile[file, center, angleStepSize, maxRadius, parameterRanges, prop, initPoints]

];


processOCTFile::wrres = "Fit could not be calculated correctly for file `1`.";
processOCTFile::winit = "Optionvalue for \"InitialPoints\" should be either Automatic or a list of points {{_,_,_,_}..} in the parameter space.";

internalProcessOCTFile[
    file_?HSFFileQ,
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

        foveaInterpol = interpolateFovea[file, prop];
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
            dataPoints = Table[getTargetFuncPoints[phi, foveaInterpol, xend, $samplingPoints], {phi, 0, 2 Pi - dphi, N[dphi]}];
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
            Join[prop,
                Association[
                    "Parameters" -> res[[All, 2, All, 2]],
                    "Errors" -> res[[All, 1]],
                    "AngleStepSize" -> N[dphi],
                    "MaxRadius" -> N[xend],
                    "ParameterRanges" -> parameterRanges
                ]
            ]
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

(* ::Section:: *)
(* Monotonic interpolants *)

Options[monotonicInterpolation] := {
    PeriodicInterpolation -> False
}

steffenEnds[{{h1_, h2_}, {d1_, d2_}}] :=
    With[{p = d1 + h1 (d1 - d2)/(h1 + h2)}, (Sign[p] + Sign[d1]) Min[
        Abs[p]/2, Abs[d1]]]


monotonicInterpolation[data_?(VectorQ[#, NumericQ] &), opts___?OptionQ] :=
    monotonicInterpolation[Transpose[{Range[Length[data]], data}],
        opts];
monotonicInterpolation[data_?MatrixQ, OptionsPattern[]] :=
    Module[{dTrans = Transpose[data], del, h, m, pp, optPeriodic,
        overhangs},
        optPeriodic = OptionValue[PeriodicInterpolation];
        h = Differences[First[dTrans]];
        del = Differences[Last[dTrans]]/h;
        overhangs = If[optPeriodic === False, {1, -1}, {-1, 1}];
        (* Note that overhangs in Partition and ListConvolve are defined differently*)
        pp = Dot @@@
            Transpose[
                MapAt[Reverse, Map[Partition[#, 2, 1, {-1, 1}] &, {h, del}], {1, All}]]/
            ListConvolve[{1, 1}, h, -1*overhangs];
        If[optPeriodic === True,
            del = ArrayPad[del, 1, "Periodic"]
        ];
        m = ListConvolve[{1, 1}, 2 UnitStep[del] - 1] *
            MapThread[Min, {Partition[Abs[del], 2, 1], Abs[pp]/2}];
        Interpolation[
            {{#1}, ##2} & @@@ Transpose[Append[dTrans,
                If[optPeriodic === True,
                    m,
                    Flatten[{
                        steffenEnds[#[[{1, 2}]] & /@ {h, del}],
                        m,
                        steffenEnds[#[[{-1, -2}]] & /@ {h, del}]
                    }]
                ]
            ]],
            PeriodicInterpolation -> optPeriodic]
    ]


End[];

EndPackage[];
