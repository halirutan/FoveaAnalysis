(* Mathematica Package *)

Package["FoveaAnalysis`"]
PackageImport["HSF`"]

(* ::Subsection:: *)
(* The mathematical model used *)

model[m_, s_, g_, a_, x_] := (x^g * m * s^2) * Exp[- x^g * m] + a (1 - Exp[- x^g * m]);

PackageExport["FoveaModel"]
FoveaModel::usage = "FoveaModel[\[Mu], \[Sigma], \[Gamma], \[Alpha], x] returns the model function for the fovea.";
FoveaModel[m_, s_, g_, a_, r_] := model[m, s, g, a, r];
FoveaModel[{m_, s_, g_, a_}, r_] := model[m, s, g, a, r];

DistributeDefinitions[model];
Block[
    {m, s, g, a, x},
    (* For fast processing I compile the model function *)
    modelC := (modelC = Compile[{{m, _Real, 0}, {s, _Real, 0}, {g, _Real,
        0}, {a, _Real, 0}, {x, _Real, 0}}, # ,
        CompilationTarget -> "C",
        Parallelization -> True,
        RuntimeAttributes -> {Listable},
        RuntimeOptions -> "Speed"
    ]&[model[m, s, g, a, x]])
];

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

getTargetFuncPoints[phi_, fovea_, rEnd_, nSamplingPoints_] := Transpose[
    Table[{r, fovea[r * Sin[phi], r * Cos[phi]]}, {r, 0, rEnd, N[rEnd / (nSamplingPoints - 1)]}]
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


PackageExport["OCTInterpolation"]
OCTInterpolation::usage = "OCTInterpolation[prop, opt] uses the information of a foveal fit property file to interpolate the retinal surface that is used in the model-fit." <>
    "For this, the properties \"VolFile\", \"Center\" and \"CentralPixelHeight\" must be present in the property file.";
Options[OCTInterpolation] = {
    "InterpolateParallel" -> False
};

OCTInterpolation[prop_Association, opts : OptionsPattern[]] := With[{file = prop["VolFile"]},
    OCTInterpolation[file, prop, opts] /; Not[MissingQ[file]]
];
OCTInterpolation::missV = "Missing center or central height value in property file.";
OCTInterpolation[file_?HSFFileQ, prop_Association, OptionsPattern[]] := Module[
    {
        center = prop["Center"],
        rescaleQ = prop["RescaleOCTMagnification"],
        centralHeight = prop["CentralPixelHeight"],
        header = HSFInfo[file],
        data = {"ILM", "BM"} /. HSFLayerSegmentation[file],
        optParallel, scaleZ, sx, sy, scanFocus
    },

    If[MissingQ[center] || MissingQ[centralHeight],
        Message[OCTInterpolation::missV];
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
        OCTInterpolation[data, center, {sy, sx}, optParallel],
        OCTInterpolation[data, center, {sy, sx}]
    ]
];

OCTInterpolation[data_, {cy_, cx_}, {sy_, sx_}] := Module[
    {
        nx, ny
    },
    {ny, nx} = Dimensions[data];
    ListInterpolation[data, {sy ({1, ny} - cy), sx ({1, nx} - cx)}, Method -> "Spline"]
];

OCTInterpolation[ddata_, {ccy_, ccx_}, {ssy_, ssx_}, sym_Symbol] := Module[
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


PackageExport["FindFoveaCenter"]
FindFoveaCenter::usage = "FindFoveaCenter[volFile] tries to calculate the center of the fovea.";
FindFoveaCenter[vol_?HSFFileQ] := FindFoveaCenter @@ Transpose[{"BM", "ILM"} /.
    HSFLayerSegmentation[vol]
];

FindFoveaCenter[rpe_?(MatrixQ[#, NumberQ] &), ilm_?(MatrixQ[#, NumberQ] &)] := Module[
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
        Throw[FindFoveaCenter],
        <|"Center" -> res, "CentralPixelHeight" -> Extract[rpe, res] - Extract[ilm, res]|>
    ]
];



(* ::Subsection:: *)
(*Processing the fitting of an a fovea OCT dataset*)

PackageExport["FindFoveaModelParameters"]
FindFoveaModelParameters::usage = "FindFoveaModelParameters[properties] finds the best fitting fovea model parameters for the dataset. " <>
    "\"properties\" is an association that must contain a valid entry for \"VolFile\" that specifies the Heidelberg OCT " <>
    "scan exported in raw HSF format that must contain the segmentation for the ILM and the RPE layer. " <>
    "Additionally, the properties need to include a \"Center\" value that specifies the position of the fovea pit.";

Options[FindFoveaModelParameters] = {
    "PropertiesPath" -> Automatic,
    "ParameterRanges" -> {{0.01, 12.0}, {0.01, 2.0}, {1.0, 10.0}, {-1, 1}},
    "AngleStepSize" -> Pi / 2,
    "MaxRadius" -> 2.0,
    "PreferPropertyFile" -> False,
    "InitialPoints" -> Automatic,
    "OCTPath" -> Automatic
};

FindFoveaModelParameters::noOCT = "Could not find OCT file";
FindFoveaModelParameters::moreOCT = "More than one matching OCT file was found. Please specify which one you want to use.";
FindFoveaModelParameters::nocentf = ",The properties are missing the center of the fovea.";
FindFoveaModelParameters::nocent = "`1` is neither a string denoting the path to the \
center file nor a tuple {cx,cy} of integers defining the center directly.";
FindFoveaModelParameters::wrongout = "`1` is not a valid output directory.";
FindFoveaModelParameters::pathNoString = "The path `1` should be a string. Using the default location.";
FindFoveaModelParameters::noCenter = "The properties file `1` does not contain a valid \"Center\" required for fitting a fovea.";
FindFoveaModelParameters::wrongParam = "Value for parameter `1` cannot be `2`";
FindFoveaModelParameters::noWrite = "Unable to write property file `1` to disk.";

FindFoveaModelParameters[prop_Association, opts : OptionsPattern[]] := Module[{file, octPath},
    octPath = OptionValue["OCTPath"];
    file = prop["VolFile"];
    If[FileExistsQ[file] && HSFFileQ[file],
        FindFoveaModelParameters[file, prop, opts],
        If[DirectoryQ[octPath],
            With[{fileSearch = FileNames[FileBaseName[file] <> ".vol", {octPath}, Infinity]},
                Switch[Length[fileSearch],
                    0, Message[FindFoveaModelParameters::nooct]; Abort[],
                    1, FindFoveaModelParameters[First[fileSearch], prop, opts],
                    _, Message[FindFoveaModelParameters::moreOCT]; Abort[]
                ]
            ]
        ]
    ]
];

FindFoveaModelParameters[
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
        Message[FindFoveaModelParameters::noCenter];
        Return[$Failed]
    ];
    maxRadius = If[preferPropQ, prop["MaxRadius"], OptionValue["MaxRadius"]];
    If[Not[NumericQ[maxRadius]] || maxRadius <= 0,
        Message[FindFoveaModelParameters::wrongParam, "MaxRadius", maxRadius];
        Return[$Failed]
    ];
    center = prop["Center"];
    If[Not[MatchQ[center, {_Integer, _Integer}]],
        Message[FindFoveaModelParameters::wrongParam, "Center", center];
        Return[$Failed]
    ];
    angleStepSize = If[preferPropQ, prop["AngleStepSize"], OptionValue["AngleStepSize"]];
    If[Not[NumericQ[angleStepSize]] || angleStepSize <= 0 || angleStepSize >= 2Pi,
        Message[FindFoveaModelParameters::wrongParam, "AngleStepSize", angleStepSize];
        Return[$Failed]
    ];
    parameterRanges = If[preferPropQ, prop["ParameterRanges"], OptionValue["ParameterRanges"]];
    If[Not[MatrixQ[parameterRanges]] || Dimensions[parameterRanges] != {4, 2},
        Message[FindFoveaModelParameters::wrongParam, "ParameterRanges", parameterRanges];
        Return[$Failed]
    ];

    iFindFoveaParameters[file, center, angleStepSize, maxRadius, parameterRanges, prop, initPoints]

];


FindFoveaModelParameters::wrres = "Fit could not be calculated correctly for file `1`.";
FindFoveaModelParameters::winit = "Optionvalue for \"InitialPoints\" should be either Automatic or a list of points {{_,_,_,_}..} in the parameter space.";

iFindFoveaParameters[
    file_?HSFFileQ,
    center : {_Integer, _Integer},
    dphi_?NumericQ,
    xe_ /; 0 < xe < 5,
    parameterRanges_,
    prop_Association,
    deInit_] := Module[
    {
        foveaInterpol,
        minm, maxm, mins, maxs, ming, maxg, mina, maxa, (* min and max values for the parameters *)
        xend = xe,
        initialPoints = prop["Parameters"],
        dataPoints, data, phi, res, method
    },

    foveaInterpol = OCTInterpolation[file, prop];
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

$foveaHeightAsPixelValue := $foveaHeightAsPixelValue = Compile[
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

$foveaHeightAsPixelValueWithCenter := $foveaHeightAsPixelValueWithCenter = Compile[
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

steffenEnds[{{h1_, h2_}, {d1_, d2_}}] := With[{p = d1 + h1 (d1 - d2) / (h1 + h2)}, (Sign[p] + Sign[d1]) Min[ Abs[p] / 2, Abs[d1]]];

PackageExport["SteffenInterpolation"]
SteffenInterpolation::usage = "SteffenInterpolation[vector, opts] interpolates the vector using Steffen monotonic interpolation.";

Options[SteffenInterpolation] := {
    PeriodicInterpolation -> False
};

SteffenInterpolation[data_?(VectorQ[#, NumericQ] &), opts___?OptionQ] := SteffenInterpolation[Transpose[{Range[Length[data]], data}], opts];
SteffenInterpolation[data_?MatrixQ, OptionsPattern[]] := Module[
    {
        dTrans = Transpose[data],
        del,
        h,
        m,
        pp,
        optPeriodic,
        overhangs
    },
    optPeriodic = OptionValue[PeriodicInterpolation];
    h = Differences[First[dTrans]];
    del = Differences[Last[dTrans]] / h;
    overhangs = If[optPeriodic === False, {1, -1}, {-1, 1}];
    (* Note that overhangs in Partition and ListConvolve are defined differently*)
    pp = Dot @@@
        Transpose[
            MapAt[Reverse, Map[Partition[#, 2, 1, {-1, 1}] &, {h, del}], {1, All}]] /
        ListConvolve[{1, 1}, h, -1 * overhangs];
    If[optPeriodic === True,
        del = ArrayPad[del, 1, "Periodic"]
    ];
    m = ListConvolve[{1, 1}, 2 UnitStep[del] - 1] *
        MapThread[Min, {Partition[Abs[del], 2, 1], Abs[pp] / 2}];
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

(* ::Section:: *)
(* Dynamic Visualizations *)

PackageExport["FoveaModelManipulate"]
FoveaModelManipulate::usage = "FoveaModelManipulate[] presents a dynamically adjustable model plot.";
FoveaModelManipulate[] := FoveaModelManipulate[{1.0, 0.5, 2.0, 0.07}];
FoveaModelManipulate[{\[Mu]Init_, \[Sigma]Init_, \[Gamma]Init_, \[Alpha]Init_}] := Block[
    {
        \[Mu], \[Sigma], \[Gamma], \[Alpha], x
    },
    With[
        {
            m = model[\[Mu], \[Sigma], \[Gamma], \[Alpha], x],
            pRanges = "ParameterRanges" /. Options[FindFoveaModelParameters]
        },
        Manipulate[
            Plot[m, {x, 0, 3},
                PlotRange -> {Automatic, {-.2, 0.3}}],
            {{\[Mu], \[Mu]Init, "\[Mu]" }, pRanges[[1, 1]], pRanges[[1, 2]]},
            {{\[Sigma], \[Sigma]Init, "\[Sigma]" }, pRanges[[2, 1]], pRanges[[2, 2]]},
            {{\[Gamma], \[Gamma]Init, "\[Gamma]" }, pRanges[[3, 1]], pRanges[[3, 2]]},
            {{\[Alpha], \[Alpha]Init, "\[Alpha]" }, pRanges[[4, 1]], pRanges[[4, 2]]}]
    ]
];

