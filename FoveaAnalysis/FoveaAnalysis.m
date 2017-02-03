(* Mathematica Package *)
(* Created by Mathematica Plugin for IntelliJ IDEA *)

(* :Title: FoveaAnalysis *)
(* :Context: FoveaAnalysis` *)
(* :Author: patrick *)
(* :Date: 2017-02-03 *)

(* :Package Version: 0.1 *)
(* :Mathematica Version: *)
(* :Copyright: (c) 2017 patrick *)
(* :Keywords: *)
(* :Discussion: *)

BeginPackage["FoveaAnalysis`"];
(* Exported symbols added here with SymbolName::usage *)
(* Mathematica Package *)

BeginPackage[ "IPCU`FoveaModelling`" , { "IPCU`IPCU`", "IPCU`Utilities`", "IPCU`Properties`" }];


(* Functions *)
IPCUFoveaModel::usage = "IPCUFoveaModel[\[Mu], \[Sigma], \[Gamma], \[Alpha], x] returns the model function for the fovea.";

IPCUInterpolateFovea::usage = "IPCUInterpolateFovea[prop, opt] uses the information of a foveal fit property file to interpolate the retinal surface that is used in the modelfit. For this, the properties \"VolFile\", \"Center\" and \"CentralPixelHeight\" must be present in the property file.";
IPCUDynamicModel::usage = "IPCUDynamicModel[] let's you interactively discover the foveal model function.";
IPCUDynamicModelResult::usage = "IPCUDynamicModelResult[volFile, properties] let's you interactively discover the quality of a specific model fit by showing the original data and the fitted model, as well as the fitting error.";
IPCUDynamicFoveaInterpolation::usage = "IPCUDynamicFoveaInterpolation[volFile, properties] displays a dynamic panel to navigate radially around the fovea center to see how the OCT foveal segmentation data looks.";

IPCUProcessOCTFile::usage = "IPCUProcessOCTFile  ";

IPCUFoveaCenter::usage = "IPCUFoveaCenter  ";

IPCUOpenZeissOCT::usage = "IPCUOpenZeissOCTCorrected[file, opts] opens a raw binary file containing an OCT scan made by an experimental Zeiss OCT.";

IPCUConvertZeissOCTToImages::usage = "IPCUConvertZeissOCTToImages[file, path] opens a file with IPCUOpenZeissOCT and stores the created images in path.";


IPCULoadFoveaProperty::usage = "IPCULoadFoveaProperty[file, props, VolFile] loads certain properties from a file. file can either be the full path to the .m file or, when the optional volFile is given (1) a directory or (2) a just a file which is searched in the directory where volFile is. props is is list of properties e.g. \"Center\", \"Parameters\", \"Errors\", \"AngleStepSize\" or \"MaxRadius\".";

IPCUFoveaProjection::usage = "IPCUFoveaProjection[prop] uses an fovea property file to create a Projection of the scan.";
IPCUPreprocessFovea::usage = "IPCUPreprocessFovea[volFile, outputDirectory] interactively checks the quality of a OCT scan and let's you find the foveal center";
IPCUFoveaCenterSegmentation::usage = "IPCUFoveaCenterSegmentation";
IPCUProcessOCTFileNonParallel::usage;
IPCUInterpolateParallel::usage;
IPCUInterpolateLayer::usage;
FoveaProperties::usage = "Stores properties as rules which can be written to files.";

IPCUFoveaCharacteristic::usage = "IPCUFoveaCharacteristic[charact_String] returns a function of foveal parameters and possibly r to calculate a specific foveal characteristic. The complete list of all available characteristics can be accessed by IPCUFovealCharacteristic[].";
IPCUFindFovealRadius::usage = "IPCUFindFovealRadius[{m, s, g, a}, opt] calculates the radius where a certain percentage of the foveal bowl area is reached. The percentage can be adjusted with the \"Percentage\" option (default is 0.95). The foveal bowl area is the area from r=0 to the rim (which is the maximum of the foveal curve).";
IPCUGetDirectionalValues::usage = "IPCUGetDirectionalParameters[prop, \"Parameters\" | \"Errors\"] extracts the 4 anatomical directions from a set radially fitted foveas.";
IPCUGetFoveaCFST::usage = "IPCUGetFoveaCFST[prop, opt] calculated the Central Retinal Thickness as mean of the thicknisses inside the 1mm central radius.";

(* Options *)

Begin[ "`Private`" ];

(* ::Subsection:: *)
(*The mathematical model used*)


(* add another 2g and make 2g:>g, m^2/s^2:>m *)
model[\[Mu]_, \[Sigma]_, \[Gamma]_, \[Alpha]_, x_] :=
        (x^\[Gamma] * \[Mu] * \[Sigma]^2) * Exp[- x^\[Gamma] * \[Mu]] + \[Alpha] (1 - Exp[- x^\[Gamma] * \[Mu]]);


IPCUFoveaModel[m_, s_, g_, a_, r_] := model[m, s, g, a, r];
IPCUFoveaModel[{m_, s_, g_, a_}, r_] := model[m, s, g, a, r];

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
(*Interpolating a fovea*)

heyexRawFileQ[file_String] := TrueQ[FileExistsQ[file] && FileExtension[file] === "vol"];
heyexRawFileQ[___] := False;

(*prepareSegmentationDataCompiled := prepareSegmentationDataCompiled = Compile[*)
(*{*)
(*{ilm, _Real, 1},*)
(*{rpe, _Real, 1},*)
(*{centralValue, _Real, 0},*)
(*{scaleZ, _Real, 0}*)
(*},*)
(*With[{vector = (rpe - ilm) - centralValue},*)
(*(If[Abs[#] > 496, 0, #] & /@ vector) * scaleZ*)
(*],*)
(*CompilationTarget -> "C",*)
(*RuntimeAttributes -> {Listable},*)
(*RuntimeOptions -> "Speed",*)
(*Parallelization -> True*)
(*];*)


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


Options[IPCUInterpolateFovea] = {
    IPCUInterpolateParallel -> False
};


IPCUInterpolateFovea[prop_IPCUProperties, opts : OptionsPattern[]] := With[{file = prop["VolFile"]},
    IPCUInterpolateFovea[file, prop, opts] /; FreeQ[file, Missing]
];
IPCUInterpolateFovea::missV = "Missing center or central height value in property file.";
IPCUInterpolateFovea[file_?heyexRawFileQ, prop_IPCUProperties, OptionsPattern[]] := Module[
    {
        center = prop["Center"],
        rescaleQ = prop["RescaleOCTMagnification"],
        centralHeight = prop["CentralPixelHeight"],
        header = Import[file, {"Heyex", "FileHeader"}],
        data = {"ILM", "RPE"} /. Import[file, {"Heyex", "SegmentationData"}],
        optParallel, scaleZ, sx, sy, scanFocus
    },

    If[Not@FreeQ[{center, centralHeight}, Missing],
        Message[IPCUInterpolateFovea::missV];
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
            {sx, sy} = 1.02302 * q * 20 / (Reverse[Dimensions[data]]-{0,1});
        ],
        {sx, sy} = { "ScaleX" , "Distance"} /. header
    ];
    optParallel = OptionValue[IPCUInterpolateParallel];
    If[ optParallel =!= False && Head[optParallel] === Symbol,
        IPCUInterpolateFovea[data, center, {sx, sy}, optParallel],
        IPCUInterpolateFovea[data, center, {sx, sy}]
    ]
];

IPCUInterpolateFovea[data_, {cx_, cy_}, {sx_, sy_}] :=
        Module[{nx, ny},
            {ny, nx} = Dimensions[data];
            ListInterpolation[data, {sy ({1, ny} - cy),
                sx ({1, nx} - cx)}, Method -> "Spline"]
        ];

IPCUInterpolateFovea[ddata_, {ccx_, ccy_}, {ssx_, ssy_}, sym_Symbol] :=
        Module[{nnx, nny},
            {nny, nnx} = Dimensions[ddata];
            Function[{nx, ny, sx, sy, cx, cy, data},
                ParallelEvaluate[
                    sym = ListInterpolation[data, {sy ({1, ny} - cy), sx ({1, nx} - cx)}, Method -> "Spline"]]
            ][nnx, nny, ssx, ssy, ccx, ccy, ddata];
            sym
        ];



(* ::Subsection:: *)
(*Processing the fitting of an a fovea OCT dataset*)

Options[IPCUProcessOCTFile] = {
    "PropertiesPath" -> Automatic,
    "ParameterRanges" -> {{0.01, 12.0}, {0.01, 2.0}, {1.0, 10.0}, {-1, 1}},
    "AngleStepSize" -> Pi / 2,
    "MaxRadius" -> 2.0,
    "PreferPropertyFile" -> False,
    "InitialPoints" -> Automatic
};


IPCUProcessOCTFile::nocentf = ",The file `1` which should contain the center of \
the fovea does not exist or could not be read.";
IPCUProcessOCTFile::nocent = "`1` is neither a string denoting the path to the \
center file nor a tuple {cx,cy} of integers defining the center directly.";
IPCUProcessOCTFile::wrongout = "`1` is not a valid output directory.";
IPCUProcessOCTFile::pathNoString = "The path `1` should be a string. Using the default location.";
IPCUProcessOCTFile::noCenter = "The properties file `1` does not contain a valid \"Center\" required for fitting a fovea.";
IPCUProcessOCTFile::wrongParam = "Value for parameter `1` cannot be `2`";
IPCUProcessOCTFile::noWrite = "Unable to write property file `1` to disk.";

IPCUProcessOCTFile[
    file_String /; FileExistsQ[file] && StringMatchQ[file, __ ~~ ".vol"],
    prop_IPCUProperties, opts : OptionsPattern[]] := Module[
    {
        maxRadius,
        center,
        angleStepSize,
        parameterRanges,
        preferPropQ = OptionValue["PreferPropertyFile"],
        initPoints
    },

    initPoints = OptionValue["InitialPoints"];

    Check[prop[Overwrite],
        Message[IPCUProcessOCTFile::noWrite, prop["FileName"]];
        Return[$Failed]
    ];

    If[prop["Center"] === Missing,
        Message[IPCUProcessOCTFile::noCenter, prop["FileName"]];
        Return[$Failed]
    ];
    maxRadius = If[preferPropQ, prop["MaxRadius"], OptionValue["MaxRadius"]];
    If[Not[NumericQ[maxRadius]] || maxRadius <= 0,
        Message[IPCUProcessOCTFile::wrongParam, "MaxRadius", maxRadius];
        Return[$Failed]
    ];
    center = prop["Center"];
    If[Not[MatchQ[center, {_Integer, _Integer}]],
        Message[IPCUProcessOCTFile::wrongParam, "Center", center];
        Return[$Failed]
    ];
    angleStepSize = If[preferPropQ, prop["AngleStepSize"], OptionValue["AngleStepSize"]];
    If[Not[NumericQ[angleStepSize]] || angleStepSize <= 0 || angleStepSize >= 2Pi,
        Message[IPCUProcessOCTFile::wrongParam, "AngleStepSize", angleStepSize];
        Return[$Failed]
    ];
    parameterRanges = If[preferPropQ, prop["ParameterRanges"], OptionValue["ParameterRanges"]];
    If[Not[MatrixQ[parameterRanges]] || Dimensions[parameterRanges] != {4, 2},
        Message[IPCUProcessOCTFile::wrongParam, "ParameterRanges", parameterRanges];
        Return[$Failed]
    ];

    internalProcessOCTFile[file, center, angleStepSize, maxRadius, parameterRanges, prop, initPoints]

];

IPCUProcessOCTFile[
    file_String /; FileExistsQ[file] && StringMatchQ[file, __ ~~ ".vol"],
    opts : OptionsPattern[]] := Module[
    {
        prop,
        propPath = OptionValue["PropertiesPath"]
    },

    If[Head[propPath] =!= String,
        If[propPath =!= Automatic,
            Message[IPCUProcessOCTFile::pathNoString, propPath]
        ];
        prop = IPCUProperties[file],
        prop = IPCUProperties[file, propPath]
    ];

    (* Check if it is a valid properties file. We need to know the Center of the fovea at this point *)
    If[prop["Center"] === Missing || Not[MatchQ[prop["Center"], {_Integer, _Integer}]],
        Message[IPCUProcessOCTFile::noCenter, prop["FileName"];
        Return[$Failed]]
    ];
    IPCUProcessOCTFile[file, prop, opts]
];

IPCUProcessOCTFile::wrres = "Fit could not be calculated correctly for file `1`.";
IPCUProcessOCTFile::winit = "Optionvalue for \"InitialPoints\" should be either Automatic or a list of points {{_,_,_,_}..} in the parameter space.";

internalProcessOCTFile[
    file_String /; FileExistsQ[file] && StringMatchQ[file, __ ~~ ".vol"],
    center : {_Integer, _Integer},
    dphi_?NumericQ,
    xe_ /; 0 < xe < 5,
    parameterRanges_,
    prop_IPCUProperties,
    deInit_] :=
        Module[
            {
                foveaInterpol,
                minm, maxm, mins, maxs, ming, maxg, mina, maxa, (* min and max values for the parameters *)
                xend = xe,
                initialPoints = prop["Parameters"],
                dataPoints, data, phi, res, method
            },

            foveaInterpol = IPCUInterpolateFovea[file, prop];
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



(* ::Subsection:: *)
(*A dynamic viewer for the model and the fitted output*)

IPCUDynamicModel[] := IPCUDynamicModel[{1.0, 0.5, 2.0, 0.07}];
IPCUDynamicModel[{\[Mu]Init_, \[Sigma]Init_, \[Gamma]Init_, \[Alpha]Init_}] := Block[
    {\[Mu], \[Sigma], \[Gamma], \[Alpha], x},
    With[
        {
            m = model[\[Mu], \[Sigma], \[Gamma], \[Alpha], x],
            pRanges = "ParameterRanges" /. Options[IPCUProcessOCTFile]
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

IPCUDynamicModelResult[p_IPCUProperties] := IPCUDynamicModelResult[p["VolFile"], p];
IPCUDynamicModelResult[volfile_?heyexRawFileQ, props_IPCUProperties ] :=
        With[
            {
                ip = IPCUInterpolateFovea[volfile, props],
                fp = (ListInterpolation[Append[#, First[#]], {0, 2 Pi}] &) /@ Transpose[props["Parameters"]],
                errorip = ListInterpolation[Append[#, First[#]], {0, 2 Pi}] &[Sqrt[props["Errors"]]],
                xend = props["MaxRadius"],
                ranges = props["ParameterRanges"],
                dphi = props["AngleStepSize"]
            },
            DynamicModule[{\[Mu] = 0, \[Sigma] = 0, \[Gamma] = 0, \[Alpha] = 0,
                cart, \[CurlyPhi] = 0, update},
                cart[{r_, phi_}] := { r * Sin[phi], r * Cos[phi]};
                update[] := {\[Mu], \[Sigma], \[Gamma], \[Alpha]} = Through[fp[\[CurlyPhi]]];
                update[];
                Panel[Column[{
                    Row[{ "\[Mu]" , Spacer[5], Slider[Dynamic[\[Mu]], ranges[[1]]],
                        Spacer[5], Dynamic[\[Mu]]}],
                    Row[{ "\[Sigma]" , Spacer[5],
                        Slider[Dynamic[\[Sigma]], ranges[[2]]], Spacer[5],
                        Dynamic[\[Sigma]]}],
                    Row[{ "\[Gamma]" , Spacer[5], Slider[Dynamic[\[Gamma]], ranges[[3]]],
                        Spacer[5], Dynamic[\[Gamma]]}],
                    Row[{ "\[Alpha]" , Spacer[5], Slider[Dynamic[\[Alpha]], ranges[[4]]],
                        Spacer[5], Dynamic[\[Alpha]]}],
                    Row[{ "\[CurlyPhi]" , Spacer[5],
                    (* We need to fix the slider range in the next line due to a bug in Mathematica 10 on Linux *)
                    (*Slider[Dynamic[\[CurlyPhi], (\[CurlyPhi] = #; update[]; #)&], {0, 2 Pi, 2Pi / 10}], Spacer[5], Dynamic[\[CurlyPhi]]}],*)
                        Slider[Dynamic[\[CurlyPhi], (\[CurlyPhi] = #; update[])&], {0.0, 2.0 Pi, dphi}], Spacer[5], Dynamic[\[CurlyPhi]]}],
                    Spacer[10],

                    Row[{
                        Dynamic[
                            Plot[{model[\[Mu], \[Sigma], \[Gamma], \[Alpha], x], ip @@ cart[{x, \[CurlyPhi]}]}, {x, 0, xend},
                                PlotRange -> {Automatic, {-0.05, 0.3}}, ImageSize -> Scaled[0.5]]],
                        Dynamic[ThermometerGauge[Quantity[10^3 * errorip[\[CurlyPhi]], "Micrometers"], {0, 10}, GaugeLabels -> Full, PerformanceGoal -> "Speed"]]
                    }]
                }]
                ]
            ]
        ];


IPCUDynamicFoveaInterpolation[file_String, props_IPCUProperties] := With[
    {
        ip = IPCUInterpolateFovea[file, props["Center"]]
    },
    DynamicModule[
        {
            cart,
            \[CurlyPhi] = 0,
            xend = Min@Abs[Flatten[First[ip]]]
        },
        cart[{r_, phi_}] := { r * Sin[phi], r * Cos[phi]};
        Panel[Column[{
            Row[{ "\[CurlyPhi]" , Spacer[5],
            (* We need to fix the slider range in the next line due to a bug in Mathematica 10 on Linux *)
            (*Slider[Dynamic[\[CurlyPhi], (\[CurlyPhi] = #; update[]; #)&], {0, 2 Pi, 2Pi / 10}], Spacer[5], Dynamic[\[CurlyPhi]]}],*)
                Slider[Dynamic[\[CurlyPhi]], {0.0, 2.0 Pi}], Spacer[5], Dynamic[\[CurlyPhi]]}],
            Spacer[10],
            Dynamic[
                Plot[ ip @@ cart[{x, \[CurlyPhi]}], {x, 0, xend},
                    PlotRange -> {Automatic, {-0.05, 0.3}}, ImageSize -> Scaled[0.5]]]
        }]
        ]
    ]
];

IPCUFoveaCenter[
    file_String /;
            FileExistsQ[file] && StringMatchQ[file, __ ~~ ".vol"]] :=
        With[{imgs = Import[file, { "Heyex" , "Images" }]},
            DynamicModule[{p, imgNumber = 1, nx, ny},
                {nx, ny} = ImageDimensions[First[imgs]];
                p = {nx, ny} / 2;
                Panel[Column[{Slider[Dynamic[imgNumber], {1, Length[imgs], 1}],
                    LocatorPane[Dynamic[p],
                        Dynamic@Show[{
                            imgs[[imgNumber]],
                            Graphics[{Red, Line[{{p[[1]], 1}, {p[[1]], ny}}]}]
                        }, ImageSize -> 500, AspectRatio -> 1 / 2],
                        Appearance -> None],
                    Dynamic[InputField[Round@{p[[1]], imgNumber}]]}]]
            ]
        ];

IPCUFoveaCenter::noin = "The file `1` does not exist or has not the file-ending \".vol\". Please choose a valid Heidelberg-OCT volume.";
IPCUFoveaCenter::noout = "The output path `1` is not a valid directory.";

IPCUFoveaCenter[file_String
        /; FileExistsQ[file] && StringMatchQ[file, __ ~~ ".vol"]] := Module[
    {outdir = FileNameJoin[{DirectoryName[file], "centers" }]},
    If[Not[DirectoryQ[outdir]],
        CreateDirectory[outdir]
    ];
    IPCUFoveaCenter[file, outdir]
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


Options[IPCUFoveaProjection] = {
    "MaxRadius" -> 2.0,
    ImageSize -> 400
};

IPCUFoveaProjection[prop_IPCUProperties, OptionsPattern[]] := With[
    {
        octFile = prop["VolFile"],
        pixelHeight = prop["CentralPixelHeight"]
    },

    Assert[FileExistsQ[octFile] && pixelHeight =!= Missing];
    Module[
        {
            hdr = Import[octFile, {"Heyex", "FileHeader"}],
            seg = {"ILM", "RPE"} /. Import[octFile, {"Heyex", "SegmentationData"}],
            segPixels,
            cx, cy, scaleX, scaleY, nx, ny,
            maxRadius = OptionValue["MaxRadius"],
            size = OptionValue[ImageSize]
        },
        {cx, cy} = prop["Center"];
        {scaleX, scaleY} = {"ScaleX", "Distance"} /. hdr;
        segPixels = ($foveaHeighAsPixelValueWithCenter[##, pixelHeight])& @@ Transpose[seg];
        {ny, nx} = Most[Dimensions[segPixels]];

        If[Not[NumericQ[maxRadius] && maxRadius > 0],
            maxRadius = 2.0
        ];

        Graphics[{Raster[segPixels, ColorFunction -> Hue], Black,
            AbsoluteThickness[1],

            Dynamic@
                    Line[{{{1, cy - .5}, {nx, cy - .5}}, {{cx, 0}, {cx, ny}}}], Red,
            Dynamic@Point[{cx, cy - .5}],

        (* max Radius Circle*)
            Dashed, Black,
            Dynamic@Circle[{cx, cy - .5}, 1.0 / {scaleX, scaleY} * maxRadius]

        }, AspectRatio -> 1,
            PlotRange -> {{0, nx}, {0, ny}},
            Frame -> True, ImageSize -> size, FrameLabel -> ToString[prop[{"Quality"}]]
        ]
    ]
];

Options[IPCUPreprocessFovea] = {
    "MaxRadius" -> 2.0,
    "FindCenter" -> False
};

(* An interactive tool for verifying the quality of the volume scan and finding the foveal center *)
IPCUPreprocessFovea[
    volFile_String, outDir_String,
    OptionsPattern[]] /; FileExistsQ[volFile] && DirectoryQ[outDir] := With[
    {
        maxRadius = OptionValue["MaxRadius"],
        hdr = Import[volFile, {"Heyex", "FileHeader"}],
        seg = {"ILM", "RPE"} /. Import[volFile, {"Heyex", "SegmentationData"}],
        imgs = Import[volFile, {"Heyex", "Images"}],
        fileBaseName = FileBaseName[volFile],
        size = 400
    },
    With[{
        segData = $foveaHeight @@ Transpose[seg],
        segPixels = $foveaHeighAsPixelValue @@ Transpose[seg],
        scaleX = ("ScaleX" /. hdr),
        scaleY = ("Distance" /. hdr)
    },
        DynamicModule[{
            cx, cy,
            nx, ny,
            updateCenter, getID,
            ip = ListInterpolation[segData, InterpolationOrder -> {2, 2}, Method -> "Spline"],
            ipF,
            minFinder,
            save
        },
            {ny, nx} = Dimensions[segData];
            {cy, cx} = Floor[{ny, nx} / 2];
            ipF[{x_?NumericQ /; 1 <= x <= nx, y_?NumericQ /; 1 <= y <= ny}] :=
                    ip[y, x];

            (* Find the center automatically with a good guess *)

            minFinder[{ix_, iy_}] := DynamicModule[{sol, x, y},
                Quiet@Check[
                    sol = FindMinimum[{ipF[{x, y}], {1 < x < nx, 1 < y < ny}}, {{x, ix}, {y, iy}}, PrecisionGoal -> 1, AccuracyGoal -> 1],
                    MessageDialog["Couldn't find center"];
                    Return[]
                ];
                updateCenter[{x, y} /. Last[sol]]
            ];

            getID[] := With[{regexp = StartOfString ~~ "B_" ~~ id : NumberString ~~ "_" ~~ ___},
                If[StringMatchQ[FileBaseName[volFile], regexp],
                    StringReplace[FileBaseName[volFile], regexp :> id],
                    Missing["NotAvailable"]
                ]
            ];

            save[] := save[True];
            save[quality_] := IPCUProperties[fileBaseName <> ".prop", outDir][
                {
                    "ID" -> getID[],
                    "ScanID" -> ("ID" /. hdr),
                    "Center" -> {cx, cy},
                    "CentralPixelHeight" -> segData[[cy, cx]],
                    "CentralHeight" -> segData[[cy, cx]] * ("ScaleZ" /. hdr),
                    "Eye" -> Switch["ScanPosition" /. hdr, "OD", Right, "OS", Left, _ , $Failed],
                    "VolFile" -> volFile,
                    "Quality" -> quality
                }
            ][Overwrite];

            updateCenter[{newX_, newY_}] := {cx, cy} = Max[1, Min[#1, #2]] & @@@ {{nx, Round[newX]}, {ny, Round[newY]}};

            If[TrueQ[OptionValue["FindCenter"]],
                minFinder[{cx, cy}];
            ];

            Deploy@Panel@
                    Column[{
                        Panel@Grid[{
                            {
                                LocatorPane[Dynamic[{cx, cy}, updateCenter[#]& ],

                                    Graphics[{Raster[segPixels, ColorFunction -> Hue], Black,
                                        AbsoluteThickness[1],

                                        Dynamic@
                                                Line[{{{1, cy - .5}, {nx, cy - .5}}, {{cx, 0}, {cx, ny}}}], Red,
                                        Dynamic@Point[{cx, cy - .5}],

                                    (* max Radius Circle*)
                                        Dashed, Black,
                                        Dynamic@Circle[{cx, cy - .5}, 1.0 / {scaleX, scaleY} * maxRadius]

                                    }, AspectRatio -> 1,
                                        PlotRange -> {{0, nx}, {0, ny}},
                                        Frame -> True, ImageSize -> size
                                    ], Appearance -> None
                                ],

                                Dynamic@
                                        ListLinePlot[segData[[Round[cy]]],
                                            PlotRange -> {Automatic, {0, 250}},
                                            DataRange -> ({-cx, (nx - cx - 1)} * scaleX),
                                            AxesOrigin -> {-cx * scaleX, 0}, Filling -> Axis,

                                            Epilog -> {AbsoluteThickness[1],
                                                Line[{{0, 0}, {0, 300}}],
                                                Dashed,
                                                Line[{{{-maxRadius, 0}, {-maxRadius, 300}}, {{maxRadius, 0}, {maxRadius, 300}}}]},
                                            ImageSize -> {Automatic, size}]

                            }
                            ,
                            {
                                Dynamic@Show[
                                    Graphics[{Raster[Reverse@ImageData[imgs[[Round[cy]]]]]}],
                                    ListLinePlot[496 - seg[[Round[cy]]], PlotRange -> {Automatic, {1, 496}}],
                                    ImageSize -> {size, Automatic}
                                ],

                                Column[{
                                    Style["Filename: " <> fileBaseName, "Text"], Spacer[5],
                                    Style["Eye Position: " <> If[("ScanPosition" /. hdr) == "OD", "Right eye", "Left eye"], "Text"],
                                    Spacer[10],

                                    Row[{Column[{
                                        Row[{
                                            Style["X-Direction", "Text"], Spacer[10],
                                            Button["\[DoubleLeftArrow]", cx = Max[1, cx - 1]],
                                            Spacer[5], Slider[Dynamic[cx], {1, nx, 1}],
                                            Spacer[5], Button["\[DoubleRightArrow]", cx = Min[cx + 1, nx]], Spacer[10],
                                            Dynamic@Style[ToString[cx], "Text"]}],

                                        Row[{
                                            Style["Y-Direction", "Text"], Spacer[10],
                                            Button["\[DoubleLeftArrow]", cy = Max[1, cy - 1]],
                                            Spacer[5], Slider[Dynamic[cy], {1, ny, 1}],
                                            Spacer[5], Button["\[DoubleRightArrow]", cy = Min[ny, cy + 1]], Spacer[10],
                                            Dynamic@Style[ToString[cy], "Text"]}]
                                    }],
                                        Spacer[5],
                                        Tooltip[Button["Find Center", minFinder[{cx, cy}], ImageSize -> {128, 64}, Background -> RGBColor[38 / 255, 139 / 255, 14 / 17]],
                                            "Note: The current position must be near the real center!"]
                                    }, Frame -> True],
                                    Spacer[10],

                                    Row[{
                                        Tooltip[
                                            Button[Style["Bad Quality", RGBColor[44 / 51, 10 / 51, 47 / 255]], save[Missing[Input["Bad quality reason:", "Cutoff"]]];
                                            MessageDialog["Saved!"], Method -> "Queued"], "Marks this file as bad quality"],

                                        Spacer[10],

                                        Tooltip[Button[Style["Save", Bold, RGBColor[38 / 255, 139 / 255, 14 / 17]], save[]; MessageDialog["The center has been saved."]],
                                            "Saves the current center!"]
                                    }]

                                }]

                            }}]

                    }]
        ]
    ]
];

IPCUFoveaCenter[file_String, out_String] := Block[{outpath},
    If[Not[FileExistsQ[file] && StringMatchQ[file, __ ~~ ".vol"]],
        Message[IPCUFoveaCenter::noin, file];
        Return[$Failed]
    ];
    If[Not[DirectoryQ[out]],
        outpath = FileNameJoin[{DirectoryName[file], out}];
        If[Not[DirectoryQ[outpath]],
            CreateDirectory[outpath]
        ],
        outpath = out;
    ];

    (*Implementation of the important stuff*)
    With[{imgs = Import[file, { "Heyex" , "Images" }]},
        DynamicModule[{p, imgNumber = 1, nx, ny, outfile},
            outfile =
                    FileNameJoin[{outpath, FileBaseName[file] <> ".m"}];
            {nx, ny} = ImageDimensions[First[imgs]];
            p = {nx, ny} / 2;
            Panel[
                Column[{Slider[Dynamic[imgNumber], {1, Length[imgs], 1}],
                    LocatorPane[Dynamic[p],
                        Dynamic@
                                Show[{imgs[[imgNumber]],
                                    Graphics[{Red, Line[{{p[[1]], 1}, {p[[1]], ny}}]}]},
                                    ImageSize -> 500, AspectRatio -> 1 / 2],
                        Appearance -> None],
                    Row[{
                        Dynamic[InputField[Round@{p[[1]], imgNumber}]],
                        Tooltip[Button[ "Store Center" ,
                            Export[outfile, Round[{p[[1]], imgNumber}]]],
                            "Stores the center into the file \"" <> outfile "\""]
                    }]
                }]
            ]
        ]
    ]
];


IPCUFoveaCenterSegmentation[file_String, out_String] := Block[{outpath},
    If[Not[FileExistsQ[file] && StringMatchQ[file, __ ~~ ".vol"]],
        Message[IPCUFoveaCenter::noin, file];
        Return[$Failed]];
    If[Not[DirectoryQ[out]],
        outpath = FileNameJoin[{DirectoryName[file], out}];
        If[Not[DirectoryQ[outpath]], CreateDirectory[outpath]],
        outpath = out;];
    (*Implementation of the important stuff*)
    With[{data =
            Map[If[Abs[#] > 1000, 0, #] &,
                Subtract @@@ ({ "RPE" , "ILM" } /.
                        Import[file, { "Heyex" , "SegmentationData" }]), {2}]},
        DynamicModule[{p, nx, ny, outfile, findCenter,
            ip},

            outfile = FileNameJoin[{outpath, FileBaseName[file] <> ".m"}];
            {ny, nx} = Dimensions[data];
            p = {nx, ny} / 2;

            ip = ListInterpolation[GaussianFilter[data, {2, 5}]];

            findCenter[{ix_, iy_}] :=
                    Module[{sol, x, y},
                        Check[
                            sol = FindMinimum[{ip[y, x], {1 < x < nx, 1 < y < ny}}, {{x,
                                ix}, {y, iy}}, PrecisionGoal -> 1, AccuracyGoal -> 1],
                            MessageDialog[ "Error in Minimize" ]
                        ];
                        {x, y} /. Last[sol]
                    ];
            p = findCenter[p];
            Export[outfile, Round[p]];
            Panel[
                Column[{Slider[Dynamic[p[[2]]], {1, ny}],
                    LocatorPane[Dynamic[{p[[1]], p[[2]] / 2}, (p = {#[[1]], #[[2]] / 2}) &],
                        Dynamic@Show[{

                            Plot[ip[p[[2]], x], {x, 1, 512},
                                PlotRange -> {Automatic, {0, 200}}]
                            , Graphics[{Red, PointSize[0.02], Point[p * {1, 2}], Line[{{p[[1]], 1}, {p[[1]], ny}}]}]},
                            ImageSize -> 500, AspectRatio -> 1 / 2], Appearance -> None],
                    Row[{Dynamic[InputField[Round@p]],
                        Tooltip[Button[ "Store Center" , Export[outfile, Round[p]]],
                            "Stores the center into the file \"" <> outfile "\""],
                        Button[ "Find Center" , p = findCenter[p]]

                    }]}]]
        ]]
];


(* ::Section:: *)
(*Extended Polynomial fit for Zemax*)

zemaxExtendePolynomial::usage = "ZemaxExtendePolynomial[x, y, d, a, c, k] creates a polynomial of degree d where the absolute part is defined as an
aspheric lens. The parameter a is the symbol which is used for the factors of the polynomial, c is the scaling of the radius and k is the \"conic section\".";
zemaxExtendePolynomial[x_Symbol, y_Symbol, degree_, a_Symbol, c_Symbol, k_Symbol] :=
        Sum[a[i, j] * x^i * y^j, {i, 0, degree}, {j, 0, degree}] /.
                a[0, 0] :> c (x^2 + y^2) / (1 + Sqrt[1 - (1 + k) c^2 (x^2 + y^2)]);



compileZemaxExtendePolynomial[x_Symbol, y_Symbol, degree_, c_Symbol, k_Symbol, a_Symbol] := Module[
    {
        parm = Flatten[{x, y, c, k, Array[Unique[a] &, (degree + 1)^2]}],
        count = 5
    },
    Function[
        Compile[#1, #2, RuntimeAttributes -> {Listable},
            Parallelization -> True, CompilationTarget :> "C"]][
        {#, _Real, 0} & /@ parm,
        Sum[parm[[count++]] * x^i * y^j, {i, 0, degree}, {j, 0, degree}] /.
                parm[[5]] :>
                        c (x^2 + y^2) / (1 + Sqrt[1 - (1 + k) c^2 (x^2 + y^2)])]
];


(* ::Section:: *)
(*Importing the dev-exported files of Zeiss OCT*)

(* "Meander", "Mean", "First" *)
Options[IPCUOpenZeissOCT] = {
    "PixelSorting" -> "Mean"};

resortC = Compile[{{v, _Integer, 1}},
    {{v[[1]], v[[4]]}, {v[[2]], v[[3]]}},
    CompilationTarget -> "C", RuntimeAttributes -> {Listable},
    Parallelization -> True];


IPCUOpenZeissOCT[file_String, OptionsPattern[]] := Module[{
    pxfunc = OptionValue[ "PixelSorting" ],
    data = BinaryReadList[file, "UnsignedInteger8" ], l},

    l = Length[data];
    If[l === 4 * 5 * 1024^2,
        Image[Reverse@(Reverse /@ #), "Byte" ] & /@
                Switch[pxfunc,
                    "Meander" ,
                    ArrayFlatten /@
                            Partition[Partition[resortC[Partition[data, 4]], 1024], 1024],
                    "Mean" ,
                    Partition[Partition[Map[Mean, Partition[N[data], 4]], 1024],
                        1024],
                    _,
                    Partition[Partition[Map[First, Partition[N[data], 4]], 1024],
                        1024]
                ],
        If[l === 5 * 1024^2,
            Image[Reverse /@ Reverse[#], "Byte" ] & /@
                    Fold[Partition[#1, #2] &, data, {1024, 1024}], $Failed]]];




IPCUConvertZeissOCTToImages[file_String, outdir_String] :=
        Module[{imgs = IPCUOpenZeissOCT[file]},
            If[MatchQ[imgs, {_Image ..}],
                CreateDirectory[FileNameJoin[{outdir, FileBaseName[file]}]];
                MapIndexed[
                    Export[FileNameJoin[{outdir,
                        FileBaseName[file]}] <> $PathnameSeparator <>
                            ToString[NumberForm[#2[[1]], 2, NumberPadding -> { "0" , "" }]] <>
                            ".tif", #1] &, imgs], $Failed]];

(* ::Section:: *)
(*Calculation of Fovea properties from parameter values*)

IPCUFoveaCharacteristic[Properties] = {
    "Slope",
    "Curvature",
    "MaximumPoint",
    "TurningPoint1",
    "TurningPoint2",
    "Steepness",
    "Area",
    "BowlArea"
};
IPCUFoveaCharacteristic[] := IPCUFoveaCharacteristic[Properties];
IPCUFoveaCharacteristic["Slope"] := Block[{m, s, g, a, r}, Function[{m, s, g, a, r}, Exp[-m r^g] * g * m * r^(g - 1) * (a + (1 - m * r^g) * s^2)]];
IPCUFoveaCharacteristic["Curvature"] := Function[{m, s, g, a, r},
    Exp[-m * r^g] * g * m * r^(g - 2) * (a * (-1 + g - g * m * r^g) + (-1 + m * r^g + g * (1 + m * r^g * (-3 + m * r^g))) * s^2)
];
IPCUFoveaCharacteristic["MaximumPoint"] := Function[{m, s, g, a},
    With[{r = (-((-a - s^2) / (m * s^2)))^(1 / g)},
        {r, IPCUFoveaModel[m, s, g, a, r]}
    ]
];
IPCUFoveaCharacteristic["TurningPoint1"] := Function[{m, s, g, a},
    With[{rturn = (((-m) * ((-a) * g + s^2 - 3 * g * s^2) - m * Sqrt[a^2 * g^2 + 2 * a * g *
            s^2 + 2 * a * g^2 * s^2 + s^4 - 2 * g * s^4 + 5 * g^2 * s^4]) / (g * m^2 * s^2))^(1 / g) / 2^g^(-1)},
        {rturn, IPCUFoveaModel[m, s, g, a, rturn]}
    ]
];
IPCUFoveaCharacteristic["TurningPoint2"] := Function[{m, s, g, a},
    With[{rturn = (((-m) * ((-a) * g + s^2 - 3 * g * s^2) + m * Sqrt[a^2 * g^2 + 2 * a * g *
            s^2 + 2 * a * g^2 * s^2 + s^4 - 2 * g * s^4 + 5 * g^2 * s^4]) / (g * m^2 * s^2))^(1 / g) / 2^g^(-1)},
        {rturn, IPCUFoveaModel[m, s, g, a, rturn]}
    ]
];
IPCUFoveaCharacteristic["Steepness"] := Function[{m, s, g, a},
    (-E^((-m) * (((a * g + (-1 + 3 * g) * s^2 -
            Sqrt[a^2 * g^2 + 2 * a * g * (1 + g) * s^2 + (1 + g * (-2 + 5 * g)) * s^4]) / (g * m *
            s^2))^(1 / g) / 2^g^(-1))^g)) * g * m * (((a * g + (-1 + 3 * g) * s^2 -
            Sqrt[a^2 * g^2 + 2 * a * g * (1 + g) * s^2 + (1 + g * (-2 + 5 * g)) * s^4]) / (g * m * s^2))^(1 /
            g) / 2^g^(-1))^(-1 + g) * (-a + s^2 * (-1 + m * (((a * g + (-1 + 3 * g) * s^2 -
            Sqrt[a^2 * g^2 + 2 * a * g * (1 + g) * s^2 + (1 + g * (-2 + 5 * g)) * s^4]) /
            (g * m * s^2))^(1 / g) / 2^g^(-1))^g))
];
IPCUFoveaCharacteristic["Area"] := Function[{m, s, g, a, r},
    (g * r * (a * g - s^2 / E^(m * r^g)) + r * (a * g - s^2) * ExpIntegralE[(-1 + g) / g, m * r^g] +
            (((-a) * g + s^2) * Gamma[1 / g]) / m^g^(-1)) / g^2
];
IPCUFoveaCharacteristic["BowlArea"] := Function[{m, s, g, a},
    (1 / g^2) * ((g * ((a + s^2) / (m * s^2))^(1 / g) * ((-a) * g + s^2 +
            g * m * s^2 * (((a + s^2) / (m * s^2))^(1 / g))^g)) /
            E^(m * (((a + s^2) / (m * s^2))^(1 / g))^g) + ((a + s^2) / (m * s^2))^(1 /
            g) * ((-a) * g + s^2) *
            ExpIntegralE[(-1 + g) / g,
                m * (((a + s^2) / (m * s^2))^(1 / g))^g] + ((a * g - s^2) * Gamma[1 / g]) /
            m^g^(-1))
];

IPCU::nobrack = "Interval boundaries don't have distinct signs. Cannot determine root.";
bisectionRootFind[rootFunc_, xStart_, xEnd_] :=
        Module[
            {
                accuracy = 10.0*^-7,
                x1, x2,
                rtb, i,
                dx, xmid,
                f, fmid},

            x1 = xStart;
            x2 = xEnd;

            f = rootFunc[x1];
            fmid = rootFunc[x2];

            If[f * fmid > 0,
                Message[IPCU::nobrack];
                Return[xStart]
            ];

            rtb = If[f < 0,
                dx = x2 - x1;
                x1,
                dx = x1 - x2;
                x2
            ];

            For[i = 0, i < $IterationLimit, i++,
                xmid = rtb + (dx *= 0.5);
                fmid = rootFunc[xmid];
                If[fmid <= 0.0, rtb = xmid];
                If[Abs[dx] < 10.0*^-7 || fmid == 0.0,
                    Return[rtb]
                ];
            ];
            Message[$IterationLimit::itlim, $IterationLimit];
        ];

Options[IPCUFindFovealRadius] = {
    "Percentage" -> 0.95
};


IPCUFindFovealRadius::nobrack = "Interval boundaries don't have distinct signs. Cannot determine root.";
IPCUFindFovealRadius::perc = "The \"Percentage\" option must be a number 0 < p < 1. Resetting it to 0.95.";
IPCUFindFovealRadius[{m_?NumericQ, s_?NumericQ, g_?NumericQ, a_?NumericQ}, OptionsPattern[]] :=
        Module[
            {
                r,
                rrim = First[IPCUFoveaCharacteristic["MaximumPoint"][m, s, g, a]],
                solution,
                percentage,
                rootFunc
            },
            percentage = OptionValue["Percentage"];
            If[Not@TrueQ[0 < percentage < 1],
                Message[IPCUFindFovealRadius::perc];
                percentage = 0.95;
            ];
            rootFunc =
                    With[
                        {p = percentage},
                        Compile[{{r, _Real, 0}},
                            (1 / g^2) * ((-E^((-m) * r^g)) * g * r * ((-a) * g + s^2 + g * m * r^g * s^2) +
                                    (g * p * ((a + s^2) / (m * s^2))^(1 / g) * ((-a) * g + s^2 +
                                            g * m * s^2 * (((a + s^2) / (m * s^2))^(1 / g))^g)) /
                                            E^(m * (((a + s^2) / (m * s^2))^(1 / g))^g) +
                                    r * (a * g - s^2) * ExpIntegralE[(-1 + g) / g, m * r^g] +
                                    p * ((a + s^2) / (m * s^2))^(1 / g) *
                                            ((-a) * g + s^2) *
                                            ExpIntegralE[(-1 + g) / g, m * (((a + s^2) / (m * s^2))^(1 / g))^g] +
                                    ((-1 + p) * (a * g - s^2) * Gamma[1 / g]) / m^g^(-1))
                        ]
                    ];
            bisectionRootFind[rootFunc, 0.01, rrim]
        ];


(* ::Section:: *)
(*Getting parameter or error values in anatomical directions*)

IPCUGetDirectionalValues::odd = "Transforming parameters into anatomical directions works only if the number of fitted directions is a multiple of 4. The provided dataset contains `1` directions.";

IPCUGetDirectionalValues[p_IPCUPropertiesList, acc___] := Thread[
    {"Nasal", "Superior", "Temporal", "Inferior"} ->
            Transpose@Part[IPCUGetDirectionalValues[#, acc] & /@ List @@ p, All, All, -1]
];
IPCUGetDirectionalValues[p_IPCUProperties, accessor_ : "Parameters"] := Module[
    {
        eye = p["Eye"],
        nDirs = Length[p[accessor]],
        values = p[accessor],
        directionNames = {"Nasal", "Superior", "Temporal", "Inferior"}
    },

    If[Mod[Length[values], 4] =!= 0 ,
        Message[IPCUGetDirectionalValues::odd, Length[values]];
        Abort[]
    ];
    If[eye === Left, values = Reverse[RotateLeft[values, Length[values] / 2 + 1]]];
    With[{dirs = Table[i * Divide[nDirs, 4], {i, 0, 3}] + 1},
        Thread[directionNames -> values[[dirs]]]
    ]
];

Options[IPCUGetFoveaCFST] = {
    "Radius" -> .5,
    "NumberOfRadialDivisions" -> 10
};

(*  This helperfunction calculates the correct radii that are required for an approximate equal distribution of
    sampling points in the circular scheme of the CSFT.
    The final formula comes from solving a simple reccurence equation:

    Block[
      {A0, Ak, Acfst},
      Acfst = Pi*rcfst^2;
      A0 = Acfst/(nd*nr + 1);
      Ak = A0*nd;
      RSolveValue[
       {A0 == Pi*rr[0], Ak == Pi*rr[n + 1] - Pi*rr[n]}, rr[n], n]
    ] // FullSimplify
*)
getRadiiForCFST[nd_, nr_, rcfst_] := Table[Sqrt[((n*nd + 1)/(nd*nr + 1))*rcfst^2], {n, 0, nr}];

IPCUGetFoveaCFST::warg = "Options are set wrong. Using radius = 0.5 and nr = 10 instead";
IPCUGetFoveaCFST[p_IPCUProperties, OptionsPattern[]] := Module[
    {
        radii,
        centerHeight,
        rcfst,
        nd,
        nr
    },
    nr = OptionValue["NumberOfRadialDivisions"];
    rcfst = OptionValue["Radius"];
    If[Not@Element[nr, Integers] || Not@NumericQ[rcfst] || nr < 1 || rcfst <= 0,
        Message[IPCUGetFoveaCFST::warg];
        nr = 10;
        rcfst = 0.5;
    ];
    {centerHeight, nd} = p[{"CentralHeight", {"Parameters", Length}}];
    radii = getRadiiForCFST[nd, nr, rcfst];
    Mean@Append[
        Flatten[Table[centerHeight + IPCUFoveaModel[##, r], {r, radii}] & @@@ p["Parameters"]],
        centerHeight
    ]
];

End[];

EndPackage[];
