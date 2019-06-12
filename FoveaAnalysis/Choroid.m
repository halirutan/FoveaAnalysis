(* Mathematica Package *)
(* Created by Mathematica Plugin for IntelliJ IDEA *)

(* :Title: ChoroidAnalysis *)
(* :Context: FoveaAnalysis`ChoroidAnalysis` *)
(* :Author: patrick *)
(* :Date: 2017-01-05 *)

(* :Package Version: 0.1 *)
(* :Mathematica Version: *)
(* :Copyright: (c) 2017 patrick *)
(* :Keywords: *)
(* :Discussion: *)

Package["FoveaAnalysis`"]
PackageImport["HSF`"]

$octInfo = {"ScanFocus", "SizeX", "NumBScans", "SizeZ", "ScaleX", "Distance",
    "ScaleZ", "ScanPosition", "ExamTime", "ScanPattern",
    "FieldSizeSlo", "SizeXSlo", "SizeYSlo", "ScaleXSlo", "ScaleYSlo"};

PackageExport["FindOCT"]
FindOCT::usage = "FindOCT[info, octPath] tries to locate info[\"VolFile\"].";
FindOCT::noFile = "Could not locate OCT file ``.";
FindOCT[info_Association, octPath_String /; FileExistsQ[octPath]] := Module[{file, files},
    If[FileExistsQ[info["VolFile"]],
        info["VolFile"],
        files = FileNames[FileBaseName[info["VolFile"] <> ".vol"], {octPath}, Infinity];
        If[files === {},
            Message[FindOCT::noFile, info["VolFile"]];
            Throw[FindOCT],
            First[files]
        ]
    ]
];


PackageScope["ExtractInfoFromOCT"]
ExtractInfoFromOCT::usage = "ExtractInfoFromOCT[file, optPath] extracts OCT information like resolution from the file.";
ExtractInfoFromOCT::several = "The following OCT file could not be found: ``";
ExtractInfoFromOCT[file_, path_] := Module[{files},
    files = FileNames[file, {path}, Infinity];
    If[Length[files] < 1,
        Message[ExtractInfoFromOCT::several, file];
        Throw[ExtractInfoFromOCT];
    ];
    {files[[1]], ExtractInfoFromOCT[files[[1]]]}
];
ExtractInfoFromOCT[file_?FileExistsQ] := Module[{elements},
    elements = HSFInfo[file][$octInfo];
    AssociationThread[$octInfo -> elements]
];

PackageScope["ExtractMarcusInfo"]
ExtractMarcusInfo::usage = "ExtractMarcusInfo[file] converts results from Marcus to make it usable for modeling.";
ExtractMarcusInfo[file_?FileExistsQ] := Module[{expr},
    expr = Get[file];
    Association @@ expr
];


PackageExport["ExtractInfo"]
ExtractInfo::usage = "ExtractInfo[marcusResultFile, octPath] collects all information for modeling.";
ExtractInfo::noVol = "No Volume file could be found in results file ``";
ExtractInfo[file_, octPath_] := Module[{marcusInfo, octInfo, octFile},
    marcusInfo = ExtractMarcusInfo[file];
    octFile = marcusInfo["VolFile"];
    If[Head[octFile] =!= String,
        Message[ExtractInfo::noVol, file];
        Throw[ExtractInfo]
    ];

    If[Catch[
        {octFile, octInfo} = ExtractInfoFromOCT[marcusInfo["VolFile"], octPath];
    ] != Null,
        Message[ExtractInfo::noVol, file];
        Throw[ExtractInfo]
    ];
    Join[marcusInfo, octInfo]
];

PackageExport["CreateFoveaProperties"]
CreateFoveaProperties::usage = "CreateFoveaProperties[infoAssociation, octPath] builds info-structure that can directly be used with processOCTFile[]";
CreateFoveaProperties[info_?AssociationQ, octPath_?FileExistsQ] := Module[{
    center = Reverse@info["Center"],
    file = AbsoluteFileName@First@FileNames[info["VolFile"], {octPath}, Infinity],
    centralPixelHeight,
    rescale = False,
    segData},

    segData = {"RPE", "ILM"} /. Import[file, {"Heyex", "SegmentationData", First[center]}];
    centralPixelHeight = Subtract @@ segData[[All, Last[center]]];
    Join[info,
        <|
            "Center" -> center,
            "CentralPixelHeight" -> centralPixelHeight,
            "VolFile" -> file,
            "RescaleOCTMagnification" -> False,
            "AngleStepSize" -> Pi,
            "ParameterRanges" -> ("ParameterRanges" /. Options[FindFoveaModelParameters]),
            "MaxRadius" -> 2
        |>
    ]
];

choroidFovealThickness[info_Association] := Module[
    {
        ip, ah, start, end, ahFov
    },

    ah = info["RPE"] - info["AH"];
    {start, end} = info["Radii"];
    ip = ListInterpolation[ah, Method -> "Spline"];
    ahFov = ah[[start ;; end]];
    <|
        "OCBIntegrated" ->
                NIntegrate[ip[x], {x, start, end}] / Abs[end - start],
        "OCBMean" -> Mean[ahFov],
        "OCBTotal" -> Total[ahFov],
        "OCBMedian" -> Median[ahFov],
        "OCBStDev" -> StandardDeviation[ahFov]
    |>
];

PackageExport["ChoroidModelFovea"]
ChoroidModelFovea::usage = "ChoroidModelFovea[resultFile, octPath] models the fovea for choroid analysis.";
ChoroidModelFovea[resultFile_?FileExistsQ, octPath_?FileExistsQ] := Module[
    {
        info,
        prop,
        result,
        radii},
    Catch[
        info = ExtractInfo[resultFile, octPath];
        prop = CreateFoveaProperties[info, octPath];
        result = FindFoveaModelParameters[prop["VolFile"], prop, "PreferPropertyFile" -> True];
        radii = Reverse[Round[(FoveaRadius /@ result["Parameters"]) * {1, -1} / prop["ScaleX"]]];
        result = Association[result, "Radii" -> Last[result["Center"]] + radii];
        Association[result, "AHMeasure" -> choroidFovealThickness[result]]
    ]
];

PackageExport["ChoroidResultPlot"]
ChoroidResultPlot::usage = "ChoroidResultPlot[info] shows the final image after modeling with segmentation of Marcus.";
ChoroidResultPlot[info_Association, octPath_String /; FileExistsQ[octPath]] := Module[
    {
        file = FindOCT[info, octPath],
        toImageCoordinates,
        img, nx, ny, cx, cy
    },
    img = HSFBScanImage[file, info["Center"] // First];
    {nx, ny} = ImageDimensions[img];
    {cy, cx} = info["Center"];
    toImageCoordinates = Function[layer,
        Line[Transpose[{Range[0, Length[layer] - 1], 496 + layer / (info["ScaleZ"] * 1000)}]]
    ];


    HighlightImage[
        HighlightImage[img,
            {
                toImageCoordinates /@ Lookup[info, {"ILM", "RPE", "AH"}],
                {Yellow,
                    Line[{{cx, 0}, {cx, ny}}],
                    Line[{{#, 0}, {#, ny}}]& /@ info["Radii"]
                }
            }
        ],
        Graphics[{Opacity[1], White, Text[Style[Column[Normal@info["AHMeasure"]],12], {nx / 5, ny - 50}]}]
    ]
];
