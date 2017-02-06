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

BeginPackage["FoveaAnalysis`Choroid`", {"FoveaAnalysis`Modeling`", "FoveaAnalysis`Characteristics`"}];

createFoveaProperties::usage = "createFoveaProperties[infoAssociation, octPath] builds info-structure that can directly \
be used with processOCTFile[]";

choroidModelFovea::usage = "choroidModelFovea[resultFile, octPath] models the fovea for choroid analysis.";
choroidResultPlot::usage = "choroidResultPlot[info] shows the final image after modeling with segmentation of Marcus.";

findOCT::usage = "findOCT[info, octPath] tries to locate info[\"VolFile\"].";

(* For testing only *)
extractInfoFromOCT::usage = "extractInfoFromOCT[file, optPath] extracts OCT information like resolution from the file.";
extractMarcusInfo::usage = "extractMarcusInfo[file] converts results from Marcus to make it usable for modeling.";
extractInfo::usage = "extractInfo[marcusResultFile, octPath] collects all information for modeling.";

Begin["`Private`"];

$octInfo = {"ScanFocus", "SizeX", "NumBScans", "SizeZ", "ScaleX", "Distance",
    "ScaleZ", "ScanPosition", "ExamTime", "ScanPattern",
    "FieldSizeSlo", "SizeXSlo", "SizeYSlo", "ScaleXSlo", "ScaleYSlo"};

findOCT::noFile = "Could not locate OCT file ``.";
findOCT[info_Association, octPath_String /; FileExistsQ[octPath]] := Module[{file, files},
    If[FileExistsQ[info["VolFile"]],
        info["VolFile"],
        files = FileNames[FileBaseName[info["VolFile"] <> ".vol"], {octPath}, Infinity];
        If[files === {},
            Message[findOCT::noFile, info["VolFile"]];
            Throw[findOCT],
            First[files]
        ]
    ]
];


extractInfoFromOCT::several = "The following OCT file could not be found: ``";
extractInfoFromOCT[file_, path_] := Module[{files},
    files = FileNames[file, {path}, Infinity];
    If[Length[files] < 1,
        Message[extractInfoFromOCT::several, file];
        Throw[extractInfoFromOCT];
    ];
    {files[[1]], extractInfoFromOCT[files[[1]]]}
];
extractInfoFromOCT[file_?FileExistsQ] := Module[{elements},
    elements = Import[file, {"Heyex", "FileHeader", $octInfo}];
    AssociationThread[$octInfo -> elements]
];

extractMarcusInfo[file_?FileExistsQ] := Module[{expr},
    expr = Get[file];
    Association @@ expr
];

extractInfo::noVol = "No Volume file could be found in results file ``";
extractInfo[file_, octPath_] := Module[{marcusInfo, octInfo, octFile},
    marcusInfo = extractMarcusInfo[file];
    octFile = marcusInfo["VolFile"];
    If[Head[octFile] =!= String,
        Message[extractInfo::noVol, file];
        Throw[extractInfo]
    ];

    If[Catch[
        {octFile, octInfo} = extractInfoFromOCT[marcusInfo["VolFile"], octPath];
    ] != Null,
        Message[extractInfo::noVol, file];
        Throw[extractInfo]
    ];
    Join[marcusInfo, octInfo]
];

createFoveaProperties[info_?AssociationQ, octPath_?FileExistsQ] := Module[{
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
            "ParameterRanges" -> ("ParameterRanges" /. Options[processOCTFile]),
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

choroidModelFovea[resultFile_?FileExistsQ, octPath_?FileExistsQ] := Module[
    {
        info,
        prop,
        result,
        radii},
    Catch[
        info = extractInfo[resultFile, octPath];
        prop = createFoveaProperties[info, octPath];
        result = processOCTFile[prop["VolFile"], prop, "PreferPropertyFile" -> True];
        radii = Reverse[Round[(foveaRadius /@ result["Parameters"]) * {1, -1} / prop["ScaleX"]]];
        result = Association[result, "Radii" -> Last[result["Center"]] + radii];
        Association[result, "AHMeasure" -> choroidFovealThickness[result]]
    ]
];

choroidResultPlot[info_Association, octPath_String /; FileExistsQ[octPath]] := Module[
    {
        file = findOCT[info, octPath],
        toImageCoordinates,
        img, nx, ny, cx, cy
    },
    img = Import[file, {"Heyex", "Images", info["Center"] // First}];
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

End[];
EndPackage[];
