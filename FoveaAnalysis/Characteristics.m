(* Mathematica Package *)
(* Created by Mathematica Plugin for IntelliJ IDEA *)

(* :Title: Characteristics *)
(* :Context: Characteristics` *)
(* :Author: patrick *)
(* :Date: 2017-02-05 *)

(* :Package Version: 0.1 *)
(* :Mathematica Version: *)
(* :Copyright: (c) 2017 patrick *)
(* :Keywords: *)
(* :Discussion: *)

BeginPackage["FoveaAnalysis`Characteristics`"];

foveaCharacteristic::usage = "foveaCharacteristic[charact_String] returns a function of foveal parameters and possibly r to calculate a specific foveal characteristic. The complete list of all available characteristics can be accessed by IPCUFovealCharacteristic[].";
foveaRadius::usage = "foveaRadius[{m, s, g, a}, opt] calculates the radius where a certain percentage of the foveal bowl area is reached. The percentage can be adjusted with the \"Percentage\" option (default is 0.95). The foveal bowl area is the area from r=0 to the rim (which is the maximum of the foveal curve).";
directionalValues::usage = "IPCUGetDirectionalParameters[prop, \"Parameters\" | \"Errors\"] extracts the 4 anatomical directions from a set radially fitted foveas.";
foveaCFST::usage = "foveaCFST[prop, opt] calculated the Central Retinal Thickness as mean of the thicknisses inside the 1mm central radius.";

Begin["`Private`"];

(* ::Section:: *)
(*Calculation of Fovea properties from parameter values*)

foveaCharacteristic[Properties] = {
    "Slope",
    "Curvature",
    "MaximumPoint",
    "TurningPoint1",
    "TurningPoint2",
    "Steepness",
    "Area",
    "BowlArea"
};
foveaCharacteristic[] := foveaCharacteristic[Properties];
foveaCharacteristic["Slope"] := Block[{m, s, g, a, r}, Function[{m, s, g, a, r}, Exp[-m r^g] * g * m * r^(g - 1) * (a + (1 - m * r^g) * s^2)]];
foveaCharacteristic["Curvature"] := Function[{m, s, g, a, r},
    Exp[-m * r^g] * g * m * r^(g - 2) * (a * (-1 + g - g * m * r^g) + (-1 + m * r^g + g * (1 + m * r^g * (-3 + m * r^g))) * s^2)
];
foveaCharacteristic["MaximumPoint"] := Function[{m, s, g, a},
    With[{r = (-((-a - s^2) / (m * s^2)))^(1 / g)},
        {r, foveaModel[m, s, g, a, r]}
    ]
];
foveaCharacteristic["TurningPoint1"] := Function[{m, s, g, a},
    With[{rturn = (((-m) * ((-a) * g + s^2 - 3 * g * s^2) - m * Sqrt[a^2 * g^2 + 2 * a * g *
        s^2 + 2 * a * g^2 * s^2 + s^4 - 2 * g * s^4 + 5 * g^2 * s^4]) / (g * m^2 * s^2))^(1 / g) / 2^g^(-1)},
        {rturn, foveaModel[m, s, g, a, rturn]}
    ]
];
foveaCharacteristic["TurningPoint2"] := Function[{m, s, g, a},
    With[{rturn = (((-m) * ((-a) * g + s^2 - 3 * g * s^2) + m * Sqrt[a^2 * g^2 + 2 * a * g *
        s^2 + 2 * a * g^2 * s^2 + s^4 - 2 * g * s^4 + 5 * g^2 * s^4]) / (g * m^2 * s^2))^(1 / g) / 2^g^(-1)},
        {rturn, foveaModel[m, s, g, a, rturn]}
    ]
];
foveaCharacteristic["Steepness"] := Function[{m, s, g, a},
    (-E^((-m) * (((a * g + (-1 + 3 * g) * s^2 -
        Sqrt[a^2 * g^2 + 2 * a * g * (1 + g) * s^2 + (1 + g * (-2 + 5 * g)) * s^4]) / (g * m *
        s^2))^(1 / g) / 2^g^(-1))^g)) * g * m * (((a * g + (-1 + 3 * g) * s^2 -
        Sqrt[a^2 * g^2 + 2 * a * g * (1 + g) * s^2 + (1 + g * (-2 + 5 * g)) * s^4]) / (g * m * s^2))^(1 /
        g) / 2^g^(-1))^(-1 + g) * (-a + s^2 * (-1 + m * (((a * g + (-1 + 3 * g) * s^2 -
        Sqrt[a^2 * g^2 + 2 * a * g * (1 + g) * s^2 + (1 + g * (-2 + 5 * g)) * s^4]) /
        (g * m * s^2))^(1 / g) / 2^g^(-1))^g))
];
foveaCharacteristic["Area"] := Function[{m, s, g, a, r},
    (g * r * (a * g - s^2 / E^(m * r^g)) + r * (a * g - s^2) * ExpIntegralE[(-1 + g) / g, m * r^g] +
        (((-a) * g + s^2) * Gamma[1 / g]) / m^g^(-1)) / g^2
];
foveaCharacteristic["BowlArea"] := Function[{m, s, g, a},
    (1 / g^2) * ((g * ((a + s^2) / (m * s^2))^(1 / g) * ((-a) * g + s^2 +
        g * m * s^2 * (((a + s^2) / (m * s^2))^(1 / g))^g)) /
        E^(m * (((a + s^2) / (m * s^2))^(1 / g))^g) + ((a + s^2) / (m * s^2))^(1 /
        g) * ((-a) * g + s^2) *
        ExpIntegralE[(-1 + g) / g,
            m * (((a + s^2) / (m * s^2))^(1 / g))^g] + ((a * g - s^2) * Gamma[1 / g]) /
        m^g^(-1))
];

bisectionRootFind::nobrack = "Interval boundaries don't have distinct signs. Cannot determine root.";
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

Options[foveaRadius] = {
    "Percentage" -> 0.95
};


foveaRadius::nobrack = "Interval boundaries don't have distinct signs. Cannot determine root.";
foveaRadius::perc = "The \"Percentage\" option must be a number 0 < p < 1. Resetting it to 0.95.";
foveaRadius[{m_?NumericQ, s_?NumericQ, g_?NumericQ, a_?NumericQ}, OptionsPattern[]] :=
    Module[
        {
            r,
            rrim = First[foveaCharacteristic["MaximumPoint"][m, s, g, a]],
            solution,
            percentage,
            rootFunc
        },
        percentage = OptionValue["Percentage"];
        If[Not@TrueQ[0 < percentage < 1],
            Message[foveaRadius::perc];
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

directionalValues::odd = "Transforming parameters into anatomical directions works only if the number of fitted directions is a multiple of 4. The provided dataset contains `1` directions.";

directionalValues[p_Association, acc___] := Thread[
    {"Nasal", "Superior", "Temporal", "Inferior"} ->
        Transpose@Part[directionalValues[#, acc] & /@ List @@ p, All, All, -1]
];
directionalValues[p_Association, accessor_ : "Parameters"] := Module[
    {
        eye = p["Eye"],
        nDirs = Length[p[accessor]],
        values = p[accessor],
        directionNames = {"Nasal", "Superior", "Temporal", "Inferior"}
    },

    If[Mod[Length[values], 4] =!= 0 ,
        Message[directionalValues::odd, Length[values]];
        Abort[]
    ];
    If[eye === Left, values = Reverse[RotateLeft[values, Length[values] / 2 + 1]]];
    With[{dirs = Table[i * Divide[nDirs, 4], {i, 0, 3}] + 1},
        Thread[directionNames -> values[[dirs]]]
    ]
];

Options[foveaCFST] = {
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

foveaCFST::warg = "Options are set wrong. Using radius = 0.5 and nr = 10 instead";
foveaCFST[p_Association, OptionsPattern[]] := Module[
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
        Message[foveaCFST::warg];
        nr = 10;
        rcfst = 0.5;
    ];
    {centerHeight, nd} = p[{"CentralHeight", {"Parameters", Length}}];
    radii = getRadiiForCFST[nd, nr, rcfst];
    Mean@Append[
        Flatten[Table[centerHeight + foveaModel[##, r], {r, radii}] & @@@ p["Parameters"]],
        centerHeight
    ]
];

End[];

EndPackage[]