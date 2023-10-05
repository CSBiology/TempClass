
#time 
#r @"C:\Users\bvenn\source\repos\TempClass\TempClass\bin\Debug\net472\FSharp.Stats.dll"
#r @"C:\Users\bvenn\source\repos\TempClass\TempClass\bin\Debug\net472\TempClass.dll"
#r "nuget: FSharp.Collections.ParallelSeq, 1.2.0"
#r "nuget: Plotly.NET, 4.2.0"

open Plotly.NET
open FSharp.Stats
open FSharp.Collections.ParallelSeq



FSharp.Stats.ServiceLocator.setEnvironmentPathVariable (@"C:\Users\bvenn\source\repos\TempClass\lib")
FSharp.Stats.Algebra.LinearAlgebra.Service()

let a = Matrix.init 10 4 (fun i j -> float (i+j))
let b = Matrix.init 4 3 (fun i j -> float (i+j))

a * b

Matrix.mul a b


open TempClass
open TemporalClassification

let sigA = (vector [1..5])
let sigB = 
    [1.;1.1;1.2;2.;2.1;1.8;1.4;1.5;1.4;4.;4.;4.1;5.;5.1;5.]
    |> Seq.chunkBySize 3 |> Seq.map Array.ofSeq |> Array.ofSeq

let sigBMeans = sigB |> Seq.map Seq.average |> vector
let sigBStdev = sigB |> Seq.map Seq.stDev |> vector


let weightingMethod = TemporalClassification.Fitting.WeightingMethod.StandardDeviationAdj

let weighting = TemporalClassification.Fitting.getWeighting sigA sigB weightingMethod

Fitting.getWeighting sigA sigB Fitting.WeightingMethod.CV 
Fitting.getWeighting sigA sigB Fitting.WeightingMethod.Equal
Fitting.getWeighting sigA sigB Fitting.WeightingMethod.StandardDeviation 
Fitting.getWeighting sigA sigB Fitting.WeightingMethod.StandardDeviationAdj 
Fitting.getWeighting sigA sigB Fitting.WeightingMethod.StandardDeviationAdjSqrt 
Fitting.getWeighting sigA sigB Fitting.WeightingMethod.StandardDeviationSqrt 
Fitting.getWeighting sigA sigB Fitting.WeightingMethod.Variance
Fitting.getWeighting sigA sigB Fitting.WeightingMethod.VarPlain
Fitting.getWeighting sigA sigB Fitting.WeightingMethod.VarRobust 

let k1 = Fitting.getBestFitOfStd sigA sigBMeans sigBStdev Fitting.WeightingMethod.StandardDeviation Fitting.Minimizer.GCV
let k2 = Fitting.getBestFit sigA sigB Fitting.WeightingMethod.StandardDeviation Fitting.Minimizer.GCV


let data =
    System.IO.File.ReadAllLines(@"C:\Users\bvenn\source\repos\TemporalClassification\assays\Testing\dataset\proteins\imputeKnnNoLow_40_HS.tsv")
    |> Array.tail
    |> Array.map (fun x -> 
        let tmp = x.Split '\t'
        tmp.[1..24] 
        |> Array.map (fun s -> if s = "" then 0. else log (float s)) 
        |> Array.chunkBySize 3
        )
    |> Array.filter (fun x -> x |> Array.concat |> Array.filter (fun t -> t = 0.) |> Array.length = 0)


let vecX = vector [1 .. 8]


let res =
    [|0..10|]
    |> PSeq.withDegreeOfParallelism 16 //16
    |> PSeq.mapi (fun i x -> 
        Fitting.getBestFit vecX data.[0] Fitting.WeightingMethod.StandardDeviation Fitting.Minimizer.AICc)
    |> PSeq.toArray

[

    [
        Chart.Line([1. .. 0.01 .. 8.] |> List.map (fun x -> x,res.[0].SplineFunction x),Name="spline")
        Chart.Point((data.[0] |> Array.mapi (fun i x -> x |> Array.map (fun s -> i+1,s)) |> Seq.concat),Name="raw")
    ]
    |> Chart.combine
    Chart.Point(res.[0].GCVArray.Value,Name="GCV")
]
|> Chart.Grid(2,1)
|> Chart.show

//todo replace gridsearch with bisection?