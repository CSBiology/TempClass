Temporal classification is an approach conceptually based on a classical time series analysis, where the dependency of the consecutive time points is exploited. 
Constrained smoothing splines with automated model selection separates between signal and noise under the assumption that high frequent changes are less likely 
to occur, simultaneously preserving information about the detected variance. 
This enables a more precise representation of the measured information and improves temporal classification in order to 
identify biologically interpretable correlations among the data. 


Usage:

```fsharp
#r "FSharp.Stats.dll"
#r "TempClass.dll"
#r "nuget: FSharp.Collections.ParallelSeq, 1.2.0"
#r "nuget: Plotly.NET, 4.2.0"

open Plotly.NET
open FSharp.Stats
open FSharp.Collections.ParallelSeq
open TempClass
open TempClass.TemporalClassification
open Plotly.NET
open Plotly.NET.StyleParam
open Plotly.NET.LayoutObjects

FSharp.Stats.ServiceLocator.setEnvironmentPathVariable (@"\..\TempClass\lib")
FSharp.Stats.Algebra.LinearAlgebra.Service()


// time points with spacing according to kinetic expectation
let timepoints = vector [|1.;2.;3.;4.;5.;6.;7.;8.|]

// three replicates where measured at each of the 8 time points
let intensitiesProteinA = 
    [|
        [|17.74781999; 17.60999355; 17.3816851|];
        [|17.44109769; 17.42662059; 17.98721015|];
        [|17.79075992; 17.6181864; 17.66741748|];
        [|17.53004396; 18.35447924; 17.84085591|];
        [|17.90062327; 17.65002708; 17.60924143|];
        [|17.77776007; 17.80117604; 17.55941645|];
        [|17.1401598; 17.73320743; 17.93044716|];
        [|18.43547806; 18.23607406; 17.99477221|]
    |]

// Time point weighting method
let weighting = Fitting.WeightingMethod.StandardDeviation

// Minimization criterion for the final model selection. Shapebased optimization is carried out using mGCV
let minimizer = Fitting.Minimizer.AICc

// smoothing spline result
let (result,modelQualityScores) = Fitting.getBestFit timepoints intensitiesProteinA weighting minimizer

// used smoothing strength
let lambda = result.Lambda //226.44802


// classification decription of the signal. If the intensities do not exceed a range of 0.05, they are classified as constant signals
// alternatively, ANOVA filtering can be applied
let classification = Classification.getClassification timepoints result.TraceA result.TraceC 0.05 1.

// function that takes a x vale and returns the predicted y value of the constrained smoothing spline
let splineFunction : float -> float = result.SplineFunction

let visualizationSpline = 
    TemporalClassification.Vis.getChart result (Some modelQualityScores)

visualizationSpline
|> Chart.show

```


