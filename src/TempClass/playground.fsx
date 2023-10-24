#time "on"
 
#r @"C:\Users\bvenn\source\repos\TempClass\TempClass\bin\Debug\net472\FSharp.Stats.dll"
#r @"C:\Users\bvenn\source\repos\TempClass\TempClass\bin\Debug\net472\TempClass.dll"
#r "nuget: FSharp.Collections.ParallelSeq, 1.2.0"
#r "nuget: Plotly.NET, 4.2.0"

open Plotly.NET
open FSharp.Stats
open FSharp.Collections.ParallelSeq



FSharp.Stats.ServiceLocator.setEnvironmentPathVariable (@"C:\Users\bvenn\source\repos\TempClass\lib")
FSharp.Stats.Algebra.LinearAlgebra.Service()

open System
open System.IO
open FSharp.Stats
open TempClass
open TempClass.TemporalClassification
open Plotly.NET
open Plotly.NET.StyleParam
open Plotly.NET.LayoutObjects


1+1





//some axis styling
module Chart = 
    let myAxis name = 
        LinearAxis.init(
            Title=Title.init name,
                Mirror=StyleParam.Mirror.All,
                Ticks=StyleParam.TickOptions.Inside,
                ShowGrid=false,
                ShowLine=true,
                TickFont=Font.init(Size=12,Family=StyleParam.FontFamily.Arial)
                //Sets tickFont
        )
    let myAxisRange name (min,max) = LinearAxis.init(Title=Title.init name,Range=Range.MinMax(min,max),Mirror=StyleParam.Mirror.All,Ticks=StyleParam.TickOptions.Inside,ShowGrid=false,ShowLine=true)
    let svgConfig =
        Config.init (
            ToImageButtonOptions = ConfigObjects.ToImageButtonOptions.init(
                Format = StyleParam.ImageFormat.SVG
            )
        )
    let withAxisTitles x y chart = 
        chart 
        |> Chart.withTemplate ChartTemplates.lightMirrored
        |> Chart.withXAxis (myAxis x) 
        |> Chart.withYAxis (myAxis y)
        |> Chart.withConfig svgConfig

    
    let withAxisTitlesXR x y (mi,ma) chart = 
        chart 
        |> Chart.withTemplate ChartTemplates.lightMirrored
        |> Chart.withXAxis (myAxisRange x (mi,ma) ) 
        |> Chart.withYAxis (myAxis y)
        |> Chart.withConfig svgConfig

    
    let withAxisTitlesYR x y (mi,ma) chart = 
        chart 
        |> Chart.withTemplate ChartTemplates.lightMirrored
        |> Chart.withXAxis (myAxis x) 
        |> Chart.withYAxis (myAxisRange y (mi,ma))
        |> Chart.withConfig svgConfig


//FSharp.Stats.ServiceLocator.setEnvironmentPathVariable (@"C:\Users\bvenn\OneDrive - Computational Systems Biology\Projects\EIT\revise0622\FSharp.StatsMSF030\lib")
//FSharp.Stats.ServiceLocator.setEnvironmentPathVariable (@"C:\Users\venn\source\repos\FSharp.StatsMSF030\bin\FSharp.Stats.MSF\lib")
//FSharp.Stats.ServiceLocator.setEnvironmentPathVariable (@"C:\Users\venn\source\repos\FSharp.StatsMSF030\lib")
//@"C:\Users\bvenn\source\repos\FSharp.Stats-fa3dce19ff859c59e8fd35640d7918fb9899d555\lib"
//FSharp.Stats.Algebra.LinearAlgebra.Service()

let a = Matrix.init 10 4 (fun i j -> float (i+j))
let b = Matrix.init 4 3 (fun i j -> float (i+j))

a * b

Matrix.mul a b

let xs = [1. .. 5.] |> vector
let ys = 
    [
        1.;1.2;1.3;
        1.5;2.5;0.6;
        0.2;0.25;0.18
        1.7;1.5;1.85;
        2.5;2.9;2.3
    ]
    |> vector

let ys_ = ys |> Seq.chunkBySize 3 |> Seq.map Seq.mean |> vector



let data = 
    System.IO.File.ReadAllLines @"C:\Users\bvenn\source\repos\TempClass\data\imputeKnnNoLow_40_HS.tsv"
    |> Array.tail
    |> Array.map (fun x -> 
        let tmp = x.Split '\t'
        let id = tmp.[0]
        let signal = tmp.[1..24] |> Array.map (fun x -> if x = "" then 0. else log (float x + 1.))
        let trivial = tmp.[25]
        let mm = tmp.[26]
        let mmdes = tmp.[27]
        let go = tmp.[28]
        let godes = tmp.[29]
        let loc = tmp.[30]
        id,signal
        )
    |> Array.filter (fun (id,signal) -> 
        let zeros = signal |> Array.filter (fun x -> x = 0.) |> Array.length
        let length = signal.Length
        let zeroamount = float zeros / float length < 0.20 
        let anova = 
            signal |> Seq.chunkBySize 3 |> Testing.Anova.oneWayAnova
            |> fun t -> t.Factor.Significance <= 0.05
        let emptyTP = signal |> Seq.chunkBySize 3 |> Seq.map Seq.mean |> Seq.tryFind (fun k -> k = 0.) |> fun t -> t.IsNone
        zeroamount && emptyTP //&& anova
        )

type Characterization = {
    ID : string
    Signal : vector 
    BestFit : Fitting.TempClassResult
    Chart :      GenericChart.GenericChart
    FitQuality: (string*float)[]
    }
    with
        static member Create id s b c fq = {ID=id; Signal = s; BestFit=b; Chart=c; FitQuality=fq}

let timepoints = [1. .. 8.] |> vector

let plot i sens = 
    let result,modelQualityScores = Fitting.getBestFit timepoints (snd data.[i] |> Array.chunkBySize 3) Fitting.WeightingMethod.StandardDeviation Fitting.Minimizer.AICc
    
    let classification = Classification.getClassification result.XValues.Value result.TraceA result.TraceC 0.05 sens
    TemporalClassification.Vis.getChart result (Some modelQualityScores)
    |> Chart.withTitle classification
    |> Chart.show



let proc (id,signal:float[]) weighting minimizer =
    let ys = signal
    let fitB,models = 
        //Fitting.getBestFit timepoints ys 3 Fitting.WeightingMethod.CV Fitting.Minimizer.GCV
        //Fitting.getBestFit timepoints ys 3 Fitting.WeightingMethod.Variance Fitting.Minimizer.GCV
        //Fitting.getBestFit timepoints ys 3 weighting minimizer //Fitting.Minimizer.GCV
        Fitting.getBestFit timepoints (ys |> Array.chunkBySize 3) weighting minimizer //Fitting.Minimizer.GCV
    
    let scatter = 
        ys 
        |> Seq.chunkBySize 3 
        |> Seq.mapi (fun i x -> x |> Array.map (fun xi -> timepoints.[i],xi)) 
        |> Seq.concat
        |> Chart.Point
        |> Chart.withTraceInfo "raw points"
    
    let linSpline = 
        let xs = timepoints |> Array.ofSeq
        let ys = ys |> Seq.chunkBySize 3 |> Seq.map Seq.mean |> Array.ofSeq
        let coefs = Interpolation.LinearSpline.initInterpolate xs ys
        
        [1. .. 0.01 .. 8.]
        |> List.map (fun x -> x,Interpolation.LinearSpline.interpolate coefs x)
        |> Chart.Line
        |> Chart.withTraceInfo "linear spline"

    let polynomial = 
        let xs = timepoints
        let ys = ys |> Seq.chunkBySize 3 |> Seq.map Seq.mean  |> vector
        let coefs = Interpolation.Polynomial.coefficients xs ys
        
        [1. .. 0.01 .. 8.]
        |> List.map (fun x -> x,Interpolation.Polynomial.fit coefs x)
        |> Chart.Line
        |> Chart.withTraceInfo "polynomial interpolation"

    let cubicSpline = 
        let xs = timepoints 
        let ys = ys |> Seq.chunkBySize 3 |> Seq.map Seq.mean |> vector
        let coefs = Interpolation.CubicSpline.Simple.coefficients Interpolation.CubicSpline.Simple.BoundaryCondition.Natural xs ys
        
        [1. .. 0.01 .. 8.]
        |> List.map (fun x -> x,Interpolation.CubicSpline.Simple.fitWithLinearPrediction coefs xs x)
        |> Chart.Line
        |> Chart.withTraceInfo "natural cubic spline"
    
    let classy = Classification.getClassification timepoints fitB.TraceA fitB.TraceC 0.05 1.

    let spline (fit:Fitting.TempClassResult) = 
        [1. .. 0.01 .. 8.]
        |> List.map (fun x -> x,fit.SplineFunction x)
        |> Chart.Line
        |> Chart.withTraceInfo "temporal classification" //(sprintf "A_%.3f  G_%.3f" fit.AICc fit.GCV)

    let description = 
        let m = 
            match minimizer with
            | Fitting.Minimizer.GCV -> "GCV"
            | Fitting.Minimizer.AICc -> "AIC"
        let w = 
            match weighting with
            | Fitting.WeightingMethod.StandardDeviation -> "StandardDeviation"
            | Fitting.WeightingMethod.StandardDeviationSqrt -> "StandardDeviationSqrt"
            | Fitting.WeightingMethod.CV -> "CV"
            | Fitting.WeightingMethod.Equal -> "Equal"
            | Fitting.WeightingMethod.StandardDeviationAdj -> "StandardDeviationAdj"
            | Fitting.WeightingMethod.StandardDeviationAdjSqrt -> "StandardDeviationAdjSqrt"
            | Fitting.WeightingMethod.Variance -> "Variance"
            | Fitting.WeightingMethod.VarPlain -> "VarPlain"
            | Fitting.WeightingMethod.VarRobust -> "VarRobust"

        sprintf "Minimizer: %s <br>Weighting: %s <br> Classification: %s" m w classy

    let combined =
        [
        scatter
        linSpline
        polynomial
        cubicSpline
        spline fitB
        ]
        |> Chart.combine
        |> Chart.withSize (900.,600.)
        |> Chart.withTemplate ChartTemplates.lightMirrored
        |> Chart.withXAxisStyle (Title=Title.init(Text="time point",Font=Font.init(FontFamily.Arial,16)),ShowGrid=false)
        |> Chart.withYAxisStyle (Title=Title.init(Text="abundance",Font=Font.init(FontFamily.Arial,16)),ShowGrid=false)
        |> Chart.withDescription [Giraffe.ViewEngine.HtmlElements.rawText description]
        |> Chart.withConfig Chart.svgConfig
        |> Chart.withLayoutStyle(Font=(Font.init(Family=FontFamily.Arial,Size=16)))

    Characterization.Create id (vector ys) fitB combined models

(proc data.[0] Fitting.WeightingMethod.StandardDeviation Fitting.Minimizer.GCV).Chart |> Chart.show

    
let resVarAic = 
    data.[0..15]
    |> PSeq.withDegreeOfParallelism 12
    |> PSeq.mapi (fun i x ->
        printfn "%i" i
        proc x Fitting.WeightingMethod.Variance Fitting.Minimizer.AICc
        )
    |> PSeq.toArray

let resVarGcv = 
    data
    |> PSeq.withDegreeOfParallelism 12
    |> PSeq.mapi (fun i x ->
        printfn "%i" i
        proc x Fitting.WeightingMethod.Variance Fitting.Minimizer.GCV
        )
    |> PSeq.toArray

let resVarPlainAic = 
    data
    |> PSeq.withDegreeOfParallelism 12
    |> PSeq.mapi (fun i x ->
        printfn "%i" i
        proc x Fitting.WeightingMethod.VarPlain Fitting.Minimizer.AICc
        )
    |> PSeq.toArray

let resVarPlainGcv = 
    data
    |> PSeq.withDegreeOfParallelism 12
    |> PSeq.mapi (fun i x ->
        printfn "%i" i
        proc x Fitting.WeightingMethod.VarPlain Fitting.Minimizer.GCV
        )
    |> PSeq.toArray

let resCVAic = 
    data
    |> PSeq.withDegreeOfParallelism 12
    |> PSeq.mapi (fun i x ->
        printfn "%i" i
        proc x Fitting.WeightingMethod.CV Fitting.Minimizer.AICc
        )
    |> PSeq.toArray

let resCVGcv = 
    data
    |> PSeq.withDegreeOfParallelism 12
    |> PSeq.mapi (fun i x ->
        printfn "%i" i
        proc x Fitting.WeightingMethod.CV Fitting.Minimizer.GCV
        )
    |> PSeq.toArray

let resstdevadjAic = 
    data
    |> PSeq.withDegreeOfParallelism 12
    |> PSeq.mapi (fun i x ->
        printfn "%i" i
        proc x Fitting.WeightingMethod.StandardDeviationAdj Fitting.Minimizer.AICc
        )
    |> PSeq.toArray

let resstdevadjGcv = 
    data
    |> PSeq.withDegreeOfParallelism 12
    |> PSeq.mapi (fun i x ->
        printfn "%i" i
        proc x Fitting.WeightingMethod.StandardDeviationAdj Fitting.Minimizer.GCV
        )
    |> PSeq.toArray

let resstdevAic = 
    data
    |> PSeq.withDegreeOfParallelism 8
    |> PSeq.mapi (fun i x ->
        if i%50=0 then printfn "%i" i
        proc x Fitting.WeightingMethod.StandardDeviation Fitting.Minimizer.AICc
        )
    |> PSeq.toArray

let resstdevGcv = 
    data
    |> PSeq.withDegreeOfParallelism 12
    |> PSeq.mapi (fun i x ->
        printfn "%i" i
        proc x Fitting.WeightingMethod.StandardDeviation Fitting.Minimizer.GCV
        )
    |> PSeq.toArray

data.[0..9]
|> PSeq.mapi (fun i x ->
    [
    (proc x Fitting.WeightingMethod.StandardDeviation Fitting.Minimizer.AICc).Chart
    (proc x Fitting.WeightingMethod.CV Fitting.Minimizer.AICc).Chart
    ]
    |> Chart.combine
    |> Chart.withTemplate ChartTemplates.lightMirrored
    |> Chart.withSize (1300.,900.)
    |> Chart.withTitle (fst x)
    |> Chart.show
    )
|> PSeq.toArray


let createHyperlink numberOfDirectoriesToStepUp relativePath nameToShow =
    let stepUp = string numberOfDirectoriesToStepUp
    let relativePath = sprintf """ "%s" """ relativePath
    let nameToShow = sprintf """ "%s" """ nameToShow
    """=HYPERLINK(MID(CELL("filename", $B$2),1,FIND("*",SUBSTITUTE(CELL("filename", $B$2),"\","*",LEN(CELL("filename", $B$2))-LEN(SUBSTITUTE(CELL("filename", $B$2),"\",""))-""" + stepUp + """))) & """ + relativePath + "," + nameToShow + ")"
   
let toString (c:Characterization) =
    let signal = c.Signal |> Seq.map string |> String.concat "\t"
    let traceA = c.BestFit.TraceA |> Seq.map string |> String.concat "\t"
    let traceC = c.BestFit.TraceC |> Seq.map string |> String.concat "\t"
    let classification = Classification.getClassification timepoints c.BestFit.TraceA c.BestFit.TraceC 0.05 1.
    let hyper = createHyperlink 0 (sprintf "plots/%s.html" c.ID.[0..30]) c.ID.[0..30]
    let id = c.ID.[0..30]
    let models = c.FitQuality |> Array.map (fun (m,qScore) -> string qScore) |> String.concat "\t"
    sprintf "%s\t%s\t%s\t%s\t%s\t%s\t%s" hyper id signal traceA traceC classification models
    
let writeResults tc path = 
    tc
    |> Array.map (fun x -> 
        [
        x.Chart
        Chart.Column(x.FitQuality,MarkerColor=Color.fromString "grey",Name="qc") |> Chart.withAxisTitles "" "quality score"
        ]
        |> Chart.Grid(2,1)
        //|> Chart.withSize (1300.,900.)
        //|> Chart.withTemplate ChartTemplates.lightMirrored
        |> Chart.saveHtml (path + @"plots\" + x.ID.[0..30] + ".html")
        toString x
        )
    |> fun data -> 
        let modelScores = "unconstrained\tIn0\tDe0\tIn1\tDe1\tIn2\tDe2\tIn3\tDe3\tIn4\tDe4"
        let header = [|"idlink\tid\tt1_1\tt1_2\tt1_3\tt2_1\tt2_2\tt2_3\tt3_1\tt3_2\tt3_3\tt4_1\tt4_2\tt4_3\tt5_1\tt5_2\tt5_3\tt6_1\tt6_2\tt6_3\tt7_1\tt7_2\tt7_3\tt8_1\tt8_2\tt8_3\ttA_1\ttA_2\ttA_3\ttA_4\ttA_5\ttA_6\ttA_7\ttA_8\ttC_1\ttC_2\ttC_3\ttC_4\ttC_5\ttC_6\ttC_7\ttC_8\tClassification\t" + modelScores|]
        System.IO.File.WriteAllLines(path + "result.tsv",Array.append header data)

//writeResults resCVAic       @"C:\Users\venn\source\repos\TemporalClassification\runs\Testing\HS40\AIC\CV\"
//writeResults resVarAic      @"C:\Users\venn\source\repos\TemporalClassification\runs\Testing\HS40\AIC\var\"
//writeResults resVarPlainAic @"C:\Users\venn\source\repos\TemporalClassification\runs\Testing\HS40\AIC\varplain\"
//writeResults resstdevadjAic @"C:\Users\venn\source\repos\TemporalClassification\runs\Testing\HS40\AIC\stdadj\"
writeResults resstdevAic    @"C:\Users\bvenn\source\repos\TemporalClassification\runs\Testing\HS40\AIC\std\"
printfn "All AIC Combinations Ready!!"

//writeResults resCVGcv       @"C:\Users\venn\source\repos\TemporalClassification\runs\Testing\HS40\GCV\CV\"
//writeResults resVarGcv      @"C:\Users\venn\source\repos\TemporalClassification\runs\Testing\HS40\AIC\var\"
//writeResults resVarPlainGcv @"C:\Users\venn\source\repos\TemporalClassification\runs\Testing\HS40\AIC\varplain\"
//writeResults resstdevadjGcv @"C:\Users\venn\source\repos\TemporalClassification\runs\Testing\HS40\GCV\stdadj\"
//writeResults resstdevGcv    @"C:\Users\venn\source\repos\TemporalClassification\runs\Testing\HS40\GCV\std\"
//printfn "All GCV Combinations Ready!!"




let showClassOccupance (input: Characterization []) aggregateSimple title descr savePath =
    input
    |> Array.map (fun c -> 
        Classification.getClassification timepoints c.BestFit.TraceA c.BestFit.TraceC 0.05 1.
        )
    |> Array.groupBy (fun x -> 
        if aggregateSimple then 
            if x.StartsWith("D") then "D" elif x.StartsWith "I" then "I" else x
        else x
        )
    |> Array.sortByDescending (snd >> Seq.length)
    |> Array.map (fun (x,i) -> x.Replace(".00",""),i.Length)
    |> Chart.Column
    |> Chart.withSize (900.,600.)
    |> Chart.withTemplate ChartTemplates.lightMirrored
    |> Chart.withXAxisStyle (Title=Title.init(Text="classes",Font=Font.init(FontFamily.Arial,16)),ShowGrid=false)
    |> Chart.withMargin(Margin.init(Bottom= 300.))
    |> Chart.withYAxisStyle (Title=Title.init(Text="count",Font=Font.init(FontFamily.Arial,16)),ShowGrid=false)
    |> Chart.withDescription [Giraffe.ViewEngine.HtmlElements.rawText descr]
    |> Chart.withConfig Chart.svgConfig
    |> Chart.withTitle (Title.init title)
    |> Chart.withLayoutStyle(Font=(Font.init(Family=FontFamily.Arial,Size=16))) 
    |> Chart.saveHtml savePath 
    
//showClassOccupance resVarGcv        false "resVarGcv"       "Minimizer: GCV<br>Weighting: var"      @"C:\Users\venn\source\repos\TemporalClassification\runs\Testing\HS40\GCV\var\ClassHistogram.html"
//showClassOccupance resVarPlainGcv   false "resVarPlainGcv"  "Minimizer: GCV<br>Weighting: varPlain" @"C:\Users\venn\source\repos\TemporalClassification\runs\Testing\HS40\GCV\varPlain\ClassHistogram.html"
//showClassOccupance resCVGcv         false "resCVGcv"        "Minimizer: GCV<br>Weighting: CV"       @"C:\Users\venn\source\repos\TemporalClassification\runs\Testing\HS40\GCV\CV\ClassHistogram.html"
//showClassOccupance resstdevadjGcv   false "resstdadjGcv"    "Minimizer: GCV<br>Weighting: stdadj"   @"C:\Users\venn\source\repos\TemporalClassification\runs\Testing\HS40\GCV\stdadj\ClassHistogram.html"
//
//showClassOccupance resVarAic        false "resVarAic"       "Minimizer: AICc<br>Weighting: var"      @"C:\Users\venn\source\repos\TemporalClassification\runs\Testing\HS40\AIC\var\ClassHistogram.html"
//showClassOccupance resVarPlainAic   false "resVarPlainAic"  "Minimizer: AICc<br>Weighting: varPlain" @"C:\Users\venn\source\repos\TemporalClassification\runs\Testing\HS40\AIC\varPlain\ClassHistogram.html"
//showClassOccupance resCVAic         false "resCVAic"        "Minimizer: AICc<br>Weighting: CV"       @"C:\Users\venn\source\repos\TemporalClassification\runs\Testing\HS40\AIC\CV\ClassHistogram.html"
//showClassOccupance resstdevadjAic   false "resstdadjAic"    "Minimizer: AICc<br>Weighting: stdadj"   @"C:\Users\venn\source\repos\TemporalClassification\runs\Testing\HS40\AIC\stdadj\ClassHistogram.html"


showClassOccupance resstdevAic   false "resstdevAic"    "Minimizer: AICc<br>Weighting: std"   @"C:\Users\bvenn\source\repos\TemporalClassification\runs\Testing\HS40\AIC\std\ClassHistogram.html"


let aicstd =  
    System.IO.File.ReadAllLines(@"C:\Users\bvenn\source\repos\TemporalClassification\runs\Testing\HS40\AIC\std\result.tsv")
    |> Array.tail
    |> Array.map (fun x -> 
        let tmp = x.Split '\t'
        tmp.[42]
        )

aicstd
|> Array.groupBy (fun x -> 
    if x.StartsWith("D") then "D" elif x.StartsWith "I" then "I" else x
    )
|> Array.sortByDescending (snd >> Seq.length)
|> Array.map (fun (x,i) -> x.Replace(".00",""),i.Length)
|> Chart.Column
|> Chart.withSize (750.,600.)
|> Chart.withTemplate ChartTemplates.lightMirrored
|> Chart.withXAxisStyle (Title=Title.init(Text="classes",Font=Font.init(FontFamily.Arial,16)),ShowGrid=false)
|> Chart.withMargin(Margin.init(Bottom= 300.))
|> Chart.withYAxisStyle (Title=Title.init(Text="count",Font=Font.init(FontFamily.Arial,16)),ShowGrid=false)
|> Chart.withConfig Chart.svgConfig
|> Chart.withTitle (Title.init "class occupancy")
|> Chart.withLayoutStyle(Font=(Font.init(Family=FontFamily.Arial,Size=16))) 
|> Chart.show

aicstd
|> Array.filter (fun x -> x.StartsWith "D" || x.StartsWith "I")
|> Array.groupBy (fun x -> x)
|> Array.sortByDescending (snd >> Seq.length)
|> Array.map (fun (x,i) -> x.Replace(".00",""),i.Length)
|> Chart.Column
|> Chart.withSize (750.,600.)
|> Chart.withTemplate ChartTemplates.lightMirrored
|> Chart.withXAxisStyle (Title=Title.init(Text="classes",Font=Font.init(FontFamily.Arial,16)),ShowGrid=false)
|> Chart.withMargin(Margin.init(Bottom= 300.))
|> Chart.withYAxisStyle (Title=Title.init(Text="count",Font=Font.init(FontFamily.Arial,16)),ShowGrid=false)
|> Chart.withConfig Chart.svgConfig
|> Chart.withTitle (Title.init "class occupancy")
|> Chart.withLayoutStyle(Font=(Font.init(Family=FontFamily.Arial,Size=16))) 
|> Chart.show

let extremaClassCount1 = 
    aicstd
    |> Array.filter (fun x -> not (x.StartsWith "D" || x.StartsWith "I" || x.StartsWith "constant"))
    |> Array.filter (fun x -> (x.Split ",") |> Array.length = 1)
    |> Array.distinct
    |> Array.length

let extremaClassCount2 = 
    aicstd
    |> Array.filter (fun x -> not (x.StartsWith "D" || x.StartsWith "I"))
    |> Array.filter (fun x -> (x.Split ",") |> Array.length = 2)
    |> Array.distinct
    |> Array.length

let extremaClassCount3 = 
    aicstd
    |> Array.filter (fun x -> not (x.StartsWith "D" || x.StartsWith "I" || x.StartsWith "constant"))
    |> Array.filter (fun x -> (x.Split ",") |> Array.length = 3)
    |> Array.distinct
    |> Array.length

let extremaClassCount4 = 
    aicstd
    |> Array.filter (fun x -> not (x.StartsWith "D" || x.StartsWith "I" || x.StartsWith "constant"))
    |> Array.filter (fun x -> (x.Split ",") |> Array.length = 4)
    |> Array.distinct
    |> Array.length

[
    Chart.StackedColumn(
        [
        "0 extrema",3
        "1 extrema",14//extremaClassCount1
        "2 extrema",38//extremaClassCount2
        "3 extrema",47//extremaClassCount3
        "4 extrema",1//extremaClassCount4
        ],Name="occupied classes",MarkerColor = Color.fromString "orange",MultiText = ["3/3";"14/16";"38/56";"47/112";"1/140"])
    Chart.StackedColumn(
        [
        "0 extrema",0
        "1 extrema",16-14//extremaClassCount1
        "2 extrema",56-38//extremaClassCount2
        "3 extrema",112-47//extremaClassCount3
        "4 extrema",140-1//extremaClassCount4
        ],Name="unoccupied classes",MarkerColor = Color.fromString "#1f77b4")
]
|> Chart.combine
|> Chart.withAxisTitles "" "shape classes"
|> Chart.withLegendStyle(Orientation=Orientation.Horizontal)
|> Chart.withSize (750.,400.)
|> Chart.withLayoutStyle(Font=(Font.init(Family=FontFamily.Arial,Size=16))) 
|> Chart.show

// make word file
// import histograms
//get example charts for some classes

(*
The 3 replicates are separated and fitted individually!!
  - Temporal classification with GCV and equal weighting (bacause n=1)
  - polynomial interpolation
  - linearspline
At 71 positions in equally spaced intervals, the variance of the three replicate fits are determined. These variances are pooled and visualised as histograms.
Expectation:
- TC (something in between)
- lin (lowest variance because all outlier even out to one "average" smoothed response)
- poly ()
- linSpline (highest variance)
*)

//let procSingle (id,signal:float[]) =
//    let ys1 = 
//        signal
//        |> Seq.chunkBySize 3 
//        |> Seq.map (fun x -> x.[0])
//        |> vector 
//    let ys2 = 
//        signal
//        |> Seq.chunkBySize 3 
//        |> Seq.map (fun x -> x.[1])
//        |> vector 
//    let ys3 = 
//        signal
//        |> Seq.chunkBySize 3 
//        |> Seq.map (fun x -> x.[2])
//        |> vector 

//    let fitB1 = Fitting.getBestFit timepoints ys1 1 Fitting.WeightingMethod.Equal Fitting.Minimizer.GCV
//    let fitB2 = Fitting.getBestFit timepoints ys2 1 Fitting.WeightingMethod.Equal Fitting.Minimizer.GCV
//    let fitB3 = Fitting.getBestFit timepoints ys3 1 Fitting.WeightingMethod.Equal Fitting.Minimizer.GCV
    
//    let poly1 = Interpolation.Polynomial.fit (Interpolation.Polynomial.coefficients timepoints ys1)
//    let poly2 = Interpolation.Polynomial.fit (Interpolation.Polynomial.coefficients timepoints ys2)
//    let poly3 = Interpolation.Polynomial.fit (Interpolation.Polynomial.coefficients timepoints ys3)
    
//    let lin1 = Interpolation.LinearSpline.initInterpolate (Array.ofSeq timepoints) (Array.ofSeq ys1) |> Interpolation.LinearSpline.interpolate
//    let lin2 = Interpolation.LinearSpline.initInterpolate (Array.ofSeq timepoints) (Array.ofSeq ys2) |> Interpolation.LinearSpline.interpolate
//    let lin3 = Interpolation.LinearSpline.initInterpolate (Array.ofSeq timepoints) (Array.ofSeq ys3) |> Interpolation.LinearSpline.interpolate

//    let scatter = 
//        signal 
//        |> Seq.chunkBySize 3 
//        |> Seq.mapi (fun i x -> x |> Array.map (fun xi -> timepoints.[i],xi)) 
//        |> Seq.concat
//        |> Chart.Point
    

//    let spline (fit:Fitting.TempClassResult) color = 
//        [1. .. 0.01 .. 8.]
//        |> List.map (fun x -> x,fit.SplineFunction x)
//        |> fun x -> Chart.Line(x,LineColor=Color.fromHex color)
//        |> Chart.withTraceInfo (sprintf "A_%.3f  G_%.3f" fit.AICc fit.GCV)

//    let polyFit (fit:float -> float) color = 
//        [1. .. 0.01 .. 8.]
//        |> List.map (fun x -> x,fit x)
//        |> fun x -> Chart.Line(x,LineColor=Color.fromHex color)

//    let linear (fit:float -> float) color = 
//        [1. .. 0.01 .. 8.]
//        |> List.map (fun x -> x,fit x)
//        |> fun x -> Chart.Line(x,LineColor=Color.fromHex color)

//    let combined =
//        [
//        scatter
//        spline fitB1 "#fc7e0f"
//        spline fitB2 "#fc7e0f"
//        spline fitB3 "#fc7e0f"
//        //polyFit poly1 "#2ca02c"
//        //polyFit poly2 "#2ca02c"
//        //polyFit poly3 "#2ca02c"
//        linear lin1 "#1f77b4"
//        linear lin2 "#1f77b4"
//        linear lin3 "#1f77b4"
//        ]
//        |> Chart.combine
//        |> Chart.withSize (1000.,800.)
//        |> Chart.withTemplate ChartTemplates.lightMirrored
//        //|> Chart.show
//    let variances = 
//        [|1. .. 0.1 .. 8.|] 
//        |> Array.map (fun x -> 
//            let varTC = 
//                [
//                    fitB1.SplineFunction x
//                    fitB2.SplineFunction x
//                    fitB3.SplineFunction x
//                ]
//                |> Seq.var
//            let varlin =
//                [
//                    lin1 x
//                    lin2 x
//                    lin3 x
//                ]
//                |> Seq.var
//            let varpoly =
//                [
//                    poly1 x
//                    poly2 x
//                    poly3 x
//                ]
//                |> Seq.var
//            [|varTC;varlin;varpoly|]
//        )
//        |> Array.transpose
//        //|> Array.map (fun x -> Seq.median x)


//    variances//[fitB1.SplineFunction;fitB2.SplineFunction;fitB3.SplineFunction;poly1;poly2;poly3]

//10 in 10 min
    
//let xx = 
//    data
//    |> Array.ofSeq
//    //|> FSharp.Stats.Array.shuffleFisherYates
//    |> fun x -> x.[0..199]
//    |> Array.map procSingle


//xx
//|> Array.transpose
//|> Array.map Array.concat
//|> Array.mapi (fun i x -> 
//    match i with
//    | 0 -> Chart.Histogram((Array.map log x),Name="TC") |> Chart.withAxisTitles "" "count"
//    | 1 -> Chart.Histogram((Array.map log x),Name="lin") |> Chart.withAxisTitles "" "count"
//    | 2 -> Chart.Histogram((Array.map log x),Name="poly") |> Chart.withAxisTitles "log variance" "count"
//    )
//|> Chart.Grid(3,1)
//|> Chart.show












////required because fitting is not possible when data is unbalanced!
//let getBestFit xVal yVal weightingMethod minimizer = 
//    let getMyWeighting (xVal:Vector<float>) (yVal:float[][]) (w_method:Fitting.WeightingMethod) =
//        match w_method with
//        | Equal ->
//            Matrix.diag(Vector.oneCreate xVal.Length)
//        //1/var centered around 1
//        | Variance ->
//            let Wtemp = Matrix.create xVal.Length xVal.Length 0.            
//            let diagVec = 
//                yVal |> Array.map (fun x -> 1. / Seq.var x)
//                //if yVal.Length%numRep = 0 then
//                //    vector [for i = 0 to xVal.Length - 1 do yield 1. / (yVal.[(i*numRep)..(i*numRep + numRep - 1)] |> Seq.var)]
//                //else failwithf "arrLength no multiple of replicate number"
//            let avg = Seq.average diagVec
//            for i = 0 to xVal.Length - 1 do Wtemp.[i,i] <- diagVec.[i] / avg
//            Wtemp
//        //just 1/var not centered
//        | VarPlain ->
//            let Wtemp = Matrix.create xVal.Length xVal.Length 0.            
//            let diagVec =                
//                yVal |> Array.map (fun x -> 1. / Seq.var x)
//                //if yVal.Length%numRep = 0 then
//                //    vector [for i = 0 to xVal.Length - 1 do yield 1. / (yVal.[(i*numRep)..(i*numRep + numRep - 1)] |> Seq.var)]
//                //else failwithf "arrLength no multiple of replicate number"
//            for i = 0 to xVal.Length - 1 do Wtemp.[i,i] <- diagVec.[i]
//            Wtemp
//        | VarRobust -> //Variance
//            let Wtemp = Matrix.create xVal.Length xVal.Length 0.            
//            let diagVec = 
//                yVal |> Array.map (fun x -> 1. / Seq.stDev x) |> vector
//                //if yVal.Length%numRep = 0 then
//                //    vector [for i = 0 to xVal.Length - 1 do yield 1. / (yVal.[(i*numRep)..(i*numRep + numRep - 1)] |> Seq.stDev)]
//                //else failwithf "arrLength no multiple of replicate number"
//            //improves robustness of standard weighting
//            let (maxWtemp,minWtemp) = 
//                let mean = diagVec |> Seq.mean
//                let std  = diagVec |> Seq.stDev
//                (mean + std),(mean - std)
                
//            let diagVecNew = 
//                diagVec |> Vector.map (fun x -> match x with
//                                                | x when x > maxWtemp -> maxWtemp
//                                                | x when x < minWtemp -> minWtemp
//                                                | _ -> x)
//            let finalDiag = 
//                let trace = (diagVecNew |> Vector.sum)
//                let length = (float Wtemp.NumCols)
//                diagVecNew |> Vector.map (fun x -> x / trace * length)
//            for i = 0 to xVal.Length - 1 do Wtemp.[i,i] <- finalDiag.[i]
//            Wtemp
//        | CV -> 
//            let Wtemp = Matrix.create xVal.Length xVal.Length 0.
//            let cvOfVec =
//                yVal 
//                |> Array.map (fun g -> 
//                    max 0. (System.Math.Log((1./(Math.Abs(Seq.cvPopulation g))),2.))
//                    )

//                //if yVal.Length%numRep = 0 then
//                //    let length = yVal.Length / numRep 
//                //    let cv = vector [for i = 0 to length-1 do yield yVal.[i*numRep..i*numRep+numRep-1] |> fun g -> max 0. (System.Math.Log((1./(Math.Abs(Seq.cvPopulation g))),2.))(*(1./((Seq.cvPopulation g) + 0.25))*)] //0.25???
//                //    cv
//                //else failwithf "arrLength no multiple of replicate number"

//            for i = 0 to xVal.Length - 1 do Wtemp.[i,i] <- cvOfVec.[i]
//            Wtemp
//        | StandardDeviation     ->
//            let Wtemp = Matrix.create xVal.Length xVal.Length 0.            
//            let diagVec = 
//                yVal |> Array.map (fun x -> 1. / Seq.stDev x) |> vector

//                //if yVal.Length%numRep = 0 then
//                //    vector [for i = 0 to xVal.Length - 1 do yield 1. / (yVal.[(i*numRep)..(i*numRep + numRep - 1)] |> Seq.stDev)]
//                //else failwithf "arrLength no multiple of replicate number"
//            let avg = Seq.average diagVec
//            for i = 0 to xVal.Length - 1 do Wtemp.[i,i] <- diagVec.[i] / avg
//            Wtemp
//        | StandardDeviationSqrt ->
//            let Wtemp = Matrix.create xVal.Length xVal.Length 0.            
//            let diagVec = 
//                yVal |> Array.map (fun x -> 1. / Seq.stDev x) |> vector
//                //if yVal.Length%numRep = 0 then
//                //    vector [for i = 0 to xVal.Length - 1 do yield Math.Sqrt (1. / (yVal.[(i*numRep)..(i*numRep + numRep - 1)] |> Seq.stDev))]
//                //else failwithf "arrLength no multiple of replicate number"
//            let avg = Seq.average diagVec
//            for i = 0 to xVal.Length - 1 do Wtemp.[i,i] <- diagVec.[i] / avg
//            Wtemp
//        | StandardDeviationAdj -> //CV_Norm          ->
//            let Wtemp = Matrix.create xVal.Length xVal.Length 0.
//            let cvOfVec =
//                yVal 
//                |> Array.map (fun x -> 
//                    sqrt ( 1. / (Seq.stDev x / Seq.average x))
//                ) 
//                |> vector
//                //if yVal.Length%numRep = 0 then
//                //    let length = yVal.Length / numRep 
//                //    let cv = vector [for i = 0 to length-1 do yield yVal.[i*numRep..i*numRep+numRep-1] |> fun x -> sqrt ( 1. / (Seq.stDev x / Seq.average x))]
//                //    cv
//                //else failwithf "arrLength no multiple of replicate number"
//            let cvAvg = Seq.average cvOfVec
//            for i = 0 to xVal.Length - 1 do Wtemp.[i,i] <- cvOfVec.[i] / cvAvg
//            Wtemp
//        | StandardDeviationAdjSqrt ->  //CV_NormSqrt      ->
//            let Wtemp = Matrix.create xVal.Length xVal.Length 0.
//            let cvOfVec =
//                yVal 
//                |> Array.map (fun x -> 
//                    sqrt(sqrt ( 1. / (Seq.stDev x / Seq.average x)))
//                ) 
//                |> vector

//                //if yVal.Length%numRep = 0 then
//                //    let length = yVal.Length / numRep 
//                //    let cv = vector [for i = 0 to length-1 do yield yVal.[i*numRep..i*numRep+numRep-1] |> fun x -> sqrt(sqrt ( 1. / (Seq.stDev x / Seq.average x)))]
//                //    cv
//                //else failwithf "arrLength no multiple of replicate number"
//            let cvAvg = Seq.average cvOfVec
//            for i = 0 to xVal.Length - 1 do Wtemp.[i,i] <- cvOfVec.[i] / cvAvg
//            Wtemp



//    let yValMeans = 
//        //Seq.getMeanOfReplicates repNumber yVal |> vector
//        yVal |> Seq.map Seq.mean
//    let weightingMatrix = (getMyWeighting xVal yVal weightingMethod)
//    let (cl,fit) = Fitting.getBestFitOfWeighting xVal yValMeans weightingMatrix minimizer
//    fit

open Fitting.LinearRegression.OrdinaryLeastSquares





type Difference = 
    | Diff2Mean
    | Diff2Orig

let getStdPlot (input:(float*float) list) colFill = 
    input
    |> List.sortBy fst
    |> List.chunkBySize (input.Length / 19)
    |> List.map (fun x -> 
        let xs,ys = List.unzip x
        //does it make difference if intensities are logged
        // for sure not global, but maybe patterns
        let xd =  List.map log xs 
        Seq.mean xd,Seq.stDev ys
        )
    |> fun t -> Chart.Line(t,LineColor=Color.fromString colFill)

//////////////////////////////////////////////////
//////////////total Timepoints are removed ///////
//////////////////////////////////////////////////
let removeTimepoints (yData: float[][]) showPlot minimizer weighting =
    let tp = [|0. .. 7.|]
    let means = yData |> Array.map Seq.mean
    let xs,ys = yData |> Array.mapi (fun i x -> x |> Array.map (fun j -> float i,j)) |> Array.concat |> Array.unzip
    
    //TempClass
    //let tc  = myGetBestFit (vector tp) yData Fitting.WeightingMethod.StandardDeviationAdj Fitting.Minimizer.GCV
    //let tc  = myGetBestFit (vector tp) yData Fitting.WeightingMethod.Variance Fitting.Minimizer.AICc
    let tc  = Fitting.getBestFit (vector tp) yData weighting minimizer |> fst
    //linear
    let linear  = Linear.Univariable.fit (Linear.Univariable.coefficient (vector xs) (vector ys))
    //linSpline
    let linSpline  = Interpolation.LinearSpline.interpolate (Interpolation.LinearSpline.initInterpolate tp means)
    //poly
    let poly  = Interpolation.Polynomial.fit (Interpolation.Polynomial.coefficients (vector tp) (vector means))
    //spline
    let spline  = Interpolation.CubicSpline.Simple.fit (Interpolation.CubicSpline.Simple.coefficients Interpolation.CubicSpline.Simple.BoundaryCondition.Natural (vector tp) (vector means)) (vector tp)

    let classif = Classification.getClassification (vector tp) tc.TraceA tc.TraceC 0.05 1.

    let errors = 
        [|1. .. 6.|]
        |> Array.map (fun tpRem ->
            let meanAtRem = means.[int tpRem]
            let varAtRem = yData.[int tpRem] |> Seq.var
        
            let tpR = tp |> Array.filter (fun x -> x <> tpRem)
            let yDataR = tpR |> Array.map (fun tp -> yData.[int tp])
            let xsR,ysR = Array.zip xs ys |> Array.filter (fun (i,x) -> i <> tpRem) |> Array.unzip
            let meansR = tpR |> Array.map (fun tp -> means.[int tp])

            //TempClass
            //let tcR = getBestFit (vector tpR) yDataR Fitting.WeightingMethod.StandardDeviationAdj Fitting.Minimizer.GCV
            let tcR = Fitting.getBestFit (vector tpR) yDataR weighting minimizer |> fst

            //linear
            let linearR = Linear.Univariable.fit (Linear.Univariable.coefficient (vector xsR) (vector ysR))

            //linSpline
            let linSplineR = Interpolation.LinearSpline.interpolate (Interpolation.LinearSpline.initInterpolate tpR meansR)

            //poly
            let polyR = Interpolation.Polynomial.fit (Interpolation.Polynomial.coefficients (vector tpR) (vector meansR))

            //cubicSpline
            let splineR = Interpolation.CubicSpline.Simple.fit (Interpolation.CubicSpline.Simple.coefficients Interpolation.CubicSpline.Simple.BoundaryCondition.Natural (vector tpR) (vector meansR)) (vector tpR)

            let difTC   = tc.SplineFunction tpRem - tcR.SplineFunction tpRem
            let difLin  = linear tpRem - linearR tpRem
            let diflSpl = linSpline tpRem - linSplineR tpRem
            let difPoly = poly tpRem - polyR tpRem
            let difCubic = spline tpRem - splineR tpRem
            
            let difMeanTC   = meanAtRem - tcR.SplineFunction tpRem
            let difMeanLin  = meanAtRem - linearR tpRem
            let difMeanlSpl = meanAtRem - linSplineR tpRem
            let difMeanPoly = meanAtRem - polyR tpRem
            let difMeanCubic = meanAtRem - splineR tpRem
            

            let chart = 
                if showPlot then 
                    [
                        //[
                        //    Chart.Point(xs,ys,MarkerColor =Color.fromString"blue") |> Chart.withLegend(false)
                        //    [0. .. 0.1 .. 7.] |> List.map (fun x -> x,linear x) |> fun t -> Chart.Line (t,LineColor =Color.fromString"blue",Name="lin")
                        //
                        //    Chart.Point(xsR,ysR,MarkerColor =Color.fromString"green") |> Chart.withLegend(false)
                        //    [0. .. 0.1 .. 7.] |> List.map (fun x -> x,linearR x) |> fun t -> Chart.Line (t,LineColor =Color.fromString"green",Name="lin")
                        //]
                        //|> Chart.combine
                        //|> Chart.withAxisTitles "" "lin"
                        [
                            Chart.Point(xs,ys,MarkerColor =Color.fromString"blue")
                            [0. .. 0.1 .. 7.] |> List.map (fun x -> x,tc.SplineFunction x) |> fun t -> Chart.Line (t,LineColor =Color.fromString"blue",Name="tc")
            
                            Chart.Point(xsR,ysR,MarkerColor =Color.fromString"green")
                            [0. .. 0.1 .. 7.] |> List.map (fun x -> x,tcR.SplineFunction x) |> fun t -> Chart.Line (t,LineColor =Color.fromString"green",Name="tc") |> Chart.withXAxisStyle "TC"
                        ]
                        |> Chart.combine
                        |> Chart.withAxisTitles "TC" ""
                        |> Chart.withXAxisStyle (ZeroLine=false)
                        [
                            Chart.Point(xs,ys,MarkerColor =Color.fromString"blue")
                            [0. .. 0.1 .. 7.] |> List.map (fun x -> x,polyR x) |> fun t -> Chart.Line (t,LineColor =Color.fromString"green",Name="poly")
            
                            Chart.Point(xsR,ysR,MarkerColor =Color.fromString"green")
                            [0. .. 0.1 .. 7.] |> List.map (fun x -> x,poly x) |> fun t -> Chart.Line (t,LineColor =Color.fromString"blue",Name="poly")
                        ]
                        |> Chart.combine
                        |> Chart.withAxisTitles "polynomial" ""
                        |> Chart.withXAxisStyle (ZeroLine=false)
                        [
                            Chart.Point(xs,ys,MarkerColor =Color.fromString"blue")
                            [0. .. 0.1 .. 7.] |> List.map (fun x -> x,linSpline x) |> fun t -> Chart.Line (t,LineColor =Color.fromString"blue",Name="linSpl")
            
                            Chart.Point(xsR,ysR,MarkerColor =Color.fromString"green")
                            [0. .. 0.1 .. 7.] |> List.map (fun x -> x,linSplineR x) |> fun t -> Chart.Line (t,LineColor =Color.fromString"green",Name="linSpl")
                        ]
                        |> Chart.combine
                        |> Chart.withAxisTitles "linear spline" ""
                        |> Chart.withXAxisStyle (ZeroLine=false)
                        [
                            Chart.Point(xs,ys,MarkerColor =Color.fromString"blue")
                            [0. .. 0.1 .. 7.] |> List.map (fun x -> x,spline x) |> fun t -> Chart.Line (t,LineColor =Color.fromString"blue",Name="spline")
            
                            Chart.Point(xsR,ysR,MarkerColor =Color.fromString"green") 
                            [0. .. 0.1 .. 7.] |> List.map (fun x -> x,splineR x) |> fun t -> Chart.Line (t,LineColor =Color.fromString"green",Name="spline")
                        ]
                        |> Chart.combine
                        |> Chart.withAxisTitles "cubic spline" ""
                        |> Chart.withXAxisStyle (ZeroLine=false)
                    ]
                    |> Chart.Grid(1,4)
                    |> Chart.withDescription 
                        [
                        Giraffe.ViewEngine.HtmlElements.rawText classif
                        Giraffe.ViewEngine.HtmlElements.rawText "<br>"
                        Giraffe.ViewEngine.HtmlElements.rawText (sprintf "TC: %f<br>lin: %f<br>linSpl: %f<br>Poly: %f<br>Spline: %f" difTC difLin diflSpl difPoly difCubic)
                        ]
                        
                    |> Chart.withSize(1200.,450.)
                    |> Chart.show
                else ()

            varAtRem,[|difMeanTC;difMeanLin;difMeanlSpl;difMeanPoly;difMeanCubic|],[|difTC;difLin;diflSpl;difPoly;difCubic|]
        )
    classif,errors


let processTimePointDelete minimizer weighting = 

    let storageMean,storageOrig =   
        let m = 
            match minimizer with
            | Fitting.Minimizer.AICc -> "_AIC"
            | Fitting.Minimizer.GCV -> "_GCV"
        let w = 
            match weighting with
            | Fitting.WeightingMethod.Variance -> "_Variance"
            | Fitting.WeightingMethod.VarRobust -> "_VarRobust"
            | Fitting.WeightingMethod.VarPlain -> "_VarPlain"
            | Fitting.WeightingMethod.CV -> "_CV"
            | Fitting.WeightingMethod.StandardDeviationAdj -> "_StandardDeviationAdj"
            | Fitting.WeightingMethod.StandardDeviationAdjSqrt -> "_StandardDeviationAdjSqrt"
            | Fitting.WeightingMethod.StandardDeviationSqrt -> "_StandardDeviationSqrt"
            | Fitting.WeightingMethod.StandardDeviation -> "_StandardDeviation"
            | Fitting.WeightingMethod.Equal -> "_Equal"


        @"C:\Users\venn\source\repos\TemporalClassification\runs\Testing\remove\removeTP\DistanceDensityToMean_" + w + m,
        @"C:\Users\venn\source\repos\TemporalClassification\runs\Testing\remove\removeTP\DistanceDensityToOrig_" + w + m

    let (tpRemoveWithVar2Mean,tpRemoveWithVar2Orig) = 
        data
        |> Array.map (fun (id,x) -> id,x ) 
        |> PSeq.mapi (fun i (id,x) -> 
            if i%10=0 then printfn "%i" i
            removeTimepoints (x |> Array.chunkBySize 3) false minimizer weighting
            )
        |> PSeq.withDegreeOfParallelism 12
        |> PSeq.toArray
        |> Array.map (fun (id,x) -> 
            let mean,orig = x |> Array.map (fun (a,b,c) -> (a,b),(a,c)) |> Array.unzip
            (id,mean),(id,orig)
            )
        |> Array.unzip

    //////////////////////////////////////////////////
    //////////////total Timepoints are removed and distance to original prediction and var are reported ///////
    //////////////////////////////////////////////////
    let mutable XcolTC   : (float*float) list= []
    let mutable XcolLin  : (float*float) list= []
    let mutable XcolSpl  : (float*float) list= []
    let mutable XcolPoly : (float*float) list= []
    let mutable XcolCubi : (float*float) list= []

    tpRemoveWithVar2Orig
    |> Array.map (fun x -> 
        snd x
        |> Array.map (fun (v,p) -> 
            XcolTC <-   (sqrt v,p.[0])::XcolTC
            XcolLin <-  (sqrt v,p.[1])::XcolLin
            XcolSpl <-  (sqrt v,p.[2])::XcolSpl
            XcolPoly <- (sqrt v,p.[3])::XcolPoly
            XcolCubi <- (sqrt v,p.[4])::XcolCubi
        )
    ) |> ignore

    ////Histos of distance to orig
    //[
    //XcolTC   |> List.map snd |> Chart.Histogram |> Chart.withTraceInfo (sprintf "colTC  : %f" (Seq.var (XcolTC   |> List.map snd))) |> Chart.withAxisTitles "Distance to OrigFit" "#count"
    //XcolLin  |> List.map snd |> Chart.Histogram |> Chart.withTraceInfo (sprintf "colLin : %f" (Seq.var (XcolLin  |> List.map snd))) |> Chart.withAxisTitles "Distance to OrigFit" "#count"
    //XcolSpl  |> List.map snd |> Chart.Histogram |> Chart.withTraceInfo (sprintf "colSpl : %f" (Seq.var (XcolSpl  |> List.map snd))) |> Chart.withAxisTitles "Distance to OrigFit" "#count"
    //XcolPoly |> List.map snd |> Chart.Histogram |> Chart.withTraceInfo (sprintf "colPoly: %f" (Seq.var (XcolPoly |> List.map snd))) |> Chart.withAxisTitles "Distance to OrigFit" "#count"
    //XcolCubi |> List.map snd |> Chart.Histogram |> Chart.withTraceInfo (sprintf "colCubic: %f" (Seq.var (XcolCubi |> List.map snd))) |> Chart.withAxisTitles "Distance to OrigFit" "#count"
    //]
    //|> Chart.Grid(5,1)
    //|> Chart.withSize(900.,1300.)
    //|> Chart.withTitle "remove complete TP, DistanceToOrigFit"
    //|> Chart.show

    //[
    //Chart.Point(XcolTC,MarkerColor=Color.fromString "#2274b7") |> Chart.withTraceInfo (sprintf "colTC  : %f" (Seq.var (XcolTC   |> List.map snd))) |> Chart.withAxisTitlesYR "Var@TP" "Distance to OrigFit (TC)" (-2.,2.)
    //Chart.Histogram(XcolTC   |> List.map snd,MarkerColor=Color.fromString "#2274b7") |> Chart.withTraceInfo (sprintf "colTC  : %f" (Seq.var (XcolTC   |> List.map snd))) |> Chart.withAxisTitlesXR "Distance to OrigFit" "#count" (-0.8,0.8)
    ////Chart.Point(XcolLin,MarkerColor=Color.fromString "#ff8208") |> Chart.withTraceInfo (sprintf "colLin : %f" (Seq.var (XcolLin  |> List.map snd))) |> Chart.withAxisTitlesYR "Var@TP" "Distance to OrigFit" (-2.,2.)
    ////Chart.Histogram(XcolLin   |> List.map snd,MarkerColor=Color.fromString "#ff8208") |> Chart.withTraceInfo (sprintf "colLin : %f" (Seq.var (XcolLin  |> List.map snd))) |> Chart.withAxisTitlesXR "Distance to OrigFit" "#count" (-0.8,0.8)
    //Chart.Point(XcolSpl,MarkerColor=Color.fromString "#2c9f2c") |> Chart.withTraceInfo (sprintf "colSpl : %f" (Seq.var (XcolSpl  |> List.map snd))) |> Chart.withAxisTitlesYR "Var@TP" "Distance to OrigFit (LinSpl)" (-2.,2.)
    //Chart.Histogram(XcolSpl   |> List.map snd,MarkerColor=Color.fromString "#2c9f2c") |> Chart.withTraceInfo (sprintf "colSpl : %f" (Seq.var (XcolSpl  |> List.map snd))) |> Chart.withAxisTitlesXR "Distance to OrigFit" "#count" (-0.8,0.8)
    //Chart.Point(XcolPoly,MarkerColor=Color.fromString "#d62728") |> Chart.withTraceInfo (sprintf "colPoly: %f" (Seq.var (XcolPoly |> List.map snd))) |> Chart.withAxisTitlesYR "Var@TP" "Distance to OrigFit (Poly)" (-2.,2.)
    //Chart.Histogram(XcolPoly   |> List.map snd,MarkerColor=Color.fromString "#d62728") |> Chart.withTraceInfo (sprintf "colPoly: %f" (Seq.var (XcolPoly |> List.map snd))) |> Chart.withAxisTitlesXR "Distance to OrigFit" "#count" (-0.8,0.8)
    //Chart.Point(XcolCubi,MarkerColor=Color.fromString "grey") |> Chart.withTraceInfo (sprintf "colCubic: %f" (Seq.var (XcolCubi |> List.map snd))) |> Chart.withAxisTitlesYR "Var@TP" "Distance to OrigFit (Cubic)" (-2.,2.)
    //Chart.Histogram(XcolCubi   |> List.map snd,MarkerColor=Color.fromString "grey") |> Chart.withTraceInfo (sprintf "colCubic: %f" (Seq.var (XcolCubi |> List.map snd))) |> Chart.withAxisTitlesXR "Distance to OrigFit" "#count" (-0.8,0.8)
    //]
    //|> Chart.Grid(4,2)
    //|> Chart.withSize(1600.,1300.)
    //|> Chart.withTitle "remove complete TP, variance@TP to DistanceToOrigFit"
    //|> Chart.saveHtml(@"C:\Users\venn\source\repos\TemporalClassification\runs\Testing\remove\removeTP\DistanceScatterToOrigFit_varaic",true)

    [
    List.unzip XcolTC |> fun (a,b) -> Chart.PointDensity(List.map log a,b,PointMarkerColor= Color.fromString "#2274b7") |> Chart.withTraceInfo (sprintf "colTC  : %f" (Seq.var (XcolTC   |> List.map snd))) |> Chart.withAxisTitlesYR "log Var@TP" "Distance to OrigFit" (-0.4,0.4)
    Chart.Histogram(XcolTC   |> List.map snd,MarkerColor=Color.fromString "#2274b7") |> Chart.withTraceInfo (sprintf "colTC  : %f" (Seq.var (XcolTC   |> List.map snd))) |> Chart.withAxisTitlesXR "Distance to OrigFit (TC)" "#count" (-0.8,0.8)
    //List.unzip XcolLin |> fun (a,b) -> Chart.PointDensity(List.map log a,b,PointMarkerColor= Color.fromString "#ff8208") |> Chart.withTraceInfo (sprintf "colLin : %f" (Seq.var (XcolLin  |> List.map snd))) |> Chart.withAxisTitlesYR "log Var@TP" "Distance to OrigFit" (-0.4,0.4)
    //Chart.Histogram(XcolLin   |> List.map snd,MarkerColor=Color.fromString "#ff8208") |> Chart.withTraceInfo (sprintf "colLin : %f" (Seq.var (XcolLin  |> List.map snd))) |> Chart.withAxisTitlesXR "Distance to OrigFit" "#count" (-0.8,0.8)
    List.unzip XcolSpl |> fun (a,b) -> Chart.PointDensity(List.map log a,b,PointMarkerColor= Color.fromString "#2c9f2c") |> Chart.withTraceInfo (sprintf "colSpl : %f" (Seq.var (XcolSpl  |> List.map snd))) |> Chart.withAxisTitlesYR "log Var@TP" "Distance to OrigFit" (-0.4,0.4)
    Chart.Histogram(XcolSpl   |> List.map snd,MarkerColor=Color.fromString "#2c9f2c") |> Chart.withTraceInfo (sprintf "colSpl : %f" (Seq.var (XcolSpl  |> List.map snd))) |> Chart.withAxisTitlesXR "Distance to OrigFit (LinSpl)" "#count" (-0.8,0.8)
    List.unzip XcolPoly |> fun (a,b) -> Chart.PointDensity(List.map log a,b,PointMarkerColor= Color.fromString "#d62728") |> Chart.withTraceInfo (sprintf "colPoly: %f" (Seq.var (XcolPoly |> List.map snd))) |> Chart.withAxisTitlesYR "log Var@TP" "Distance to OrigFit" (-0.4,0.4)
    Chart.Histogram(XcolPoly   |> List.map snd,MarkerColor=Color.fromString "#d62728") |> Chart.withTraceInfo (sprintf "colPoly: %f" (Seq.var (XcolPoly |> List.map snd))) |> Chart.withAxisTitlesXR "Distance to OrigFit (Poly)" "#count" (-0.8,0.8)
    List.unzip XcolCubi |> fun (a,b) -> Chart.PointDensity(List.map log a,b,PointMarkerColor= Color.fromString "grey") |> Chart.withTraceInfo (sprintf "colCubic: %f" (Seq.var (XcolCubi |> List.map snd))) |> Chart.withAxisTitlesYR "log Var@TP" "Distance to OrigFit" (-0.4,0.4)
    Chart.Histogram(XcolCubi   |> List.map snd,MarkerColor=Color.fromString "grey") |> Chart.withTraceInfo (sprintf "colCubic: %f" (Seq.var (XcolCubi |> List.map snd))) |> Chart.withAxisTitlesXR "Distance to OrigFit (Cubic)" "#count" (-0.8,0.8)
    ]
    |> Chart.Grid(4,2)
    |> Chart.withSize(1600.,1300.)
    |> Chart.withTitle "remove complete TP, variance@TP to DistanceToOrigFit"
    |> Chart.saveHtml(storageOrig,true)


    //////////////////////////////////////////////////
    //////////////total Timepoints are removed and distance to mean and var are reported ///////
    //////////////////////////////////////////////////
    let mutable X2MeancolTC   : (float*float) list= []
    let mutable X2MeancolLin  : (float*float) list= []
    let mutable X2MeancolSpl  : (float*float) list= []
    let mutable X2MeancolPoly : (float*float) list= []
    let mutable X2MeancolCubic : (float*float) list= []

    tpRemoveWithVar2Mean
    |> Array.map (fun x -> 
        snd x
        |> Array.map (fun (v,p) -> 
            X2MeancolTC   <-   (sqrt v,p.[0])::X2MeancolTC  
            X2MeancolLin  <-   (sqrt v,p.[1])::X2MeancolLin 
            X2MeancolSpl  <-   (sqrt v,p.[2])::X2MeancolSpl 
            X2MeancolPoly <-   (sqrt v,p.[3])::X2MeancolPoly
            X2MeancolCubic <-   (sqrt v,p.[4])::X2MeancolCubic
        )
    ) |> ignore

    //Scatter Var to distance 2Mean
    //[
    //X2MeancolTC   |> Chart.Point |> Chart.withTraceInfo (sprintf "colTC  : %f" (Seq.var (X2MeancolTC   |> List.map snd))) |> Chart.withAxisTitles "Distance to orig fit" "#count"
    ////X2MeancolLin  |> Chart.Point |> Chart.withTraceInfo (sprintf "colLin : %f" (Seq.var (X2MeancolLin  |> List.map snd))) |> Chart.withAxisTitles "Distance to orig fit" "#count"
    //X2MeancolSpl  |> Chart.Point |> Chart.withTraceInfo (sprintf "colSpl : %f" (Seq.var (X2MeancolSpl  |> List.map snd))) |> Chart.withAxisTitles "Distance to orig fit" "#count"
    //X2MeancolPoly |> Chart.Point |> Chart.withTraceInfo (sprintf "colPoly: %f" (Seq.var (X2MeancolPoly |> List.map snd))) |> Chart.withAxisTitles "Distance to orig fit" "#count"
    //X2MeancolCubic |> Chart.Point |> Chart.withTraceInfo (sprintf "colCubic: %f" (Seq.var (X2MeancolCubic |> List.map snd))) |> Chart.withAxisTitles "Distance to orig fit" "#count"
    //]
    //|> Chart.Grid(4,1)
    //|> Chart.withSize(900.,1300.)
    //|> Chart.withTitle "remove complete TP, variance@TP to DistanceToTPmean"
    //|> Chart.show

    //[
    //X2MeancolTC   |> List.map snd |> Chart.Histogram |> Chart.withTraceInfo (sprintf "colTC  : %f" (Seq.var (X2MeancolTC   |> List.map snd))) |> Chart.withAxisTitles "" "count"
    ////X2MeancolLin  |> List.map snd |> Chart.Histogram |> Chart.withTraceInfo (sprintf "colLin : %f" (Seq.var (X2MeancolLin  |> List.map snd))) |> Chart.withAxisTitles "" "count"
    //X2MeancolSpl  |> List.map snd |> Chart.Histogram |> Chart.withTraceInfo (sprintf "colSpl : %f" (Seq.var (X2MeancolSpl  |> List.map snd))) |> Chart.withAxisTitles "" "count"
    //X2MeancolPoly |> List.map snd |> Chart.Histogram |> Chart.withTraceInfo (sprintf "colPoly: %f" (Seq.var (X2MeancolPoly |> List.map snd))) |> Chart.withAxisTitles "" "count"
    //X2MeancolCubic |> List.map snd |> Chart.Histogram |> Chart.withTraceInfo (sprintf "colCubic: %f" (Seq.var (X2MeancolCubic |> List.map snd))) |> Chart.withAxisTitles "Distance to mean" "count"
    //]
    //|> Chart.Grid(4,1)
    //|> Chart.withTitle "remove complete TP, DistanceToTPmean"
    //|> Chart.saveHtml(@"C:\Users\venn\source\repos\TemporalClassification\runs\Testing\remove\removeTP\DistanceToMean_varaic",true)
    //
    //[
    //Chart.Point(X2MeancolTC,MarkerColor=Color.fromString "#2274b7") |> Chart.withTraceInfo (sprintf "colTC  : %f" (Seq.var (X2MeancolTC   |> List.map snd))) |> Chart.withAxisTitlesYR "Var@TP" "Distance to TPMean (TC)" (-2.,2.)
    //Chart.Histogram(X2MeancolTC   |> List.map snd,MarkerColor=Color.fromString "#2274b7") |> Chart.withTraceInfo (sprintf "colTC  : %f" (Seq.var (X2MeancolTC   |> List.map snd))) |> Chart.withAxisTitlesXR "Distance to Mean" "#count" (-0.8,0.8)
    ////Chart.Point(X2MeancolLin,MarkerColor=Color.fromString "#ff8208") |> Chart.withTraceInfo (sprintf "colLin : %f" (Seq.var (X2MeancolLin  |> List.map snd))) |> Chart.withAxisTitlesYR "Var@TP" "Distance to TPMean" (-2.,2.)
    ////Chart.Histogram(X2MeancolLin   |> List.map snd,MarkerColor=Color.fromString "#ff8208") |> Chart.withTraceInfo (sprintf "colLin : %f" (Seq.var (X2MeancolLin  |> List.map snd))) |> Chart.withAxisTitlesXR "Distance to Mean" "#count" (-0.8,0.8)
    //Chart.Point(X2MeancolSpl,MarkerColor=Color.fromString "#2c9f2c") |> Chart.withTraceInfo (sprintf "colSpl : %f" (Seq.var (X2MeancolSpl  |> List.map snd))) |> Chart.withAxisTitlesYR "Var@TP" "Distance to TPMean (LinSpline)" (-2.,2.)
    //Chart.Histogram(X2MeancolSpl   |> List.map snd,MarkerColor=Color.fromString "#2c9f2c") |> Chart.withTraceInfo (sprintf "colSpl : %f" (Seq.var (X2MeancolSpl  |> List.map snd))) |> Chart.withAxisTitlesXR "Distance to Mean" "#count" (-0.8,0.8)
    //Chart.Point(X2MeancolPoly,MarkerColor=Color.fromString "#d62728") |> Chart.withTraceInfo (sprintf "colPoly: %f" (Seq.var (X2MeancolPoly |> List.map snd))) |> Chart.withAxisTitlesYR "Var@TP" "Distance to TPMean (Polynomial)" (-2.,2.)
    //Chart.Histogram(X2MeancolPoly   |> List.map snd,MarkerColor=Color.fromString "#d62728") |> Chart.withTraceInfo (sprintf "colPoly: %f" (Seq.var (X2MeancolPoly |> List.map snd))) |> Chart.withAxisTitlesXR "Distance to Mean" "#count" (-0.8,0.8)
    //Chart.Point(X2MeancolCubic,MarkerColor=Color.fromString "grey") |> Chart.withTraceInfo (sprintf "colCubic: %f" (Seq.var (X2MeancolCubic |> List.map snd))) |> Chart.withAxisTitlesYR "Var@TP" "Distance to TPMean (Cubic)" (-2.,2.)
    //Chart.Histogram(X2MeancolCubic   |> List.map snd,MarkerColor=Color.fromString "grey") |> Chart.withTraceInfo (sprintf "colCubic: %f" (Seq.var (X2MeancolCubic |> List.map snd))) |> Chart.withAxisTitlesXR "Distance to Mean" "#count" (-0.8,0.8)
    //]
    //|> Chart.Grid(4,2)
    //|> Chart.withSize(1600.,1300.)
    //|> Chart.withTitle "remove complete TP, variance@TP to DistanceToTPmean"
    //|> Chart.saveHtml(@"C:\Users\venn\source\repos\TemporalClassification\runs\Testing\remove\removeTP\DistanceScatterToMean_varaic",true)

    [
    List.unzip X2MeancolTC |> fun (a,b) -> Chart.PointDensity(List.map log a,b,PointMarkerColor= Color.fromString "#2274b7") |> Chart.withTraceInfo (sprintf "colTC  : %f" (Seq.var (X2MeancolTC   |> List.map snd))) |> Chart.withAxisTitlesYR "log Var@TP" "Distance to TPMean" (-0.4,0.4)
    Chart.Histogram(X2MeancolTC   |> List.map snd,MarkerColor=Color.fromString "#2274b7") |> Chart.withTraceInfo (sprintf "colTC  : %f" (Seq.var (X2MeancolTC   |> List.map snd))) |> Chart.withAxisTitlesXR "Distance to Mean (TC)" "#count" (-0.8,0.8)
    //List.unzip X2MeancolLin |> fun (a,b) -> Chart.PointDensity(List.map log a,b,PointMarkerColor= Color.fromString "#ff8208") |> Chart.withTraceInfo (sprintf "colLin : %f" (Seq.var (X2MeancolLin  |> List.map snd))) |> Chart.withAxisTitlesYR "log Var@TP" "Distance to TPMean" (-0.4,0.4)
    //Chart.Histogram(X2MeancolLin   |> List.map snd,MarkerColor=Color.fromString "#ff8208") |> Chart.withTraceInfo (sprintf "colLin : %f" (Seq.var (X2MeancolLin  |> List.map snd))) |> Chart.withAxisTitlesXR "Distance to Mean" "#count" (-0.8,0.8)
    List.unzip X2MeancolSpl |> fun (a,b) -> Chart.PointDensity(List.map log a,b,PointMarkerColor= Color.fromString "#2c9f2c") |> Chart.withTraceInfo (sprintf "colSpl : %f" (Seq.var (X2MeancolSpl  |> List.map snd))) |> Chart.withAxisTitlesYR "log Var@TP" "Distance to TPMean" (-0.4,0.4)
    Chart.Histogram(X2MeancolSpl   |> List.map snd,MarkerColor=Color.fromString "#2c9f2c") |> Chart.withTraceInfo (sprintf "colSpl : %f" (Seq.var (X2MeancolSpl  |> List.map snd))) |> Chart.withAxisTitlesXR "Distance to Mean (LinSpline)" "#count" (-0.8,0.8)
    List.unzip X2MeancolPoly |> fun (a,b) -> Chart.PointDensity(List.map log a,b,PointMarkerColor= Color.fromString "#d62728") |> Chart.withTraceInfo (sprintf "colPoly: %f" (Seq.var (X2MeancolPoly |> List.map snd))) |> Chart.withAxisTitlesYR "log Var@TP" "Distance to TPMean" (-0.4,0.4)
    Chart.Histogram(X2MeancolPoly   |> List.map snd,MarkerColor=Color.fromString "#d62728") |> Chart.withTraceInfo (sprintf "colPoly: %f" (Seq.var (X2MeancolPoly |> List.map snd))) |> Chart.withAxisTitlesXR "Distance to Mean (Polynomial)" "#count" (-0.8,0.8)
    List.unzip X2MeancolCubic |> fun (a,b) -> Chart.PointDensity(List.map log a,b,PointMarkerColor= Color.fromString "grey") |> Chart.withTraceInfo (sprintf "colPoly: %f" (Seq.var (X2MeancolCubic |> List.map snd))) |> Chart.withAxisTitlesYR "log Var@TP" "Distance to TPMean" (-0.4,0.4)
    Chart.Histogram(X2MeancolCubic   |> List.map snd,MarkerColor=Color.fromString "grey") |> Chart.withTraceInfo (sprintf "colPoly: %f" (Seq.var (X2MeancolCubic |> List.map snd))) |> Chart.withAxisTitlesXR "Distance to Mean (Cubic)" "#count" (-0.8,0.8)
    ]
    |> Chart.Grid(4,2)
    |> Chart.withSize(1600.,1300.)
    |> Chart.withTitle "remove complete TP, variance@TP to DistanceToTPmean"
    |> Chart.saveHtml(storageMean,true)

    


    [
    getStdPlot XcolTC   "#2274b7" |> Chart.withTraceInfo "TC"
    //getStdPlot XcolLin  "#ff7f0e" |> Chart.withTraceInfo "Lin"
    getStdPlot XcolSpl  "#2ca02c" |> Chart.withTraceInfo "Spl"
    getStdPlot XcolPoly "#d62728" |> Chart.withTraceInfo "Poly"
    getStdPlot XcolCubi "grey" |> Chart.withTraceInfo "Cubic"
    ]
    |> Chart.combine
    |> Chart.withAxisTitles "log var@TP" "stDev of Distances to Orig fit" 
    |> Chart.saveHtml (storageOrig + "standardDeviationOfDistances",true)

    [
    getStdPlot X2MeancolTC   "#2274b7" |> Chart.withTraceInfo "TC"
    //getStdPlot X2MeancolLin  "#ff7f0e" |> Chart.withTraceInfo "Lin"
    getStdPlot X2MeancolSpl  "#2ca02c" |> Chart.withTraceInfo "Spl"
    getStdPlot X2MeancolPoly "#d62728" |> Chart.withTraceInfo "Poly"
    getStdPlot X2MeancolCubic "grey" |> Chart.withTraceInfo "Cubic"
    ]
    |> Chart.combine
    |> Chart.withAxisTitles "log var@TP" "stDev of Distances to mean"
    |> Chart.saveHtml (storageMean + "standardDeviationOfDistances",true)

open TemporalClassification.Fitting
    
processTimePointDelete Minimizer.GCV  WeightingMethod.VarPlain
processTimePointDelete Minimizer.GCV  WeightingMethod.Variance
processTimePointDelete Minimizer.GCV  WeightingMethod.StandardDeviationAdj
processTimePointDelete Minimizer.GCV  WeightingMethod.StandardDeviation
processTimePointDelete Minimizer.GCV  WeightingMethod.CV
printfn "All GCV Combinations of TP deletion ready!"

processTimePointDelete Minimizer.AICc WeightingMethod.VarPlain
processTimePointDelete Minimizer.AICc WeightingMethod.Variance
processTimePointDelete Minimizer.AICc WeightingMethod.StandardDeviationAdj
processTimePointDelete Minimizer.AICc WeightingMethod.StandardDeviation
processTimePointDelete Minimizer.AICc WeightingMethod.CV
printfn "All AIC Combinations of TP deletion ready!"

//cursor 230815 14:11 Here!!


(*
Histograme von
	- Distanz zu Orig wenn TP fehlt - ok
	- Distanz zu Orig wenn Rep fehlt - ok
	- Distanz zu Mean wenn TP fehlt - 
	- Distanz zu Orig wenn Rep fehlt - 
ScatterPlots von
	- Var vs Distance to Orig wenn TP fehlt
	- Var vs Distance to Mean wenn TP fehlt
	- Var vs Distance to Orig wenn Rep fehlt
	- Var vs Distance to Mean wenn Rep fehlt

Was wurde bei complexomics gemacht?
- Silouette index auf Replikaten=?!

*)



































//////////////////////////////////////////////////
//////////////single relicate points are removed ///////
//////////////////////////////////////////////////

let removeSingleReplicates (yData: float[][]) minimizer weighting =
    let tp = [|0. .. 7.|]
    let means = yData |> Array.map Seq.mean
    let xs,ys = yData |> Array.mapi (fun i x -> x |> Array.map (fun j -> float i,j)) |> Array.concat |> Array.unzip
    
    //TempClass
    //let tc  = getBestFit (vector tp) yData Fitting.WeightingMethod.StandardDeviationAdj Fitting.Minimizer.GCV
    //let tc  = myGetBestFit (vector tp) yData Fitting.WeightingMethod.Variance Fitting.Minimizer.AICc |> fst
    let tc  = Fitting.getBestFit (vector tp) yData weighting minimizer |> fst
    //linear
    let linear  = Linear.Univariable.fit (Linear.Univariable.coefficient (vector xs) (vector ys))
    //linSpline
    let linSpline  = Interpolation.LinearSpline.interpolate (Interpolation.LinearSpline.initInterpolate tp means)
    //poly
    let poly  = Interpolation.Polynomial.fit (Interpolation.Polynomial.coefficients (vector tp) (vector means))
    //cubic spline
    let cubic  = Interpolation.CubicSpline.Simple.fit (Interpolation.CubicSpline.Simple.coefficients Interpolation.CubicSpline.Simple.BoundaryCondition.Natural (vector tp) (vector means)) (vector tp)

    let classif = Classification.getClassification (vector tp) tc.TraceA tc.TraceC 0.05 1.

    let errors = 
        [|1. .. 6.|]
        |> Array.map (fun tpRem ->
            [|0..2|]
            |> Array.map (fun repRem -> 
                let mutable mean = 0.
                let mutable var = 0.
                let yDataR = 
                    yData 
                    |> Array.mapi (fun tp dat -> 
                        if tp = int tpRem then
                            let trData = 
                                Array.except [|repRem|] [|0 .. 2|]
                                |> Array.map (fun z -> dat.[z ])
                            mean <- Seq.mean dat
                            var <- Seq.var dat
                            trData
                        else
                            dat
                    )
                let xsR,ysR = yDataR |> Array.mapi (fun i x -> x |> Array.map (fun t -> float i,t)) |> Array.concat |> Array.unzip
                let meansR = yDataR |> Array.map Seq.mean

                //TempClass
                //let tcR = getBestFit (vector tp) yDataR Fitting.WeightingMethod.StandardDeviationAdj Fitting.Minimizer.GCV
                let tcR = Fitting.getBestFit (vector tp) yDataR weighting minimizer  |> fst

                //linear
                let linearR = Linear.Univariable.fit (Linear.Univariable.coefficient (vector xsR) (vector ysR))

                //linSpline
                let linSplineR = Interpolation.LinearSpline.interpolate (Interpolation.LinearSpline.initInterpolate tp meansR)

                //poly
                let polyR = Interpolation.Polynomial.fit (Interpolation.Polynomial.coefficients (vector tp) (vector meansR))

                //cubic
                let cubicR  = Interpolation.CubicSpline.Simple.fit (Interpolation.CubicSpline.Simple.coefficients Interpolation.CubicSpline.Simple.BoundaryCondition.Natural (vector tp) (vector meansR)) (vector tp)

    
                let difTC   = mean - tcR.SplineFunction tpRem,tc.SplineFunction tpRem - tcR.SplineFunction tpRem
                let difLin  = mean - linearR tpRem           ,                linear tpRem- linearR tpRem
                let diflSpl = mean - linSplineR tpRem        ,                linSpline tpRem- linSplineR tpRem
                let difPoly = mean - polyR tpRem             ,                poly tpRem - polyR tpRem
                let difCubi = mean - cubicR tpRem             ,               cubic tpRem - cubicR tpRem

                var,[|difTC;difLin;diflSpl;difPoly;difCubi|]
            )
        )
    classif,errors


let summarizeSingleRepRemove minimizer weighting = 

    let storageMean,storageOrig =   
        let m = 
            match minimizer with
            | Fitting.Minimizer.AICc -> "_AIC"
            | Fitting.Minimizer.GCV -> "_GCV"
        let w = 
            match weighting with
            | Fitting.WeightingMethod.Variance -> "_Variance"
            | Fitting.WeightingMethod.VarRobust -> "_VarRobust"
            | Fitting.WeightingMethod.VarPlain -> "_VarPlain"
            | Fitting.WeightingMethod.CV -> "_CV"
            | Fitting.WeightingMethod.StandardDeviationAdj -> "_StandardDeviationAdj"
            | Fitting.WeightingMethod.StandardDeviationAdjSqrt -> "_StandardDeviationAdjSqrt"
            | Fitting.WeightingMethod.StandardDeviationSqrt -> "_StandardDeviationSqrt"
            | Fitting.WeightingMethod.StandardDeviation -> "_StandardDeviation"
            | Fitting.WeightingMethod.Equal -> "_Equal"

        @"C:\Users\venn\source\repos\TemporalClassification\runs\Testing\remove\removeSingleRepatTP\DistanceDensityToMean" + w + m,
        @"C:\Users\venn\source\repos\TemporalClassification\runs\Testing\remove\removeSingleRepatTP\SummaryToOrigFit" + w + m


    let singleRepRemoveToMean = 
        data
        |> Array.map (fun (id,x) -> id,x |> Array.map Math.Exp) // unlog
        |> PSeq.mapi (fun i (id,x) -> 
            if i%10=0 then printfn "%i" i
            removeSingleReplicates (x |> Array.chunkBySize 3) minimizer weighting
            )
        |> PSeq.withDegreeOfParallelism 16
        |> PSeq.toArray


    let mutable tomeanTC   : (float*float) list= []
    let mutable tomeanLin  : (float*float) list= []
    let mutable tomeanSpl  : (float*float) list= []
    let mutable tomeanPoly : (float*float) list= []
    let mutable tomeanCubic : (float*float) list= []

    let mutable toOrigTC   : (float*float) list= []
    let mutable toOrigLin  : (float*float) list= []
    let mutable toOrigSpl  : (float*float) list= []
    let mutable toOrigPoly : (float*float) list= []
    let mutable toOrigCubic : (float*float) list= []


    singleRepRemoveToMean
    |> Array.iter (fun x -> 
        x
        |> snd
        |> Array.concat
        |> Array.iter (fun (vari,a) -> 
            let tcDiffToMean,tcDiffToOrig = a.[0]
            let liDiffToMean,liDiffToOrig = a.[1]
            let spDiffToMean,spDiffToOrig = a.[2]
            let poDiffToMean,poDiffToOrig = a.[3]
            let cuDiffToMean,cuDiffToOrig = a.[4]
            tomeanTC   <- (vari,tcDiffToMean)::tomeanTC  
            tomeanLin  <- (vari,liDiffToMean)::tomeanLin 
            tomeanSpl  <- (vari,spDiffToMean)::tomeanSpl 
            tomeanPoly <- (vari,poDiffToMean)::tomeanPoly
            tomeanCubic <- (vari,cuDiffToMean)::tomeanCubic

            toOrigTC   <- (vari,tcDiffToOrig)::toOrigTC  
            toOrigLin  <- (vari,liDiffToOrig)::toOrigLin 
            toOrigSpl  <- (vari,spDiffToOrig)::toOrigSpl 
            toOrigPoly <- (vari,poDiffToOrig)::toOrigPoly
            toOrigCubic <- (vari,cuDiffToOrig)::toOrigCubic
        )
    )

    [
    List.unzip tomeanTC |> fun (a,b) -> Chart.PointDensity(List.map log a,b,PointMarkerColor= Color.fromString "#2274b7") |> Chart.withTraceInfo (sprintf "colTC  : %f " (Seq.stDev (tomeanTC   |> List.map snd))) |> Chart.withAxisTitlesYR "log Var@TP" "Distance to mean@TP" (-0.4,0.4)
    Chart.Histogram(tomeanTC   |> List.map snd,MarkerColor=Color.fromString "#2274b7") |> Chart.withTraceInfo "colTC" |> Chart.withAxisTitlesXR "Distance to mean@TP" "#count" (-0.4,0.4)
    //List.unzip tomeanLin |> fun (a,b) -> Chart.PointDensity(List.map log a,b,PointMarkerColor= Color.fromString "#ff8208") |> Chart.withTraceInfo (sprintf "colLin : %f " (Seq.stDev (tomeanLin  |> List.map snd))) |> Chart.withAxisTitlesYR "log Var@TP" "Distance to mean@TP" (-0.4,0.4)
    //Chart.Histogram(tomeanLin   |> List.map snd,MarkerColor=Color.fromString "#ff8208") |> Chart.withTraceInfo "colLin" |> Chart.withAxisTitlesXR "Distance to mean@TP" "#count" (-0.4,0.4)
    List.unzip tomeanSpl |> fun (a,b) -> Chart.PointDensity(List.map log a,b,PointMarkerColor= Color.fromString "#2c9f2c") |> Chart.withTraceInfo (sprintf "colSpl : %f " (Seq.stDev (tomeanSpl  |> List.map snd))) |> Chart.withAxisTitlesYR "log Var@TP" "Distance to mean@TP" (-0.4,0.4)
    Chart.Histogram(tomeanSpl   |> List.map snd,MarkerColor=Color.fromString "#2c9f2c") |> Chart.withTraceInfo "colSpl" |> Chart.withAxisTitlesXR "Distance to mean@TP" "#count" (-0.4,0.4)
    List.unzip tomeanPoly |> fun (a,b) -> Chart.PointDensity(List.map log a,b,PointMarkerColor= Color.fromString "#d62728") |> Chart.withTraceInfo (sprintf "colPoly: %f " (Seq.stDev (tomeanPoly |> List.map snd))) |> Chart.withAxisTitlesYR "log Var@TP" "Distance to mean@TP" (-0.4,0.4)
    Chart.Histogram(tomeanPoly   |> List.map snd,MarkerColor=Color.fromString "#d62728") |> Chart.withTraceInfo "colPoly" |> Chart.withAxisTitlesXR "Distance to mean@TP" "#count" (-0.4,0.4)
    List.unzip tomeanCubic |> fun (a,b) -> Chart.PointDensity(List.map log a,b,PointMarkerColor= Color.fromString "grey") |> Chart.withTraceInfo (sprintf "colCubic: %f " (Seq.stDev (tomeanCubic |> List.map snd))) |> Chart.withAxisTitlesYR "log Var@TP" "Distance to mean@TP" (-0.4,0.4)
    Chart.Histogram(tomeanCubic   |> List.map snd,MarkerColor=Color.fromString "grey") |> Chart.withTraceInfo "colCubic" |> Chart.withAxisTitlesXR "Distance to mean@TP" "#count" (-0.4,0.4)
    ]
    |> Chart.Grid(4,2)
    |> Chart.withSize(1600.,1300.)
    |> Chart.withTitle "remove single Rep@TP, variance@TP to DistanceToMean"
    |> Chart.saveHtml (storageMean,true)

    [
    List.unzip toOrigTC |> fun (a,b) -> Chart.PointDensity(List.map log a,b,PointMarkerColor= Color.fromString "#2274b7") |> Chart.withTraceInfo (sprintf "colTC  : %f" (Seq.stDev (toOrigTC   |> List.map snd))) |> Chart.withAxisTitlesYR "log Var@TP" "Distance to OrigFit" (-0.4,0.4)
    Chart.Histogram(toOrigTC   |> List.map snd,MarkerColor=Color.fromString "#2274b7") |> Chart.withTraceInfo "colTC" |> Chart.withAxisTitlesXR "Distance to OrigFit" "#count" (-0.4,0.4)
    //List.unzip toOrigLin |> fun (a,b) -> Chart.PointDensity(List.map log a,b,PointMarkerColor= Color.fromString "#ff8208") |> Chart.withTraceInfo (sprintf "colLin : %f" (Seq.stDev (toOrigLin  |> List.map snd))) |> Chart.withAxisTitlesYR "log Var@TP" "Distance to OrigFit" (-0.4,0.4)
    //Chart.Histogram(toOrigLin   |> List.map snd,MarkerColor=Color.fromString "#ff8208") |> Chart.withTraceInfo "colLin" |> Chart.withAxisTitlesXR "Distance to OrigFit" "#count" (-0.4,0.4)
    List.unzip toOrigSpl |> fun (a,b) -> Chart.PointDensity(List.map log a,b,PointMarkerColor= Color.fromString "#2c9f2c") |> Chart.withTraceInfo (sprintf "colSpl : %f" (Seq.stDev (toOrigSpl  |> List.map snd))) |> Chart.withAxisTitlesYR "log Var@TP" "Distance to OrigFit" (-0.4,0.4)
    Chart.Histogram(toOrigSpl   |> List.map snd,MarkerColor=Color.fromString "#2c9f2c") |> Chart.withTraceInfo "colSpl" |> Chart.withAxisTitlesXR "Distance to OrigFit" "#count" (-0.4,0.4)
    List.unzip toOrigPoly |> fun (a,b) -> Chart.PointDensity(List.map log a,b,PointMarkerColor= Color.fromString "#d62728") |> Chart.withTraceInfo (sprintf "colPoly: %f" (Seq.stDev (toOrigPoly |> List.map snd))) |> Chart.withAxisTitlesYR "log Var@TP" "Distance to OrigFit" (-0.4,0.4)
    Chart.Histogram(toOrigPoly   |> List.map snd,MarkerColor=Color.fromString "#d62728") |> Chart.withTraceInfo "colPoly" |> Chart.withAxisTitlesXR "Distance to OrigFit" "#count" (-0.4,0.4)
    List.unzip toOrigCubic |> fun (a,b) -> Chart.PointDensity(List.map log a,b,PointMarkerColor= Color.fromString "grey") |> Chart.withTraceInfo (sprintf "colCubic: %f" (Seq.stDev (toOrigCubic |> List.map snd))) |> Chart.withAxisTitlesYR "log Var@TP" "Distance to OrigFit" (-0.4,0.4)
    Chart.Histogram(toOrigCubic   |> List.map snd,MarkerColor=Color.fromString "grey") |> Chart.withTraceInfo "colCubic" |> Chart.withAxisTitlesXR "Distance to OrigFit" "#count" (-0.4,0.4)
    ]
    |> Chart.Grid(4,2)
    |> Chart.withSize(1600.,1300.)
    |> Chart.withTitle "remove single Rep@TP, variance@TP to DistanceToOrig"
    |> Chart.saveHtml (storageOrig,true)
    
    [
    getStdPlot tomeanTC   "#2274b7" |> Chart.withTraceInfo "TC"
    //getStdPlot tomeanLin  "#ff7f0e" |> Chart.withTraceInfo "Lin"
    getStdPlot tomeanSpl  "#2ca02c" |> Chart.withTraceInfo "Spl"
    getStdPlot tomeanPoly "#d62728" |> Chart.withTraceInfo "Poly"
    getStdPlot tomeanCubic "grey" |> Chart.withTraceInfo "Cubic"
    ]
    |> Chart.combine
    |> Chart.withAxisTitles "log var@TP" "stDev of Distances to Mean"
    |> Chart.saveHtml (storageMean + "standardDeviationOfDistances",true)



    [
    getStdPlot toOrigTC   "#2274b7" |> Chart.withTraceInfo "TC"
    //getStdPlot toOrigLin  "#ff7f0e" |> Chart.withTraceInfo "Lin"
    getStdPlot toOrigSpl  "#2ca02c" |> Chart.withTraceInfo "Spl"
    getStdPlot toOrigPoly "#d62728" |> Chart.withTraceInfo "Poly"
    getStdPlot toOrigCubic "grey" |> Chart.withTraceInfo "Cubic"
    ]
    |> Chart.combine
    |> Chart.withAxisTitles "log var@TP" "stDev of Distances to Orig fit"
    |> Chart.saveHtml (storageOrig + "standardDeviationOfDistances",true)


    [|tomeanTC;tomeanLin;tomeanSpl;tomeanPoly;tomeanCubic|],[|toOrigTC;toOrigLin;toOrigSpl;toOrigPoly;toOrigCubic|]

printfn "start VarAic %A" (System.DateTime.Now)
let singleRepVarAic         = summarizeSingleRepRemove Minimizer.AICc WeightingMethod.Variance
printfn "start VarPlainAic %A" (System.DateTime.Now)
let singleRepVarPlainAic    = summarizeSingleRepRemove Minimizer.AICc WeightingMethod.VarPlain
printfn "start CVAic %A" (System.DateTime.Now)
let singleRepCVAic          = summarizeSingleRepRemove Minimizer.AICc WeightingMethod.CV
printfn "start StDevAdjAic %A" (System.DateTime.Now)
let StandardDeviation    = summarizeSingleRepRemove Minimizer.AICc WeightingMethod.StandardDeviation
printfn "start StDevAic %A" (System.DateTime.Now)
let singleRepStDevAdjAic    = summarizeSingleRepRemove Minimizer.AICc WeightingMethod.StandardDeviationAdj
printfn "All combintations AIC of singleRepDeletion ready!"

printfn "start VarGcv %A" (System.DateTime.Now)
let singleRepVarGcv         = summarizeSingleRepRemove Minimizer.GCV WeightingMethod.Variance
printfn "start VarPlainGcv %A" (System.DateTime.Now)
let singleRepVarPlainGcv    = summarizeSingleRepRemove Minimizer.GCV WeightingMethod.VarPlain
printfn "start CVGcv %A" (System.DateTime.Now)
let singleRepCVGcv          = summarizeSingleRepRemove Minimizer.GCV WeightingMethod.CV
printfn "start StDevAdjGcv %A" (System.DateTime.Now)
let singleRepStDevAdjGcv    = summarizeSingleRepRemove Minimizer.GCV WeightingMethod.StandardDeviationAdj
printfn "start StDevGcv %A" (System.DateTime.Now)
let singleRepStDevGcv    = summarizeSingleRepRemove Minimizer.GCV WeightingMethod.StandardDeviation
printfn "All combintations GCV of singleRepDeletion ready!"



//let getBoxPlot (input:(float*float) list) colOutline colFill title yAxisTitle = 
//    input
//    |> List.sortBy fst
//    |> List.chunkBySize (input.Length / 9)
//    |> List.map (fun x -> 
//        let xs,ys = List.unzip x
//        let xd =  List.map log xs 
//        let xV = Array.init x.Length (fun _ -> Seq.mean xd)
//        Chart.BoxPlot(X=xV,Y=ys,Name=sprintf "%s %.3f" title (Seq.stDev ys),FillColor=Color.fromString colFill,OutlineColor=Color.fromString colOutline)
//        )
//    |> Chart.combine
//    |> Chart.withAxisTitlesYR "log var@TP" yAxisTitle (-0.2,0.2)
//[
//getBoxPlot tomeanTC  "#2274b7" "#91b6dc" "TC" "Distance to Mean"
//getBoxPlot tomeanLin "#ff7f0e" "#ffbe86" "Lin" "Distance to Mean"
//getBoxPlot tomeanSpl "#2ca02c" "#95cf95" "Spl" "Distance to Mean"
//getBoxPlot tomeanPoly "#d62728" "#ea9293" "Poly" "Distance to Mean"
//]
//|> Chart.Grid(4,1)
//|> Chart.withSize(900.,1300.)
//|> Chart.withTitle "remove single Rep@TP, variance@TP to DistanceToMean"
//|> Chart.saveHtml (@"C:\Users\venn\source\repos\TemporalClassification\runs\Testing\remove\removeSingleRepatTP\box_toMean",true)
//[
//getBoxPlot toOrigTC  "#2274b7" "#91b6dc" "TC" "Distance to OrigFit"
//getBoxPlot toOrigLin "#ff7f0e" "#ffbe86" "Lin" "Distance to OrigFit"
//getBoxPlot toOrigSpl "#2ca02c" "#95cf95" "Spl" "Distance to OrigFit"
//getBoxPlot toOrigPoly "#d62728" "#ea9293" "Poly" "Distance to OrigFit"
//]
//|> Chart.Grid(4,1)
//|> Chart.withSize(900.,1300.)
//|> Chart.withTitle "remove single Rep@TP, variance@TP to DistanceToOrigFit"
//|> Chart.saveHtml (@"C:\Users\venn\source\repos\TemporalClassification\runs\Testing\remove\removeSingleRepatTP\box_toOrigFit",true)


//move to main function, as well above in remove timepoint!!















































//////testTimePoint 
let remove4TP =
    data
    |> Array.take 10
    |> PSeq.mapi (fun i (id,values) -> 
        if i%100=0 then printfn "%i/%i" i data.Length
        let timepoints = vector [0.;1.;2.;3.;4.;5.;6.;7.]
        let original = values |> Array.chunkBySize 3
        let originalFit,models = 
            Fitting.getBestFit timepoints original Fitting.WeightingMethod.StandardDeviation Fitting.Minimizer.AICc

        let removeTimepoints = vector [0.;1.;2.;3.;(*4.;*)5.;6.;7.]
        let removeDataPoint = removeTimepoints |> Seq.map (fun x -> original.[int x]) |> Array.ofSeq
        let removeFit,modelsRemove = 
            Fitting.getBestFit removeTimepoints removeDataPoint Fitting.WeightingMethod.StandardDeviation Fitting.Minimizer.AICc
    
        (original,originalFit),(removeDataPoint,removeFit)
    )
    |> PSeq.withDegreeOfParallelism 8
    |> PSeq.toArray



[|remove4TP.[2];remove4TP.[4];remove4TP.[7]|]
|> Array.mapi (fun i ((orig,origFit),(rem,remFit)) -> 
    let origx,origy = orig |> Array.mapi (fun i x -> x |> Array.map (fun k -> float i,k)) |> Array.concat |> Array.unzip
    let oMx,oMy = orig |> Array.mapi (fun i x -> float i,Seq.mean x) |> Array.unzip
    let orFit = [0. .. 0.1 .. 7.] |> List.map (fun x -> x,origFit.SplineFunction x) 

    let origlinSpline = Interpolation.LinearSpline.interpolate (Interpolation.LinearSpline.initInterpolate oMx oMy)
    let origcubic = Interpolation.CubicSpline.Simple.fit (Interpolation.CubicSpline.Simple.coefficients Interpolation.CubicSpline.Simple.BoundaryCondition.Natural (vector oMx) (vector oMy)) (vector oMx)
    let origpolyInt = Interpolation.Polynomial.fit (Interpolation.Polynomial.coefficients (vector oMx) (vector oMy))

    let t = [0;1;2;3;5;6;7]
    let remx,remy  = rem |> Array.mapi (fun i x -> x |> Array.map (fun k->  float t.[i],k)) |> Array.concat |> Array.unzip
    let rMx,rMy  = rem |> Array.mapi (fun i x -> float t.[i],Seq.mean x) |> Array.unzip
    let reFit = [0. .. 0.1 .. 7.] |> List.map (fun x -> x,remFit.SplineFunction x)
    
    let remlinSpline = Interpolation.LinearSpline.interpolate (Interpolation.LinearSpline.initInterpolate rMx rMy) 
    let remcubic = Interpolation.CubicSpline.Simple.fit (Interpolation.CubicSpline.Simple.coefficients Interpolation.CubicSpline.Simple.BoundaryCondition.Natural (vector rMx) (vector rMy)) (vector rMx)
    let rempolyInt = Interpolation.Polynomial.fit (Interpolation.Polynomial.coefficients (vector rMx) (vector rMy))

    let difSpline = Math.Abs (remlinSpline 4. - origlinSpline 4.)
    let difLinear = Math.Abs (remcubic  4. - origcubic 4.)
    let difPoly = Math.Abs (rempolyInt  4. - origpolyInt 4.)
    let difTC = Math.Abs ((remFit.SplineFunction 4.) - (origFit.SplineFunction 4.))

    [
        Chart.Point((origx |> Array.map (fun x -> x + 1.)),origy,MarkerColor =Color.fromString "red") |> Chart.withMarkerStyle(Size=5,Color=Color.fromString "red")
        Chart.Line((orFit|> Seq.map (fun (x,y) -> x + 1.,y)),LineColor =Color.fromString "#1f77b4") |> Chart.withTraceInfo "origTC"
        [0. .. 0.1 .. 7.] |> List.map (fun x -> x + 1.,origlinSpline x) |> fun t -> Chart.Line (t,LineColor =Color.fromString "#1f77b4",Name="linSpl")
        [0. .. 0.1 .. 7.] |> List.map (fun x -> x + 1.,origcubic x) |> fun t -> Chart.Line (t,LineColor =Color.fromString "#1f77b4",Name="cubic")
        [0. .. 0.1 .. 7.] |> List.map (fun x -> x + 1.,origpolyInt x) |> fun t -> Chart.Line (t,LineColor =Color.fromString "#1f77b4",Name="poly")

        Chart.Point((remx |> Array.map (fun x -> x + 1.)),remy,MarkerColor =Color.fromString"black") |> Chart.withMarkerStyle(Size=5,Color=Color.fromString "black")
        Chart.Line((reFit|> Seq.map (fun (x,y) -> x + 1.,y)),LineColor =Color.fromString "#ff7f0e") |> Chart.withTraceInfo "remTC"
        [0. .. 0.1 .. 7.] |> List.map (fun x -> x + 1.,remlinSpline x) |> fun t -> Chart.Line (t,LineColor =Color.fromString "#ff7f0e",Name="linSpl")
        [0. .. 0.1 .. 7.] |> List.map (fun x -> x + 1.,remcubic x) |> fun t -> Chart.Line (t,LineColor =Color.fromString "#ff7f0e",Name="cubic")
        [0. .. 0.1 .. 7.] |> List.map (fun x -> x + 1.,rempolyInt x) |> fun t -> Chart.Line (t,LineColor =Color.fromString "#ff7f0e",Name="poly")
        

    ]
    |> Chart.combine
    |> Chart.withAxisTitles "" ""

    |> Chart.withDescription 
        [
        Giraffe.ViewEngine.HtmlElements.rawText (sprintf "poly: %f<br>linSpl: %f<br>lin: %f<br>TC: %f" difPoly difSpline difLinear difTC)
        
        ]
    |> Chart.withTitle (string i )
    |> Chart.withLayoutStyle(Font=(Font.init(Family=FontFamily.Arial,Size=14)))
    |> Chart.withSize(500.,500)
    |> Chart.show
    )



[|remove4TP.[7]|]
|> Array.mapi (fun i ((orig,origFit),(rem,remFit)) -> 
    let origx,origy = orig |> Array.mapi (fun i x -> x |> Array.map (fun k -> float i,k)) |> Array.concat |> Array.unzip
    let oMx,oMy = orig |> Array.mapi (fun i x -> float i,Seq.mean x) |> Array.unzip
    let orFit = [0. .. 0.1 .. 7.] |> List.map (fun x -> x,origFit.SplineFunction x) 

    let origlinSpline = Interpolation.LinearSpline.interpolate (Interpolation.LinearSpline.initInterpolate oMx oMy)
    let origcubic = Interpolation.CubicSpline.Simple.fit (Interpolation.CubicSpline.Simple.coefficients Interpolation.CubicSpline.Simple.BoundaryCondition.Natural (vector oMx) (vector oMy)) (vector oMx)
    let origpolyInt = Interpolation.Polynomial.fit (Interpolation.Polynomial.coefficients (vector oMx) (vector oMy))

    let t = [0;1;2;3;5;6;7]
    let remx,remy  = rem |> Array.mapi (fun i x -> x |> Array.map (fun k->  float t.[i],k)) |> Array.concat |> Array.unzip
    let rMx,rMy  = rem |> Array.mapi (fun i x -> float t.[i],Seq.mean x) |> Array.unzip
    let reFit = [0. .. 0.1 .. 7.] |> List.map (fun x -> x,remFit.SplineFunction x)
    
    let remlinSpline = Interpolation.LinearSpline.interpolate (Interpolation.LinearSpline.initInterpolate rMx rMy) 
    let remcubic = Interpolation.CubicSpline.Simple.fit (Interpolation.CubicSpline.Simple.coefficients Interpolation.CubicSpline.Simple.BoundaryCondition.Natural (vector rMx) (vector rMy)) (vector rMx)
    let rempolyInt = Interpolation.Polynomial.fit (Interpolation.Polynomial.coefficients (vector rMx) (vector rMy))

    let difSpline = Math.Abs (remlinSpline 4. - origlinSpline 4.)
    let difLinear = Math.Abs (remcubic  4. - origcubic 4.)
    let difPoly = Math.Abs (rempolyInt  4. - origpolyInt 4.)
    let difTC = Math.Abs ((remFit.SplineFunction 4.) - (origFit.SplineFunction 4.))
    
    let chartTC =
        [
            Chart.Point((origx |> Array.map (fun x -> x + 1.)),origy,MarkerColor =Color.fromString "red") |> Chart.withMarkerStyle(Size=5,Color=Color.fromString "red")
            Chart.Point((remx |> Array.map (fun x -> x + 1.)),remy,MarkerColor =Color.fromString"black") |> Chart.withMarkerStyle(Size=5,Color=Color.fromString "black")
            Chart.Line((orFit|> Seq.map (fun (x,y) -> x + 1.,y)),LineColor =Color.fromString "#1f77b4") |> Chart.withTraceInfo "origTC"
            Chart.Line((reFit|> Seq.map (fun (x,y) -> x + 1.,y)),LineColor =Color.fromString "#ff7f0e") |> Chart.withTraceInfo "remTC"
        ] |> Chart.combine |> Chart.withAxisTitles "" ""
        
    let chartLS =
        [
            Chart.Point((origx |> Array.map (fun x -> x + 1.)),origy,MarkerColor =Color.fromString "red") |> Chart.withMarkerStyle(Size=5,Color=Color.fromString "red")
            Chart.Point((remx |> Array.map (fun x -> x + 1.)),remy,MarkerColor =Color.fromString"black") |> Chart.withMarkerStyle(Size=5,Color=Color.fromString "black")
            [0. .. 0.1 .. 7.] |> List.map (fun x -> x + 1.,origlinSpline x) |> fun t -> Chart.Line (t,LineColor =Color.fromString "#1f77b4",Name="linSpl")
            [0. .. 0.1 .. 7.] |> List.map (fun x -> x + 1.,remlinSpline x) |> fun t -> Chart.Line (t,LineColor =Color.fromString "#ff7f0e",Name="linSpl")
        ] |> Chart.combine |> Chart.withAxisTitles "" ""
        
    let chartCS =
        [
            Chart.Point((origx |> Array.map (fun x -> x + 1.)),origy,MarkerColor =Color.fromString "red") |> Chart.withMarkerStyle(Size=5,Color=Color.fromString "red")
            Chart.Point((remx |> Array.map (fun x -> x + 1.)),remy,MarkerColor =Color.fromString"black") |> Chart.withMarkerStyle(Size=5,Color=Color.fromString "black")
            [0. .. 0.1 .. 7.] |> List.map (fun x -> x + 1.,origcubic x) |> fun t -> Chart.Line (t,LineColor =Color.fromString "#1f77b4",Name="cubic")
            [0. .. 0.1 .. 7.] |> List.map (fun x -> x + 1.,remcubic x) |> fun t -> Chart.Line (t,LineColor =Color.fromString "#ff7f0e",Name="cubic")
        ] |> Chart.combine |> Chart.withAxisTitles "" ""
        
    let chartPO =
        [
            Chart.Point((origx |> Array.map (fun x -> x + 1.)),origy,MarkerColor =Color.fromString "red") |> Chart.withMarkerStyle(Size=5,Color=Color.fromString "red")
            Chart.Point((remx |> Array.map (fun x -> x + 1.)),remy,MarkerColor =Color.fromString"black") |> Chart.withMarkerStyle(Size=5,Color=Color.fromString "black")
            [0. .. 0.1 .. 7.] |> List.map (fun x -> x + 1.,origpolyInt x) |> fun t -> Chart.Line (t,LineColor =Color.fromString "#1f77b4",Name="poly")
            [0. .. 0.1 .. 7.] |> List.map (fun x -> x + 1.,rempolyInt x) |> fun t -> Chart.Line (t,LineColor =Color.fromString "#ff7f0e",Name="poly")
        ] |> Chart.combine |> Chart.withAxisTitles "" ""

    [
        chartTC |> Chart.withXAxis(LinearAxis.init(ShowTickLabels=false))
        chartPO |> Chart.withYAxis(LinearAxis.init(ShowTickLabels=false)) |> Chart.withXAxis(LinearAxis.init(ShowTickLabels=false))
        chartLS
        chartCS |> Chart.withYAxis(LinearAxis.init(ShowTickLabels=false))
    ]
    |> Chart.Grid(2,2,YGap=(0.0),XGap=(0.0))
    |> Chart.withLayoutStyle(Font=(Font.init(Family=FontFamily.Arial,Size=14)))
    |> Chart.withSize(700.,500)
    |> Chart.show
    )



remove4TP
|> Array.map (fun ((orig,origFit),(rem,remFit)) -> 
    let origx,origy = orig |> Array.mapi (fun i x -> x |> Array.map (fun k -> float i,k)) |> Array.concat |> Array.unzip
    let oMx,oMy = orig |> Array.mapi (fun i x -> float i,Seq.mean x) |> Array.unzip
    let orFit = [0. .. 0.1 .. 7.] |> List.map (fun x -> x,origFit.SplineFunction x) 

    let origlinSpline = Interpolation.LinearSpline.interpolate (Interpolation.LinearSpline.initInterpolate oMx oMy)
    let origlinear = Linear.Univariable.fit (Linear.Univariable.coefficient (vector oMx) (vector oMy))
    let origpolyInt = Interpolation.Polynomial.fit (Interpolation.Polynomial.coefficients (vector oMx) (vector oMy))

    let t = [0;1;2;3;5;6;7]
    let remx,remy  = rem |> Array.mapi (fun i x -> x |> Array.map (fun k->  float t.[i],k)) |> Array.concat |> Array.unzip
    let rMx,rMy  = rem |> Array.mapi (fun i x -> float t.[i],Seq.mean x) |> Array.unzip
    let reFit = [0. .. 0.1 .. 7.] |> List.map (fun x -> x,remFit.SplineFunction x)
    
    let remlinSpline = Interpolation.LinearSpline.interpolate (Interpolation.LinearSpline.initInterpolate rMx rMy) 
    let remlinear = Linear.Univariable.fit (Linear.Univariable.coefficient (vector rMx) (vector rMy)) 
    let rempolyInt = Interpolation.Polynomial.fit (Interpolation.Polynomial.coefficients (vector rMx) (vector rMy))

    let difSpline = Math.Abs (remlinSpline 4. - origlinSpline 4.)
    let difLinear = Math.Abs (remlinear  4. - origlinear 4.)
    let difPoly = Math.Abs (rempolyInt  4. - origpolyInt 4.)
    let difTC = Math.Abs ((remFit.SplineFunction 4.) - (origFit.SplineFunction 4.))

    [|difPoly;difSpline;difLinear;difTC|]
    )
|> JaggedArray.transpose
|> Array.map Chart.Histogram
|> Chart.Grid(4,1)
|> Chart.show




let plotFourFitsInParallel i xaxisText = 
    
    remove4TP.[i]
    |> (fun ((orig,origFit),(rem,remFit)) -> 
        let origx,origy = orig |> Array.mapi (fun i x -> x |> Array.map (fun k -> float i,k)) |> Array.concat |> Array.unzip
        let oMx,oMy = orig |> Array.mapi (fun i x -> float i,Seq.mean x) |> Array.unzip
        let orFit = [0. .. 0.1 .. 7.] |> List.map (fun x -> x,origFit.SplineFunction x) 

        let origlinSpline = Interpolation.LinearSpline.interpolate (Interpolation.LinearSpline.initInterpolate oMx oMy)
        let origcubic = Interpolation.CubicSpline.Simple.fit (Interpolation.CubicSpline.Simple.coefficients Interpolation.CubicSpline.Simple.BoundaryCondition.Natural (vector oMx) (vector oMy)) (vector oMx)
        let origpolyInt = Interpolation.Polynomial.fit (Interpolation.Polynomial.coefficients (vector oMx) (vector oMy))

        let t = [0;1;2;3;5;6;7]
        let remx,remy  = rem |> Array.mapi (fun i x -> x |> Array.map (fun k->  float t.[i],k)) |> Array.concat |> Array.unzip
        let rMx,rMy  = rem |> Array.mapi (fun i x -> float t.[i],Seq.mean x) |> Array.unzip
        let reFit = [0. .. 0.1 .. 7.] |> List.map (fun x -> x,remFit.SplineFunction x)
    
        let remlinSpline = Interpolation.LinearSpline.interpolate (Interpolation.LinearSpline.initInterpolate rMx rMy) 
        let remcubic = Interpolation.CubicSpline.Simple.fit (Interpolation.CubicSpline.Simple.coefficients Interpolation.CubicSpline.Simple.BoundaryCondition.Natural (vector rMx) (vector rMy)) (vector rMx)
        let rempolyInt = Interpolation.Polynomial.fit (Interpolation.Polynomial.coefficients (vector rMx) (vector rMy))

        let difSpline = Math.Abs (remlinSpline 4. - origlinSpline 4.)
        let difLinear = Math.Abs (remcubic  4. - origcubic 4.)
        let difPoly = Math.Abs (rempolyInt  4. - origpolyInt 4.)
        let difTC = Math.Abs ((remFit.SplineFunction 4.) - (origFit.SplineFunction 4.))

        let markersize = 4
        let myAxisWoD name = 
            LinearAxis.init(
                Title=Title.init name,
                Mirror=StyleParam.Mirror.All,
                Ticks=StyleParam.TickOptions.Inside,
                ShowGrid=false,
                ShowLine=true,
                TickVals = ("" |> Set.ofSeq),
                TickFont=Font.init(Size=12,Family=StyleParam.FontFamily.Arial)
            )

        let min,max = 
            origy 
            |> fun x -> 
                Seq.min x * 0.99,Seq.max x *1.01
            
        [
            [
            Chart.Point((origx |> Array.map (fun x -> x + 1.)),origy,ShowLegend=false) |> Chart.withMarkerStyle(Color=Color.fromString "black",Size=markersize)
            Chart.Line((orFit |> List.map (fun (x,y) -> x + 1.,y)),LineColor =Color.fromString "blue",ShowLegend=false) |> Chart.withTraceInfo "constrained spline"
            ] |> Chart.combine |> Chart.withAxisTitles (if xaxisText then "constrained spline" else "") "protein intensity" |> Chart.withYAxisStyle(MinMax=(min,max))
            [
            Chart.Point((origx |> Array.map (fun x -> x + 1.)),origy,ShowLegend=false) |> Chart.withMarkerStyle(Color=Color.fromString "black",Size=markersize)
            [0. .. 0.1 .. 7.] |> List.map (fun x -> x + 1.,origpolyInt x) |> fun t -> Chart.Line (t,LineColor =Color.fromString"blue",ShowLegend=false,Name = "polynomial interpolant")
            ] |> Chart.combine |> Chart.withAxisTitles (if xaxisText then "polynomial interpolant" else "") "" |> Chart.withYAxis (myAxisWoD "") |> Chart.withYAxisStyle(MinMax=(min,max))
            [
            Chart.Point((origx |> Array.map (fun x -> x + 1.)),origy,ShowLegend=false) |> Chart.withMarkerStyle(Color=Color.fromString "black",Size=markersize)
            [0. .. 0.1 .. 7.] |> List.map (fun x -> x + 1.,origlinSpline x) |> fun t -> Chart.Line (t,LineColor =Color.fromString"blue",ShowLegend=false,Name = "linear spline")
            ] |> Chart.combine |> Chart.withAxisTitles (if xaxisText then "linear spline" else "") "" |> Chart.withYAxis (myAxisWoD "")|> Chart.withYAxisStyle(MinMax=(min,max))
            [
            Chart.Point((origx |> Array.map (fun x -> x + 1.)),origy,ShowLegend=false) |> Chart.withMarkerStyle(Color=Color.fromString "black",Size=markersize)
            [0. .. 0.1 .. 7.] |> List.map (fun x -> x + 1.,origcubic x) |> fun t -> Chart.Line (t,LineColor =Color.fromString"blue",ShowLegend=false,Name = "cubic spline")
            ] |> Chart.combine |> Chart.withAxisTitles (if xaxisText then "cubic spline" else "") "" |> Chart.withYAxis (myAxisWoD "")|> Chart.withYAxisStyle(MinMax=(min,max))
        ] 
        )





[
plotFourFitsInParallel 1 false
plotFourFitsInParallel 21 false
plotFourFitsInParallel 2 false
plotFourFitsInParallel 11 false
plotFourFitsInParallel 20 true
]
|> List.concat
|> Chart.Grid(5,4,YGap = 0.1,XGap = 0.1)
|> Chart.withSize(800,1000)
|> Chart.withLayoutStyle(Font=(Font.init(Family=FontFamily.Arial,Size=14)))
|> Chart.show

(remove4TP.[1] |> fst |> snd)

[1;21;2;11;20] |> List.map (fun i -> fst data.[i])