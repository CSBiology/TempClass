namespace TempClass

open Microsoft.SolverFoundation
open System
open System.Collections.Generic
open TempClass.Optimization
open FSharp.Stats

module TemporalClassification =    

    open FSharp.Stats
    open FSharp.Stats.Optimization

    module Aux =

        /// <summary>Takes a collection of quadrupels and sorts them in to four individual arrays</summary>
        /// <param name="input">collection of values to be separated</param>
        /// <returns>A tuple ('a [] * 'b [] * 'c [] * 'd [])</returns>
        /// <example> 
        /// <code> 
        /// 
        /// let xData = [|(1,"a",1.,'a');(2,"b",2.,'b');(3,"c",3.,'c');|]
        /// 
        /// unzip4 xData
        /// </code> 
        /// </example>
        let unzip4 input =
            let length = Seq.length input
            let la = Array.zeroCreate length
            let lb = Array.zeroCreate length
            let lc = Array.zeroCreate length
            let ld = Array.zeroCreate length
            input 
            |> Seq.iteri (fun i (a,b,c,d) ->   
                (la.[i] <- a)
                (lb.[i] <- b)
                (lc.[i] <- c)
                (ld.[i] <- d)
                )
            (la,lb,lc,ld)

        /// <summary>Takes a collection of quadrupels and sorts them in to four individual arrays</summary>
        /// <param name="xVal">a single x value that should be rounded to the neighboring value within xVals</param>
        /// <param name="xVals">collection of x values</param>
        /// <returns>The nearest x value</returns>
        /// <example> 
        /// <code> 
        /// 
        /// let xData = [1.;2.;5.;6.;]
        /// 
        /// roundToNext  3.4 xData
        /// // returns 5.
        /// </code> 
        /// </example>
        let roundToNext (xVal:float) (xVals:seq<float>) =
            xVals 
            |> Seq.minBy (fun xI -> Math.Abs (xI - xVal))

        /// <summary>Gets indices of minima and maxima by brute force search.</summary>
        /// <param name="pct">Allowed percentage deviation (+/-) for a value to be considered as equal</param>
        /// <param name="arr">collection of values to investigate</param>
        /// <returns>Lists of minima indices and maxima indices.</returns>
        /// <example> 
        /// <code> 
        /// 
        /// let xData = [|1.;2.;5.;3.;|]
        /// 
        /// investigateTrace 0.025 xData
        /// // returns [],[2], indicating a maximum at index 2 
        /// </code> 
        /// </example>
        let investigateTrace pct (arr: float[]) =
            let almostEqual (pct:float) f1 f2 =
                let (min,max) = f2 - Math.Abs(f2 * pct),f2 + Math.Abs(f2 * pct)
                f1 >= min && f1 <= max
   
            let defGreaterThan (pct:float) f1 f2 =
                if f1 > 0. then 
                    let max = f2 * (1.+pct)
                    f1 > max
                else 
                    let max = f2 * (1.-pct) 
                    f1 > max

            let defLowerThan (pct:float) f1 f2 =
                if f1 > 0. then 
                    let min = f2 * (1.-pct)  //2. * GradientDescent.eps
                    f1 < min
                else 
                    let min = f2 * (1.+pct) 
                    f1 < min
            let length = arr.Length - 1
            let rec loop x accMin accMax=
                if x < length && x > 0 then
                    let tmp = arr.[x]
                    let prev = arr.[x-1] 
                    let next = arr.[x+1]
                    if not (almostEqual pct tmp prev) && not (almostEqual pct tmp prev) then 
                        let globalprevMax = arr.[0..x-1]|> Array.max |> defLowerThan pct tmp
                        let globalprevMin = arr.[0..x-1]|> Array.min |> defGreaterThan pct tmp 
                        let globalnextMax = arr.[x+1..] |> Array.max |> defLowerThan pct tmp 
                        let globalnextMin = arr.[x+1..] |> Array.min |> defGreaterThan pct tmp 
                        ////simple max
                        if defGreaterThan pct tmp prev && defGreaterThan pct tmp next then 
                            loop (x+1) accMin (x::accMax)
                        else
                            //simple min
                            if defLowerThan pct tmp prev && defLowerThan pct tmp next then 
                                loop (x+1) (x::accMin) accMax
                            else
                                //comp min
                                if defLowerThan pct tmp prev && almostEqual pct tmp next && globalnextMax && globalprevMax && not globalnextMin  then 
                                    loop (x+1) (x::accMin) accMax
                                else
                                //comp max
                                    if defGreaterThan pct tmp prev && almostEqual pct tmp next && globalnextMin && globalprevMin && not globalnextMax then 
                                        loop (x+1) accMin (x::accMax) 
                                    else loop (x+1) accMin accMax
                    else loop (x+1) accMin accMax
                else 
                    accMin,accMax
            loop 1 [] []

    module Fitting =
        
        type Extremum = {
            /// 1:Maximum; -1:Minimum
            Indicator : int
            /// xValue of extremum
            Xvalue    : float
            } with
                static member Create indicator xValue = {Indicator = indicator; Xvalue = xValue}

        /// <summary>Converts a collection of `Extremum` to a String in the form of Min2.00,Max1.00.</summary>
        /// <param name="extrema">Extrema to be converted</param>
        /// <returns>A single string describing the extrema configuration</returns>
        let extremaToString (extrema: Extremum list) =
            let indicatorToString ind =
                if ind = 1 then "Max" else "Min"
            extrema
            |> List.map (fun x -> sprintf "%s%.2f" (indicatorToString x.Indicator) x.Xvalue)
            |> String.concat ","

        /// <summary>Parent shape descriptor (e.g. In2 for 2 extrema starting increasing (max,min))</summary>
        type Condition =
            /// monotonically increasing
            | In0
            /// 1 maximum
            | In1
            /// maximum, then minimum
            | In2
            /// maximum, then minimum, then maximum
            | In3
            /// maximum, then minimum, then maximum, then minimum
            | In4
            /// monotonically decreasing
            | De0
            /// 1 minimum
            | De1
            /// minimum, then maximum
            | De2
            /// minimum, then maximum, then minimum
            | De3
            /// minimum, then maximum, then minimum, then maximum
            | De4
            /// contains more than 4 extrema
            | Complex //more complex
                       
        /// <summary>Parent shape descriptor (e.g. In2 for 2 extrema starting increasing (max,min))</summary>
        type WeightingMethod =
            /// weight: 1
            | Equal
            /// weight: minmax((1/stDev),Wmin)
            | VarRobust//Variance
            /// 1/var centered around 1
            | Variance//Variance
            /// 1/var not centered
            | VarPlain
            /// max 0. (log(1./|Seq.cvPopulation g|),2.))
            | CV
            /// weight: (1/stDev)
            | StandardDeviation
            /// weight: sqrt(1/stDev)
            | StandardDeviationSqrt
            /// weight: sqrt(1/(stDev/mean))
            | StandardDeviationAdj
            /// weight: sqrt(sqrt(1/(stDev/mean)))
            | StandardDeviationAdjSqrt
            
        /// <summary>Parent shape descriptor (e.g. In2 for 2 extrema starting increasing (max,min))</summary>
        type Minimizer = 
            /// modified generalized cross validation (rho 1.3)
            | GCV
            /// Akaike information criterion corrected for small sample sizes
            | AICc

        ///
        type TempClassResult = {
            XValues : vector option
            YValues : float [][] option
            WeightingMethod : WeightingMethod option
            Minimizer : Minimizer option
            ///spline function values at knots
            TraceA   : vector
            ///spline second derivative values at knots
            TraceC   : vector
            ///weighted error of a and y
            Error    : float
            ///mGCV
            GCV      : float
            ///used smoothing factor lambda
            Lambda   : float
            ///used constraint matrix for shape determination
            Ctemp    : float [,]
            ///AIC determined by SSE and the number of used extrema
            AICc     : float
            ///function, that returns the splines' value at a given x_value
            SplineFunction  : (float -> float)
            
            GCVArray: (float*float)[] option
            } with 
                static member Create x y w m a c e gcv lambda ctemp aic splinefunction gcvArray= {
                    XValues         = x
                    YValues         = y
                    WeightingMethod = w
                    Minimizer       = m
                    TraceA          = a
                    TraceC          = c
                    Error           = e
                    GCV             = gcv
                    Lambda          = lambda
                    Ctemp           = ctemp
                    AICc            = aic
                    SplineFunction  = splinefunction
                    GCVArray        = gcvArray
                    } 

        /// record type for easy handling inside the shape classification process
        type InnerResult = {
            TraceA  : Vector<float>
            GCV     : float
            Var     : float
            Error   : float
            CTemp   : float [,]
            Lambda  : float
            } with
                static member Create a g c e v l = {
                    TraceA  = a
                    GCV     = g
                    Var     = v
                    Error   = e
                    CTemp   = c
                    Lambda  = l}    


        /// <summary>Takes x values, y values and curvature of function at x values and returns a function that maps x values to its corresponding y value.</summary>
        /// <param name="x">vector of x values</param>
        /// <param name="a">vector of y values at the position of the x values</param>
        /// <param name="c">vector of curvature (second derivative) at the position of the x values</param>
        /// <returns>gives smoothing spline function</returns>
        /// <example> 
        /// <code> 
        /// 
        /// let func = initEvalAt (vector [1.;2.;3.;4.]) (vector [1.;2.;3.;2.]) (vector [0.1;0.02;0.11;0.2])
        /// 
        /// func 1.3
        /// // returns the position of the defined spline at position 1.3
        /// </code> 
        /// </example>
        let initEvalAt (x:Vector<float>) (a:Vector<float>) (c:Vector<float>) =
            let n = x.Length

            //define interval stepwidths
            let h = Vector.init (n-1) (fun i -> x.[i+1] - x.[i] )
            // helper functions f(i,1-4)(t) for smoothing spline
            let calcF1 (h:Vector<float>) (x:Vector<float>) (i:int) (t:float) =
                (x.[i+1] - t) / h.[i]

            let calcF2 (h:Vector<float>) (x:Vector<float>) (i:int) (t:float) =
                (t - x.[i]) / h.[i]

            let calcF3 (h:Vector<float>) (x:Vector<float>) (i:int) (t:float) =
                ((calcF1 h x i t)**3.0 - (calcF1 h x i t) ) * (h.[i]**2.0) / 6.0

            let calcF4 (h:Vector<float>) (x:Vector<float>) (i:int) (t:float) =
                ((calcF2 h x i t)**3.0 - (calcF2 h x i t) ) * (h.[i]**2.0) / 6.0

            // helper function s(i) for smoothing spline 
            let calcS (h:Vector<float>) (x:Vector<float>) (a:Vector<float>) (c:Vector<float>) (i:int) (t:float) =
                a.[i] * (calcF1 h x i t) + a.[i+1] * (calcF2 h x i t) + c.[i] * (calcF3 h x i t) + c.[i+1] * (calcF4 h x i t)

            //function: a * ((y - t) / h) + b * (t - x) / h + c * ((((y - t) / h)^3.0 - ((y - t) / h) ) * (h^2.0) / 6.0  )+ d * ((((t - x) / h)^3.0 - ((t - x) / h) ) * (h^2.0) / 6.0) where b = a.[i+1] ...
            let evalAt =  calcS h x a c 
            
            (fun t ->
                let i = 
                    match Array.tryFindIndexBack (fun xs -> xs <= t) (x.InternalValues) with 
                    | Some z -> if t = Seq.last x then h.Length-1 else z //the last supported x Value can be processed by the f_i-1
                    | None   -> failwith "The x-value is not part of the section used to estimate the spline coefficients, thus a monotone function progression can not be guaranteed"
                evalAt i t
                )

        /// <summary>Takes x values, y values and curvature of function at x values and returns a function that maps x values to its corresponding y value slope.</summary>
        /// <param name="x">vector of x values</param>
        /// <param name="a">vector of y values at the position of the x values</param>
        /// <param name="c">vector of curvature (second derivative) at the position of the x values</param>
        /// <returns>gives function of function slope</returns>
        /// <example> 
        /// <code> 
        /// 
        /// let func = initEvalFstDerivativeAt (vector [1.;2.;3.;4.]) (vector [1.;2.;3.;2.]) (vector [0.1;0.02;0.11;0.2])
        /// 
        /// func 1.3
        /// // returns the slope of the defined spline at position 1.3
        /// </code> 
        /// </example>
        let initEvalFstDerivativeAt (x:Vector<float>) (a:Vector<float>) (c:Vector<float>) =
            let n = x.Length

            //define interval stepwidths
            let h = Vector.init (n-1) (fun i -> x.[i+1] - x.[i] )

            let calcF1 (h:Vector<float>) (x:Vector<float>) (i:int) (t:float) =
                (pown h.[i] 2) * (3. * (pown (t-x.[i]) 2) / (pown h.[i] 3) - (1. / h.[i]))
                |> fun x -> x / 6.

            let calcF2 (h:Vector<float>) (x:Vector<float>) (i:int) (t:float) =
                (pown h.[i] 2) * ((1./h.[i]) - ((3. * (pown (x.[i+1] - t) 2)/pown h.[i] 3))) 
                |> fun x -> x / 6.

            let calcF3 (h:Vector<float>) (x:Vector<float>) (i:int) (t:float) =
                1. / h.[i]

            let calcF4 (h:Vector<float>) (x:Vector<float>) (i:int) (t:float) =
                1. / h.[i]

            // helper function for fst derivative
            let calcS (h:Vector<float>) (x:Vector<float>) (a:Vector<float>) (c:Vector<float>) (i:int) (t:float) =
                c.[i+1] * (calcF1 h x i t) + c.[i] * (calcF2 h x i t) + a.[i+1] * (calcF3 h x i t) - a.[i] * (calcF4 h x i t)

            let evalAt =  calcS h x a c 
        
            (fun t ->
                let i = 
                    match Array.tryFindIndexBack (fun xs -> xs <= t) (x.InternalValues) with 
                    | Some z -> if t = Seq.last x then h.Length-1 else z 
                    | None   -> failwith "The x-value is not part of the section used to estimate the spline coefficients, thus a monotone function progression can not be guaranteed"
                evalAt i t
                )

        /// <summary>Takes x values, y values and curvature of function at x values and returns a function that maps x values to its corresponding y value curvature.</summary>
        /// <param name="x">vector of x values</param>
        /// <param name="a">vector of y values at the position of the x values</param>
        /// <param name="c">vector of curvature (second derivative) at the position of the x values</param>
        /// <returns>gives function of function curvature</returns>
        /// <example> 
        /// <code> 
        /// 
        /// let func = initEvalSndDerivativeAt (vector [1.;2.;3.;4.]) (vector [1.;2.;3.;2.]) (vector [0.1;0.02;0.11;0.2])
        /// 
        /// func 1.3
        /// // returns the curvature of the defined spline at position 1.3
        /// </code> 
        /// </example>
        let initEvalSndDerivativeAt (x:Vector<float>) (a:Vector<float>) (c:Vector<float>) =
            let n = x.Length

            //define interval stepwidths
            let h = Vector.init (n-1) (fun i -> x.[i+1] - x.[i] )

            let calcF1 (h:Vector<float>) (x:Vector<float>) (i:int) (t:float) =
                (t - x.[i]) / h.[i]

            let calcF2 (h:Vector<float>) (x:Vector<float>) (i:int) (t:float) =
                (x.[i+1] - t) / h.[i]

            // helper function s(i) for smoothing spline 
            let calcS (h:Vector<float>) (x:Vector<float>) (a:Vector<float>) (c:Vector<float>) (i:int) (t:float) =
                c.[i+1] * (calcF1 h x i t) + c.[i] * (calcF2 h x i t)

            let evalAt =  calcS h x a c 

            (fun t ->
                let i = 
                    match Array.tryFindIndexBack (fun xs -> xs <= t) (x.InternalValues) with 
                    | Some z -> if t = Seq.last x then h.Length-1 else z 
                    | None   -> failwith "The x-value is not part of the section used to estimate the spline coefficients, thus a monotone function progression can not be guaranteed"
                evalAt i t
                )

        /// <summary>Takes x values, y values and curvature of function at x values and returns a function that maps x values to its corresponding third derivative.</summary>
        /// <param name="x">vector of x values</param>
        /// <param name="a">vector of y values at the position of the x values</param>
        /// <param name="c">vector of curvature (second derivative) at the position of the x values</param>
        /// <returns>gives function of third derivative</returns>
        /// <example> 
        /// <code> 
        /// 
        /// let func = initEvalTrdDerivativeAt (vector [1.;2.;3.;4.]) (vector [1.;2.;3.;2.]) (vector [0.1;0.02;0.11;0.2])
        /// 
        /// func 1.3
        /// // returns third derivative of the defined spline at position 1.3
        /// </code> 
        /// </example>
        let initEvalTrdDerivativeAt (x:Vector<float>) (a:Vector<float>) (c:Vector<float>) =
            let n = x.Length

            //define interval stepwidths
            let h = Vector.init (n-1) (fun i -> x.[i+1] - x.[i] )

            // helper function s(i) for smoothing spline 
            let calcS (h:Vector<float>) (x:Vector<float>) (a:Vector<float>) (c:Vector<float>) (i:int) (t:float) =
                (c.[i+1] - c.[i]) / h.[i]

            let evalAt =  calcS h x a c 
            
            (fun t ->
                let i = 
                    match Array.tryFindIndexBack (fun xs -> xs <= t) (x.InternalValues) with 
                    | Some z -> if t = Seq.last x then h.Length-1 else z 
                    | None   -> failwith "The x-value is not part of the section used to estimate the spline coefficients, thus a monotone function progression can not be guaranteed"
                evalAt i t
                )

        /// <summary>Calculates extrema present in the given smoothing spline. Based on first and second
        /// derivative, the curve is analyzed for extrema.
        /// </summary>
        /// <param name="x">vector of x values</param>
        /// <param name="a">vector of y values at the position of the x values</param>
        /// <param name="c">vector of curvature (second derivative) at the position of the x values</param>
        /// <returns>Lists of extreme points</returns>
        let getExtrema (x:Vector<float>) (a:Vector<float>) (c:Vector<float>) =
            let n = x.Length 
            //define interval stepwidths
            let h = Vector.init (n-1) (fun i -> x.[i+1] - x.[i] )  
            // helper function for identification of roots in first derivative
            let calcS (h:Vector<float>) (x:Vector<float>) (a:Vector<float>) (c:Vector<float>) (i:int) (t:float) =      
                let f1 = (c.[i] * x.[i+1] - c.[i+1] * x.[i]) / (c.[i] - c.[i+1])
                let f2 = sqrt 2. / (c.[i] - c.[i+1])
                let root = 
                    let s1 = 0.5*c.[i]**2.*x.[i+1]**2.
                    let s2 = -c.[i]*c.[i+1]*x.[i+1]*x.[i]
                    let s3 = -1./6. * c.[i] * (6.*a.[i]-6.*a.[i+1]-c.[i]*h.[i]**2. + 3. *c.[i]*x.[i+1]**2. + c.[i+1]*h.[i]**2. - 3.* c.[i+1]*x.[i]**2.)
                    let s4 = 0.5 * c.[i+1]**2. * x.[i]**2.
                    let s5 = 1./6. * c.[i+1] * (6.*a.[i]-6.*a.[i+1]-c.[i]*h.[i]**2. + 3. *c.[i]*x.[i+1]**2. + c.[i+1]*h.[i]**2. - 3.* c.[i+1]*x.[i]**2.)
                    (s1 + s2 + s3 + s4 + s5)
                if root < 0. then 
                    []
                else 
                    [
                    f1 - f2 * sqrt root
                    f1 + f2 * sqrt root
                    ]
            let evalAt = calcS h x a c 
            // calculate amplitude of second derivative in interval i at xValue t
            let calcSndDeriv (i:int) (t:float) =
                c.[i+1] * (t - x.[i]) / h.[i] + c.[i] * (x.[i+1] - t) / h.[i]
            let sigBigger a b = a > (b+0.0000001)
            let splineFunction = initEvalAt x a c
            x.[0.. Seq.length x - 2]
            |> Seq.map 
                    (fun t ->
                        let i = 
                            match Array.tryFindIndexBack (fun xs -> xs <= t) (x.InternalValues) with 
                            | Some x -> x 
                            | None   -> failwith "The x-value is not part of the section used to estimate the spline coefficients, thus a monotone function progression can not be guaranteed"
                        //get an open interval of the current signal interval
                        let interval = Intervals.create (x.[i]-0.0000001) (x.[i+1]+0.0000001)
                        let secondDerivative t = calcSndDeriv i t
                        //calculate roots of fst derivative and check snd derivative
                        let extrema = 
                            evalAt i t
                            |> List.filter (fun x -> Intervals.liesInInterval x interval)
                            |> List.choose (fun xAtS0 -> 
                                let s1 = secondDerivative xAtS0
                                match s1 with
                                | w when round 2 s1 = 0. -> 
                                    if round 1 (xAtS0-0.05) <= Seq.head x  || round 1 (xAtS0+0.05) >= Seq.last x then 
                                        None
                                    else 
                                        let yAtX = splineFunction xAtS0
                                        let yAtpX = splineFunction (xAtS0 - 0.05)
                                        let yAtaX = splineFunction (xAtS0 + 0.05)
                                        if   sigBigger yAtX yAtpX && sigBigger yAtX yAtaX then Some (Extremum.Create  1 xAtS0) //(1.,xAtS0)
                                        elif sigBigger yAtpX yAtX && sigBigger yAtaX yAtX then Some (Extremum.Create -1 xAtS0) //(-1.,xAtS0)
                                        else None//real saddle point
                                | w when s1 < 0. -> Some (Extremum.Create 1 xAtS0) //(1.,xAtS0)
                                | _ -> Some (Extremum.Create -1 xAtS0)//(-1.,xAtS0)
                                )
                        extrema
                        //|> List.filter (fun (i,xv) -> i <> 0.)
                        )
                    |> List.concat
                    |> List.sortBy (fun ex -> ex.Xvalue)
                    |> List.map (fun ex -> {ex with Xvalue = Aux.roundToNext ex.Xvalue x})
                    |> List.distinct
        
        /// <summary>Determines extrema present in the given smoothing spline. Subsequently, the
        /// first interval behaviour (in- or decreasing) and the number of extrema is encoded (e.g. In2).
        /// </summary>
        /// <param name="x">vector of x values</param>
        /// <param name="a">vector of y values at the position of the x values</param>
        /// <param name="c">vector of curvature (second derivative) at the position of the x values</param>
        /// <returns>Condition, that encodes the parent shape configuration</returns>
        let getParentShape (x:Vector<float>) (a:Vector<float>) (c:Vector<float>) =
            let extrema = getExtrema x a c
            let extremaCount = extrema.Length
            let n = a.Length
            match extremaCount with
            | extremaCount when extremaCount = 0 -> 
                if a.[0] < a.[n-1] then In0
                else De0
            | extremaCount when extremaCount = 1 -> 
                if extrema.[0].Indicator = 1 then In1
                else De1
            | extremaCount when extremaCount = 2 -> 
                if extrema.[0].Indicator = 1 && extrema.[1].Indicator = -1 then In2
                else De2
            | extremaCount when extremaCount = 3 ->
                if extrema.[0].Indicator = 1 && extrema.[1].Indicator = -1 then In3
                else De3
            | extremaCount when extremaCount = 4 ->
                if extrema.[0].Indicator = 1 && extrema.[1].Indicator = -1 then In4
                else De4
            | _ -> Complex


        [<Obsolete("")>]
        //check the spline for the predefined condition
        let private checkshape (x:Vector<float>) (a:Vector<float>) (c:Vector<float>) (con:Condition)=
            let extrema = getExtrema x a c
            let extremaCount = extrema.Length
            let n = a.Length
            match con with
            | In0       -> extremaCount = 0 && a.[0] <= a.[n-1]
            | In1       -> extremaCount = 1 && extrema.[0].Indicator = 1
            | In2       -> extremaCount = 2 && extrema.[0].Indicator = 1 && extrema.[1].Indicator = -1
            | In3       -> extremaCount = 3 && extrema.[0].Indicator = 1
            | In4       -> extremaCount = 4 && extrema.[0].Indicator = 1 && extrema.[1].Indicator = -1
            | De0       -> extremaCount = 0 && a.[0] >= a.[n-1]
            | De1       -> extremaCount = 1 && extrema.[0].Indicator = -1
            | De2       -> extremaCount = 2 && extrema.[0].Indicator = -1 && extrema.[1].Indicator = 1
            | De3       -> extremaCount = 3 && extrema.[0].Indicator = -1
            | De4       -> extremaCount = 4 && extrema.[0].Indicator = -1 && extrema.[1].Indicator = 1
            | Complex   -> true

        [<Obsolete("")>]
        let private mapCondition (operator:Matrix<float> -> matrix) con =
            let mat = matrix [[1.]]
            match con with
            | x when con = 0 -> if mat = operator mat then Condition.De0 else Condition.In0
            | x when con = 1 -> if mat = operator mat then Condition.De1 else Condition.In1
            | x when con = 2 -> if mat = operator mat then Condition.De2 else Condition.In2
            | x when con = 3 -> if mat = operator mat then Condition.De3 else Condition.In3
            | x when con = 4 -> if mat = operator mat then Condition.De4 else Condition.In4
            | _ -> failwithf "con>5 not supported"

        /// <summary>Generates a weighting matrix out of the x-, and y-Values, and 
        /// the given weighting method.
        /// </summary>
        /// <param name="xVal">vector of x values</param>
        /// <param name="yVal">collection of y value samples</param>
        /// <param name="w_method">Method to determine the weighting (e.g. 1/var)</param>
        /// <returns>Matrix with weights encoded in its diagonal</returns>
        let getWeighting (xVal:Vector<float>) (yVal: float [] []) (w_method:WeightingMethod) =
            match w_method with
            | Equal ->
                Matrix.diag(Vector.oneCreate xVal.Length)
            //1/var centered around 1
            | Variance ->
                let Wtemp = Matrix.create xVal.Length xVal.Length 0.            
                let diagVec = 
                    //vector [for i = 0 to xVal.Length - 1 do yield 1. / (yVal.[(i*numRep)..(i*numRep + numRep - 1)] |> Seq.var)]
                    yVal 
                    |> Array.map (fun x -> 1. / max 0.00001 (Seq.var x))
                    |> vector

                let avg = Seq.average diagVec
                for i = 0 to xVal.Length - 1 do Wtemp.[i,i] <- diagVec.[i] / avg
                Wtemp
            //just 1/var not centered
            | VarPlain ->
                let Wtemp = Matrix.create xVal.Length xVal.Length 0.            
                let diagVec = 
                    //vector [for i = 0 to xVal.Length - 1 do yield 1. / (yVal.[(i*numRep)..(i*numRep + numRep - 1)] |> Seq.var)]
                    yVal 
                    |> Array.map (fun x -> 1. / max 0.00001 (Seq.var x))
                    |> vector
                for i = 0 to xVal.Length - 1 do Wtemp.[i,i] <- diagVec.[i]
                Wtemp
            | VarRobust -> //Variance
                let Wtemp = Matrix.create xVal.Length xVal.Length 0.            
                let diagVec = 
                    //if yVal.Length%numRep = 0 then
                    //    vector [for i = 0 to xVal.Length - 1 do yield 1. / (yVal.[(i*numRep)..(i*numRep + numRep - 1)] |> Seq.stDev)]
                    //else failwithf "arrLength no multiple of replicate number"
                    yVal 
                    |> Array.map (fun x -> 1. / max 0.00001 (Seq.stDev x))
                    |> vector
                //improves robustness of standard weighting
                let (maxWtemp,minWtemp) = 
                    let mean = diagVec |> Seq.mean
                    let std  = diagVec |> Seq.stDev
                    (mean + std),(mean - std)
                
                let diagVecNew = 
                    diagVec |> Vector.map (fun x -> match x with
                                                    | x when x > maxWtemp -> maxWtemp
                                                    | x when x < minWtemp -> minWtemp
                                                    | _ -> x)
                let finalDiag = 
                    let trace = (diagVecNew |> Vector.sum)
                    let length = (float Wtemp.NumCols)
                    diagVecNew |> Vector.map (fun x -> x / trace * length)
                for i = 0 to xVal.Length - 1 do Wtemp.[i,i] <- finalDiag.[i]
                Wtemp
            | CV -> 
                let Wtemp = Matrix.create xVal.Length xVal.Length 0.
                let cvOfVec =
                    yVal 
                    |> Array.map (fun x -> 
                        max 0. (System.Math.Log((1./(System.Math.Abs(Seq.cvPopulation x))),2.)))
                    |> vector
                    |> Vector.map (fun w -> if nan.Equals w then 0. else w)
                    
                    //if yVal.Length%numRep = 0 then
                    //    let length = yVal.Length / numRep 
                    //    let cv = vector [for i = 0 to length-1 do yield yVal.[i*numRep..i*numRep+numRep-1] |> fun g -> max 0. (System.Math.Log((1./(Math.Abs(Seq.cvPopulation g))),2.))(*(1./((Seq.cvPopulation g) + 0.25))*)] //0.25???
                    //    cv
                    //else failwithf "arrLength no multiple of replicate number"

                for i = 0 to xVal.Length - 1 do Wtemp.[i,i] <- cvOfVec.[i]
                Wtemp
            | StandardDeviation     ->
                let Wtemp = Matrix.create xVal.Length xVal.Length 0.            
                let diagVec = 
                    //if yVal.Length%numRep = 0 then
                    //    vector [for i = 0 to xVal.Length - 1 do yield 1. / (yVal.[(i*numRep)..(i*numRep + numRep - 1)] |> Seq.stDev)]
                    //else failwithf "arrLength no multiple of replicate number"
                    yVal 
                    |> Array.map (fun x -> 1. / max 0.00001 (Seq.stDev x))
                    |> vector
                let avg = Seq.average diagVec
                for i = 0 to xVal.Length - 1 do Wtemp.[i,i] <- diagVec.[i] / avg
                Wtemp
            | StandardDeviationSqrt ->
                let Wtemp = Matrix.create xVal.Length xVal.Length 0.            
                let diagVec = 
                    //if yVal.Length%numRep = 0 then
                    //    vector [for i = 0 to xVal.Length - 1 do yield Math.Sqrt (1. / (yVal.[(i*numRep)..(i*numRep + numRep - 1)] |> Seq.stDev))]
                    //else failwithf "arrLength no multiple of replicate number"
                    yVal 
                    |> Array.map (fun x -> Math.Sqrt (1. / max 0.00001 (Seq.stDev x)))
                    |> vector
                let avg = Seq.average diagVec
                for i = 0 to xVal.Length - 1 do Wtemp.[i,i] <- diagVec.[i] / avg
                Wtemp
            | StandardDeviationAdj -> //CV_Norm          ->
                let Wtemp = Matrix.create xVal.Length xVal.Length 0.
                let cvOfVec =
                    //if yVal.Length%numRep = 0 then
                    //    let length = yVal.Length / numRep 
                    //    let cv = vector [for i = 0 to length-1 do yield yVal.[i*numRep..i*numRep+numRep-1] |> fun x -> sqrt ( 1. / (Seq.stDev x / Seq.average x))]
                    //    cv
                    //else failwithf "arrLength no multiple of replicate number"
                    
                    yVal 
                    |> Array.map (fun x -> sqrt ( 1. / (max 0.00001 (Seq.stDev x) / Seq.average x)))
                    |> vector
                let cvAvg = Seq.average cvOfVec
                for i = 0 to xVal.Length - 1 do Wtemp.[i,i] <- cvOfVec.[i] / cvAvg
                Wtemp
            | StandardDeviationAdjSqrt ->  //CV_NormSqrt      ->
                let Wtemp = Matrix.create xVal.Length xVal.Length 0.
                let cvOfVec =
                    //if yVal.Length%numRep = 0 then
                    //    let length = yVal.Length / numRep 
                    //    let cv = vector [for i = 0 to length-1 do yield yVal.[i*numRep..i*numRep+numRep-1] |> fun x -> sqrt(sqrt ( 1. / (Seq.stDev x / Seq.average x)))]
                    //    cv
                    //else failwithf "arrLength no multiple of replicate number"
                    yVal 
                    |> Array.map (fun x -> sqrt (sqrt ( 1. / (max 0.00001 (Seq.stDev x) / Seq.average x))))
                    |> vector
                let cvAvg = Seq.average cvOfVec
                for i = 0 to xVal.Length - 1 do Wtemp.[i,i] <- cvOfVec.[i] / cvAvg
                Wtemp
                   
        /// <summary>Calculates the Akaike Information criterion from actual values to a prediction
        /// </summary>
        /// <param name="d">complexity measure(count of extrema)</param>
        /// <param name="a">vector of predicted y values</param>
        /// <param name="y">vector of actual y values</param>
        /// <param name="W">Matrix with weightings encoded as diagonal</param>
        /// <returns>Akaike Information Criterion corrected for small sample sizes.</returns>
        let getAIC d (a:Vector<float>) (y:Vector<float>) (W:Matrix<float>) =
            let sse =
                let weightedResiduals =
                    W * (y - a)
                let nSSE =
                    pown (weightedResiduals |> Vector.norm) 2
                nSSE
            let n = a.Length |> float
            n * Math.Log(sse/n) + 2. * (float d) + ((2. * float d * (float d + 1.))/(float n - d - 1.))

        /// <summary>calculates an smoothing spline with the given weighting without shape constraints</summary>
        /// <param name="D">Matrix D as defined in Wood 1994</param>
        /// <param name="H">Matrix H as defined in Wood 1994</param>
        /// <param name="W">Matrix with weightings encoded as diagonal</param>
        /// <param name="y">Vector of y value estimates</param>
        /// <param name="xVal">Vector of x values</param>
        /// <returns>Object for fitting the smoothing splines, together with quality measures.</returns>
        let getInitialEstimate (D:Matrix<float>) (H:Matrix<float>) (W:Matrix<float>) y (xVal:Vector<float>) = 
            let n = xVal.Length
            let I = Matrix.identity n
            let Z = Matrix.identity n
            let calcGLambda lambd = 
                let dInverse = Algebra.LinearAlgebra.Inverse D
                (2.0 * ( H.Transpose * dInverse * H + (lambd / float n) * (W.Transpose * W)))
                
            let A_tild (lambd:float)    = 
                let zproInverse = Algebra.LinearAlgebra.Inverse (Z.Transpose * (calcGLambda lambd) * Z)
                2. * (lambd / (float n)) * Z * zproInverse * Z.Transpose * W.Transpose * W
        
            let a_star (lambd:float)= 
                (A_tild lambd) * y
                
            let V_tild (lambd:float) =     
                //robustGCV for small sample size
                let rho = 1.3

                let no = 
                    (W * ((a_star lambd) - y))
                    |> Vector.fold (fun acc x -> acc + (pown x 2)) 0.   
                
                let tr = 
                    I - (rho * (A_tild lambd)) 
                    |> Matrix.trace
                (float n) * no / (tr * tr)
        
            let Var_tild lambd =
                let aMat = (a_star lambd) |> Matrix.ofVector
                let yMat = y |> Matrix.ofVector
                let no = 
                    (W * (aMat - yMat))
                    |> fun m -> Matrix.getCol m 0
                    |> Vector.norm
                    |> fun x -> x*x
                let rho = 1.3
                let tr = I - (rho * (A_tild lambd)) |> Matrix.getDiag |> Vector.sum //
                (float n) * no / (tr)  
        
            //NelderMeadSolver for global minimum
            let lambda_hat = 
                //let lamlist = [for i = 0 to 64 do yield 0.01*(1.2**(float i))]
                //lamlist
                //|> List.map V_tild
                //|> List.min
                NelderMead.minimizeSingleWith V_tild 3. 0.0001 2000.    
            
            let variance_unconstrained = Var_tild lambda_hat
            let GCV_unconstrained = V_tild lambda_hat
            let a_unconstrained = a_star lambda_hat
            let c =         
                let tmpH = H
                let tmpD = D
                let tmp = 
                    let tmp' = Algebra.LinearAlgebra.SolveLinearSystems tmpD tmpH
                    tmp' * a_unconstrained
                Vector.init (tmp.Length+2) (fun i -> if i>0 && i<=tmp.Length then tmp.[i-1] else 0.0)

            let getShape = getExtrema xVal a_unconstrained c

            let aic = getAIC (float getShape.Length) a_unconstrained y W

            ///function that calculates the y_Value of the spline corresponding to the given x_Value
            let splineFunction = initEvalAt xVal a_unconstrained c

            TempClassResult.Create (Some xVal) None None None a_unconstrained c 0. GCV_unconstrained 0. (Array2D.zeroCreate 0 0) aic splineFunction None
        
        /// <summary>calculates an smoothing spline with the given weighting without shape constraints</summary>
        /// <param name="W">Matrix with weightings encoded as diagonal</param>
        /// <param name="y">Vector of y value estimates</param>
        /// <param name="xVal">Vector of x values</param>
        /// <returns>Object for fitting the smoothing splines, together with quality measures.</returns>
        let getInitialEstimateOfWeighting W (y:Vector<float>) (xVal:Vector<float>) =
            let n = y |> Seq.length
        
            let h = Array.init (n-1) (fun i -> xVal.[i+1] - xVal.[i] )
            let H = Array2D.zeroCreate (n-2) n
            let D = Array2D.zeroCreate (n-2) (n-2)
            for i = 1 to (n-2) do
                H.[i-1,i-1] <-  1.0/h.[i-1]
                H.[i-1,i]   <- -( (1.0/h.[i-1]) + (1.0/h.[i]) )
                H.[i-1,i+1] <-  1.0/h.[i]
                D.[i-1,i-1] <-  (h.[i-1] + h.[i]) / 3.0
            for i = 1 to (n-3) do
                D.[i-1,i]   <- h.[i]/6.0
                D.[i,i-1] <- h.[i]/6.0        
        
            getInitialEstimate (Matrix.ofArray2D D) (Matrix.ofArray2D H) W y xVal

        /// <summary>calculates an unconstrained smoothing spline</summary>
        /// <param name="xVal">Vector of x values</param>
        /// <param name="yVal">Vector of y value estimates</param>
        /// <returns>Object for fitting the smoothing splines, together with quality measures.</returns>
        let getUnconstrainedSpline (xVal:Vector<float>) (yVal:Vector<float>) =
            let W = Matrix.diag (Vector.oneCreate xVal.Length)
            let n = yVal |> Seq.length
    
            let h = Array.init (n-1) (fun i -> xVal.[i+1] - xVal.[i] )
            let H = Array2D.zeroCreate (n-2) n
            let D = Array2D.zeroCreate (n-2) (n-2)
            for i = 1 to (n-2) do
                H.[i-1,i-1] <-  1.0/h.[i-1]
                H.[i-1,i]   <- -( (1.0/h.[i-1]) + (1.0/h.[i]) )
                H.[i-1,i+1] <-  1.0/h.[i]
                D.[i-1,i-1] <-  (h.[i-1] + h.[i]) / 3.0
            for i = 1 to (n-3) do
                D.[i-1,i]   <- h.[i]/6.0
                D.[i,i-1] <- h.[i]/6.0        
    
            getInitialEstimate (Matrix.ofArray2D D) (Matrix.ofArray2D H) W yVal xVal

        
        /// <summary>calculates the generalized cross validation(GCV), 
        /// and the modified GCV for small sample sizes(mGCV)
        /// </summary>
        /// <param name="C">constraint matrix</param>
        /// <param name="D">Matrix D as defined in Wood 1994</param>
        /// <param name="H">Matrix H as defined in Wood 1994</param>
        /// <param name="W">Matrix with weightings encoded as diagonal</param>
        /// <param name="y">vector of actual y values</param>
        /// <param name="a">vector of predicted y values</param>
        /// <param name="lambda">regularization parameter lambda</param>
        /// <returns>Tuple in the form of (mGCV,GCV)</returns>
        let getGCV (C:Matrix<float>) (D:Matrix<float>) (H:Matrix<float>) (W:Matrix<float>) (y:vector) (a:vector) (lambda:float) =
            let n = a.Length
            //identify inequality constraints where C*a<=[0]
            let constr = C * (a |> Matrix.ofVector) 
            let tol = 0.001 
            //identify 'Active set' where C*a=[0]
            let Ca =    
                let rec loop i caPrev =
                    if i < C.NumRows then
                        if System.Math.Abs(constr.[i,0]) <= tol then 
                            loop (i+1) (((C.Row i) |> Seq.toArray)::caPrev)
                        else loop (i+1) caPrev
                    else 
                        if caPrev |> List.isEmpty then
                            Matrix.create 0 0 0.
                        else 
                        matrix (caPrev |> List.rev)
                loop 0 []
            ///calculate null space of Ca (based on SVD)
            let Z =
                if Ca.NumRows * Ca.NumCols = 0 then 
                    Matrix.identity n            
                else
                    let nullspace = 
                        let (eps,U,V) =
                            Algebra.LinearAlgebra.SVD Ca
                        //let (U,eps,V) =
                        //    TC_SVD.compute (Matrix.toArray2D Ca)
                        //    |> fun (u,e,v) -> u,e,(Matrix.ofArray2D v)
                        let rank = eps |> Seq.filter (fun x -> x >= 1e-5) |> Seq.length //Threshold if a singular value is considered as 0. //mat.NumRows - eps.Length
                        let count = V.NumRows - rank
                        Matrix.getRows V (rank) (count)
                        |> Matrix.transpose
                    nullspace
           
            let I = Matrix.identity n

            let G_lambda (lambd:float)  = 
                let dInverse = Algebra.LinearAlgebra.Inverse D
                2. * ( H.Transpose * dInverse * H + lambd / (float n)  * W.Transpose * W)
        
            let A_tild (lambd:float)    = 
                let zproInverse = Algebra.LinearAlgebra.Inverse (Z.Transpose * (G_lambda lambd) * Z)
                2. * (lambd / (float n)) * Z * zproInverse * Z.Transpose * W.Transpose * W

            let V_tild rho (lambd:float) =     
                let no =
                    W * (a - y)
                    |> Vector.fold (fun acc x -> acc + (pown x 2)) 0.   
  
                let tr = 
                    I - rho * (A_tild lambd)
                    |> Matrix.trace

                (float n) * no / (tr * tr)

            let var_ny (lambd:float) =
                let a = ((W * (a-y)) |>Vector.norm ) ** 2.
                let b = I - (A_tild lambd) |> Matrix.getDiag |> Vector.sum
                a/b

            //standard gcv
            let gcv  = V_tild 1.0 lambda
            //modified gcv
            let mGcv = V_tild 1.3 lambda
            //robust gcv
            let rGcv() =
                let A = A_tild lambda
                let gamma = 0.2
                let mu = ((A * A) |> Matrix.trace) / (float n)
                (gamma + (1. - gamma) * mu) * gcv


            let variance = var_ny lambda
            (mGcv,gcv)

        /// <summary>calculates a spline with the shape defined by the operator(increasing or decreasing) 
        /// and the condition(0 - 4 extrema)</summary>
        /// <param name="operator">id (~+) or negative transforming (~-)</param>
        /// <param name="x">vector of x values</param>
        /// <param name="y">vector of actual y values</param>
        /// <param name="W">Matrix with weightings encoded as diagonal</param>
        /// <param name="con">number of extrema that shall be enforced</param>
        /// <returns>Tuple of constraint matrices, that fulfill the conditions (operator and con)
        /// and the final spline solution that minimizes the mGCV.</returns>
        let internal spline operator (x:Vector<float>) (y:Vector<float>) (W:Matrix<float>) con =

            /// Sets the diagonal to value inplace
            let setDiagonalInplace value (m:Matrix<_>) =
                let min = min m.NumRows m.NumCols
                for i=0 to min-1 do
                    m.[i,i] <- value
                m

            let calcGlambda  (D:Matrix<float>) (H:Matrix<float>) (W:Matrix<float>) (n:int) lambd = 
                let dInverse = Algebra.LinearAlgebra.Inverse D
                2.0 * ((H.Transpose * dInverse * H)   +  (lambd / float n) * W.Transpose * W)         

            let calcclambda (W:Matrix<float>) (n:int) (y:Vector<float>) lambd =
                //-2.0 * ((float lambd) / (float n)) * (W.Transpose * W) * y |> Vector.toArray//(rowvec * matrix) <> (matrix * vec)
                -2.0 * (float lambd) / (float n) * y.Transpose * W.Transpose * W
                
            ///calculates the inequality constraint matrix C for the quadratic programming approach
            let getCtemP n (B:Matrix<float>) (h:Vector<float>) =
                let mutable list :Matrix<float> list= []
                let Ctemp = Matrix.zero (4*(n-1)+1) n

                if con=0 then
                    for i =1 to (n-1) do    
                        let temp = 
                            (List.init (i-1) (fun x -> 0.))@((-3./h.[i-1])::(3./h.[i-1])::(List.init (n-i-1) (fun x -> 0.)))
                            |> vector
                        Matrix.setRow Ctemp (4*i-4) temp
                        Matrix.setRow Ctemp (4*i-3) (B.Row (i-1) |> vector)
                        Matrix.setRow Ctemp (4*i-2) (temp - (B.Row (i-1) |> vector))
                        Matrix.setRow Ctemp (4*i-1) (temp - (B.Row (i) |> vector))
                    Matrix.setRow Ctemp (4*(n-1)) (B.Row (n - 1) |> vector)
                    [|Ctemp|]

                elif (con=1) then
                    //for j = 3 to (n - 1) do 
                    for j = 2 to (n - 1) do 
                        for i = 1 to (j-1) do
                            let temp = 
                                (List.init (i-1) (fun x -> 0.))@((-3./h.[i-1])::(3./h.[i-1])::(List.init (n-i-1) (fun x -> 0.)))
                                |> vector
                            Matrix.setRow Ctemp (4*i-4) temp
                            Matrix.setRow Ctemp (4*i-3) (B.Row (i-1) |> vector)
                            Matrix.setRow Ctemp (4*i-2) (temp - (B.Row (i-1) |> vector))
                            Matrix.setRow Ctemp (4*i-1) (temp - (B.Row (i) |> vector))
                            if i = j - 1 then
                                Matrix.setRow Ctemp (4*i-4) (Vector.zeroCreate Ctemp.NumCols)
                                Matrix.setRow Ctemp (4*i-3) (B.Row (i-1) |> vector) //0
                                Matrix.setRow Ctemp (4*i-2) (Vector.zeroCreate Ctemp.NumCols)
                                Matrix.setRow Ctemp (4*i-1) (Vector.zeroCreate Ctemp.NumCols)

                        for i = j to (n-1) do 
                            let temp = 
                                (List.init (i-1) (fun x -> 0.))@((-3./h.[i-1])::(3./h.[i-1])::(List.init (n-i-1) (fun x -> 0.)))
                                |> vector
                            Matrix.setRow Ctemp (4*i-4) (- temp)
                            Matrix.setRow Ctemp (4*i-3) (- (B.Row (i-1) |> vector))
                            Matrix.setRow Ctemp (4*i-2) (- (temp - (B.Row (i-1) |> vector)))
                            Matrix.setRow Ctemp (4*i-1) (- (temp - (B.Row (i) |> vector)))
                        Matrix.setRow Ctemp (4*(n-1)) (- (B.Row (n-1) |> vector))
                        list <- (Matrix.copy Ctemp)::list
                    list
                    |> Array.ofList 
                    
                elif con = 2 then
                    //for m = 3 to (n-1) do 
                    for m = 2 to (n-1) do 
                        for j = (m+1) to (n-1) do //from 7 to 6 do
                            for i=1 to (m-1) do 
                                let temp = 
                                    (List.init (i-1) (fun x -> 0.))@((-3./h.[i-1])::(3./h.[i-1])::(List.init (n-i-1) (fun x -> 0.)))
                                    |> vector
                                Matrix.setRow Ctemp (4*i-4) temp
                                Matrix.setRow Ctemp (4*i-3) (B.Row (i-1) |> vector)
                                Matrix.setRow Ctemp (4*i-2) (temp - (B.Row (i-1) |> vector))
                                Matrix.setRow Ctemp (4*i-1) (temp - (B.Row (i) |> vector))
                                if i = m - 1 then
                                    Matrix.setRow Ctemp (4*i-4) (Vector.zeroCreate Ctemp.NumCols)
                                    Matrix.setRow Ctemp (4*i-3) (B.Row (i-1) |> vector)
                                    Matrix.setRow Ctemp (4*i-2) (Vector.zeroCreate Ctemp.NumCols)
                                    Matrix.setRow Ctemp (4*i-1) (Vector.zeroCreate Ctemp.NumCols)
                            for i = m to j-1 do
                                let temp =                
                                    (List.init (i-1) (fun x -> 0.))@((-3./h.[i-1])::(3./h.[i-1])::(List.init (n-i-1) (fun x -> 0.)))
                                    |> vector
                                Matrix.setRow Ctemp (4*i-4) (- temp)
                                Matrix.setRow Ctemp (4*i-3) (- (B.Row (i-1) |> vector))
                                Matrix.setRow Ctemp (4*i-2) (- (temp - (B.Row (i-1) |> vector)))
                                Matrix.setRow Ctemp (4*i-1) (- (temp - (B.Row (i) |> vector)))
                                if i = j - 1 then
                                    Matrix.setRow Ctemp (4*i-4) (Vector.zeroCreate Ctemp.NumCols)
                                    Matrix.setRow Ctemp (4*i-3) -(B.Row (i-1) |> vector)
                                    Matrix.setRow Ctemp (4*i-2) (Vector.zeroCreate Ctemp.NumCols)
                                    Matrix.setRow Ctemp (4*i-1) (Vector.zeroCreate Ctemp.NumCols)
                            for i = j to n-1 do
                                let temp = 
                                    (List.init (i-1) (fun x -> 0.))@((-3./h.[i-1])::(3./h.[i-1])::(List.init (n-i-1) (fun x -> 0.)))
                                    |> vector
                                Matrix.setRow Ctemp (4*i-4) (temp)
                                Matrix.setRow Ctemp (4*i-3) ((B.Row (i-1) |> vector))
                                Matrix.setRow Ctemp (4*i-2) ((temp - (B.Row (i-1) |> vector)))
                                Matrix.setRow Ctemp (4*i-1) ((temp - (B.Row (i) |> vector)))
                            Matrix.setRow Ctemp (4*(n-1)) (B.Row (n-1) |> vector)
                            list <- (Matrix.copy Ctemp)::list
                    list
                    |> Array.ofList 

                elif con = 3 then
                    //for m = 3 to (n-1) do 
                    for m = 2 to (n-1) do 
                        for j = (m+1) to (n-1) do
                            for k = (j+1) to (n-1) do
                                for i=1 to (m-1) do 
                                    let temp = 
                                        (List.init (i-1) (fun x -> 0.))@((-3./h.[i-1])::(3./h.[i-1])::(List.init (n-i-1) (fun x -> 0.)))
                                        |> vector
                                    Matrix.setRow Ctemp (4*i-4) temp
                                    Matrix.setRow Ctemp (4*i-3) (B.Row (i-1) |> vector)
                                    Matrix.setRow Ctemp (4*i-2) (temp - (B.Row (i-1) |> vector))
                                    Matrix.setRow Ctemp (4*i-1) (temp - (B.Row (i) |> vector))
                                    if i = m - 1 then
                                        Matrix.setRow Ctemp (4*i-4) (Vector.zeroCreate Ctemp.NumCols)
                                        Matrix.setRow Ctemp (4*i-3) (B.Row (i-1) |> vector)
                                        Matrix.setRow Ctemp (4*i-2) (Vector.zeroCreate Ctemp.NumCols)
                                        Matrix.setRow Ctemp (4*i-1) (Vector.zeroCreate Ctemp.NumCols)
                                for i = m to j-1 do
                                    let temp =                
                                        (List.init (i-1) (fun x -> 0.))@((-3./h.[i-1])::(3./h.[i-1])::(List.init (n-i-1) (fun x -> 0.)))
                                        |> vector
                                    Matrix.setRow Ctemp (4*i-4) (- temp)
                                    Matrix.setRow Ctemp (4*i-3) (- (B.Row (i-1) |> vector))
                                    Matrix.setRow Ctemp (4*i-2) (- (temp - (B.Row (i-1) |> vector)))
                                    Matrix.setRow Ctemp (4*i-1) (- (temp - (B.Row (i) |> vector)))
                                    if i = j - 1 then
                                        Matrix.setRow Ctemp (4*i-4) (Vector.zeroCreate Ctemp.NumCols)
                                        Matrix.setRow Ctemp (4*i-3) -(B.Row (i-1) |> vector)
                                        Matrix.setRow Ctemp (4*i-2) (Vector.zeroCreate Ctemp.NumCols)
                                        Matrix.setRow Ctemp (4*i-1) (Vector.zeroCreate Ctemp.NumCols)
                                for i = j to k-1 do
                                    let temp =                
                                        (List.init (i-1) (fun x -> 0.))@((-3./h.[i-1])::(3./h.[i-1])::(List.init (n-i-1) (fun x -> 0.)))
                                        |> vector
                                    Matrix.setRow Ctemp (4*i-4) (temp)
                                    Matrix.setRow Ctemp (4*i-3) (B.Row (i-1) |> vector)
                                    Matrix.setRow Ctemp (4*i-2) (temp - (B.Row (i-1) |> vector))
                                    Matrix.setRow Ctemp (4*i-1) (temp - (B.Row (i) |> vector))
                                    if i = k - 1 then //if i = j - 1 then
                                        Matrix.setRow Ctemp (4*i-4) (Vector.zeroCreate Ctemp.NumCols)
                                        Matrix.setRow Ctemp (4*i-3) (B.Row (i-1) |> vector) 
                                        Matrix.setRow Ctemp (4*i-2) (Vector.zeroCreate Ctemp.NumCols)
                                        Matrix.setRow Ctemp (4*i-1) (Vector.zeroCreate Ctemp.NumCols)
                                for i = k to n-1 do
                                    let temp = 
                                        (List.init (i-1) (fun x -> 0.))@((-3./h.[i-1])::(3./h.[i-1])::(List.init (n-i-1) (fun x -> 0.)))
                                        |> vector
                                    Matrix.setRow Ctemp (4*i-4) (- temp)
                                    Matrix.setRow Ctemp (4*i-3) (-(B.Row (i-1) |> vector))
                                    Matrix.setRow Ctemp (4*i-2) (-(temp - (B.Row (i-1) |> vector)))
                                    Matrix.setRow Ctemp (4*i-1) (-(temp - (B.Row (i) |> vector)))
                                Matrix.setRow Ctemp (4*(n-1)) (- B.Row (n-1) |> vector)
                                list <- (Matrix.copy Ctemp)::list
                    list
                    |> Array.ofList 
                    
                elif con = 4 then
                    //for m = 3 to (n-1) do 
                    for m = 2 to (n-1) do 
                        for j = (m+1) to (n-1) do
                            for k = (j+1) to (n-1) do
                                for l = (k+1) to (n-1) do //
                                    for i=1 to (m-1) do 
                                        let temp = 
                                            (List.init (i-1) (fun x -> 0.))@((-3./h.[i-1])::(3./h.[i-1])::(List.init (n-i-1) (fun x -> 0.)))
                                            |> vector
                                        Matrix.setRow Ctemp (4*i-4) temp
                                        Matrix.setRow Ctemp (4*i-3) (B.Row (i-1) |> vector)
                                        Matrix.setRow Ctemp (4*i-2) (temp - (B.Row (i-1) |> vector))
                                        Matrix.setRow Ctemp (4*i-1) (temp - (B.Row (i) |> vector))
                                        if i = m - 1 then
                                            Matrix.setRow Ctemp (4*i-4) (Vector.zeroCreate Ctemp.NumCols)
                                            Matrix.setRow Ctemp (4*i-3) (B.Row (i-1) |> vector)
                                            Matrix.setRow Ctemp (4*i-2) (Vector.zeroCreate Ctemp.NumCols)
                                            Matrix.setRow Ctemp (4*i-1) (Vector.zeroCreate Ctemp.NumCols)
                                    for i = m to j-1 do
                                        let temp =                
                                            (List.init (i-1) (fun x -> 0.))@((-3./h.[i-1])::(3./h.[i-1])::(List.init (n-i-1) (fun x -> 0.)))
                                            |> vector
                                        Matrix.setRow Ctemp (4*i-4) (- temp)
                                        Matrix.setRow Ctemp (4*i-3) (- (B.Row (i-1) |> vector))
                                        Matrix.setRow Ctemp (4*i-2) (- (temp - (B.Row (i-1) |> vector)))
                                        Matrix.setRow Ctemp (4*i-1) (- (temp - (B.Row (i) |> vector)))
                                        if i = j - 1 then
                                            Matrix.setRow Ctemp (4*i-4) (Vector.zeroCreate Ctemp.NumCols)
                                            Matrix.setRow Ctemp (4*i-3) -(B.Row (i-1) |> vector)
                                            Matrix.setRow Ctemp (4*i-2) (Vector.zeroCreate Ctemp.NumCols)
                                            Matrix.setRow Ctemp (4*i-1) (Vector.zeroCreate Ctemp.NumCols)
                                    for i = j to k-1 do
                                        let temp =                
                                            (List.init (i-1) (fun x -> 0.))@((-3./h.[i-1])::(3./h.[i-1])::(List.init (n-i-1) (fun x -> 0.)))
                                            |> vector
                                        Matrix.setRow Ctemp (4*i-4) (temp)
                                        Matrix.setRow Ctemp (4*i-3) (B.Row (i-1) |> vector)
                                        Matrix.setRow Ctemp (4*i-2) (temp - (B.Row (i-1) |> vector))
                                        Matrix.setRow Ctemp (4*i-1) (temp - (B.Row (i) |> vector))
                                        if i = k - 1 then //if i = j - 1 then
                                            Matrix.setRow Ctemp (4*i-4) (Vector.zeroCreate Ctemp.NumCols)
                                            Matrix.setRow Ctemp (4*i-3) (B.Row (i-1) |> vector) 
                                            Matrix.setRow Ctemp (4*i-2) (Vector.zeroCreate Ctemp.NumCols)
                                            Matrix.setRow Ctemp (4*i-1) (Vector.zeroCreate Ctemp.NumCols)
                                    for i = k to l-1 do
                                        let temp =                
                                            (List.init (i-1) (fun x -> 0.))@((-3./h.[i-1])::(3./h.[i-1])::(List.init (n-i-1) (fun x -> 0.)))
                                            |> vector
                                        Matrix.setRow Ctemp (4*i-4) (temp)
                                        Matrix.setRow Ctemp (4*i-3) (B.Row (i-1) |> vector)
                                        Matrix.setRow Ctemp (4*i-2) (temp - (B.Row (i-1) |> vector))
                                        Matrix.setRow Ctemp (4*i-1) (temp - (B.Row (i) |> vector))
                                        if i = l - 1 then //if i = j - 1 then
                                            Matrix.setRow Ctemp (4*i-4) (Vector.zeroCreate Ctemp.NumCols)
                                            Matrix.setRow Ctemp (4*i-3) (B.Row (i-1) |> vector) 
                                            Matrix.setRow Ctemp (4*i-2) (Vector.zeroCreate Ctemp.NumCols)
                                            Matrix.setRow Ctemp (4*i-1) (Vector.zeroCreate Ctemp.NumCols)
                                    for i = l to n-1 do
                                        let temp = 
                                            (List.init (i-1) (fun x -> 0.))@((-3./h.[i-1])::(3./h.[i-1])::(List.init (n-i-1) (fun x -> 0.)))
                                            |> vector
                                        Matrix.setRow Ctemp (4*i-4) (- temp)
                                        Matrix.setRow Ctemp (4*i-3) (-(B.Row (i-1) |> vector))
                                        Matrix.setRow Ctemp (4*i-2) (-(temp - (B.Row (i-1) |> vector)))
                                        Matrix.setRow Ctemp (4*i-1) (-(temp - (B.Row (i) |> vector)))
                                    Matrix.setRow Ctemp (4*(n-1)) (- B.Row (n-1) |> vector)
                                    list <- (Matrix.copy Ctemp)::list
                    list
                    |> Array.ofList  
                else failwith "Condition max 4"

            ///calculates the weighted squared error of the spline functions and the observation values at the knots
            let getError (y:Vector<float>) (a:Vector<float>) (W:Matrix<float>) =
                let tmp =  W * (y - a)
                tmp |> Seq.averageBy (fun x -> x**2.) |> sqrt
            
            //datalength
            let n = x.Length

            //Define intervals stepwidth
            let h = Array.init (n-1) (fun i -> x.[i+1] - x.[i] )
            
            //generation of matrices D and H (Wood94)
            let H = Array2D.zeroCreate (n-2) n
            let D = Array2D.zeroCreate (n-2) (n-2)

            for i = 1 to (n-2) do
                H.[i-1,i-1] <-  1.0/h.[i-1]
                H.[i-1,i]   <- -( (1.0/h.[i-1]) + (1.0/h.[i]) )
                H.[i-1,i+1] <-  1.0/h.[i]
                D.[i-1,i-1] <-  (h.[i-1] + h.[i]) / 3.0

            for i = 1 to (n-3) do
                D.[i-1,i]   <- h.[i]/6.0
                D.[i,i-1] <- h.[i]/6.0

            //generation of matrices P, U and B (Wood94) --- Matrix P corrected!
            let P = Array2D.zeroCreate n n
            let U = Array2D.zeroCreate n n

            for i = 1 to n do
                P.[i-1,i-1] <- 2.0

            for i = 2 to (n-1) do
                P.[i-1,i-2] <- h.[i-1] / (h.[i-1] + h.[i-2])//uncorrected: h.[i-1] / (h.[i-1] + h.[i]) --> h.[i] not existent when i = n-1
                P.[i-1,i] <- 1.0 - P.[i-1,i-2]
                U.[i-1,i-2] <- -3.0 * (P.[i-1,i-2] / h.[i-2])
                U.[i-1,i] <- 3.0 * (P.[i-1,i] / h.[i-1])
                U.[i-1,i-1] <- -(U.[i-1,i] + U.[i-1,i-2])

            P.[0,1] <- 1.0
            U.[0,0] <- -3.0 / h.[0]
            U.[0,1] <- - U.[0,0]
            P.[n-1,n-2] <- 1.0
            U.[n-1,n-2] <- -3.0 / h.[n-2]
            U.[n-1,n-1] <- - U.[n-1,n-2]

            ///Calculation of B by P*B = U
            let PMat = Matrix.ofArray2D P
            let UMat = Matrix.ofArray2D U
            let B = Algebra.LinearAlgebra.SolveLinearSystems PMat UMat

            ///generation of the constraint matrix
            let ctList = getCtemP n B (Vector.ofArray h)
            
            ///the constraint-matrices for a given shape are transformed into 2D-Arrays (solver compatibility) and the operator defines the monotonicity direction
            let A = 
                ctList 
                |> Array.map (fun x -> 
                    operator(x) 
                    |> Matrix.toArray2D)
            
            let lamlist = [for i = 0 to 64 do yield 0.01*(1.2**(float i))]
            //let lamlist = [6.]//here attention
            
            let lambdaStrength = 
                A
                |> Array.map (fun aMat ->
                    lamlist
                    |> List.mapi (fun z lambd ->
                        let Q = calcGlambda (Matrix.ofArray2D D) (Matrix.ofArray2D H) W n lambd //outsource (performance)
                        let c = (calcclambda W n y lambd) |> RowVector.toArray                  //outsource (performance)
                        //borders to solve the inequality constraints 
                        let b = Array.create ((A.[0] |> Array2D.length1)) (-infinity,0.)
                        //solver specific transformation (because symmetry of the matrix)
                        let a' = 
                            let tmp = Matrix.create Q.NumRows Q.NumCols 1. |> setDiagonalInplace 0.5
                            let Q' = (Q .* tmp) |> Matrix.toArray2D
                            QP.minimize aMat b Q' c |> Vector.ofArray
                        let e' = getError y a' W
                        //let c' =         
                        //    let tmpH = Matrix.ofArray2D H
                        //    let tmpD = Matrix.ofArray2D D
                        //    let tmp = 
                        //        let tmp' = Algebra.LinearAlgebra.SolveLinearSystems tmpD tmpH //tested
                        //        tmp' * a'                                                     //tested
                        //    Vector.init (tmp.Length+2) (fun i -> if i>0 && i<=tmp.Length then tmp.[i-1] else 0.0)

                        let (mgcv,gcv) = getGCV (aMat |> Matrix.ofArray2D) (Matrix.ofArray2D D) (Matrix.ofArray2D H) W y a' lambd
                        InnerResult.Create a' mgcv aMat e' gcv lambd
                        )
                //minimize the gcv in respect to the used smoothing factor lambda
                
                )
                //|> fun (a,e,gcvvar,c) -> 
                //    let gcv = gcvvar |> Array.map fst

                //    gcv 
                //    |> Array.indexed 
                //    |> Array.map (fun (i,x) -> sprintf "%i\t%f" i x)
                //    |> fun asd -> System.IO.File.WriteAllLines(@"C:\Users\bvenn\Documents\Projects\EIT\June\TestGCV\" + (string aind) + ".txt",asd)

                //    let var = gcvvar |> Array.map snd
                //    let gcvminIndex = gcv |> Array.findIndex (fun x -> x = (gcv |> Array.min))
                //    let varmin = var.[gcvminIndex]
                //    let gcvmin = gcv.[gcvminIndex]
                //    let efin   = e.[gcvminIndex]
                //    let afin   = a.[gcvminIndex]
                //    let cfin   = c.[gcvminIndex]
                //    let lambdafin = 0.01*(1.2**(float gcvminIndex))
                //    ((efin,gcvmin),(afin,varmin),lambdafin,cfin)
                //        )
            lambdaStrength
            |> fun models -> 
                let minimal = 
                    models
                    |> Array.map (List.minBy (fun innerResult -> innerResult.GCV))
                minimal
                |> fun xt -> 
                    //additional case, when no spline satisfies the given conditions (with additional constraints that checks extrema count)
                    if xt = [||] then
                        //for shape restriction if coefficient shrinkage is to intense (uncomment above)
                        A,TempClassResult.Create None None None None (Vector.init n (fun _ -> 0.)) (Vector.init n (fun _ -> 0.)) infinity infinity infinity (Array2D.zeroCreate 1 1) infinity id None
                    else 
                        xt
                        //among all shape possibilities under a given parent shape ((min, then max);(max);...) minimize the mGCV to obtain the most promising spline candidate
                        |> Array.indexed
                        |> Array.minBy (fun (_,innerResult) -> innerResult.GCV)
                        |> fun (resultindex,result) -> 
                            let traceC = 
                                let tmpH = Matrix.ofArray2D H
                                let tmpD = Matrix.ofArray2D D
                                let tmp = 
                                    let tmp' = Algebra.LinearAlgebra.SolveLinearSystems tmpD tmpH 
                                    tmp' * result.TraceA                                          
                                Vector.init (tmp.Length+2) (fun i -> if i>0 && i<=tmp.Length then tmp.[i-1] else 0.0)
                
                            let aic = getAIC (float con) result.TraceA y W
                            let splineFunction = initEvalAt x result.TraceA traceC
                            A,TempClassResult.Create None None None None result.TraceA traceC result.Error result.GCV result.Lambda result.CTemp aic splineFunction (Some (models.[resultindex] |> List.map (fun g -> g.Lambda,g.GCV)|> Array.ofList) )
            //|> fun (x,afinvar,lambdafin,c) -> ((x |> Array.unzip),afinvar,lambdafin,c)
            //|> fun tmp -> 
            //    let e =         tmp |> fun (a,b,c,d) -> fst a
            //    let gcvmin =    tmp |> fun (a,b,c,d) -> snd a
            //    let afin =      tmp |> fun (a,b,c,d) -> b |> Array.map fst
            //    let varfin =    tmp |> fun (a,b,c,d) -> b |> Array.map snd
            //    let lambdafin = tmp |> fun (a,b,c,d) -> c
            //    let cfin =      tmp |> fun (a,b,c,d) -> d
            
            //    let eminIndex = e |> Array.findIndex (fun x -> x = (e |> Array.min))
            //    let efinal = e.[eminIndex]
            //    let afinal = afin.[eminIndex]
            //    let gcvfinal = gcvmin.[eminIndex]
            //    let lambdafinal = lambdafin.[eminIndex]
            //    let cfinal = cfin.[eminIndex]
            //    let varfinal = varfin.[eminIndex]
            //    createTempClassResultSimple afinal cfinal efinal gcvfinal lambdafinal (A.[eminIndex]) varfinal

        /// <summary>calculate the constrained smoothing spline with given constraintmatrix and given lambda</summary>
        /// <param name="ctemp">constraint matrix to enforce specific curve configuration (shape)</param>
        /// <param name="x">vector of x values</param>
        /// <param name="y">vector of actual y values</param>
        /// <param name="W">Matrix with weightings encoded as diagonal</param>
        /// <param name="lambd">spline smoothing strength</param>
        /// <returns>spline solution that minimizes the mGCV</returns>
        let splineManual ctemp (x:Vector<float>) (y:Vector<float>) (W:Matrix<float>) lambd =

            /// Sets the diagonal to value inplace
            let setDiagonalInplace value (m:Matrix<_>) =
                let min = min m.NumRows m.NumCols
                for i=0 to min-1 do
                    m.[i,i] <- value
                m

            let calcGlambda  (D:Matrix<float>) (H:Matrix<float>) (W:Matrix<float>) (n:int) lambd = 
                let dInverse = Algebra.LinearAlgebra.Inverse D
                (2.0 * ( (H.Transpose * dInverse * H)   +  ((lambd / float n)) * (W.Transpose * W)    ) )           

            let calcclambda (W:Matrix<float>) (n:int) (y:Vector<float>) lambd =
                //-2.0 * ((float lambd) / (float n)) * (W.Transpose * W) * y
                -2.0 * ((float lambd) / (float n)) * y.Transpose * (W.Transpose * W) 


            let getError (y:Vector<float>) (a:Vector<float>) (W:Matrix<float>) =
                let tmp =  W * (y - a)
                tmp |> Seq.averageBy (fun x -> x**2.) |> sqrt

        
            let n = x.Length

            //define intervals stepwidth
            let h = Array.init (n-1) (fun i -> x.[i+1] - x.[i] )
            
            //generation of matrices D and H (Wood)
            let H = Array2D.zeroCreate (n-2) n
            let D = Array2D.zeroCreate (n-2) (n-2)

            for i = 1 to (n-2) do
                H.[i-1,i-1] <-  1.0/h.[i-1]
                H.[i-1,i]   <- -( (1.0/h.[i-1]) + (1.0/h.[i]) )
                H.[i-1,i+1] <-  1.0/h.[i]
                D.[i-1,i-1] <-  (h.[i-1] + h.[i]) / 3.0

            for i = 1 to (n-3) do
                D.[i-1,i]   <- h.[i]/6.0
                D.[i,i-1] <- h.[i]/6.0

            let A = [|ctemp|]
                
            A
            |> Array.map (fun aMat ->
                [lambd]
                |> List.mapi (fun z lambd ->
                    let Q = calcGlambda (Matrix.ofArray2D D) (Matrix.ofArray2D H) W n lambd 
                    let c = calcclambda W n y lambd |> RowVector.toArray
                    let b = Array.create ((A.[0] |> Array2D.length1)) (-infinity,0.)
                    let a' = 
                        let tmp = Matrix.create Q.NumRows Q.NumCols 1. |> setDiagonalInplace 0.5
                        let Q' = (Q .* tmp) |> Matrix.toArray2D
                        QP.minimize aMat b Q' c |> Vector.ofArray
                    let e' = getError y a' W
                    let c' =         
                        let tmpH = Matrix.ofArray2D H
                        let tmpD = Matrix.ofArray2D D
                        let tmp = 
                            let tmp' = Algebra.LinearAlgebra.SolveLinearSystems tmpD tmpH
                            tmp' * a'
                        Vector.init (tmp.Length+2) (fun i -> if i>0 && i<=tmp.Length then tmp.[i-1] else 0.0)
                    let (mgcv,gcv) = getGCV (aMat |> Matrix.ofArray2D) (Matrix.ofArray2D D) (Matrix.ofArray2D H) W y a' lambd
                    a',e',mgcv,c' 
                            )
                |> Aux.unzip4
                |> fun (a,e,gcv,c) -> 
                    let gcvminIndex = gcv |> Array.findIndex (fun x -> x = (gcv |> Array.min))
                    let gcvmin = gcv.[gcvminIndex]
                    let efin = e.[gcvminIndex]
                    let afin = a.[gcvminIndex]
                    let cfin = c.[gcvminIndex]
                    let lambdafin = lambd
                    ((efin,gcvmin),afin,lambdafin,cfin)
                        )
            |> Aux.unzip4
            |> fun (x,afin,lambdafin,c) -> ((x |> Array.unzip),afin,lambdafin,c)
            |> fun tmp -> 
                let e = tmp |> fun (a,b,c,d) -> fst a
                let gcvmin = tmp |> fun (a,b,c,d) -> snd a
                let afin = tmp |> fun (a,b,c,d) -> b
                let lambdafin = tmp |> fun (a,b,c,d) -> c
                let cfin = tmp |> fun (a,b,c,d) -> d
            
                let eminIndex = e |> Array.findIndex (fun x -> x = (e |> Array.min))
                let efinal = e.[eminIndex]
                let afinal = afin.[eminIndex]
                let gcvfinal = gcvmin.[eminIndex]
                let lambdafinal = lambdafin.[eminIndex]
                let cfinal = cfin.[eminIndex]
                let splineFunction = initEvalAt x afinal cfinal
                TempClassResult.Create (Some x) None None None  afinal cfinal gcvfinal efinal lambdafinal ctemp infinity splineFunction None

        /// <summary>calculates a constrained, initially increasing spline</summary>
        /// <param name="x">vector of x values</param>
        /// <param name="y">vector of actual y values</param>
        /// <param name="W">Matrix with weightings encoded as diagonal</param>
        /// <param name="con">number of extrema to be enforced in the spline</param>
        /// <returns>spline solution that optimizes the smoothing strength lambda by minimizing mGCV</returns>
        let splineIncreasing (x:Vector<float>) (y:Vector<float>) (W:Matrix<float>) con = 
            spline (~-) x y W con
        
        /// <summary>calculates a constrained, initially decreasing spline</summary>
        /// <param name="x">vector of x values</param>
        /// <param name="y">vector of actual y values</param>
        /// <param name="W">Matrix with weightings encoded as diagonal</param>
        /// <param name="con">number of extrema to be enforced in the spline</param>
        /// <returns>spline solution that optimizes the smoothing strength lambda by minimizing mGCV</returns>
        let splineDecreasing (x:Vector<float>) (y:Vector<float>) (W:Matrix<float>) con = 
            spline (~+) x y W con

        /// <summary>Returns the best fit, selected from 0 - 4 present extrema</summary>
        /// <param name="x_values">vector of x values</param>
        /// <param name="y_values">vector of actual y values</param>
        /// <param name="wMat">Matrix with weightings encoded as diagonal</param>
        /// <param name="minimizer">either Minimizer.AICc or Minimizer.GCV</param>
        /// <returns>spline solution that optimizes the smoothing strength lambda by minimizing mGCV
        /// and optimizes the curve configuration (shape) by minimizing the defined minimizer.</returns>
        let getBestFitOfWeighting x_values y_values wMat minimizer = 
            let getinitialestimate = getInitialEstimateOfWeighting wMat (vector y_values) x_values
            //monotone splines
            let fst1 =(In0,splineIncreasing x_values (vector y_values) wMat 0) |> fun (shapeClass,(constraintMatrices,result)) -> shapeClass,result 
            let fst2 =(De0,splineDecreasing x_values (vector y_values) wMat 0) |> fun (shapeClass,(constraintMatrices,result)) -> shapeClass,result 

            //splines with one extremum
            let snd1 =(In1,splineIncreasing x_values (vector y_values) wMat 1) |> fun (shapeClass,(constraintMatrices,result)) -> shapeClass,result 
            let snd2 =(De1,splineDecreasing x_values (vector y_values) wMat 1) |> fun (shapeClass,(constraintMatrices,result)) -> shapeClass,result 

            //splines with two extrema                                                                                                              
            let trd1 =(In2,splineIncreasing x_values (vector y_values) wMat 2) |> fun (shapeClass,(constraintMatrices,result)) -> shapeClass,result 
            let trd2 =(De2,splineDecreasing x_values (vector y_values) wMat 2) |> fun (shapeClass,(constraintMatrices,result)) -> shapeClass,result 

            //splines with three extremum
            let qua1 =(In3,splineIncreasing x_values (vector y_values) wMat 3) |> fun (shapeClass,(constraintMatrices,result)) -> shapeClass,result 
            let qua2 =(De3,splineDecreasing x_values (vector y_values) wMat 3) |> fun (shapeClass,(constraintMatrices,result)) -> shapeClass,result 

            //splines with four extrema                                                                                                             
            let qui1 =(In4,splineIncreasing x_values (vector y_values) wMat 4) |> fun (shapeClass,(constraintMatrices,result)) -> shapeClass,result 
            let qui2 =(De4,splineDecreasing x_values (vector y_values) wMat 4) |> fun (shapeClass,(constraintMatrices,result)) -> shapeClass,result 

            match minimizer with
            | Minimizer.AICc -> 
                let initialSelectionCriterion = getinitialestimate.AICc
            
                let bestfst = [fst1;fst2] |> List.minBy (fun (cl,x) -> x.AICc)                                                                
                let bestsnd = [snd1;snd2] |> List.minBy (fun (cl,x) -> x.AICc)                                                               
                let besttrd = [trd1;trd2] |> List.minBy (fun (cl,x) -> x.AICc)
                let bestqua = [qua1;qua2] |> List.minBy (fun (cl,x) -> x.AICc)                                                               
                let bestqui = [qui1;qui2] |> List.minBy (fun (cl,x) -> x.AICc)
                let models = 
                    let extract ((c,r) : (Condition*TempClassResult)) = 
                        c.ToString(),r.AICc
                    [|("unconstrained",initialSelectionCriterion);extract fst1;extract fst2;extract snd1;extract snd2;extract trd1;extract trd2;extract qua1;extract qua2;extract qui1;extract qui2|]
            
                //selection of optimal shape possibility by model selection via AICc 
                [bestfst;bestsnd;besttrd;bestqua;bestqui]
                |> List.indexed 
                |> List.minBy (fun (i,(cl,result)) -> (*1.05**(float i) **) result.AICc)
                |> fun (i,(cl,result)) -> 
                    if result.AICc < initialSelectionCriterion then 
                        (cl,result,models) 
                    else Complex,result,models
            | Minimizer.GCV -> 
                let initialSelectionCriterion = getinitialestimate.GCV             
                let bestfst = [fst1;fst2] |> List.minBy (fun (cl,x) -> x.GCV)                                                                
                let bestsnd = [snd1;snd2] |> List.minBy (fun (cl,x) -> x.GCV)                                                               
                let besttrd = [trd1;trd2] |> List.minBy (fun (cl,x) -> x.GCV)
                let bestqua = [qua1;qua2] |> List.minBy (fun (cl,x) -> x.GCV)                                                               
                let bestqui = [qui1;qui2] |> List.minBy (fun (cl,x) -> x.GCV)
                let models = 
                    let extract ((c,r) : (Condition*TempClassResult)) = 
                        c.ToString(),r.GCV
                    [|("unconstrained",initialSelectionCriterion);extract fst1;extract fst2;extract snd1;extract snd2;extract trd1;extract trd2;extract qua1;extract qua2;extract qui1;extract qui2|]
                //selection of optimal shape possibility by model selection via AICc 
                [bestfst;bestsnd;besttrd;bestqua;bestqui]
                |> List.indexed 
                |> List.minBy (fun (i,(cl,result)) -> result.GCV)
                |> fun (i,(cl,result)) -> 
                    if result.GCV < initialSelectionCriterion then 
                        (cl,result,models) 
                    else Complex,result,models          

        [<Obsolete("Only applicable at equal x spacing. Use initEvalAt instead")>]
        //same as initEvalAt, but with recalculated polynomial coefficients
        let initFunctionWithCoefficients (x:Vector<float>) (a:Vector<float>) (c:Vector<float>) =  
            let n = x.Length
            //define interval stepwidth
            let h = Vector.init (n-1) (fun i -> x.[i+1] - x.[i] )

            let deriv i (a:Vector<float>) (c:Vector<float>) xV=
                let tmpT = (xV - x.[i]) / h.[i]
                let pa = -1./6. * c.[0] + 1./6. * c.[1]
                let pb = 0.5 * c.[0] 
                let pc = -1./3. * c.[0] - 1./6. * c.[1] - a.[0] + a.[1]
                let pd = a.[0]
                pa * tmpT**3. + pb * tmpT**2. + pc * tmpT + pd

            (fun t ->
                let i =
                    if t = Seq.last x then 
                        Seq.length x - 2
                    else
                        x
                        |> Seq.findIndex(fun x_Knot -> (x_Knot - t) > 0.)
                        |> fun nextInterval -> nextInterval - 1
                deriv i a.[i .. i+1] c.[i .. i+1] t
                )
        
        (*
        How to determine coefficients?
        f(x)    =ax^3+bx^2+cx+d //traceA
        f'(x)   =3ax^2+2bx+c
        f''(x)  =6ax+2b         //traceC
        - at knots x is 0 or 1 (begin and end of interval)
        - c0 = 2b
        - c1 = 6a+2b
        - a0 = d
        - a1 = a+b+c+d
        solve system by substitution and you get the coefficients for this interval
        *)

        /// <summary>Takes raw data, weighting, and minimizer and returns the best fit, selected from 0 - 4 present extrema</summary>
        /// <param name="xVal">vector of x values</param>
        /// <param name="yVal">vector of actual y value samples</param>
        /// <param name="weightingMethod">y value weighting method (e.g. 1/var)</param>
        /// <param name="minimizer">either Minimizer.AICc or Minimizer.GCV</param>
        /// <returns>spline solution that optimizes the smoothing strength lambda by minimizing mGCV
        /// and optimizes the curve configuration (shape) by minimizing the defined minimizer.</returns>
        let getBestFit xVal (yVal: float [][]) weightingMethod minimizer = 
            let yValMeans = yVal |> Seq.map Seq.mean |> vector
            let weightingMatrix = (getWeighting xVal yVal weightingMethod)
            let (cl,fit,models) = getBestFitOfWeighting xVal yValMeans weightingMatrix minimizer
            {fit with XValues = Some xVal;YValues = Some yVal;WeightingMethod = Some weightingMethod;Minimizer = Some minimizer},models

        /// <summary>Takes raw data, weighting, and minimizer and returns the best fit, selected from 0 - 4 present extrema</summary>
        /// <param name="xVal">vector of x values</param>
        /// <param name="yVal">vector of actual y values</param>
        /// <param name="stdVal">vector of y value standard deviations</param>
        /// <param name="weightingMethod">y value weighting method (e.g. 1/var)</param>
        /// <param name="minimizer">either Minimizer.AICc or Minimizer.GCV</param>
        /// <returns>spline solution that optimizes the smoothing strength lambda by minimizing mGCV
        /// and optimizes the curve configuration (shape) by minimizing the defined minimizer.</returns>
        let getBestFitOfStd (xVal:Vector<float>) (yVal:Vector<float>) (stdVal:Vector<float>) weightingMethod minimizer= 
            let weightingMatrix = 
                //(getWeighting xVal (yVal |> vector) weightingMethod repNumber)
                let weigths =   
                    match weightingMethod with
                    | WeightingMethod.Equal ->
                        Vector.init xVal.Length (fun x -> 1.)
                    | WeightingMethod.StandardDeviation -> 
                        Vector.init xVal.Length (fun i -> 1. / stdVal.[i])
                    | WeightingMethod.StandardDeviationSqrt -> 
                        Vector.init xVal.Length (fun i -> sqrt (1. / stdVal.[i]))                   
                    | WeightingMethod.StandardDeviationAdj -> 
                        let w = Vector.init xVal.Length (fun i -> sqrt ( 1. / (stdVal.[i] / yVal.[i]) ))
                        let mean = Seq.mean w
                        w |> Vector.map (fun e -> e / mean)
                    | WeightingMethod.StandardDeviationAdjSqrt -> 
                        let w = Vector.init xVal.Length (fun i -> sqrt( sqrt ( 1. / (stdVal.[i] / yVal.[i]))))
                        let mean = Seq.mean w
                        w |> Vector.map (fun e -> e / mean)
                    | _ -> failwithf "not implemented yet"
                Matrix.init xVal.Length xVal.Length (fun i j -> if i=j then weigths.[i] else 0.)
            let (cl,fit,models) = getBestFitOfWeighting xVal yVal weightingMatrix minimizer
            fit,models

        ///helper function for initWithLinearInterpol
        let private leftSegmentIdx arr value = 
            let idx = 
                let tmp = Array.BinarySearch(arr, value)
                let idx = if tmp < 0 then ~~~tmp-1 else tmp
                idx
            idx 

        /// <summary>Takes spline data and returns prediciton function, that predicts with the last known 
        /// slope and curvature of 0 in ranges outside the x value coverage</summary>
        /// <param name="x">vector of x values</param>
        /// <param name="a">vector of predicted y values</param>
        /// <param name="c">vector of predicted curvature at positions of x</param>
        /// <returns>prediction function that takes x and returns y_pred</returns>
        let initEvalAtWithLinearInterpol (x:Vector<float>) (a:Vector<float>) (c:Vector<float>) =
            let splineFunction  = initEvalAt x a c
            let splineFunction' = initEvalFstDerivativeAt x a c
            let linearInterPolLowerBorder t =
                let yAtLeftMostKnot = splineFunction x.[0]
                let m = splineFunction' x.[0]
                let b = yAtLeftMostKnot - m * x.InternalValues.[0]
                m * t + b
            let linearInterPolUpperBorder t =
                let yAtRightMostKnot = splineFunction (Seq.last x)
                let m = splineFunction' (Seq.last x)
                let b = yAtRightMostKnot - (m * x.InternalValues.[x.InternalValues.Length-1])
                m * t + b           
            (fun t ->
                    match leftSegmentIdx x.InternalValues t with 
                    | idx when idx < 0 -> 
                        linearInterPolLowerBorder t 
                    | idx when idx > (x.InternalValues.Length-2) -> 
                        linearInterPolUpperBorder t
                    | idx -> initEvalAt x a c t 
                )   
                
    /// <summary>Module to identify classes based on curve configuration</summary>
    module Classification =

        /// <summary>Takes spline data and returns prediciton function, that predicts with the last known 
        /// slope and curvature of 0 in ranges outside the x value coverage</summary>
        /// <param name="xValues">vector of x values</param>
        /// <param name="traceA">vector of predicted y values</param>
        /// <param name="traceC">vector of predicted curvature at positions of x</param>
        /// <param name="linearThreshold">if the range of the predicted values at the knots is lower than this
        /// value, the curge is expected to be constant.</param>
        /// <param name="sensitivity">
        /// if no extrema exist, the classification is based on curvature properties. 
        /// 1: all curvature changes are interpreted for classification
        /// </param>
        /// <returns>Classification string in the form of e.g. Max2.00,Min3.00,Max5.00</returns>
        let getClassification xValues traceA traceC linearThreshold sensitivity = 
            let extrema = Fitting.getExtrema xValues traceA traceC
        
            // if the absolute maximal change of the signal is below the threshold, the signal is assumed constant
            if Seq.max traceA - Seq.min traceA < linearThreshold then //0.05
                "constant"
            else
                let signalFu   = Fitting.initEvalAt xValues traceA traceC 
                let fstDerivFu = Fitting.initEvalFstDerivativeAt xValues traceA traceC 
                let sndDerivFu = Fitting.initEvalSndDerivativeAt xValues traceA traceC
            
                // if extrema are present, the classification is only based on extrema, if no are present, perform a subclassification based on maximal/minimal curvature
                if extrema = [] then 
                
                    // check if there is a plateau by dividing the minimal absolute slope by its standard deviation.
                    // if the result is lower than 2. a plateau is assumed
                    let plateau = 
                        let fstDer =
                            [(Seq.head xValues) + 0.01 .. 0.01 .. (Seq.last xValues) - 0.01]
                            |> List.map fstDerivFu
                        let stdOfFstDerivative = Seq.stDevPopulation fstDer
                        let closestToZero = 
                            fstDer 
                            |> Seq.minBy (fun fd -> Math.Abs fd ) 
                            |> Math.Abs
                        closestToZero / stdOfFstDerivative < 2.
                    // if there is no plateau, no further subclassification is possible 
                    if not plateau then
                        if traceA.[0] < (Seq.last traceA) then "I" else "D"
                    else
                        //get Array of secondderivatives with reasonably spacing (xValues and halfsteps)
                        let sndDerivativeValues =
                            [0..xValues.Length-2]
                            |> List.map (fun x -> 
                                let xValA = xValues.[x]
                                let xValB = (xValues.[x] + xValues.[x+1]) / 2.
                                [
                                xValA,sndDerivFu xValA
                                xValB,sndDerivFu xValB
                                ]
                                )
                            |> List.concat
                        let pseudoExtrema =
                            let processPseudoExtrema (pseudoExtrema :(float*float) list) =
                                pseudoExtrema
                                |> List.rev 
                                |> List.filter (fun (a,b) -> a <> 0.) 
                                // round to nearest xValue and afterwards round the xValue to 2 decimals for printing
                                |> List.map (fun (xVal,yVal) -> round 2 (Aux.roundToNext xVal xValues))
                            //let startFstDerivative = Core.Operators.sign (fstDerivFu xValues.[1])
                            let startSndDerivative = Core.Operators.sign (sndDerivFu xValues.[1])
                            
                            // walk through second derivative and identify all local extrema
                            let rec loop sign pMa pMi accMa accMi i = 
                                if i = sndDerivativeValues.Length - 1 then 
                                    // add last potential PseudoExtremum to respective list
                                    // first true PseudoExtremum indicates direction
                                    let (maFin,miFin) =
                                        if sign = 1 then  
                                            let pseudoMaxima = pMa::accMa |> processPseudoExtrema
                                            let pseudoMinima = accMi      |> processPseudoExtrema
                                            pseudoMaxima,pseudoMinima
                                        else                                                                                                                
                                            let pseudoMaxima = accMa      |> processPseudoExtrema
                                            let pseudoMinima = pMi::accMi |> processPseudoExtrema
                                            pseudoMaxima,pseudoMinima
                                    
                                    // if every pseudo extremum is below the sensitivity threshold classify by "I" and "D"
                                    if maFin = [] && miFin = [] then
                                        ""
                                    //if sensitivity is maximal, take all local extrema for sub classification
                                    elif sensitivity = 1. then
                                        sprintf "%A,%A" maFin miFin
                                    else 
                                        // only one kind of pseudo maximum is reasonable in monoton signals. 
                                        let firstPseudoMaximum = if maFin = [] then infinity else maFin.[0]
                                    
                                        let firstPseudoMinimum = if miFin = [] then infinity else miFin.[0]
                                        
                                        let firstPseudoExtremum = min firstPseudoMaximum firstPseudoMinimum
                                        
                                        let slopeAtFstPseudoExtremum = fstDerivFu firstPseudoExtremum

                                        if firstPseudoMaximum < firstPseudoMinimum then 
                                            //sprintf "%A,%A" maFin []  
                                            if Core.Operators.sign slopeAtFstPseudoExtremum = -1 then 
                                                sprintf "%A,%A" maFin []    
                                            else sprintf "%A,%A" [] miFin 
                                        else 
                                            if Core.Operators.sign slopeAtFstPseudoExtremum = 1 then 
                                                sprintf "%A,%A" [] miFin   
                                            else sprintf "%A,%A" maFin [] 
                                else 
                                    let (_,prevY)   = sndDerivativeValues.[i-1] 
                                    let (tmpX,tmpY) = sndDerivativeValues.[i] 
                                    let (_,nextY)   = sndDerivativeValues.[i+1]
                                    let newSign    = Core.Operators.sign tmpY
                                    // was there a xAxis crossing?
                                    let signChange = newSign <> sign

                                    // if there is a sign change in the second derivative add the potential extremum to the varified extrema
                                    let newAccMa,newAccMi = 
                                        if signChange then 
                                            if sign = 1 then pMa::accMa,accMi //sign = 0?
                                            else accMa,pMi::accMi
                                        else accMa,accMi

                                    // update potential minimum
                                    let newPMi = 
                                        let maxIntensityByMaxCurvature =
                                            let maximalAmplitudeDifference = (traceA |> Seq.max) - (traceA |> Seq.min)
                                            let maximalCurvature = traceC |> Seq.map Math.Abs |> Seq.max
                                            maximalAmplitudeDifference/maximalCurvature
                                        let curvatureBySlope = Math.Abs(tmpY / fstDerivFu tmpX) 
                                        //if curv = 0 and slope = 0 then ratio is high and pseudomaxima are valid, even if there is none
                                        let reasonableCurvature = Math.Abs(tmpY) > 0.00001
                                        // amplitude minimum, sign remains negative, pMi amplitude should be greater than current amplitude, absolute amplitude must by greater than threshold (0.05)
                                        if reasonableCurvature && prevY > tmpY && tmpY < nextY && newSign = -1 && (snd pMi) > tmpY && curvatureBySlope > 10.*(1.-sensitivity) //Math.Abs(tmpY) > (10.*(1.-sensitivity)**10.) //0.05 //sensitivity
                                            then (tmpX,tmpY) 
                                        else 
                                            if signChange then (0.,0.) else pMi
                                
                                    // update potential maximum
                                    let newPMa = 
                                        
                                        let curvatureBySlope = Math.Abs(tmpY / fstDerivFu tmpX) 
                                        //if curv = 0 and slope = 0 then ratio is big and pseudomaxima are valid, even if there is none
                                        let reasonableCurvature = Math.Abs(tmpY) > 0.00001
                                        // amplitude maximum, sign remains positive, pMa amplitude should be lower than current amplitude, absolute amplitude must by greater than threshold (0.05)
                                        if reasonableCurvature && prevY < tmpY && tmpY > nextY && newSign =  1 && (snd pMa) < tmpY && curvatureBySlope > 10.*(1.-sensitivity) //Math.Abs(tmpY) > (10.*(1.-sensitivity)**10.) //0.05 //sensitivity
                                            then (tmpX,tmpY) 
                                        else
                                            if signChange then (0.,0.) else pMa

                                    loop newSign newPMa newPMi newAccMa newAccMi (i+1)
                            loop startSndDerivative (0.,0.) (0.,0.) [] [] 1
                        

                        // add the general descriptor, if trace is in- or decreasing
                        let pseudoExtremaDescriptor = 
                            if pseudoExtrema = "[],[]" then "" else pseudoExtrema
                        if traceA.[0] > Seq.last traceA then "D" + pseudoExtremaDescriptor
                        else "I" + pseudoExtremaDescriptor
                else 
                    extrema 
                    |> List.map (fun ex -> {ex with Xvalue = Aux.roundToNext ex.Xvalue xValues})
                    |> Fitting.extremaToString
               
    /// <summary>Module to visually inspect the fitting and classification result</summary>
    module Vis =

        open Plotly.NET

        /// <summary>Takes a spline and returns a visualization of the spline and quality parameters</summary>
        /// <param name="tc">Fitting result</param>
        /// <param name="modelQualityScores">collection of minimizer values (e.g. ["In0",12.;"In1",2.])</param>
        /// <returns>Plotly chart that summarizes the fit properties</returns>
        let getChart (tc: Fitting.TempClassResult) (modelQualityScores: ((string*float)[]) option) = 
            let min,max =
                Seq.min tc.XValues.Value,
                Seq.max tc.XValues.Value
                
            let transformClassNames str =
                match str with 
                | "unconstrained" -> "unconstrained"
                | "In0" -> "monotonically increasing"
                | "De0" -> "monotonically decreasing"
                | "In1" -> "1 maximum"
                | "De1" -> "1 minimum"
                | "In2" -> "max + min"
                | "De2" -> "min + max"
                | "In3" -> "max + min + max"
                | "De3" -> "min + max + min"
                | "In4" -> "max + min + max + min"
                | "De4" -> "min + max + min + max"
                | _ -> "misc"

            let classification = Classification.getClassification tc.XValues.Value tc.TraceA tc.TraceC 0.05 1.

            let visualizationSpline = 
                [min .. 0.05 .. max]
                |> List.map (fun x -> 
                    x,tc.SplineFunction x
                    )
                |> Chart.Line
                |> Chart.withTraceInfo "constrained spline"

            let visualizationRawSignal = 
                tc.YValues.Value 
                |> Array.mapi (fun i x -> 
                    x |> Array.map (fun a -> tc.XValues.Value.[i],a)
                    )
                |> Array.concat
                |> Chart.Point
                |> Chart.withTraceInfo "raw signal"

            let signalChart =
                [
                    visualizationRawSignal
                    visualizationSpline
                ]
                |> Chart.combine
                |> Chart.withXAxisStyle "time point" 
                |> Chart.withYAxisStyle "intensity"

            let qualityParameter = 
                if modelQualityScores.IsSome then 
                    let columnvalues = 
                        modelQualityScores.Value
                        |> Array.map (fun (a,b) -> transformClassNames a, b)
                    
                    Chart.Column(columnvalues,MarkerColor=Color.fromString "grey",Name="qc") |> Chart.withXAxisStyle "" |> Chart.withYAxisStyle "quality score"
                else Chart.Invisible()

            let weightingChart = 
                let weightings = 
                    Fitting.getWeighting tc.XValues.Value tc.YValues.Value tc.WeightingMethod.Value
                    |> Matrix.getDiag
                    |> Vector.toArray
                    |> Array.zip (Vector.toArray tc.XValues.Value)
                Chart.Column(weightings,Name = "point weightings")
                |> Chart.withXAxisStyle "time point" 
                |> Chart.withYAxisStyle "weighting"

            [signalChart;weightingChart;qualityParameter]
            |> Chart.Grid(3,1)
            |> Chart.withTemplate ChartTemplates.lightMirrored
            |> Chart.withTitle classification
            |> Chart.withSize(600.,900.)
            |> Chart.withConfig(Config.init (ToImageButtonOptions = ConfigObjects.ToImageButtonOptions.init(Format = StyleParam.ImageFormat.SVG)))
            |> Chart.withMarginSize(Bottom=200)




