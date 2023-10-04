#r @"C:\Users\bvenn\source\repos\TemporalClassification_test\TemporalClassification_test\bin\Debug\net472\FSharp.Stats.dll"
#r @"C:\Users\bvenn\source\repos\TemporalClassification_test\TemporalClassification_test\bin\Debug\net472\TemporalClassification_test.dll"

open FSharp.Stats

//FSharp.Stats.ServiceLocator.setEnvironmentPathVariable (@"C:\Users\venn\source\repos\FSharp.StatsMSF030\lib")
//FSharp.Stats.Algebra.LinearAlgebra.Service()

let a = Matrix.init 10 4 (fun i j -> float (i+j))
let b = Matrix.init 4 3 (fun i j -> float (i+j))

a * b

Matrix.mul a b


open TemporalClassification_test
open TemporalClassification

let sigA = (vector [1..5])
let sigB = (vector [1.;1.1;1.2;2.;2.1;1.8;1.4;1.5;1.4;4.;4.;4.1;5.;5.1;5.])
let sigBMeans = sigB |> Seq.chunkBySize 3 |> Seq.map Seq.average |> vector

let weightingMethod = TemporalClassification.Fitting.WeightingMethod.StandardDeviationAdj

let weighting = TemporalClassification.Fitting.getWeighting sigA sigB weightingMethod 3


TemporalClassification.Fitting.getInitialEstimateOfWeighting weighting sigBMeans sigA
TemporalClassification.Fitting.splineDecreasing sigA sigBMeans weighting 1
TemporalClassification.Fitting.splineDecreasing sigA sigBMeans weighting 2
TemporalClassification.Fitting.splineDecreasing sigA sigBMeans weighting 3
TemporalClassification.Fitting.splineDecreasing sigA sigBMeans weighting 1
TemporalClassification.Fitting.splineDecreasing sigA sigBMeans weighting 2
TemporalClassification.Fitting.splineDecreasing sigA sigBMeans weighting 3

Fitting.getBestFit sigA sigB 3 Fitting.WeightingMethod.Equal Fitting.Minimizer.GCV



let getinitialestimate = Fitting.getInitialEstimateOfWeighting weighting (vector sigBMeans) sigA
let fst1 =(Fitting.In0,Fitting.splineIncreasing sigA (vector sigBMeans) weighting 0) |> fun (shapeClass,(constraintMatrices,result)) -> shapeClass,result 
let fst2 =(Fitting.De0,Fitting.splineDecreasing sigA (vector sigBMeans) weighting 0) |> fun (shapeClass,(constraintMatrices,result)) -> shapeClass,result 
let snd1 =(Fitting.In1,Fitting.splineIncreasing sigA (vector sigBMeans) weighting 1) |> fun (shapeClass,(constraintMatrices,result)) -> shapeClass,result 
let snd2 =(Fitting.De1,Fitting.splineDecreasing sigA (vector sigBMeans) weighting 1) |> fun (shapeClass,(constraintMatrices,result)) -> shapeClass,result 



Fitting.getBestFitOfWeighting sigA sigBMeans weighting Fitting.Minimizer.GCV

//record type for easy handling inside the shape classification process
type InnerResult = {
    TraceA  : Vector<float>
    GCV     : float
    Var     : float
    Error   : float
    CTemp   : float [,]
    Lambda  : float}

let createInnerResult a g c e v l = {
    TraceA  = a
    GCV     = g
    Var     = v
    Error   = e
    CTemp   = c
    Lambda  = l}

let spline operator (x:Vector<float>) (y:Vector<float>) (W:Matrix<float>) con =

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
    printfn "Halloho"
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
    
    printfn "Halloho2"
    ///the constraint-matrices for a given shape are transformed into 2D-Arrays (solver compatibility) and the operator defines the monotonicity direction
    let A = 
        ctList 
        |> Array.map (fun x -> 
            operator(x) 
            |> Matrix.toArray2D)
            
    let lamlist = [for i = 0 to 64 do yield 0.01*(1.2**(float i))]

    A
    |> Array.mapi (fun k aMat ->

        lamlist
        |> List.mapi (fun z lambd ->
            let Q = calcGlambda (Matrix.ofArray2D D) (Matrix.ofArray2D H) W n lambd //outsource (performance)
            let c = (calcclambda W n y lambd) |> RowVector.toArray                  //outsource (performance)
            //borders to solve the inequality constraints 
            let b = Array.create ((A.[0] |> Array2D.length1)) (-infinity,0.)
            //solver specific transformation (because symmetry of the matrix)
            
            printfn "Halloho3"
            let a' = 
                let tmp = Matrix.create Q.NumRows Q.NumCols 1. |> setDiagonalInplace 0.5
                let Q' = (Q .* tmp) |> Matrix.toArray2D
                Optimization.QP.minimize aMat b Q' c |> Vector.ofArray
            
            printfn "Halloho4"
            let e' = getError y a' W
            
            printfn "Halloho5 AMatrix: %i | Lambda:%i | lambd: %f | a': %A | D: %A |H: %A" k z lambd a' D H
            printfn "Fitting.getGCV (%A |> Matrix.ofArray2D) (Matrix.ofArray2D %A) (Matrix.ofArray2D %A) %A %i %A %A %f" aMat D H W y.Length y a' lambd
            let (mgcv,gcv) = Fitting.getGCV (aMat |> Matrix.ofArray2D) (Matrix.ofArray2D D) (Matrix.ofArray2D H) W y.Length y a' lambd
            
            printfn "Halloho6"
            createInnerResult a' mgcv aMat e' gcv lambd
                )
        //minimize the gcv in respect to the used smoothing factor lambda
        |> List.minBy (fun innerResult -> innerResult.GCV)
        )
    
    |> fun xt -> 
        printfn "Halloho7"
        //additional case, when no spline satisfies the given conditions (with additional constraints that checks extrema count)
        if xt = [||] then
            //for shape restriction if coefficient shrinkage is to intense (uncomment above)
            A,Fitting.createTempClassResult (Vector.init n (fun _ -> 0.)) (Vector.init n (fun _ -> 0.)) infinity infinity infinity (Array2D.zeroCreate 1 1) infinity id
        else 
            xt
            //among all shape possibilities under a given parent shape ((min, then max);(max);...) minimize the mGCV to obtain the most promising spline candidate
            |> Array.minBy (fun innerResult -> innerResult.GCV)
            |> fun result -> 
                let traceC = 
                    let tmpH = Matrix.ofArray2D H
                    let tmpD = Matrix.ofArray2D D
                    let tmp = 
                        let tmp' = Algebra.LinearAlgebra.SolveLinearSystems tmpD tmpH 
                        tmp' * result.TraceA                                          
                    Vector.init (tmp.Length+2) (fun i -> if i>0 && i<=tmp.Length then tmp.[i-1] else 0.0)
                
                let aic = Fitting.getAIC (float con) result.TraceA y W
                let splineFunction = Fitting.initEvalAt x result.TraceA traceC
                A,Fitting.createTempClassResult result.TraceA traceC result.Error result.GCV result.Lambda result.CTemp aic splineFunction


spline (~-) sigA (vector sigBMeans) weighting 0








let myGetGCV (C:Matrix<float>) (D:Matrix<float>) (H:Matrix<float>) (W:Matrix<float>) (n:int) (y:vector) (a:vector) (lambda:float) =
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

    printfn "Halloe 1.0 %A" Ca
    ///calculate null space of Ca (based on SVD)
    let Z =
        if Ca.NumRows * Ca.NumCols = 0 then 
            Matrix.identity n            
        else
            let nullspace = 
                let (eps,U,V) =
                    Algebra.LinearAlgebra.SVD Ca
                let rank = eps |> Seq.filter (fun x -> x >= 1e-5) |> Seq.length //Threshold if a singular value is considered as 0. //mat.NumRows - eps.Length
                let count = V.NumRows - rank
                Matrix.getRows V (rank) (count)
                |> Matrix.transpose
            nullspace
    
    printfn "Halloe 2.0"
    let I = Matrix.identity n
    
    printfn "Halloe 3.0"
    let G_lambda (lambd:float)  = 
        let dInverse = Algebra.LinearAlgebra.Inverse D
        2. * ( H.Transpose * dInverse * H + lambd / (float n)  * W.Transpose * W)
    
    printfn "Halloe 4.0"
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
        
    printfn "Halloe 5.0"
    let var_ny (lambd:float) =
        let a = ((W * (a-y)) |>Vector.norm ) ** 2.
        let b = I - (A_tild lambd) |> Matrix.getDiag |> Vector.sum
        a/b
    
    printfn "Halloe 6.0"
    //standard gcv
    let gcv  = V_tild 1.0 lambda
    printfn "Halloe 7.0"
    //modified gcv
    let mGcv = V_tild 1.3 lambda
    printfn "Halloe 8.0"
    //robust gcv
    let rGcv() =
        let A = A_tild lambda
        let gamma = 0.2
        let mu = ((A * A) |> Matrix.trace) / (float n)
        (gamma + (1. - gamma) * mu) * gcv

    
    printfn "Halloe 9.0"
    let variance = var_ny lambda
    printfn "Halloe 10.0"
    //(gcv,variance)
    (mGcv,gcv)

myGetGCV
    (matrix
        [[3.0; -3.0; -0.0; -0.0; -0.0]
         [1.267857143; -1.607142857; 0.4285714286; -0.1071428571; 0.01785714286]
         [1.732142857; -1.392857143; -0.4285714286; 0.1071428571; -0.01785714286]
         [2.535714286; -3.214285714; 0.8571428571; -0.2142857143; 0.03571428571]
         [-0.0; 3.0; -3.0; -0.0; -0.0]
         [0.4642857143; 0.2142857143; -0.8571428571; 0.2142857143; -0.03571428571]
         [-0.4642857143; 2.785714286; -2.142857143; -0.2142857143; 0.03571428571]
         [0.125; 2.25; -3.0; 0.75; -0.125]
         [-0.0; -0.0; 3.0; -3.0; -0.0]
         [-0.125; 0.75; -0.0; -0.75; 0.125]
         [0.125; -0.75; 3.0; -2.25; -0.125]
         [-0.03571428571; 0.2142857143; 2.142857143; -2.785714286; 0.4642857143]
         [-0.0; -0.0; -0.0; 3.0; -3.0]
         [0.03571428571; -0.2142857143; 0.8571428571; -0.2142857143; -0.4642857143]
         [-0.03571428571; 0.2142857143; -0.8571428571; 3.214285714; -2.535714286]
         [0.01785714286; -0.1071428571; 0.4285714286; 1.392857143; -1.732142857]
         [-0.01785714286; 0.1071428571; -0.4285714286; 1.607142857; -1.267857143]])
     
     (matrix
        [[0.6666666667; 0.1666666667; 0.0]
         [0.1666666667; 0.6666666667; 0.1666666667]
         [0.0; 0.1666666667; 0.6666666667]]) 
    (matrix 
        [[1.0; -2.0; 1.0; 0.0; 0.0]
         [0.0; 1.0; -2.0; 1.0; 0.0]
         [0.0; 0.0; 1.0; -2.0; 1.0]]) 
    (matrix 
        [[0.560570857; 0.0; 0.0; 0.0; 0.0]
         [0.0; 0.60646446; 0.0; 0.0; 0.0]
         [0.0; 0.0; 0.842147352; 0.0; 0.0]
         [0.0; 0.0; 0.0; 1.412688939; 0.0]
         [0.0; 0.0; 0.0; 0.0; 1.578128392]]) 
    5 
    (vector [|1.1; 1.966666667; 1.433333333; 4.033333333; 5.033333333|]) 
    (vector [|1.162920666; 1.485139503; 1.926132803; 3.899231341; 5.063632497|]) 
    131.046309




matrix  [[0.125; 2.25; -3.0; 0.75; -0.125]] * (vector sigBMeans)



Algebra.LinearAlgebra.SVD (matrix [[0.125; 2.25; -3.0; 0.75; -0.125]])




