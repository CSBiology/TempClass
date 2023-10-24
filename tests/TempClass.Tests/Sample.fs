module Tests

open Expecto
open TempClass
open FSharp.Stats

[<Tests>]
let tests =
    testList "tests" [
        test "can use lapack" {
            FSharp.Stats.ServiceLocator.setEnvironmentPathVariable (@"./")
            FSharp.Stats.Algebra.LinearAlgebra.Service() |> ignore
            Expect.isTrue true "this message will not show on fail because above line failed xd"
        }
    ]
