(* ::Package:: *)

  (* Input parameters                                  *)
  k = 8          (* order of approximating polynomial   *)  
  q = 2 k        (* precision q                         *)
  steps = 4 k    (* steps = number of data points       *) 
  
   (* data generation *)
  f[y_] := N[Gamma[1+y],q];
  
   (* model components as,ys *)
  as = (Symbol["a"<>ToString[#]]&)/@ Range[k];
  ys = (y^#&) /@ Range[k];
  model[y_] = 1 + Total[as ys];
  
   (* solve linear system of equations for parmameters  *)
  getCoeffs[] :=  
    Module[{xs, array, m, b, sol, parms},
    xs = Range[1/k, 1, 1/k]; (* list of equidistant xs  *)
    array = model[#]& /@ xs;
    m = Coefficient[#,as]& /@ array;  (* build matrix m *)     
    b = f[#]-1& /@ xs;          (* a0 = 1 => right side *)
    sol = LinearSolve[m, b];    (* solve the system     *) 
    parms = Thread[as->sol];   (* list of  parameters   *) 
    Return[parms];
  ];
  
   (*get error function *)
  getError[] :=
    Module[{maxErr, ers},
    coeffs = getCoeffs[];      (* now global variable  *)
    ers= {#,Log10[Abs[1-model[#]/f[#]]]} /. coeffs & /@ 
        Range[0,1,1/steps];    (* list of errors       *) 
    Return[ers];
  ];
  
   (* start code execution                       *)
  errs = getError[];
  ListPlot[ errs, Joined->True ]
  Print[coeffs //TableForm]; 
