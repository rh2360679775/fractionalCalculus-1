(* ::Package:: *)

(* Input parameters                                  *)
k = 19         (*  order of approximating polynomial *)  
q = 2 k        (*  precision q                       *)
steps = 4 k    (*  steps = number of data points     *) 
d = 3/10       (*  order of approximating polynomial *)
s0=18/10       (*  start value s                     *) 
n=3            (* initial number of function calls   *)
nmax =4       (* maximum function calls             *) 

(*data generation*)
f[y_]:=N[Gamma[1+y],q];

(*model components as,ys*)
as=(Symbol["a"<>ToString[#]]&)/@Range[k];
ys=(y^#&)/@Range[k];
model[y_]=1+Total[as ys];

(*shift positions x*)
shift[x_,s_]:=1/2^(1-s)If[x<1/2,x^s,2^(1-s)-(1-x)^s];

(*parabola and minimum*)
parabola[y_]:=c0+c1 y+c2 y^2;
minP:=Rationalize[-c1/(2 c2),10^-q];

getminP[pnts_]:=
Module[{x1m2,x1m3,x2m3,x1p2,x1p3,x2p3,k1,k2,k3,minP},
x1m2 = pnts[[1,1]]-pnts[[2,1]];
x1m3 = pnts[[1,1]]-pnts[[3,1]];
x2m3 = pnts[[2,1]]-pnts[[3,1]];
x1p2 = pnts[[1,1]]+pnts[[2,1]];
x1p3 = pnts[[1,1]]+pnts[[3,1]];
x2p3 = pnts[[2,1]]+pnts[[3,1]];
k1 = pnts[[1,2]] x2m3;
k2 = pnts[[2,2]] x1m3;
k3 = pnts[[3,2]] x1m2;
minP = 1/2(k1 x2p3 - k2 x1p3 + k3 x1p2)/(k1 -k2 +k3);
minP = Rationalize[minP,10^-q];
Return[minP];];

(*solve linear system of equations for coeffs*)
getCoeffs[s_]:=
Module[{xshifted,array,mat,b,sol,parms},
xshifted=shift[#,s] & /@ Range[1/k,1,1/k];
array=model[#] & /@ xshifted;
mat=Coefficient[#,as] & /@ array;
b=f[#]-1 & /@ xshifted;     (*a0=1*)
sol=LinearSolve[mat,b];
parms = Thread[as->sol];
Return[parms];];

(*get maximum error*)
getMaxError[s_]:=
Module[{maxErr,ers},
coeffs=getCoeffs[s];
ers=Log10[Abs[1-model[#]/f[#]]]/. coeffs&/@Range[0,1,1/steps];
maxErr=TakeLargest[ers,1][[1]];
Return[maxErr];];

(* start code execution                       *)
start = AbsoluteTime[];
errsall = Table[{s,getMaxError[s]},{s,s0-d,s0+d,d}];
errs3best = errsall;
Print["errs3best =  ", N[errs3best]];
(* loop                       *)	 
While[ n < nmax,
	pb = FindFit[errs3best,parabola[x],{c0,c1,c2},x];
	snew = minP /. pb;
	snew = Rationalize[snew, 10^-q]; (* for high precision*) 
	Print[N[snew], " f ", N[getminP[errs3best]]];
	point = {snew, getMaxError[snew]};
	errsall=SortBy[Append[errsall,point],First];
	pos=First@Ordering[errsall[[All,2]]];
	errs3best = Take[errsall,{pos-1,pos+1}];
    Print["errs3best =  ", N[errs3best]];
n++];
Print["ersMax =  ", N[10^errs3best[[2,2]]]];

stop =  TimeUsed[];
Print["coeffs =  ", coeffs];

minVal = -15.4
ersMin = ({#[[1]], Max[minVal,#[[2]]]} &) /@ ersX;
xport=StringJoin["d",ToString[k],"s.txt"]
Export[xport,N[ersMin],"Table"]
ListPlot[ ersMin, Joined->True ]






