(* ::Package:: *)

IRKTruncationError[A_,b_,c_,p_]:=Module[{k,s,y0,y1,exactSeries},
s=Length[A];
y0=1;
exactSeries=Series[Exp[dt],{dt,0,p+1}];
k=LinearSolve[(IdentityMatrix[s]-dt*A),Table[1,s]];
y1=y0+dt*k.b;
(y1-exactSeries)/dt//FullSimplify]
