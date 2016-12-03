(* ::Package:: *)

BeginPackage["RadauIIA`"]


A::usage="A[s] returns the symbolic Radau IIA Butcher array with s stages";
b::usage="b[s] returns the symbolic Radau IIA weights for s stages (given by the bottom row of A[s])"
c::usage="c[s] returns the symbolic Radau IIA abscissa for s stages";
An::usage="An[s] returns the numeric Radau IIA Butcher array with s stages";
bn::usage="bn[s] returns the numeric Radau IIA weights for s stages (given by the bottom row of An[s])"
cn::usage="cn[s] returns the numeric Radau IIA abscissa for s stages";


Begin["Private`"];


(* ::Text:: *)
(*Using the derivation of Owe Axelson, "A class of A-stable methods," BIT (1969)*)


(* ::Text:: *)
(*We define the following orthogonal polynomial (corresponding to a=-1,b=0 in Axelson's paper)*)


q[n_,t_]=LegendreP[n,2t-1]-LegendreP[n-1,2t-1];


(* ::Text:: *)
(*The roots of the polynomial are the quadrature points*)


c[n_]:=(t/.Solve[q[n,t]==0,t])//(Sort[#,Less]&)


(* ::Text:: *)
(*We define the following function (note that there appears to be a minor error in Axelson's paper, where Subscript[x, k]appears instead of Subscript[u, k])*)


l[n_,k_,t_]:=q[n,t]/((t-c[n][[k]])(D[q[n,t],t]/.t->c[n][[k]]))


(* ::Text:: *)
(*and then integrate in order to compute the quadrature weights.*)


a[n_,i_,k_]:=Integrate[l[n,k,t],{t,0,c[n][[i]]}]


(* ::Text:: *)
(*The resulting matrix is the Butcher tableau:*)


A[n_]:=Table[a[n,i,k],{i,1,n},{k,1,n}]


(* ::Text:: *)
(*The weights are given by the bottom row of the Butcher tableau:*)


b[n_]:=A[n][[-1,All]]


(* ::Section:: *)
(*Numerical Computations*)


(* ::Text:: *)
(*For more than 3 stages, the symbolic calculations can become infeasible, so we compute the coefficients numerically.*)


cn[n_]:=(t/.NSolve[q[n,t]==0,t])//Re//(Sort[#,Less]&)
ln[n_,k_,t_]:=q[n,t]/((t-cn[n][[k]])(D[q[n,t],t]/.t->cn[n][[k]]))
an[n_,i_,k_]:=NIntegrate[ln[n,k,t],{t,0,cn[n][[i]]}]
An[n_]:=Table[an[n,i,k],{i,1,n},{k,1,n}]
bn[n_]:=An[n][[-1,All]]
End[];


EndPackage[]
