[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/kelfriend/1-knots/master)# The below functions are included in [HAP](http://hamilton.nuigalway.ie/Hap/www/index.html) as of version 1.22.
### Complements/knotcomp.g
**KnotComplement** will input an arc presentation representing some knot and output a regular CW-decomposition of its complement.
The below code computes the complement of the trefoil, displays the number of cells in the complex and obtains its 0<sup>th</sup>,
1<sup>st</sup> & 2<sup>nd</sup> homology groups.
```
gap> arc:=ArcPresentation(PureCubicalKnot(3,1));
[ [ 2, 5 ], [ 1, 3 ], [ 2, 4 ], [ 3, 5 ], [ 1, 4 ] ]
gap> comp:=KnotComplement(arc);
Regular CW-complex of dimension 3

gap> Size(comp);
395
gap> for i in [0..2] do Print(Homology(comp,i),"\n"); od;
[ 0 ]
[ 0 ]
[ 0 ]
```
An optional argument of **"rand"** will enable random generation of 2-cells. This change results in the same complex, but the
alternate ordering may offer a different (but equivalent) presentation of the fundamental group.
```
gap> pi_1:=FundamentalGroup(comp,1);
#I  there are 2 generators and 1 relator of total length 6
<fp group of size infinity on the generators [ f1, f2 ]>
gap> RelatorsOfFpGroup(pi_1);
[ f1*f2^-1*f1^-1*f2*f1^-1*f2^-1 ]
gap> comp:=KnotComplement(arc,"rand");
Random 2-cell selection is enabled.
Regular CW-complex of dimension 3

gap> pi_1:=FundamentalGroup(comp,1);; RelatorsOfFpGroup(pi_1);
#I  there are 2 generators and 1 relator of total length 6
[ f2^-1*f1*f2^-1*f1^-1*f2*f1^-1 ]
```
### Complements/knotcompbound.g
**KnotComplementWithBoundary** takes an arc presentation as above and returns an inclusion map of regular CW-complexes from the
boundary of a knot to its complement. The below code calculates this inclusion from the 3rd prime knot on 11
crossings. The induced homomorphism from their fundamental groups is then computed.
```
gap> arc:=ArcPresentation(PureCubicalKnot(11,3));
[ [ 3, 9 ], [ 2, 5 ], [ 1, 3 ], [ 7, 10 ], [ 4, 8 ], [ 6, 9 ], [ 7, 11 ], 
  [ 5, 10 ], [ 2, 6 ], [ 4, 11 ], [ 1, 8 ] ]
gap> i:=KnotComplementWithBoundary(arc);
Map of regular CW-complexes

gap> FundamentalGroup(i,1);
[ f1, f2 ] -> 
[ 
  f2*f1*f2^-1*f1*f3^-1*f2^-1*f3*f2*f3*f1^-1*f3*f2*f3^-1*f2^-1*f3^-1*f2*f3*f2^-1\
*f1*f3^-1*f2^-1*f3*f2*f3*f1^-1*f2*(f1*f2^-1*f3^-1)^3*f1*f3^-1*f2^-1*f3^-1*f2*f\
3*f1^-1*f3^-1*f2^-1*(f3*f2)^2*f3*f1^-1, f1*f3^-1*f2^-1*f3^-1*f2*f3*f1^-1 ]
```

#### These algorithms, as well as previously implemented methods for lifting chain maps, allow for computation of (co)homology with local coeffeicients. Please see [this preprint](http://hamilton.nuigalway.ie/preprints/LocalCoho.pdf) (submitted to Mathematics of Computation) for a more detailed treatment of this (and the above) processes, as well as a description of a certain knot invariant which can be obtained from a finite-index covering of a knot complement.
