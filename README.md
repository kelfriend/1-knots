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
gap> for i in [0..3] do Print(Homology(comp,i),"\n"); od;
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
### Complements/chainmapbound.g
**ChainMapOfKnotBoundaryToComplement** inputs the inclusion map from **KnotComplementWithBoundary** and outputs a chain map. This function uses existing HAP methods to compute a chain map from the chain complex of the universal covers of each of the two spaces. The final output is a chain map between the two complexes after having been tensored with the integers over some finite index subgroup of the knot group. Note that a specific subgroup may be used, but this is not implemented.
This allows for computation of (co)homology with local coeffeicients. Please see [this preprint](http://hamilton.nuigalway.ie/preprints/LocalCoho.pdf) for a more detailed treatment of this (and the above) processes. Below, I will compute this chain map for a 4-fold cover of the granny knot and then I will obtain group homomorphisms associated to the various homology groups.
```
gap> tre:=PureCubicalKnot(3,1);
prime knot 1 with 3 crossings

gap> granny:=ArcPresentation(KnotSum(tre,tre));
[ [ 7, 10 ], [ 6, 8 ], [ 7, 9 ], [ 8, 10 ], [ 2, 6 ], [ 5, 9 ], [ 1, 3 ], 
  [ 2, 4 ], [ 3, 5 ], [ 1, 4 ] ]
gap> i:=KnotComplementWithBoundary(granny);
Map of regular CW-complexes

gap> iota:=ChainMapOfKnotBoundaryToComplement(i,4);
Chain Map between complexes of length 3 . 

gap> Homology(iota,0);
[ g1, g2, g3 ] -> [ g1, g1, g1 ]
gap> Homology(iota,1);
[ g1, g2, g3, g4, g5, g6 ] -> [ g1^2*g2^2, g1^2*g2^-1*g3^2, g2^2*g3^2, g1^3*g2^2, g1^2*g3^2, g2^2*g3^3 ]
gap> Homology(iota,2);
[ g1, g2, g3 ] -> [ g1*g3^-1*g5^-1, g1^-1*g2*g4^-1*g6^-1, g2^-1 ]
```
