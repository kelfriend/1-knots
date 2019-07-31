#Read("spin_mod.g");

K:=PureCubicalKnot(3,1);;
A:=K!.binaryArray;;
A[1][2][2]:=1;
A[1][2][11]:=1;
for i in [3..10] do
A[2][2][i]:=0;
od;

B:=0*A;
B[1]:=1*A[1];

Y:=PureCubicalComplex(A);;
U:=PureCubicalComplex(B);;
inc:=HAP_PureCubicalPairToCWMap(Y,U);
S:=Spin(inc);

cY:=PureComplexComplement(Y);
cY!.binaryArray[1][2][2]:=0;
cY!.binaryArray[1][2][11]:=0;
B:=0*A;
B[1]:=B[1]+1;
B[1][2][2]:=0;
B[1][2][11]:=0;
cU:=PureCubicalComplex(B);;
cinc:=HAP_PureCubicalPairToCWMap(cY,cU);
cS:=Spin(cinc);

