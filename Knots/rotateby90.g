LoadPackage("HAP");

RotateBy90:=function(B)
    local n, R, i, j, current;
    
    n:=Length(B);
    R:=ShallowCopy(B);

    for i in [1..Int(n/2)]
        do
        for j in [i..n-i]
            do
            current:=R[i][j];
            R[i][j]:=R[j][n-i+1];
            R[j][n-i+1]:=R[n-i+1][n-j+1];
            R[n-i+1][n-j+1]:=R[n-j+1][i];
            R[n-j+1][i]:=current;
        od;
    od;
    return R;
end;

A:=ReadImageAsPureCubicalComplex("~/Pictures/a.png",1);;
a:=ShallowCopy(A!.binaryArray);;
r:=RotateBy90(a);;
R:=PureCubicalComplex(r);;