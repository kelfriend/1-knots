LoadPackage("HAP");

# inputs some 3-dimensional cubical complex/knot, embeds it into 
# 4-space and then spins it about an axis, forming a 
# 4-dimensional complex/knot
####################################################################
SpinCubicalComplex:=function(K)
    
    local 
          array, embedding, rotation_matrix, empty, centre,
          n, theta, i, j, k, l, v, union;

    array:=K!.binaryArray;
    embedding:=PureCubicalComplex([0*array,array,0*array]);
    array:=embedding!.binaryArray;

    rotation_matrix:=function(theta)
        local R;

        R:=[
        [Cos(theta),-Sin(theta),0,0],
        [Sin(theta),Cos(theta),0,0],
        [0,0,1,0],
        [0,0,0,1]
        ];
        return R;
    end;

    empty:=0*array;
    
    centre:=[
    [Int(Length(empty)/2)],
    [Int(Length(empty[1])/2)],
    [Int(Length(empty[1][1])/2)],
    [Int(Length(empty[1][1][1])/2)]
    ];

    union:=[];

    for n in [1..360]
        do
        theta:=2*(3.142/n);
        for i in [1..Length(array)]
            do
            for j in [1..Length(array[1])]
                do
                for k in [1..Length(array[1][1])]
                    do
                    for l in [1..Length(array[1][1][1])]
                        do
                        if array[i][j][k][l]=1
                            then
                            v:=rotation_matrix(theta)*
                            ([[i],[j],[k],[l]]-centre) + centre;

                            v:=[Int(v[1][1]),Int(v[2][1]),
                                Int(v[3][1]),Int(v[4][1])];

                            if v[1]>0 and v[1]<=Length(empty) and 
                            v[2]>0 and v[2]<Length(empty[1]) and
                            v[3]>0 and v[3]<Length(empty[1][1]) and
                            v[4]>0 and v[4]<Length(empty[1][1][1])
                                then
                                empty[v[1]][v[2]][v[3]][v[4]]:=1;
                            fi;
                        fi;
                    od;
                od;
            od;
            if theta=2*3.142
                then
                union:=PureCubicalComplex(empty);
            else
                union:=PureCubicalComplexUnion(
                       PureCubicalComplex(empty),union
                       );
            fi;
        od;
    od;
    return union;
end;
####################################################################


# takes a 2-dimensional array and rotates it by 90 degrees
####################################################################
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
####################################################################
