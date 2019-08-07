#Read("spin_mod.g");
#########################################################################################
######################## The Complement of Spun Knots ###################################
#########################################################################################
################# Input: - a list of integers [n,p] corresponding #######################
######################## to the pth prime knot on n crossings. ##########################
######################## - alternatively, a pure cubical knot. ##########################
#########################################################################################
################ Output: a 5-dimensional regular CW-complex S(K*) #######################
######################## homotopy equivalent to the complement ##########################
######################## of the spinning of K* (a knotted arc formed #################### 
######################## from K by removing an unknotted segment) #######################
######################## about a plane. #################################################
#########################################################################################
SpunKnotComplement:=function(k)

    local
          d, K, DeletedPair, pair, C, U, C2, omicron;

    d:=[0,0,1,1,2,3,7,21,49,165,552];

    if IsList(k)
        then
        if not (IsInt(k[1]) and IsInt(k[2]) and SignInt(k[1])=1 and SignInt(k[2])=1)
            then
            Error("the input must be a pair of positive integers.\n");
        elif k[1]>11
            then
            Error("only knots with less than 12 crossings are stored in HAP.\n");
        elif k[2]>d[k[1]]
            then
            Error("no such prime knot exists.\n");
        fi;
        K:=PureCubicalKnot(k[1],k[2]);
    elif IsPureComplex(k)
        then
        K:=ShallowCopy(k);
    else
        Error("input must be an integer pair or a cubical knot.\n");
    fi;

    DeletedPair:=function(K)
        # inputs a cubical knot
        # outputs a list of cubical complexes [C,U]
        # C being the complement of K with a certain line removed
        # (except for its end-points)
        # U is a plane intersecting C at just these endpoints
        local
              C, array, 0array, 0row,
              i, j, pos, a, b, subarray;

        C:=PureComplexComplement(K);
        array:=ShallowCopy(C!.binaryArray);
        0array:=0; # location of 1st non-empty array
        0row:=0; # location of the line to be altered

        for i in [1..Length(array)]
            do
            if 0row=0
                then
                for j in [1..Length(array[i])]
                    do
                    if 0 in array[i][j]
                        then
                        0array:=ShallowCopy(i);
                        0row:=0row+j;
                        break;
                    fi;
                od;
            else
                break;
            fi;
        od;

        if 0array=1
            then
            array:=FrameArray(array);
            0array:=0array+1;
        fi;

        pos:=Positions(array[0array][0row],0);
        a:=pos[1]; # the endpoints of 0row
        b:=pos[Length(pos)];

        array[0array][0row]:=0*ShallowCopy(array[0array][0row])+1;
        array[0array][0row][a]:=0;
        array[0array][0row][b]:=0;
        array[0array-1][0row]:=ShallowCopy(array[0array][0row]);

        subarray:=0*ShallowCopy(array); # the plane about which the knot
        subarray[0array-1]:=ShallowCopy(array[0array-1]); # complement will be spun

        return [PureCubicalComplex(array), PureCubicalComplex(subarray)];

    end;

    pair:=DeletedPair(K);
    C:=ShallowCopy(pair[1]);
    U:=ShallowCopy(pair[2]);
    C2:=ContractedComplex(C,U);

    omicron:=RegularCWMap(C2,U);

    return Spin(omicron);

end;
