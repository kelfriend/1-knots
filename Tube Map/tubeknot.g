# Input: List of signed integer pairs.
# Output: A 2-dimensional regular CW-complex.

TubeKnot:=function(l)
    local
        len, signless, bound,
        IsIntersection, grid, i,
        j, 0c, no0s, no0, o, b, 1c,
        l1, l2, l3, l4, filter, 3c;

    len:=Length(l);
    signless:=List(l,x->[AbsInt(x[1]),AbsInt(x[2])]);
    bound:=[[],[],[],[]];

    IsIntersection:=function(i,j)

        if grid[i][j]=0
            then
            if 1 in grid[i]{[1..j]}
                then
                if 1 in grid[i]{[j..len]}
                    then
                    if 1 in List([1..i],x->grid[x][j])
                        then
                        if 1 in List([i..len],x->grid[x][j])
                            then
                            return true;
                        fi;
                    fi;
                fi;
            fi;
        fi;

        return false;
    end;

    grid:=List([1..len],x->List([1..len],x->0));
    for i in [1..len] do
        for j in [1..len] do
            if j=Reversed(signless)[i][1] then
                grid[i][j]:=1;
            elif j=Reversed(signless)[i][2] then
                grid[i][j]:=1;
            fi;
        od;
    od;
    for i in [2..len-1] do
        for j in [2..len-1] do
            if IsIntersection(i,j) then
                grid[i][j]:=2;
            fi;
        od;
    od;

    0c:=1;
    for i in [1..len] do
        for j in [1..len] do
            if grid[i][j]=1 then
                grid[i][j]:=[0c,0c+1];
                0c:=0c+2;
            elif grid[i][j]=2 then
                grid[i][j]:=[[0c,0c+1],[0c+2,0c+3]];
                0c:=0c+4;
            fi;
        od;
    od;

    for i in [1..0c-1] do
        Add(bound[1],[1,0]); 
    od;

    for i in [1..len] do
        no0s:=Filtered(grid[i],x->x<>0);
        no0:=[];
        for o in no0s do
            if not IsInt(o[1]) then
                Add(no0,o[1]);
                Add(no0,o[2]);
            else
                Add(no0,o);
            fi;
        od;
        b:=-4;
        1c:=Length(bound[2])+1;
        for j in [1..Length(no0)-1] do
            b:=b+4;

            if j=1 then
                Add(bound[2],[2,no0[j][1],no0[j][2]]);
                Add(bound[2],[2,no0[j][1],no0[j][2]]);
            fi;

            Add(bound[2],[2,no0[j][1],no0[j+1][1]]);
            Add(bound[2],[2,no0[j][2],no0[j+1][2]]);

            Add(bound[2],[2,no0[j+1][1],no0[j+1][2]]);
            Add(bound[2],[2,no0[j+1][1],no0[j+1][2]]);
            # here is where the 2-cells for the horizontal tube are formed
            Add(bound[3],[4,1c,1c+2,1c+4,1c+3]+[0,b,b,b,b]);
            Add(bound[3],[4,1c+1,1c+2,1c+5,1c+3]+[0,b,b,b,b]);
        od;
    od;

    for i in [1..len] do
        no0:=Filtered(List([1..len],x->grid[x][i]),x->x<>0);
        for j in [1..Length(no0)-1] do
            if not IsList(no0[j][1]) then
                1c:=Length(bound[2]);
                l1:=Position(bound[2],[2,no0[j][1],no0[j][2]]);
                if not IsList(no0[j+1][1]) then
                    l2:=Position(bound[2],[2,no0[j+1][1],no0[j+1][2]]);

                    Add(bound[2],[2,no0[j][1],no0[j+1][1]]); # connect
                    Add(bound[2],[2,no0[j][2],no0[j+1][2]]); # vertically

                    Add(bound[3],[4,1c+1,1c+2,l1+1,l2+1]); # add the 2-cells
                    Add(bound[3],[4,1c+1,1c+2,l1,l2]);
                else
                    Add(bound[2],[2,no0[j+1][1][1],no0[j+1][2][1]]);
                    Add(bound[2],[2,no0[j+1][1][2],no0[j+1][2][2]]);

                    l2:=Position(bound[2],[2,no0[j+1][1][1],no0[j+1][1][2]]);
                    l3:=Position(bound[2],[2,no0[j+1][2][1],no0[j+1][2][2]]);

                    Add(bound[2],[2,no0[j][1],no0[j+1][1][1]]);
                    Add(bound[2],[2,no0[j][2],no0[j+1][1][2]]);

                    Add(bound[3],[4,1c+3,1c+4,l1+1,l2]);
                    Add(bound[3],[6,1c+3,1c+4,l1,1c+1,1c+2,l3]);
                    Add(bound[3],[2,l2,l2+1]); # patches up the holes
                    Add(bound[3],[2,l3,l3+1]); # in the vertical tubes
                fi;
            else
                1c:=Length(bound[2]);
                l1:=Position(bound[2],[2,no0[j][1][1],no0[j][1][2]]);
                l2:=Position(bound[2],[2,no0[j][2][1],no0[j][2][2]]);
                if not IsList(no0[j+1][1]) then
                    l3:=Position(bound[2],[2,no0[j+1][1],no0[j+1][2]]);

                    Add(bound[2],[2,no0[j][1][1],no0[j+1][1]]);
                    Add(bound[2],[2,no0[j][1][2],no0[j+1][2]]);

                    Add(bound[3],[4,1c+1,1c+2,l1+1,l3+1]);
                    Add(bound[3],[6,1c+1,1c+2,l2+1,l3,1c-3,1c-2]);
                else
                    Add(bound[2],[2,no0[j+1][1][1],no0[j+1][2][1]]);
                    Add(bound[2],[2,no0[j+1][1][2],no0[j+1][2][2]]);

                    l3:=Position(bound[2],[2,no0[j+1][1][1],no0[j+1][1][2]]);
                    l4:=Position(bound[2],[2,no0[j+1][2][1],no0[j+1][2][2]]);

                    Add(bound[2],[2,no0[j][1][1],no0[j+1][1][1]]);
                    Add(bound[2],[2,no0[j][1][2],no0[j+1][1][2]]);

                    Add(bound[3],[4,1c+3,1c+4,l1+1,l3]);
                    Add(bound[3],[8,1c+1,1c+2,1c-2,1c-3,1c+3,1c+4,l2+1,l4]);
                    Add(bound[3],[2,l3,l3+1]);
                    Add(bound[3],[2,l4,l4+1]);
                fi;
            fi;
        od;
    od;

    #filter:=[];
    #3c:=[];
    #for i in Filtered(bound[3],x->x[1]=2) do
    #    Add(filter,i[2]);
    #    Add(filter,i[3]);
    #od;

    #for i in [1..Length(bound[3])] do
    #    if not Length(Intersection(
    #            bound[3][i]{[2..Length(bound[3][i])]},
    #            filter))=2 then
    #        Add(3c,i);
    #    fi;
    #od;

    #Add(3c,Length(3c),1);
    #Add(bound[4],3c);
    #Add(bound,[]);

    return RegularCWComplex(bound);
end;
K:=[[2,4],[1,3],[2,4],[1,3]];;
