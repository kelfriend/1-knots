# Input: List of signed integer pairs.
# Output: A 2-dimensional regular CW-complex.

TubeKnot:=function(l);
    local
        len, signless, bound,
        IsIntersection, grid, i,
        j, 0c, no0s, no0, o, j, k;

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
            if Size(o[1])=2 then
                Add(no0,o[1]);
                Add(no0,o[2]);
            else
                Add(no0,o);
            fi;
        od;
        for j in [1..Length(no0)-1] do
            for k in [1,2] do
                Add(bound[2],[2,no0[j][1],no0[j][2]]);
                Add(bound[2],[2,no0[j][k],no0[j+1][k]]);
                Add(bound[2],[2,no0[j+1][1],no0[j+1][2]]);
            od;
        od;
    od;

    for i in [1..len] do
        no0:=Filtered(List([1..len],x->grid[x][i]),x->x<>0);
        for j in [1..Length(no0)-1] do
            if Size(no0[j][1])=2 then
                if Size(no0[j+1][1])=2 then
                    for k in [1,2] do
                        Add(bound[2],[2,no0[j][k][1],no0[j+1][k][1]]);
                        Add(bound[2],[2,no0[j][k][2],no0[j+1][k][2]]);
                    od;
                else
                    for k in [1,2] do
                        Add(bound[2],[2,no0[j][1][k],no0[j+1][k]]);
                        Add(bound[2],[2,no0[j][2][k],no0[j+1][k]]);
                    od;
                fi;
            else
                if Size(no0[j+1][1])=2 then
                    for k in [1,2] do
                        Add(bound[2],[2,no0[j][1],no0[j+1][k][1]]);
                        Add(bound[2],[2,no0[j][2],no0[j+1][k][2]]);
                    od;
                else
                    Add(bound[2],[2,no0[j][1],no0[j+1][1]]);
                    Add(bound[2],[2,no0[j][2],no0[j+1][2]]);
                fi;
            fi;
        od;
    od;

    return bound;
end;
K:=[[2,5],[1,3],[2,4],[3,5],[1,4]];;