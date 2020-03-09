#  Input: an arc presentation.
# Output: inclusion of regular CW-complexes from a 1-dimesional link
#         to a 3-dimensional ball.
KnotComplementOn1Skeleton:=function(arc)
    local
        gn, grid, i, kbnd, bnd, map, imap,
        IsIntersection, ints, j, B2Decomposition;

    gn:=Length(arc); # the grid number
    
    grid:=List([1..3*gn],x->[1..gn]*0); # form a (3*gn) x gn
    for i in [0..gn-1] do # matrix from the arc presentation
        grid[3*(gn-i)][arc[i+1][1]]:=1;
        grid[3*(gn-i)][arc[i+1][2]]:=1;
    od;

    kbnd:=List([1..4],x->[]); # boundary list of the knot
    bnd:=List([1..4],x->[]); # boundary list of the complement of the knot
    map:=List([1..3],x->[]); # inclusion map from kbnd to bnd
    imap:=List([1..3],x->[]); # inverse image of the above inclusion map

    IsIntersection:=function(i,j) # finds where crossings occur in the grid
        if grid[i][j]=0 then
            if 1 in grid[i]{[1..j]} then
                if 1 in grid[i]{[j..gn]} then
                    if 1 in List([1..i],x->grid[x][j]) then
                        if 1 in List([i..3*gn],x->grid[x][j]) then
                            return true;
                        fi;
                    fi;
                fi;
            fi;
        fi;

        return false;
    end;

    ints:=[]; # records the coordinates of each crossing
    for i in [1..3*gn] do
        for j in [1..gn] do
            if IsIntersection(i,j) then
                Add(ints,[i,j]);
                grid[i-1][j]:=1;
                grid[i][j]:=1;
                grid[i+1][j]:=1;
            fi;
        od;
    od;

    B2Decomposition:=function();
# takes what we have so far and uses it to form a regular CW-decomposition of
# the 2-ball with the appropriate inclusion map
        local 0c, i, j, hslice;

        0c:=2; # label the entries of grid so that it models the 0-skeleton
        for i in [1..3*gn] do
            for j in [1..gn] do
                if grid[i][j]=1 then
                    grid[i][j]:=0c;
                    0c:=0c+1;
                fi;
            od;
        od;
        0c:=0c+1;

        kbnd[1]:=[1..0c-2]; bnd[1]:=[1..0c];
        map[1]:=kbnd[1]+1, imap[1]:=Concatenation(['*'],map[1],['*']);

        Add(bnd[2],[2,1,2]);
        Add(bnd[2],[2,1,0c]); Add(bnd[2],[2,1,0c]);
        for i in [1..3*gn] do
            hslice:=[];
            for j in [1..gn] do
                if grid[i][j]<>0 then
                    Add(hslice,grid[i][j]);
                fi;
            od;
            for j in [1..Length(hslice)-1] do
                Add(kbnd[2],[2,hslice[1]-1,hslice[i+1]-1]);
                Add(bnd[2],[2,hslice[1],hslice[i+1]]);
                Add(map[2],Length(bnd[2]));
            od;
        od;
        for j in [1..gn] do
            vslice:=[];
            for i in [1..3*gn] do
                if grid[i][j]<>0 then
                    Add(vslice,grid[i][j]);
                fi;
            od;
            for j in [1..Length(vslice)-1] do
                Add(bnd[2],[2,j,j+1]); ####################
                if [i,j] in ints then
                fi;
            od;
        od;

                 
        Add(bnd[2],[2,0c-1,0c]); # 1-cell to maintain regularity


    end;

    return imap;
end;