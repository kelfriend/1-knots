#  Input: an arc presentation.
# Output: inclusion of regular CW-complexes from a 1-dimesional link
#         to a 3-dimensional ball.
KnotComplementOn1Skeleton:=function(arc)
    local
        gn, grid, i, kbnd, bnd, map, imap,
        IsIntersection, ints, ext_ints, j, 0c, B2Decomposition;

    gn:=Length(arc); # the grid number
    
    grid:=List([1..5*gn],x->[1..5*gn]*0); # form a (3*gn) x gn
    for i in [0..gn-1] do # matrix from the arc presentation
        grid[5*(gn-i)-2][5*arc[i+1][1]-2]:=1;
        grid[5*(gn-i)-2][5*arc[i+1][2]-2]:=1;
    od;

    kbnd:=List([1..5],x->[]); # boundary list of the knot
    bnd:=List([1..5],x->[]); # boundary list of the complement of the knot
    map:=List([1..4],x->[]); # inclusion map from kbnd to bnd
    imap:=List([1..4],x->[]); # inverse image of the above inclusion map

    IsIntersection:=function(i,j) # finds where crossings occur in the grid
        if grid[i][j]=0 then
            if 1 in grid[i]{[1..j]} then
                if 1 in grid[i]{[j..5*gn]} then
                    if 1 in List([1..i],x->grid[x][j]) then
                        if 1 in List([i..5*gn],x->grid[x][j]) then
                            return true;
                        fi;
                    fi;
                fi;
            fi;
        fi;

        return false;
    end;

    ints:=[]; # records the coordinates of each crossing
    ext_ints:=[]; # records those 0-cells which arise in intersection points
    for i in [1..5*gn] do
        for j in [1..5*gn] do
            if IsIntersection(i,j) then
                Add(ext_ints,[i-1,j]); Add(ext_ints,[i+1,j]);
                Add(ints,[i,j]);
                grid[i-1][j]:='*';
                grid[i][j]:='*';
                grid[i+1][j]:='*';
            fi;
        od;
    od;

    0c:=4; # label the entries of grid so that it models the 0-skeleton
    for i in [1..5*gn] do
        for j in [1..5*gn] do
            if grid[i][j]<>0 then
                grid[i][j]:=0c;
                0c:=0c+1;
            fi;
        od;
    od;
    0c:=0c+2;

    ints:=List(ints,x->grid[x[1],x[2]]);

    B2Decomposition:=function()
# takes what we have so far and uses it to form a regular CW-decomposition of
# the 2-ball with the appropriate inclusion map
        local i, j, hslice, vslice, 2CellTracer;

        kbnd[1]:=[1..0c-6]; bnd[1]:=[1..0c];
        map[1]:=kbnd[1]+3; imap[1]:=Concatenation([0,0,0],map[1],[0,0,0]);

        Add(bnd[2],[2,1,2]); Add(bnd[2],[2,1,3]); Add(bnd[2],[2,2,0c-2]);
        Add(bnd[2],[2,3,0c-1]); Add(bnd[2],[2,0c-2,0c]); Add(bnd[2],[2,0c-1,0c]);
        Add(bnd[2],[2,3,4]); Add(bnd[2],[2,0c-3,0c-2]);
        # add 1-cells to just the complement to stay regular

        for i in [1..8] do
            Add(imap[2],0);
        od;

        for i in [1..5*gn] do
            hslice:=[];
            for j in [1..5*gn] do
                if grid[i][j]<>0 and not [i,j] in ext_ints then
                    Add(hslice,grid[i][j]);
                fi;
            od;
            for j in [1..Length(hslice)-1] do
                Add(kbnd[2],[2,hslice[j]-3,hslice[j+1]-3]);
                Add(bnd[2],[2,hslice[j],hslice[j+1]]);
                Add(map[2],Length(bnd[2]));
                Add(imap[2],Length(kbnd[2]));
            od;
        od;
        for j in [1..5*gn] do
            vslice:=[];
            for i in [1..5*gn] do
                if grid[i][j]<>0 then
                    Add(vslice,grid[i][j]);
                fi;
            od;
            for i in [1..Length(vslice)-1] do
                Add(bnd[2],[2,vslice[i],vslice[i+1]]);
                if not(vslice[i] in ints or vslice[i+1] in ints) then
                    Add(kbnd[2],[2,vslice[i]-3,vslice[i+1]-3]);
                    Add(map[2],Length(bnd[2]));
                    Add(imap[2],Length(kbnd[2]));
                else
                    Add(imap[2],0);
                fi;
            od;
        od;

        2CellTracer:=function()
            local
                ori, i, j, cell, top, rgt,
                btm, lft, s, dir, path, unselected,
                tick, t, edge, p_edge, s2;

            grid:=FrameArray(grid);
            grid[1][1]:=1; grid[1][5*gn+2]:=2; grid[4][1]:=3;
            grid[5*gn-1][5*gn+2]:=0c-2; grid[5*gn+2][1]:=0c-1;
            grid[5*gn+2][5*gn+2]:=0c;

            ori:=List([1..0c],x->[1..4]*0);
            # each 0-cell will have the 0-cells N/E/S/W of it
            # listed in that order 

            for i in [1..5*gn+2] do
                for j in [1..5*gn+2] do
                    cell:=grid[i][j];
                    if cell<>0 then
                        top:=List([1..i-1],x->grid[x][j]);
                        top:=Filtered(top,x->x<>0);
                        if top<>[] then
                            ori[cell][1]:=top[Length(top)];
                        fi;
                        rgt:=grid[i]{[j+1..5*gn+2]};
                        rgt:=Filtered(rgt,x->x<>0);
                        if rgt<>[] then
                            ori[cell][2]:=rgt[1];
                        fi;
                        btm:=List([i+1..5*gn+2],x->grid[x][j]);
                        btm:=Filtered(btm,x->x<>0);
                        if btm<>[] then
                            ori[cell][3]:=btm[1];
                        fi;
                        lft:=grid[i]{[1..j-1]};
                        lft:=Filtered(lft,x->x<>0);
                        if lft<>[] then
                            ori[cell][4]:=lft[Length(lft)];
                        fi;
                        if [i,j] in ext_ints then # 0-cells in ext_ints
                            ori[cell][2]:=0; # never have edges from the left
                            ori[cell][4]:=0; # or right of them
                        fi;
                    fi;
                od;
            od;

# having assigned orientations to all 0-cells, we now perform an anticlockwise walk
# along the 1-skeleton so as to trace all of the 2-cells of bnd
            p_edge:=1; # init: select the first edge of bnd
            path:=[]; # 1-cells present in the boundary of the 2-cell
            unselected:=Concatenation
                (
                    [1..Length(bnd[2])],
                    [7..Length(bnd[2])] # the first 6 1-cells occur just once
                );
            tick:={n}->(n mod 4)+1; # do the weird anticlock walk that flips directions
            while unselected<>[] do # then delete cells from ori per edge use ###maybe not though
                edge:=bnd[2][p_edge]; # add a function to switch source and target so that things dont get stuck
                s:=edge[2]; # source vertex
                t:=edge[3]; # target vertex
                dir:=Position(ori[s],t);
                if not p_edge in path then
                    ####ori[s][dir]:=0;
                    Unbind(unselected[Position(unselected,p_edge)]);
                    Add(path,p_edge);

                    s2:=t*1; # the new source vertex
                    dir:=tick(dir);
                    t:=ori[s2][dir];
                    while t=0 or t=s do
                        dir:=tick(dir);
                        t:=ori[s2][dir];
                    od;
                    p_edge:=Position(bnd[2],Concatenation([2],Set([s2,t])));
                else
                    path:=Set(path);
                    Add(path,Length(path),1);
                    Add(bnd[3],path);
                    Add(imap[3],0);
                    path:=[];
                    p_edge:=Filtered(unselected,IsInt)[1];
                fi;
            od;
        end;
        2CellTracer();
    end;
    B2Decomposition();

    return [grid,kbnd,bnd,map,imap];
end;
