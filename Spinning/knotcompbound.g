KnotComplementWithBoundary:=function(arc)
    local
        D, len, signless, PuncturedDisk;

    D:=arc;
    len:=Length(D);
    signless:=List(D,x->[AbsInt(x[1]),AbsInt(x[2])]);

    PuncturedDisk:=function(D)
        local
            grid, i, IsIntersection,
            CornerConfiguration, bound,
            bigGrid, GridFill, j, tick,
            hslice1, hslice2, l1, vslice1,
            vslice2;

        grid:=List([1..len],x->List([1..len],y->0));
        for i in [1..len]
            do
            grid[len-i+1][D[i][1]]:=1;
            grid[len-i+1][D[i][2]]:=1;
        od;

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

        CornerConfiguration:=function(i,j);

            if grid[i][j]=1
                then
                if Size(Positions(grid[i]{[j..len]},1))=2
                    then
                    if Size(Positions(List([i..len],x->grid[x][j]),1))=2
                        then # Corner type 1, i.e : __
                        return 1; #                |
                    elif Size(Positions(List([1..i],x->grid[x][j]),1))=2
                        then # Corner type 3, i.e :
                        return 3; #                |__
                    fi;
                elif Size(Positions(grid[i]{[1..j]},1))=2
                    then
                    if Size(Positions(List([i..len],x->grid[x][j]),1))=2
                        then # Corner type 2, i.e : __
                        return 2; #                   |
                    elif Size(Positions(List([1..i],x->grid[x][j]),1))=2
                        then # Corner type 4, i.e :
                        return 4; #                 __|
                    fi;
                fi;
            fi;

            return 0;
        end;

        bound:=[[],[],[],[],[]];
        bigGrid:=List([1..2*len],x->List([1..2*len],y->0));

        GridFill:=function(c,i,j);
            if c=1 or c=4
                then
                bigGrid[(2*i)-1][(2*j)-1]:='*';
                bigGrid[2*i][2*j]:='*';
            elif c=2 or c=3
                then
                bigGrid[(2*i)-1][2*j]:='*';
                bigGrid[2*i][(2*j)-1]:='*';
            fi;
        end;

        for i in [1..len]
            do
            for j in [1..len]
                do
                if IsIntersection(i,j)
                    then
                    bigGrid[(2*i)-1][(2*j)-1]:='*';
                    bigGrid[(2*i)-1][2*j]:='*';
                    bigGrid[2*i][(2*j)-1]:='*';
                    bigGrid[2*i][2*j]:='*';
                elif grid[i][j]=1
                    then
                    GridFill(CornerConfiguration(i,j),i,j);
                fi;
            od;
        od;

        tick:=1;
        for i in [1..2*len]
            do
            for j in [1..2*len]
                do
                if bigGrid[i][j]='*'
                    then
                    bigGrid[i][j]:=tick;
                    tick:=tick+1;
                fi;
            od;
        od;

        ### add the 0, 1 & 2-cells ###
        ########## to bound ##########
        for i in [1..Maximum(bigGrid[Length(bigGrid)])] do
            Add(bound[1],[1,0]);
        od;

        for i in [1..len] do
            hslice1:=Filtered(bigGrid[2*i-1],x->x<>0);
            hslice2:=Filtered(bigGrid[2*i],x->x<>0);
            for j in [1..Length(hslice1)-1] do
                Add(
                    bound[2],
                    [2,hslice1[j],hslice2[j]]
                );
                if j<>1 then
                    l1:=Length(bound[2]);
                    Add(
                        bound[3],
                        [4,l1-3,l1-2,l1-1,l1]
                    );
                fi;
                Add(
                    bound[2],
                    Concatenation([2],hslice1{[j,j+1]})
                );
                Add(
                    bound[2],
                    Concatenation([2],hslice2{[j,j+1]})
                );
                if j=Length(hslice1)-1 then
                    Add(
                        bound[2],
                        [2,hslice1[j+1],hslice2[j+1]]
                    );
                    l1:=Length(bound[2]);
                    Add(
                        bound[3],
                        [4,l1-3,l1-2,l1-1,l1]
                    );
                fi;
            od;
        od;

        for i in [1..len] do
            vslice1:=[];
            vslice2:=[];
        od;
        ##############################

        return bound;
    end;
    return PuncturedDisk(D);
end;
trefoil:=[[2,5],[1,3],[2,4],[3,5],[1,4]];;
