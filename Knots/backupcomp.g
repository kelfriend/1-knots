#########################################################################################
######################## Regular CW-complex of knot complements #########################
#########################################################################################
################# Input: a list of (signed) integer pairs corresponding #################
######################## to the lengths of horizontal segments of a knot/link ###########
######################## and the depth of each segment in 4-space. ###################### 
#########################################################################################
################ Output: a regular CW-complex representing the complement of ############
######################## the input knot/link. ###########################################
##################################################################################### k.k
KnotComplement:=function(D)
    local
        len, signless, PuncturedDisk,
        P, DxI, Patches;

    len:=Length(D);
    signless:=List(D,x->[AbsInt(x[1]),AbsInt(x[2])]);

    PuncturedDisk:=function(D)
        local
            grid, i, IsIntersection,
            CornerConfiguration, bound,
            horizontal, tick, HorizontalFill,
            j, hslice, vslice, k, Edgewalker;

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
        horizontal:=List([1..2*len],x->List([1..2*len],y->0));
        tick:=1;
        Add(bound[1],[1,0]); # circumferential 0-cell no. 1

        HorizontalFill:=function(c,i,j);
# creates the loops present at each horizontal endpoint
            tick:=tick+1;
            if c=1 or c=4
                then
                horizontal[(2*i)-1][(2*j)-1]:=tick;
                tick:=tick+1;
                horizontal[2*i][2*j]:=tick;

                Add(bound[2],[2,tick-1,tick]);
                Add(bound[2],[2,tick-1,tick]);
            elif c=2 or c=3
                then
                horizontal[(2*i)-1][2*j]:=tick;
                tick:=tick+1;
                horizontal[2*i][(2*j)-1]:=tick;

                Add(bound[2],[2,tick-1,tick]);
                Add(bound[2],[2,tick-1,tick]);
            fi;

        end;

        for i in [1..len]
            do # add the 0-cells
            for j in [1..len]
                do
                if IsIntersection(i,j)
                    then # four 0-cells at an intersection
                    tick:=tick+1; horizontal[(2*i)-1][(2*j)-1]:=tick;
                    tick:=tick+1; horizontal[(2*i)-1][2*j]:=tick;
                    tick:=tick+1; horizontal[2*i][(2*j)-1]:=tick;
                    tick:=tick+1; horizontal[2*i][2*j]:=tick;
                    Add(bound[1],[1,0]); Add(bound[1],[1,0]);
                    Add(bound[1],[1,0]); Add(bound[1],[1,0]);
                elif grid[i][j]=1
                    then # two 0-cells at the endpoints of each horizontal bar
                    HorizontalFill(CornerConfiguration(i,j),i,j);
                    Add(bound[1],[1,0]); Add(bound[1],[1,0]);
                fi;
            od;
        od;
        Add(bound[1],[1,0]); # circumferential 0-cell no. 2

        for i in [1..2*len]
            do # connect all 0-cells that lie in the same
            hslice:=[]; # horizontal/vertical 'slice'
            vslice:=[];
            for j in [1..2*len]
                do
                if not horizontal[i][j]=0
                    then
                    Add(hslice,horizontal[i][j]);
                fi;
                if not horizontal[j][i]=0
                    then
                    Add(vslice,horizontal[j][i]);
                fi;
            od;
            for k in [1..Length(hslice)]
                do
                if Length(hslice)>k
                    then
                    Add(bound[2],[2,hslice[k],hslice[k+1]]);
                fi;
            od;
            for k in [1..Length(vslice)]
                do
                if Length(vslice)>k
                    then
                    Add(bound[2],[2,vslice[k],vslice[k+1]]);
                fi;
            od;
        od;

        Add(bound[2],[2,1,2]); # connect the central knot to the circumference
        Add(bound[2],[2,Length(bound[1]-1),Length(bound[1])]); # of the disk
        Add(bound[2],[2,1,Length(bound[1])]);
        Add(bound[2],[2,1,Length(bound[1])]);

        Edgewalker:=function(bound)
            local
                unchosen, neighbours, i, j;

            unchosen:=List(ShallowCopy(bound[2]),x->[x[2],x[3]]);
            neighbours:=List(ShallowCopy(bound[1]),x->[]);

            for i in [1..Length(bound[1])]
                do
                for j in [1..Length(unchosen)]
                    do
                    if i in unchosen[j]
                        then
                        Add(neighbours[i],j);
                    fi;
                od;
            od;

            return neighbours;
        end;

        return horizontal;
    end;

    P:=PuncturedDisk(D);

    return P;
end;