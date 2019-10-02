Read("~/proj/Knots/knotcomp.g");
KnotComplementWithBoundary:=function(arc)
    local
        RegularCWKnot, knot, comp, inclusion, iota;

    RegularCWKnot:=function(arc)
        local
            D, len, signless, HollowTubes, max,
            bigGrid, correction, TubeJoiner;

        D:=arc;
        len:=Length(D);
        signless:=List(D,x->[AbsInt(x[1]),AbsInt(x[2])]);

        HollowTubes:=function(D)
            local
                grid, i, IsIntersection,
                CornerConfiguration, bound,
                bigGrid, GridFill, j, tick, correction,
                hslice1, hslice2, l1, l2, l3,
                max, vslice1, vslice2;

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

    # UPDATE: needed to account for configuration of corners at the end-step
    # (i.e. when matching the loops of one layer to the other).
    # There are sometimes disparities in the ordering on 1-cells from
    # left-to-right vs. when ordering from top-to-bottom. Not realising this was
    # causing 111 of the pre-stored knots in HAP to yield incorrect
    # CW-decomposition.
            correction:=[];
            for i in [1..len] do
                for j in [1..len] do
                    if CornerConfiguration(j,i)<>0 then
                        if CornerConfiguration(j,i) in [1,4] then
                            Add(correction,1);
                            Add(correction,-1);
                        else
                            Add(correction,0);
                            Add(correction,0);
                        fi;
                    fi;
                od;
            od;

            ### add the 0, 1 & 2-cells ###
            ########## to bound ##########
            for i in [1..2*Maximum(bigGrid[Length(bigGrid)])] do
                Add(bound[1],[1,0]);
            od;

            for i in [1..len] do # add the 'horizontal' 2-cells
                hslice1:=Filtered(bigGrid[2*i-1],x->x<>0);
                hslice2:=Filtered(bigGrid[2*i],x->x<>0);
                l2:=Length(bound[2]);
                for j in [1..Length(hslice1)-1] do
                    Add(
                        bound[2],
                        [2,hslice1[j],hslice2[j]]
                    );
                    if j=1 then
                        Add(
                            bound[2],
                            [2,hslice1[j],hslice2[j]]
                        );
                    fi;
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
                            bound[2],
                            [2,hslice1[j+1],hslice2[j+1]]
                        );
                        Add(
                            bound[3],
                            [4,l1-3,l1-2,l1-1,l1]
                        );
                    fi;
                od;
                l3:=Concatenation(
                    [l2+1],
                    Filtered(
                        [l2+3..Length(bound[2])-2],
                        x->AbsInt(bound[2][x][2]-bound[2][x][3])=1
                    ),
                    [Length(bound[2])]
                );
                Add(l3,Length(l3),1);
                Add(bound[3],l3);
            od;

            max:=Maximum(bigGrid[Length(bigGrid)]);
            for i in [1..2*len] do
                for j in [1..2*len] do
                    if bigGrid[i][j]<>0 then
                        bigGrid[i][j]:=bigGrid[i][j]+max;
                    fi;
                od;
            od;

            for i in [1..len] do # add the 'vertical' 2-cells
                vslice1:=Filtered(List([1..2*len],x->bigGrid[x][2*i-1]),x->x<>0);
                vslice2:=Filtered(List([1..2*len],x->bigGrid[x][2*i]),x->x<>0);
                l2:=Length(bound[2]);
                for j in [1..Length(vslice1)-1] do
                    Add(
                        bound[2],
                        Concatenation(
                            [2],
                            Set(
                                [vslice1[j],
                                vslice2[j]]
                            )
                        )
                    );
                    if j=1 then
                        Add(
                            bound[2],
                            Concatenation(
                                [2],
                                Set(
                                    [vslice1[j],
                                    vslice2[j]]
                                )
                            )
                        );
                    fi;
                    if j<>1 then
                        l1:=Length(bound[2]);
                        Add(
                            bound[3],
                            [4,l1-3,l1-2,l1-1,l1]
                        );
                    fi;
                    Add(
                        bound[2],
                        Concatenation([2],Set(vslice1{[j,j+1]}))
                    );
                    Add(
                        bound[2],
                        Concatenation([2],Set(vslice2{[j,j+1]}))
                    );
                    if j=Length(vslice1)-1 then
                        Add(
                            bound[2],
                            Concatenation(
                                [2],
                                Set(
                                    [vslice1[j+1],
                                    vslice2[j+1]]
                                )
                            )
                        );
                        l1:=Length(bound[2]);
                        Add(
                            bound[2],
                            Concatenation(
                                [2],
                                Set(
                                    [vslice1[j+1],
                                    vslice2[j+1]]
                                )
                            )
                        );
                        Add(
                            bound[3],
                            [4,l1-3,l1-2,l1-1,l1]
                        );
                    fi;
                od;
                l3:=Concatenation(
                    [l2+1],
                    Filtered(
                        [l2+3..Length(bound[2])-2],
                        x->AbsInt(bound[2][x][2]-bound[2][x][3])<>1
                    ),
                    [Length(bound[2])]
                );
                Add(l3,Length(l3),1);
                Add(bound[3],l3);
            od;
            ##############################

            return [max,bound,bigGrid,correction];
        end;

        D:=HollowTubes(D);
        max:=D[1];
        bigGrid:=D[3];
        correction:=D[4];
        D:=D[2];

        TubeJoiner:=function(D)
            local
                loops, size, hloops, vloops,
                i, l;

            loops:=Filtered([1..Length(D[2])],x->Length(Positions(D[2],D[2][x]))=2);
            size:=Length(loops)/2;
            hloops:=loops{[1..size]};
            vloops:=List([1..size],x->Position(D[2],D[2][hloops[x]]+[0,max,max]));
            for i in [1..size] do
                if i mod 2 = 0 then
                    vloops[i]:=vloops[i]+1;
                fi;
            od;
            vloops:=vloops+correction;

            for i in [1..size/2] do
                Add(D[2],[2,D[2][hloops[2*i]][2],D[2][vloops[2*i]][2]]);
                Add(D[2],[2,D[2][hloops[2*i]][3],D[2][vloops[2*i]][3]]);
                l:=Length(D[2]);

                Add(D[3],[4,hloops[2*i-1],vloops[2*i-1],l-1,l]);
                Add(D[3],[4,hloops[2*i],vloops[2*i],l-1,l]);
            od;

            return D;
        end;

        D:=TubeJoiner(D);
        D:=RegularCWComplex(D);
        D!.grid:=bigGrid;
        D!.arcPresentation:=arc;

        return D;
    end;

    knot:=RegularCWKnot(arc);
    comp:=KnotComplement(arc);

    inclusion:=function(bound)
        local
            bound1, bound2, len,
            1c1, 1c2, 2c2, inc;

        bound1:=bound[1];
        bound2:=bound[2];

        len:=Length(bound1[1])/2;

        1c1:=bound1[2]*1;
        1c1:=List(
            1c1,
            x->Concatenation(
                [2],
                [x[2]+1+2*Int((x[2]-1)/len),
                x[3]+1+2*Int((x[3]-1)/len)]
            )
        );

        1c2:=bound2[2]*1;
        1c2:=List(
            1c2,
            x->Concatenation(
                [2],
                Set(
                    [x[2],x[3]]
                )
            )
        );

        2c2:=List(bound2[3],x->Concatenation([x[1]],Set(x{[2..Length(x)]})));

        inc:=function(n,k)
            local
                ind, 2cell;

            if n=0 then
                return k+1+2*Int((k-1)/len);
            elif n=1 then
                ind:=1;
                if k>1 and 1c1[k-1]=1c1[k] then
                    ind:=2;
                fi;
                return Positions(1c2,1c1[k])[ind];
            else
                2cell:=List(bound1[3][k]{[2..4]},x->inc(1,x));
                2cell:=Concatenation([4],Set(2cell));
                return Position(2c2,2cell);
            fi;
        end;

        return inc;
    end;

    iota:=inclusion([knot!.boundaries,comp!.boundaries]);

    return Objectify(
        HapRegularCWMap,
        rec(
            source:=knot,
            target:=comp,
            mapping:=iota
        ));
end;
