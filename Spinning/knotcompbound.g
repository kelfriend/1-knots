KnotComplementWithBoundary:=function(knot)
    local
        D, len, signless, PuncturedDisk,
        P, grid, ContestedEdge, Interior,
        PuncturedTube;

    D:=knot;
    len:=Length(D);
    signless:=List(D,x->[AbsInt(x[1]),AbsInt(x[2])]);

    PuncturedDisk:=function(D)
        local
            grid, i, IsIntersection,
            CornerConfiguration, bound,
            bigGrid, GridFill, j, tick,
            hslice, vslice, k, 0c, Orient,
            path, FaceTrace, cgrid;

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
# places an * at each point where a 0-cell is to be added to bigGrid
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
            do # loop through bigGrid and add temporary *s
            for j in [1..len]
                do
                if IsIntersection(i,j)
                    then # four 0-cells at an intersection
                    bigGrid[(2*i)-1][(2*j)-1]:='*';
                    bigGrid[(2*i)-1][2*j]:='*';
                    bigGrid[2*i][(2*j)-1]:='*';
                    bigGrid[2*i][2*j]:='*';
                elif grid[i][j]=1
                    then # two 0-cells at the endpoints of each horizontal bar
                    GridFill(CornerConfiguration(i,j),i,j);
                fi;
            od;
        od;

        tick:=2;
        for i in [1..2*len]
            do # number the 0-cells row-by-row
            for j in [1..2*len]
                do
                if bigGrid[i][j]='*'
                    then
                    bigGrid[i][j]:=tick;
                    tick:=tick+1;
                fi;
            od;
        od;

        for i in [1..2*len]
            do # connect all 0-cells that lie in the same
            hslice:=[]; # horizontal/vertical 'slice'
            vslice:=[];
            for j in [1..2*len]
                do
                if not bigGrid[i][j]=0
                    then
                    Add(hslice,bigGrid[i][j]);
                fi;
                if not bigGrid[j][i]=0
                    then
                    Add(vslice,bigGrid[j][i]);
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

        for i in [1..len]
            do # add the looping 1-cells to the 1-skeleton
            for j in [1..len]
                do
                if CornerConfiguration(i,j) in [1,4]
                    then
                    Add(bound[2],[
                        2,
                        bigGrid[(2*i)-1][(2*j)-1],
                        bigGrid[2*i][2*j]
                        ]
                    );
                    Add(bound[2],[
                        2,
                        bigGrid[(2*i)-1][(2*j)-1],
                        bigGrid[2*i][2*j]
                        ]
                    );
                elif CornerConfiguration(i,j) in [2,3]
                    then
                    Add(bound[2],[
                        2,
                        bigGrid[(2*i)-1][2*j],
                        bigGrid[2*i][(2*j)-1]
                        ]
                    );
                    Add(bound[2],[
                        2,
                        bigGrid[(2*i)-1][2*j],
                        bigGrid[2*i][(2*j)-1]
                        ]
                    );
                fi;
            od;
        od;

        0c:=Maximum(List(bigGrid,x->Maximum(x)));
        for i in [1..0c+1]
            do
            Add(bound[1],[1,0]);
        od;

        Add(bound[2],[2,1,2]); # connect the central component to the
        Add(bound[2],[2,Length(bound[1])-1,Length(bound[1])]); # circumference
        Add(bound[2],[2,1,Length(bound[1])]); # of the disk
        Add(bound[2],[2,1,Length(bound[1])]);

        bigGrid:=FrameArray(bigGrid);
        bigGrid[1][2]:=1; # Adds the first and last 0-cells to bigGrid
        bigGrid[Length(bigGrid)][Length(bigGrid[1])-1]:=0c+1;

        Orient:=function(bound)
# traces the 1-skeleton in a clockwise walk to yield the 2-cells
            local 
                unchosen, neighbours, i, j,
                Clockwise;

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

            Clockwise:=function(neighbours)
                local # orders clockwise the neighbours of each 0-cell
                    oriented, first0, last0,
                    i, j, x, k, l, posi, posx;

                oriented:=List(neighbours,x->List([1..12],y->"pass"));
                first0:=SortedList(neighbours[1]);
                last0:=SortedList(neighbours[Length(neighbours)]);

                oriented[1][7]:=first0[1];
                oriented[1][6]:=first0[3];
                oriented[1][8]:=first0[2];
            # these two orderings are always fixed;
            # they correspond to the circumferential edges
                oriented[Length(oriented)][1]:=last0[1]; 
                oriented[Length(oriented)][2]:=last0[3];
                oriented[Length(oriented)][12]:=last0[2];

                for i in [2..Length(neighbours)-1]
                    do # excludes the 1st and last 0-cells
                    for j in [1..Length(neighbours[i])]
                        do # x is a neighbouring 0-cell to i
                        x:=bound[2][neighbours[i][j]];
                        x:=Filtered(x{[2,3]},y->y<>i)[1];
                        for k in [1..Length(bigGrid)]
                            do
                            for l in [1..Length(bigGrid[1])]
                                do
                                if i=bigGrid[k][l]
                                    then
                                    posi:=[k,l];
                                fi;
                                if x=bigGrid[k][l]
                                    then
                                    posx:=[k,l];
                                fi;
                            od;
                        od;
                        # below are the checks for orientation,
                        # there are 12 in total (two for each diagonal):
                        # _\\|//_
                        #  //|\\
                        # ! ugly code warning !
                        if posi[1]>posx[1]
                            then
                            if posi[2]=posx[2]
                                then
                                oriented[i][1]:=neighbours[i][j];
                            elif posi[2]<posx[2]
                                then
                                if oriented[i][2]="pass"
                                    then # *always assigns the upper loop first*
                                    oriented[i][2]:=neighbours[i][j];
                                else
                                    oriented[i][3]:=neighbours[i][j];
                                fi;
                            elif posi[2]>posx[2]
                                then
                                if oriented[i][12]="pass"
                                    then
                                    oriented[i][12]:=neighbours[i][j];
                                else
                                    oriented[i][11]:=neighbours[i][j];
                                fi;
                            fi;
                        elif posi[1]=posx[1]
                            then
                            if posi[2]<posx[2]
                                then
                                oriented[i][4]:=neighbours[i][j];
                            elif posi[2]>posx[2]
                                then
                                oriented[i][10]:=neighbours[i][j];
                            fi;
                        elif posi[1]<posx[1]
                            then
                            if posi[2]=posx[2]
                                then
                                oriented[i][7]:=neighbours[i][j];
                            elif posi[2]<posx[2]
                                then
                                if oriented[i][5]="pass"
                                    then
                                    oriented[i][5]:=neighbours[i][j];
                                else
                                    oriented[i][6]:=neighbours[i][j];
                                fi;
                            elif posi[2]>posx[2]
                                then
                                if oriented[i][9]="pass"
                                    then
                                    oriented[i][9]:=neighbours[i][j];
                                else
                                    oriented[i][8]:=neighbours[i][j];
                                fi;
                            fi;
                        fi;
                    od;
                od;
                
                return oriented;
            end;

            return Clockwise(neighbours);
        end;

        path:=Orient(bound);
        # this is an ordered list of the neighbours of each 1-cell

        FaceTrace:=function(path)
            local
                unselectedEdges, sourceORtarget,
                x, ClockwiseTurn, 2cell,
                sORt, ori, e1, e0, i;

            unselectedEdges:=List([1..Length(bound[2])-2]);
            unselectedEdges:=Concatenation(unselectedEdges,unselectedEdges);
            Add(unselectedEdges,Length(bound[2])-1);
            Add(unselectedEdges,Length(bound[2]));
# list of two of each edge except for the circumferential edges

            ClockwiseTurn:=function(p,e)
# inputs the orientation list of a node and the number of an edge in that list,
# outputs the next edge after a clockwise turn
                local
                    f;
                
                f:=(Position(p,e) mod 12)+1;
                while p[f]="pass"
                    do
                    f:=(f mod 12)+1;
                od;
                
                return p[f];
            end;

            sourceORtarget:=List([1..Length(bound[2])],y->[3,2]);
            x:=1;
            while unselectedEdges<>[]
                do # main loop, locates all 2-cells
                while (not x in unselectedEdges) and (not e1 in unselectedEdges)
                    do
                    x:=x+1;
                od;
                2cell:=[x]; # the 2-cell begins with just x in its boundary
                sORt:=sourceORtarget[x][Length(sourceORtarget[x])];
                Unbind(sourceORtarget[x][Length(sourceORtarget[x])]);
                ori:=path[bound[2][x][sORt]]; # the orientation of x's target
                e0:=bound[2][x][sORt];
                e1:=ClockwiseTurn(ori,x); # next edge to travel along
                while e1<>x
                    do
                    Add(2cell,e1);
                    e0:=Filtered(bound[2][e1]{[2,3]},y->y<>e0)[1]; # e1's target
                    ori:=path[e0];
                    e1:=ClockwiseTurn(ori,e1);
                od;
                Add(2cell,Length(2cell),1);
                if (not Set(2cell) in List(bound[3],x->Set(x)))
                    then
                    for i in Filtered(2cell{[2..Length(2cell)]},
                    y->y in unselectedEdges)
                        do
                        Unbind(unselectedEdges[Position(unselectedEdges,i)]);
                    od;
                    Add(bound[3],2cell);
                fi;
            od;

            bound[3]:=Filtered(bound[3],y->y[1]<>2);
            return bound;
        end;

        cgrid:=grid*0; # this is needed at the very end when
        for i in [1..Length(grid)] # patching the tubes together
            do
            for j in [1..Length(grid)]
                do
                cgrid[i][j]:=CornerConfiguration(i,j);
            od;
        od;

        return [FaceTrace(path),cgrid];
    end;

    P:=PuncturedDisk(D);
    grid:=P[2];
    P:=P[1];

    ##### filtering process #####
    ContestedEdge:=function(Interior)
        local
            faces, faceoff, i, j,
            interiorfaces, k,
            exteriorfaces;

        faces:=Interior[3];
        faceoff:=[];
        for i in faces do
            for j in Filtered(faces,x->x<>i) do
                if Length(Intersection(i,j))>1 then
                    Add(faceoff,i);
                    Add(faceoff,j);
                fi;
            od;
        od;

        interiorfaces:=[];
        for i in [1..Length(faceoff)] do
            for j in faceoff[i] do
                for k in [2..Length(j)] do
                    if Length(Positions(Interior[2],Interior[2][k]))>1 then
                        if not j in interiorfaces then
                            Add(interiorfaces,j);
                        fi;
                    fi;
                od;
            od;
        od;

        exteriorfaces:=Difference(Union(faceoff),interiorfaces);
        
        return interiorfaces;#Filtered(Interior[3],x->not x in exteriorfaces);
    end;

    Interior:=[
        List([1..Length(P[1])-2],x->[1,0]),
        List(P[2]{[1..Length(P[2])-4]},x->x+[0,-1,-1]),
        Filtered(P[3],x->x[1]=4),
        [],
        []
    ];

    P:=ContestedEdge(Interior);
    #############################
    return Interior;
end;
trefoil:=[[2,5],[1,3],[2,4],[3,5],[1,4]];;