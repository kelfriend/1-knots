#########################################################################################
######################## Regular CW-complex of knot complement ##########################
#########################################################################################
################# Input: a list of (signed) integer pairs corresponding #################
######################## to the lengths of horizontal segments of a knot/link ###########
######################## and the depth of each segment in 4-space. ###################### 
#########################################################################################
################ Output: a regular CW-complex representing the complement of ############
######################## the input knot/link. ###########################################
###################################################################################### kk
KnotComplement:=function(l)
    local 
        signless, len, Entrench,
        boundaries, Ductwork;

    signless:=List(l, x->[SignInt(x[1])*x[1], SignInt(x[2])*x[2]]);
    len:=Length(l);

# creates the upper and lower 'hemispheres' of the complement
# input is how many horizontal/vertical segments there are in the knot
    Entrench:=function(len)
        local 
            bound, i, x, bool,
            recursiveCell, j;
        
        bound:=[[],[],[],[],[]];
        
        for i in [1..2*((4*len)+2)]
            do # add the 0-cells
            Add(bound[1],[1,0]);
        od;

        for i in [1..2*((4*len)+2)]
            do # add all 1-cells to the lower hemisphere
            if not i mod ((4*len)+2)=0
                then
                if i<(4*len)+3
                    then
                    x:=((i-1) mod 4)+1;
                    bool:=true;
                else
                    x:=((i+1) mod 4)+1;
                    bool:=false;
                fi;
                if i=1 or i=(4*len)+3
                    then
                    Add(bound[2],[2,i,i+(4*len)+1]); # circumferential 1-cell
                    Add(bound[2],[2,i,i+(4*len)+1]); # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^
                    if bool
                        then
                        Add(bound[2],[2,i,i+1]); # add connecting
                    fi; # 1-cells to the lower hemisphere -only-
                elif x=2
                    then
                    Add(bound[2],[2,i,i+1]);
                    Add(bound[2],[2,i,i+2]);
                    Add(bound[2],[2,i,i+2]);
                elif x=3
                    then
                    Add(bound[2],[2,i,i+2]);
                    Add(bound[2],[2,i,i+2]);
                elif x=4
                    then
                    Add(bound[2],[2,i,i+1]);
                    if bool
                        then
                        Add(bound[2],[2,i,i+2]);
                    fi;
                fi;
            fi;
        od;

        Add(bound[3],[1,3]); # beginnings of the large 2-cells to the
        Add(bound[3],[2,3]); # left and right of each bar

        Add(bound[3],[1,2]); # underside of the lower hemisphere
        Add(bound[3],[(7*len)+4,(7*len)+5]); # topside of the upper hemisphere

        recursiveCell:=function(n,x,y)
            Add(bound[n+1][x],bound[n+1][x][Length(bound[n+1][x])]+y);
        end;

        for i in [1..len]
            do # add the exterior 2-cells
            recursiveCell(2,1,2); recursiveCell(2,1,5);

            recursiveCell(2,2,1); recursiveCell(2,2,3); 
            recursiveCell(2,2,2); recursiveCell(2,2,1);
        od; 

        for i in [0,(7*len)+2]
            do
            for j in [0..len-1]
                do # add the 2-cells to each trench
                if i>0
                    then
                    i:=i-j; # to account for the lack
                fi; # of connecting 1-cells above
                Add(bound[3],[4,5,7,9]+(7*j)+i);
                Add(bound[3],[4,6,8,9]+(7*j)+i);
            od;
        od;

        Add(bound[4],[1,2,3]);
        for i in [1..len]
            do # add the lower 3-cell
            recursiveCell(3,1,2);
        od;

        # the following just adds the length to the beginning of each boundary
        for i in bound{[3,4]}
            do
            for j in i
                do
                Add(j,Length(j),1);
            od;
        od;

        return bound;
    end;

    boundaries:=Entrench(len);

# creates the tubes connecting the hemispheres
# according to the input
    Ductwork:=function(bound)
        local
            pre, horizontal, vertical,
            i, x, y, c1, c2, side,
            temp, u, v, w;

        pre:=Length(bound[2])+0;

        Add(bound[2],[2,1,(4*len)+3]); # connects the hemispheres

        horizontal:=List([1..len],x->[]);
        Apply(horizontal,x->[2,4,3,5]+(Position(horizontal,x)-1)*4);
        vertical:=ShallowCopy(horizontal)+(4*len)+2;

        for i in [1..len]
            do # add the 1-skeleton of the tubes as well as the upper 
            x:=signless[i][1]; # hemisphere's connecting 1-cells
            y:=signless[i][2];

            Add(bound[2],[2,horizontal[i][1],vertical[x][1]]);
            Add(bound[2],[2,horizontal[i][2],vertical[x][2]]);
            c1:=ShallowCopy(vertical[x][1]);
            vertical[x]:=Difference(vertical[x],vertical[x]{[1,2]});

            Add(bound[2],[2,horizontal[i][3],vertical[y][1]]);
            Add(bound[2],[2,horizontal[i][4],vertical[y][2]]);
            c2:=ShallowCopy(vertical[y][1]);
            vertical[y]:=Difference(vertical[y],vertical[y]{[1,2]});

        od;
        
        Add(bound[2],[2,(4*len)+2,2*((4*len)+2)]); # connects the hemispheres

        for i in [0..(2*len)-1]
            do # add the 2-cells enclosing the tubes
            for side in [0,1]
                do
                Add(
                    bound[3],[
                        4,
                        Position(bound[2], # lower-hemispherical edge
                        [2,bound[2][pre+2+(2*i)][2],bound[2][pre+3+(2*i)][2]])
                        +side, # same vertices but multiple edges
                        Position(bound[2], # upper-hemispherical edge
                        [2,bound[2][pre+2+(2*i)][3],bound[2][pre+3+(2*i)][3]])
                        +side,
                        pre+3+(2*i), # first tubular edge
                        pre+4+(2*i) # second tubular edge
                        ]
                    );
            od;
        od;

        temp:=[3];
        Add(bound[2],[2,(4*len)+3,bound[2][pre+2][3]]);
        Add(bound[3],[4,3,pre+1,pre+2,Length(bound[2])]);
        for i in [1..len]
            do
            u:=temp[Length(temp)]+1;
            Add(temp,u);
            Add(bound[2],[2,bound[2][pre+2+4*(i-1)][3],bound[2][pre+4+4*(i-1)][3]]);
            Add(
                bound[3],[
                    4,
                    u,
                    pre+2+4*(i-1),
                    pre+4+4*(i-1),
                    Length(bound[2])
                    ]
                );

            v:=temp[Length(temp)]+5;
            Add(temp,v);
            Add(bound[2],[2,bound[2][pre+3+4*(i-1)][3],bound[2][pre+5+4*(i-1)][3]]);
            Add(
                bound[3],[
                    4,
                    v,
                    pre+3+4*(i-1),
                    pre+5+4*(i-1),
                    Length(bound[2])
                    ]
                );

            w:=temp[Length(temp)]+1;
            Add(temp,w);
            Add(bound[2],[2,bound[2][pre+3+4*(i-1)][3],bound[2][pre+6+4*(i-1)][3]]);
            Add(
                bound[3],[
                    4,
                    w,
                    pre+3+4*(i-1),
                    pre+6+4*(i-1),
                    Length(bound[2])
                    ]
                );
        od;

        return bound;
    end;

    boundaries:=Ductwork(boundaries);

    return boundaries; #RegularCWComplex(complement);
end;

# ABANDONED : (
# this method for constructing the complement of a knot is flawed.
# there doesn't seem to be a possible regular CW-structure to put on
# even the simplest example (e.g. the unknot [ [ 1, 2 ], [ 1, 2 ] ])
# see the below for a sample output:
#
# gap> K:=KnotComplement([[1,2],[1,2]]);
#[ [ [ 1, 0 ], [ 1, 0 ], [ 1, 0 ], [ 1, 0 ], [ 1, 0 ], [ 1, 0 ], [ 1, 0 ], 
#      [ 1, 0 ], [ 1, 0 ], [ 1, 0 ], [ 1, 0 ], [ 1, 0 ], [ 1, 0 ], [ 1, 0 ], 
#      [ 1, 0 ], [ 1, 0 ], [ 1, 0 ], [ 1, 0 ], [ 1, 0 ], [ 1, 0 ] ], 
#  [ [ 2, 1, 10 ], [ 2, 1, 10 ], [ 2, 1, 2 ], [ 2, 2, 3 ], [ 2, 2, 4 ], 
#      [ 2, 2, 4 ], [ 2, 3, 5 ], [ 2, 3, 5 ], [ 2, 4, 5 ], [ 2, 4, 6 ], 
#      [ 2, 6, 7 ], [ 2, 6, 8 ], [ 2, 6, 8 ], [ 2, 7, 9 ], [ 2, 7, 9 ], 
#      [ 2, 8, 9 ], [ 2, 8, 10 ], [ 2, 11, 20 ], [ 2, 11, 20 ], [ 2, 12, 13 ], 
#      [ 2, 12, 14 ], [ 2, 12, 14 ], [ 2, 13, 15 ], [ 2, 13, 15 ], 
#      [ 2, 14, 15 ], [ 2, 16, 17 ], [ 2, 16, 18 ], [ 2, 16, 18 ], 
#      [ 2, 17, 19 ], [ 2, 17, 19 ], [ 2, 18, 19 ], [ 2, 1, 11 ], 
#      [ 2, 2, 12 ], [ 2, 4, 14 ], [ 2, 3, 16 ], [ 2, 5, 18 ], [ 2, 6, 13 ], 
#      [ 2, 8, 15 ], [ 2, 7, 17 ], [ 2, 9, 19 ], [ 2, 10, 20 ], [ 2, 11, 12 ], 
#      [ 2, 12, 16 ], [ 2, 14, 18 ], [ 2, 14, 13 ], [ 2, 13, 17 ], 
#      [ 2, 15, 19 ], [ 2, 15, 20 ] ], 
#  [ [ 6, 1, 3, 5, 10, 12, 17 ], [ 10, 2, 3, 4, 7, 9, 10, 11, 14, 16, 17 ], 
#      [ 2, 1, 2 ], [ 2, 18, 19 ], [ 4, 4, 5, 7, 9 ], [ 4, 4, 6, 8, 9 ], 
#      [ 4, 11, 12, 14, 16 ], [ 4, 11, 13, 15, 16 ], [ 4, 20, 21, 23, 25 ], 
#      [ 4, 20, 22, 24, 25 ], [ 4, 26, 27, 29, 31 ], [ 4, 26, 28, 30, 31 ], 
#      [ 4, 5, 21, 34, 35 ], [ 4, 6, 22, 34, 35 ], [ 4, 7, 27, 36, 37 ], 
#      [ 4, 8, 28, 36, 37 ], [ 4, 12, 23, 38, 39 ], [ 4, 13, 24, 38, 39 ], 
#      [ 4, 14, 29, 40, 41 ], [ 4, 15, 30, 40, 41 ], [ 4, 3, 32, 33, 42 ], 
#      [ 4, 4, 33, 35, 43 ], [ 4, 9, 34, 36, 44 ], [ 4, 10, 34, 37, 45 ], 
#      [ 4, 11, 37, 39, 46 ], [ 4, 16, 38, 40, 47 ], [ 4, 17, 38, 41, 48 ] ], 
#  [ [ 5, 1, 2, 3, 5, 7 ] ], [  ] ]