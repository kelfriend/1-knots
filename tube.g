#########################################################################################
######################## The Image of Shin Satoh's Tube Map #############################
#########################################################################################
################# Input: a list of signed integers corresponding ########################
######################## to a link where the signs indicate the 'depth' #################
######################## of each horizontal bar of the link in the 4th ##################
######################## dimension. #####################################################
#########################################################################################
################ Output: a pure cubical complex equivalent to ###########################
######################## Tube(the specified link). ######################################
##################################################################################### K.K
Tube:=function(D)

    local
        signless, gridSize, width, height, null,
        newX, newY, i, j, k, vert, midpoints;

# obtains the binary array dimensions
    signless:=List(D,x->[AbsoluteValue(x[1]),AbsoluteValue(x[2])]);
    gridSize:=Maximum(List(signless,x->Maximum(x)));

    width:=((gridSize-1)*8)+7;
    height:=((gridSize-1)*5)+3;

# constructs an empty array of appropriate size
# notes: - 3-d depth will always be 7 and 4-d depth will always be 3
#        - a buffer layer of 0s is added to each dimension
    null:=List([1..5],v->(
        List([1..9],w->(
            List([1..height+2],x->(
                List([1..width+2],y->0)))))));

# lists for changing grid co-ordinates
    newX:=List([1..gridSize],x->0);
    newY:=List([1..gridSize],x->0);
    newX[1]:=5;
    newY[1]:=3;

    for i in [2..gridSize]
        do
        newX[i]:=newX[i-1]+8;
        newY[i]:=newY[i-1]+5;
    od;

# adds the endpoints of each horizontal tube to the array
# if those points are on the topmost or bottommost tubes, they are connected
    for i in [1..Length(signless)]
        do
        for j in signless[i]
            do
            null[3][5][newY[Length(newY)+1-i]][newX[j]]:=1;
        od;
        if i=1 or i=gridSize
            then
            for k in [1..Length(null[3][5][newY[Length(newY)+1-i]])]
                do
                if k>newX[signless[i][1]] and k<newX[signless[i][2]]
                    then
                    null[3][5][newY[Length(newY)+1-i]][k]:=1;
                fi;
            od;
        fi;
    od;

# creates an empty copy of the array for the vertical tubes to be added to.
# adds an array entry at the midpoint of those endpoints with the same
# height coordinate.
    vert:=ShallowCopy(null)*0;
    midpoints:=[];

    for i in [1..Length(signless)]
        do
		dist:=AbsoluteValue(newY[signless[i][1]]-newY[signless[i][2]]);
        Add(midpoints,
        [newX[Length(newX)+1-i],AbsoluteValue(newY[signless[i][1]]-newY[signless[i][2]])],
        i);
    od;

	for i in midpoints
		do
	od;

	return midpoints;
end;