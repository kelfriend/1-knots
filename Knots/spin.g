####################################################################
################### Spinning CW-complexes ##########################
####################################################################
############ Input: inclusion map from a subcomplex U to ###########
################### its parent (n-dimensional) regular #############
################### CW-complex Y ###################################
########### Output: (n+1)-dimensional CW-complex S(Y) ##############
################### corresponding to the `spinning' ################
################### of the original complex about U ################
####################################################################
Spin:=function(iota)

    local 
          Y, U, Imi, n, k, bndY, S, bndS, bnd,
          quad2pair, count, i, j, ij, y, s, BND, 
          a, b, ia, ib, deleted, deletedbnd, m, 
          diffs, o, p, q, r, l;

    if not IsHapRegularCWMap(iota)
        then
        Error(
        "the input must be a map between regular CW-complexes.\n"
        );
    fi;

    Y:=iota!.target; # U < Y
    U:=iota!.source; # iota: U -> Y

    Imi:=List([1..Dimension(U)+1],x->[]); # Image of iota
    for n in [0..Dimension(U)] # in terms of Y's indexing
        do
        for k in [1..U!.nrCells(n)]
            do
            Imi[n+1][k]:=iota!.mapping(n,k);
        od;
    od;

    bndY:=Y!.boundaries;

    S:=[[[1,0],[1,0]],[[2,1,2],[2,2,1]],[[2,1,2]],[]];
    S:=RegularCWComplex(S); # 2-dimensional disk
    bndS:=S!.boundaries;

    bnd:=List([0..Dimension(Y)+Dimension(S)],x->[]);
    # eventual boundary list

    quad2pair:=[]; 
    # (i,y) x (j,s) |-> quad2pair[i][j][y][s]
    count:=List([1..1+Dimension(Y)+Dimension(S)],x->0);

    for i in [1..1+Dimension(Y)]
        do
        quad2pair[i]:=[];
        for j in [1..1+Dimension(S)]
            do
            quad2pair[i][j]:=[];
            ij:=i-1+j;
            for y in [1..Length(bndY[i])]
                do
                quad2pair[i][j][y]:=[];
                for s in [1..Length(bndS[j])]
                    do
                    count[ij]:=count[ij]+1;
                    quad2pair[i][j][y][s]:=[ij,count[ij]];
                od;
            od;
        od;
    od;

    for i in [1..1+Dimension(Y)]
        do
        for j in [1..1+Dimension(S)]
            do
            for y in [1..Length(bndY[i])]
                do
                if j<3 or (j=3 and not y in Imi[i]) # filters out
                    then # those cells which are the product of
                    # S's 2-cell and something in the image of U

                    for s in [1..Length(bndS[j])]
                        do
                        BND:=[0];
                        if i>1
                            then
                            a:=bndY[i][y];
                            for ia in a{[2..Length(a)]}
                                do
                                Add(BND,
                                quad2pair[i-1][j][ia][s][2]);
                                BND[1]:=BND[1]+1;
                            od;
                        fi;
                        if j>1
                            then 
                            b:=bndS[j][s];
                            for ib in b{[2..Length(b)]}
                                do
                                Add(BND,
                                quad2pair[i][j-1][y][ib][2]);
                                BND[1]:=BND[1]+1;
                            od;
                        fi;
                        bnd[i-1+j][quad2pair[i][j][y][s][2]]:=BND;
                    od;
                fi;
            od;
        od;
    od;

    bnd[1]:=List(bnd[1],i->[1,0]);

    deleted:=List([1..Length(bnd)],x->[]); # tracks removed cells

    deletedbnd:=List([1..Length(bnd)],x->[]); # removes empty spaces
    for n in [1..Length(bnd)] # left in bnd after filtering
        do
        for m in [1..Length(bnd[n])]
            do
            if IsBound(bnd[n][m])
                then
                Add(deletedbnd[n], bnd[n][m]);
            else
                Add(deleted[n],m);
            fi;
        od;
    od;
    
    # lastly, perform the reindexing to allow for construction of a
    # regular CW-complex
    for p in [4..Length(deletedbnd)]
        do
        for q in [1..Length(deletedbnd[p])]
            do
            l:=0;
            for r in [2..Length(deletedbnd[p][q])]
                do
                l:=Length(Filtered(deleted[p-1],
                                   x->x<deletedbnd[p][q][r]
                                   ));
                deletedbnd[p][q][r]:=deletedbnd[p][q][r]-l;
            od;
        od;
    od;

    if not deletedbnd[Length(deletedbnd)]=[]
        then
        Add(deletedbnd,[]);
    fi;

    return deletedbnd; #RegularCWComplex(deletedbnd);

end;
####################################################################
########### Test complex (i): [ ][ ][ ]   [x][x][x] ################
############################# [ ][ ][ ] < [x][x][x] ################
############################# [x][x][x]   [x][x][x] ################
####################################################################
####################### (ii): Complement of trefoil ################
############################# and a 2-dimensional ##################
############################# subspace #############################
####################################################################
# (i)
Y:=[[1,1,1],[1,1,1],[1,1,1]];;
Y:=PureCubicalComplex(Y);;

U:=[[0,0,0],[0,0,0],[1,1,1]];;
U:=PureCubicalComplex(U);;

psi:=RegularCWMap(Y,U);;

# (ii)
K:=PureCubicalKnot(3,1);;
K:=PureComplexComplement(K);;

k:=0*K!.binaryArray;;
for cell in [1..Length(K!.binaryArray[1][1])]
    do
    if K!.binaryArray[1][1][cell]=1
        then
        k[1][1][cell]:=1;
    fi;
od;

k:=PureCubicalComplex(k);;
k:=PureComplexComplement(k);;

omicron:=RegularCWMap(K,k);;