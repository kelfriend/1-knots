LiftedChainMap:=function(Y,H)
    local
        U, G, L, n;

        U:=UniversalCover(Y);
        G:=U!.group;
        
        if IsInt(H) then
            L:=LowIndexSubgroupsFpGroup(G,n);
            n:=Position(List(L,x->Index(G,x)),H);
            H:=L[n];
        fi;

        L:=EquivariantCWComplexToRegularCWComplex(U,H);
        L:=BoundaryMap(L);
        L:=ChainMap(L);

    return L;
end; 
BadLiftedChainMap:=function(arg...)
    local 
        iota, subgroup, inv_mapping,
        inclusion_preimage, C_star, C_star_copy,
        Bound, D_star, H, D_star_tensor,
        C_star_tensor, D_i2p, C_p2i, alpha;

    iota:=arg[1];

    if not Length(arg)>1 then
        Error("please input an inclusion map and a finite index subgroup");
    else
        subgroup:=arg[2];
    fi;

    if EvaluateProperty(iota,"image")=fail then
        inv_mapping:=[[],[],[]];
        inv_mapping[1]:=List([1..iota!.source!.nrCells(0)],x->iota!.mapping(0,x));
        inv_mapping[2]:=List([1..iota!.source!.nrCells(1)],x->iota!.mapping(1,x));
        inv_mapping[3]:=List([1..iota!.source!.nrCells(2)],x->iota!.mapping(2,x));
        Add(iota!.properties,["image",inv_mapping]);
    fi;

    inclusion_preimage:={n,k}->Position(EvaluateProperty(iota,"image")[n+1],k);
    
    C_star:=ChainComplexOfUniversalCover(iota!.target,false);
    C_star_copy:=ChainComplexOfUniversalCover(iota!.target,false);

    Bound:={n,k}->List(
        C_star_copy!.boundary(n,iota!.mapping(n,k)),
        x->[SignInt(x[1])*inclusion_preimage(n-1,AbsInt(x[1])),x[2]]
    );

    D_star:=Objectify(
        HapEquivariantChainComplex,
        rec(
            dimension:=iota!.source!.nrCells,
            boundary:=Bound,
            elts:=ShallowCopy(C_star!.elts),
            group:=C_star!.group,
            properties:=[
                ["dimension",2],
                ["characteristic",0],
                ["length",EvaluateProperty(C_star,"length")]
            ]
        )
    );

    if IsInt(subgroup) then # obtains the first subgroup of that index in the
        H:=LowIndexSubgroupsFpGroup(C_star!.group,index); # below list
        H:=H[Position(List(H,x->Index(C_star!.group,x)),index)];
    else
        H:=index;
    fi;

    C_star_tensor:=TensorWithIntegersOverSubgroup(C_star,H);
    D_star_tensor:=TensorWithIntegersOverSubgroup(D_star,H);
    D_i2p:=EvaluateProperty(D_star_tensor,"int2pair");
    C_p2i:=EvaluateProperty(C_star_tensor,"pair2int");

    alpha:=function(v,n)
        local
            vec, i, pair, int, inc;

        vec:=List([1..C_star_tensor!.dimension(n)],x->0);
        for i in [1..Length(v)] do
            if v[i]<>0 then
                pair:=D_i2p(i);
                inc:=SignInt(pair[1])*iota!.mapping(n,AbsInt(pair[1]));
                int:=C_p2i([inc,pair[2]]);
                vec[int]:=v[i];
            fi;
        od;
        
        return vec;
    end;
#BoundaryMatrix(C_star_tensor,0); BoundaryMatrix(C_star_tensor,1);
#BoundaryMatrix(C_star_tensor,2); BoundaryMatrix(C_star_tensor,3);
#BoundaryMatrix(C_star_tensor,4);
    return Objectify(
        HapChainMap,
        rec(
            source:=D_star_tensor,
            target:=C_star_tensor,
            mapping:=alpha,
            properties:=[
                ["type","chainMap"],
                ["characteristic",0]
            ]
       )
    );
end;
################################################################################
############################## CHANGE LOG ######################################
################################################################################
##### File: ####### RegularCWComplexes/universalcover.gi (line 224) ############
################################################################################ 
##### Before: ##### X:=ContractedComplex(arg[1]); ##############################
################################################################################
##### Now: ######## X:=arg[1]; #################################################
################################################################################
################################################################################
################################################################################
##### File: ####### Resolutions/resInfSubgroup.gi (line 131) ###################
################################################################################
##### Before: ##### ["reduced",false] ])); #####################################
################################################################################
##### Now: ######## ["reduced",false], #########################################
################### ["pair2int",Pair2Int], #####################################
################### ["int2pair",Int2Pair] ])); #################################
################################################################################
################################################################################
################################################################################
##### File: ####### Functors/tensorWithZ.gi (line 66) ##########################
################################################################################
##### Before: ##### EvaluateProperty(R,"characteristic")] ])); #################
################################################################################
##### Now: ######## EvaluateProperty(R,"characteristic")], #####################
################### ["pair2int",EvaluateProperty(R,"pair2int")], ###############
################### ["int2pair",EvaluateProperty(R,"int2pair")] ])); ###########
################################################################################
