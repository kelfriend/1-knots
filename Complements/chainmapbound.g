Read("~/proj/Knots/knotcompbound.g");
ChainMapOfKnotBoundaryToComplement:=function(arg...)
    local 
        iota, index, inv_mapping,
        inclusion_preimage, C_star, C_star_copy,
        Bound, D_star, H, D_star_tensor,
        C_star_tensor, D_i2p, C_p2i, alpha;

    iota:=arg[1];

    if not Length(arg)>1 then
        Error("please input an inclusion map and a positive integer");
    else
        index:=arg[2];
    fi;

    inv_mapping:=[[],[],[]];
    inv_mapping[1]:=List([1..iota!.source!.nrCells(0)],x->iota!.mapping(0,x));
    inv_mapping[2]:=List([1..iota!.source!.nrCells(1)],x->iota!.mapping(1,x));
    inv_mapping[3]:=List([1..iota!.source!.nrCells(2)],x->iota!.mapping(2,x));

    inclusion_preimage:={n,k}->Position(inv_mapping[n+1],k);
    
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

    H:=LowIndexSubgroupsFpGroup(C_star!.group,index);
    H:=H[Position(List(H,x->Index(C_star!.group,x)),index)];
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
