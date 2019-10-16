ChainMapOfKnotBoundaryToComplement:=function(arg...)
    local 
        iota, index, mapping,
        inclusion_preimage, C_star,
        Bound, D_star, H, D_star_tensor,
        C_star_tensor;

    iota:=arg[1];

    if Length(arg)>1 then
        index:=arg[2];
    fi;

    mapping:=[[],[],[]];
    mapping[1]:=List([1..iota!.source!.nrCells(0)],x->iota!.mapping(0,x));
    mapping[2]:=List([1..iota!.source!.nrCells(1)],x->iota!.mapping(1,x));
    mapping[3]:=List([1..iota!.source!.nrCells(2)],x->iota!.mapping(2,x));

    inclusion_preimage:={n,k}->Position(mapping[n+1],k);
    
    C_star:=ChainComplexOfUniversalCover(iota!.target,false); # disable the
    # discrete vector field for now
    # Check universalcover.gi and disable the contraction !!! (done, line 224)

    Bound:={n,k}->List(
        C_star!.boundary(n,iota!.mapping(n,k)),
        x->[SignInt(x[1])*inclusion_preimage(n-1,AbsInt(x[1])),x[2]]
    );
    
    D_star:=Objectify(
        HapEquivariantChainComplex,
        rec(
            dimension:=iota!.source!.nrCells,
            boundary:=Bound,
            elts:=C_star!.elts,
            group:=C_star!.group,
            properties:=[
                ["dimension",2],
                ["characteristic",0],
                ["length",EvaluateProperty(C_star,"length")]
            ]
        )
    );

    H:=LowIndexSubgroupsFpGroup(D_star!.group,index)[index];
    D_star_tensor:=TensorWithIntegersOverSubgroup(D_star,H);
    C_star_tensor:=TensorWithIntegersOverSubgroup(C_star,H);

    # output: chain map from D_star_tensor -> C_star_tensor 
    return [D_star_tensor,C_star_tensor];
end;
