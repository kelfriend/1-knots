ChainMapOfKnotBoundaryToComplement:=function(iota)
    local 
        mapping, inclusion_preimage, C_star, Bound, D_star;

    mapping:=[[],[],[]];
    mapping[1]:=List([1..iota!.source!.nrCells(0)],x->iota!.mapping(0,x));
    mapping[2]:=List([1..iota!.source!.nrCells(1)],x->iota!.mapping(1,x));
    mapping[3]:=List([1..iota!.source!.nrCells(2)],x->iota!.mapping(2,x));

    inclusion_preimage:={n,k}->Position(mapping[n+1],k);
    
    C_star:=ChainComplexOfUniversalCover(iota!.target,false); # disable the
    # discrete vector field for now

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

    return [D_star,C_star,inclusion_preimage];
end;
