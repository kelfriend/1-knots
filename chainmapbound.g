ChainMapOfKnotBoundaryToComplement:=function(iota)
    local 
        inclusion_preimage, map, C_star, D_star, Bound;

    inclusion_preimage:=function(bound)
        local
            bound1, bound2, len,
            1c1, 1c2, 2c2, inc_pre;

        bound1:=bound[1];
        bound2:=bound[2];

        len:=Length(bound1[1])/2;

        1c1:=bound1[2]*1;
        1c1:=List(
            1c1,
            x->Concatenation(
                [2],
                Set(
                    [x[2],x[3]]
                )
            )
        );

        2c2:=List(bound2[3],x->Concatenation([x[1]],Set(x{[2..Length(x)]})));

        inc_pre:=function(n,k)
            local
                ind, pos, 2cell;

            if n=0 then
                return Position(List([1..2*len],x->x+1+2*Int((x-1)/len)),k);
            elif n=1 then
                ind:=1;
                if k>1 and 1c2[k-1]=1c2[k] then
                    ind:=2;
                fi;
                pos:=Positions(1c1,1c2[k]);
                if pos<>[] then
                    return Positions(1c1,1c2[k])[ind];
                else
                    return [];
                fi;
            elif n=2 then
                2cell:=List(bound2[3][k]{[2..4]},x->inc_pre(1,x));
                2cell:=Concatenation([4],Set(2cell));
                return Position(bound1[3],2cell);
            else
                return [];
            fi;
        end;

        1c2:=bound2[2]*1;
        1c2:=List(
            1c2,
            x->Concatenation(
                [2],
                [inc_pre(0,x[2]),
                inc_pre(0,x[3])]
            )
        );

        return inc_pre;
    end;

    map:=inclusion_preimage(
        [iota!.source!.boundaries,iota!.target!.boundaries]
    );
    
    C_star:=ChainComplexOfUniversalCover(iota!.target,false); # disable the
    # discrete vector field for now

    Bound:={n,k}->List(
        C_star!.boundary(n,iota!.mapping(n,k)),
        x->[SignInt(x[1])*map(n-1,AbsInt(x[1])),x[2]]
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

    return [D_star,C_star,map];
end;