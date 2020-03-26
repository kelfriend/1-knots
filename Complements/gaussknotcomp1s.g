#  Input: a Gauss code.
# Output: inclusion of regular CW-complexes from a 1-dimensional link
#         to a 3-dimensional ball.
EmbeddedKnot:=function(gauss)
    local bnd, cbnd, map, max, Disk;

    bnd:=[[],[],[]];
    cbnd:=[[],[],[],[],[]];
    map:=[[],[]];
    max:=Maximum(gauss);

    Disk:=function()
        local verts;

        # add the 0-cells of the disk
        bnd[1]:=List([1..3*max],x->[1,0]);
        cbnd[1]:=List([1..3*max+2],x->[1,0]); # the last two 0-cells
        map[1]:=List([1..3*max]); # are in cbnd only

        verts:=List([1..max],x->[1,2,3]+3*x-3);
        for i in [1..2*max-1] do
            if gauss[i]<0 then
                if gauss[i+1]<0
                    Add(
                        cbnd[2],
                        [
                            2,
                            verts[-gauss[i]][2],
                            verts[-gauss[i+1]][2]
                        ]
                    );
                    Add(
                        bnd[2],
                        [
                            2,
                            verts[-gauss[i]][2],
                            verts[-gauss[i+1]][2]
                        ]
                    );
                    Add(map[2],Length(cbnd[2]));
                else
                    Add(
                        cbnd[2],
                        [
                            2,
                            verts[-gauss[i]][2],
                            verts[gauss[i+1]][1]
                        ]
                    );
                    Add(
                        bnd[2],
                        [
                            2,
                            verts[-gauss[i]][2],
                            verts[gauss[i+1]][1]
                        ]
                    );
                    Add(map[2],Length(cbnd[2]));
                    Add(
                        cbnd[2],
                        [
                            2,
                            verts[gauss[i+1]][1],
                            verts[gauss[i+1]][2]
                        ]
                    );
                    Add(
                        cbnd[2],
                        [
                            2,
                            verts[gauss[i+1]][2],
                            verts[gauss[i+1]][3]
                        ]
                    );
                fi;
            else
                if gauss[i+1]<0 then

                else

                fi;
            fi;
    end;
    Disk();

    return bnd;
end;
k:=[ -1, 3, -2, 1, -3, 2 ];