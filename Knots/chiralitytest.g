###############################################################################
########## Produces a list of all prime knots K with mirror image L ###########
########## where ##############################################################
###############################################################################
########## pi_1( K + K ) = pi_1( K + L ) ######################################
########## and ################################################################
########## K + K != K + L #####################################################
###############################################################################

#ChiralityTest:=function()

#    local
#        d, chiral, i, j, K, L, 
#        P, Q, YP, YQ, CP, CQ, LP,
#        LQ, invP, invQ;

    d:=[1,1,2,3,7,21,49,165,552];
    chiral:=[];

    for i in [1..Length(d)]
    do
        for j in [1..d[i]]
        do
            K:=PureCubicalKnot(i+2,j);
            L:=ReflectedCubicalKnot(K);

            P:=KnotSum(K,L);
            Q:=KnotSum(K,K);

            YP:=PureComplexComplement(P);
            YQ:=PureComplexComplement(Q);
            YP:=ContractedComplex(YP);
            YQ:=ContractedComplex(YQ);
            YP:=RegularCWComplex(YP);
            YQ:=RegularCWComplex(YQ);

            CP:=ChainComplexOfUniversalCover(YP);
            CQ:=ChainComplexOfUniversalCover(YQ);

            LP:=LowIndexSubgroupsFpGroup(CP!.group,5);
            LQ:=LowIndexSubgroupsFpGroup(CQ!.group,5);

            invP:=SortedList(List(LP,g->Homology(
                TensorWithIntegersOverSubgroup(CP,g),3)));

            invQ:=SortedList(List(LQ,g->Homology(
                TensorWithIntegersOverSubgroup(CQ,g),3)));
                    
            if not invP=invQ
                then
                Add(chiral,[i+2,j]);
            fi;

            Print([i+2,j]," ",invP=invQ,"\n");
        od;
    od;

    #return chiral;
    
#end;