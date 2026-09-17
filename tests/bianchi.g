# The reported BianchiGcomplex(-1) regression. Reduce BEFORE computing
# homology: HAP's resolution code extends the input element tables.
if LoadPackage("hap") = fail then
    Error("HAP and its dependencies must be installed for the Bianchi test.");
fi;
Read("fred2.g");
Read("tests/helpers.g");

Fred2TestBianchi := function()
    local K, KK, L, HomologyList, HK, HKK, HL;
    K:=BianchiGcomplex(-1);
    KK:=BarycentricSubdivision(K);
    L:=ReduceGComplex(KK);
    HomologyList:=function(X)
        local C;
        C:=TensorWithIntegers(FreeGResolution(X,4));
        return List([0..3],i->Homology(C,i));
    end;
    HK:=HomologyList(K);
    HKK:=HomologyList(KK);
    HL:=HomologyList(L);
    Print("Original Bianchi homology: ",HK,"\n");
    Print("Subdivided Bianchi homology: ",HKK,"\n");
    Print("Reduced Bianchi homology: ",HL,"\n");
    Fred2Check(HK,[[0],[2,2],[6],[4,24]],"Bianchi reference homology");
    Fred2Check(HKK,HK,"Bianchi subdivision preserves homology");
    Fred2Check(HL,HK,"Bianchi reduction preserves homology in degrees 0 through 3");
    Print("PASS: Bianchi homology regression\n");
end;

Fred2TestBianchi();
