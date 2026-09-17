# The reported regression: direct reduction of SL3Zs must preserve homology.
# This test is slower than the small integration fixtures; run explicitly.
if LoadPackage("hap") = fail then
    Error("HAP and its dependencies must be installed for the SL3Zs test.");
fi;
Read("fred2.g");
Read("tests/helpers.g");

Fred2TestSL3Zs := function()
    local K, L, d, CK, CL, HK, HL;
    K:=ContractibleGcomplex("SL3Zs");
    d:=3;
    CK:=TensorWithIntegers(FreeGResolution(K,d+2));
    HK:=List([0..d+1],i->Homology(CK,i));
    Print("Original SL3Zs homology: ",HK,"\n");
    L:=ReduceGComplex(K);
    CL:=TensorWithIntegers(FreeGResolution(L,d+2));
    HL:=List([0..d+1],i->Homology(CL,i));
    Print("Reduced SL3Zs homology: ",HL,"\n");
    Fred2Check(HL,HK,"SL3Zs reduction preserves homology in degrees 0 through 4");
    Print("PASS: SL3Zs homology regression\n");
end;

Fred2TestSL3Zs();
