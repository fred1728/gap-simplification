Read("fred2.g");

#####################################################################
#####################################################################
TestGReduce:=function(K)
local d,KK,L,CK,CKK,CL,HK,HKK,HL;

######Find dimension d of K#######
d:=0;
while K!.dimension(d)<>0 do
    d:=d+1;
od;
d:=d-1;
##################################

KK:=BarycentricSubdivision(K);
L:=ReduceGComplex(KK,d);

CK:= TensorWithIntegers(FreeGResolution(K,d+2));
HK:=List([0..d+1],i->Homology(CK,i)); Print(HK,"\n");
CKK:= TensorWithIntegers(FreeGResolution(KK,d+2));
HKK:=List([0..d+1],i->Homology(CKK,i)); Print(HKK,"\n");
CL:= TensorWithIntegers(FreeGResolution(L,d+2));
HL:=List([0..d+1],i->Homology(CL,i));
Print(HL,"\n");
return HK=HL;
end;
#####################################################################
#####################################################################

