DeclareGlobalFunction("SimplifiedRegularGCWComplex_OneStep");

############################################################################
##
InstallOtherMethod(SimplifiedComplex,
"Simplified non-free resolution",
[IsHapNonFreeResolution],
function(K)
local L,new_sz,sz;

L:=StructuralCopy(K);
sz:=0;
new_sz:=-1;

while new_sz < sz do
    sz:=Sum(List([0..Length(L)],i->L!.dimension(i)));
    L:=SimplifiedRegularGCWComplex_OneStep(L);
    new_sz:=Sum(List([0..Length(L)],i->L!.dimension(i)));
od;

return L;

end);
##
############################################################################

############################################################################
##
InstallOtherMethod(Length,
"The dimension of GCW complex",
[IsHapNonFreeResolution],
function(K)
local d;
d:=-1;
while true do
   if K!.dimension(d+1)>0 then 
      d:=d+1;
   else
      break;
   fi;
od;
return d;
end);
##
############################################################################

############################################################################
##
InstallGlobalFunction(SimplifiedRegularGCWComplex_OneStep,
function(K)
local L, kStar, IsGBalanced, SimplifiedUnion, RemoveUnion;

L:=StructuralCopy(K);

#############################
##
## Inputs the degree n of a cell, 
## the position e of the cell in the list of n-cells, 
## and an integer k>0 such that n+k<=Length(K).
## returns a list of all cells of degrees in the range n..n+k that touch e^n.
kStar:=function(e,n,k)
local star, dim, f, bnd, cell;
star := List([1..k+1], i -> []);
star[1]:=[e,K!.stabilizer(n,e)];

for dim in [1..k] do
    for f in [1..L!.dimension(n+dim)] do
        bnd := L!.boundary(n+dim, f);
        for cell in bnd do
            if AbsoluteValue(cell[1]) in (star[dim]) then
                Add(star[dim+1], f);
                break;
            fi;
        od;
    od;
od;
return star;. 
end;
##
#############################

#############################
##
## Inputs a list U of cells and returns true if U is G-balanced, 
## false otherwise.
IsGBalanced:=function(U)

return 42;
end;
##
#############################

#############################
##
## Inputs a list U of cells and simplifies then to a list of cells
## V. It returns V.
SimplifiedUnion:=function(U)

return 42;
end;
##
#############################

#############################
##
## Inputs a list of cells U and cell f to which U is homeomorphic.
## The function does not return anything. 
## Rather, it modifies the cell structure of L, replacing the cells in U 
## by a single cell f, sets the stabilizer group of f, and sets the +/-1 
## action of G on f.
RemoveUnion:=function(U,f)


end;
##
#############################

## Here we loop through all orbit representatives e in K and for each 
## e, of degree n say, and each suitable k constructs U:= kStar(e,n,k),
## and tests if U is G-balanced and homeomorphic to a single cell f. 
## If such a U and f is found then RemoveUnion(U,f) is applied, after 
## which we break out of the k loop and break out of the e loop. 

return L;
end);
##
############################################################################
