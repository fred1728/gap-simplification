#Generates the k-star of the e-th cell in dimension n in the equivariant complex L 
#Output: The i-th entry of the output is a list of cells in dimension k+i of the star. The entries have the form [cell orbit,]
kStar:=function(L,e,n,k)
local star, stardup, dim, f, bnd, cell, Li, g, g1, h, id, cell2, Stabf, Stabcell, K1, K2, l,T,Stabe;

stardup := List([1..k+1], i -> []);
Stabe:=L!.stabilizer(n,e);
id := Identity(L!.elts[1]);
stardup[1]:=[[e,id,Stabe]];
dim:=1;

for f in [1..L!.dimension(n+dim)] do
    bnd := List(L!.boundary(n+dim, f),x-> [AbsoluteValue(x[1]),x[2]]);
    for cell in bnd do
        Li:= Filtered(stardup[dim],x -> x[1]=cell[1]);
        g:=L!.elts[cell[2]];
        for cell2 in Li do
            Stabf:=L!.stabilizer(n+dim,f);
            Stabcell:=L!.stabilizer(n+dim-1,cell[1]);
            K1:=Intersection(Stabf,Stabcell);
            K2:=RightTransversal(Stabcell,K1);
            for l in K2 do
                g1:= Inverse(l)*Inverse(g);
                Add(stardup[dim+1],[f,g1,Stabf]);
            od;
        od;
    od;
od;

for dim in [2..k] do
    for f in [1..L!.dimension(n+dim)] do
        bnd := List(L!.boundary(n+dim, f),x-> [AbsoluteValue(x[1]),x[2]]);
        for cell in bnd do
            Li:= Filtered(stardup[dim],x -> x[1]=cell[1]);
            g:=L!.elts[cell[2]];
            for cell2 in Li do
                h:=cell2[2];
                Stabf:=L!.stabilizer(n+dim,f);
                Stabcell:=L!.stabilizer(n+dim-1,cell[1]);
                K1:=Intersection(Stabf,Stabcell);
                K2:=RightTransversal(Stabcell,K1);
                for l in K2 do
                    g1:= h*Inverse(l)*Inverse(g);
                    Add(stardup[dim+1],[f,g1,Stabf]);
                od;
            od;
        od;
    od;
od;

star := List([1..k+1], i -> []);
star[1]:=[[e,id,Stabe]];
for dim in [2..k+1] do
    for cell in stardup[dim] do
        T:=true;
        for cell2 in Filtered(star[dim],x->x[1]=cell[1]) do
            Stabcell:=L!.stabilizer(n+dim-1,cell[1]);
            if Inverse(cell[2])*cell2[2] in Stabcell then
                T:=false;
                break;
            fi;
        od;
        if T then
            Add(star[dim],cell);
        fi;
    od;
od;

return star;
end;
##########################################################################

ForgetG:=function(L,U,n)
local FL, Y, k, d, cell, boundarycell, dcell, orbits, grouped, LO, i, g, g1, h, gh, t, r, z, Le, Stab, glist, j;

k:=Length(U)-1;
FL:=[];
FL[1]:=[[1,0]];

for d in [2..k+1] do
    Le:=Length(U[d]);
    FL[d]:=List([1..Le], i -> []);
    for r in [1..Le] do
        cell:=U[d][r];
        t:=0; #t counts how many cells are in the boundary
        FL[d][r][1]:=1; #Platzhalter
        boundarycell:=L!.boundary(d+n-1,cell[1]);
        orbits:=Set(List(U[d-1],x->x[1]));
        grouped:=List(orbits, y -> Filtered(U[d-1], x -> x[1] = y));
        for dcell in boundarycell do
            Stab:=L!.stabilizer(d+n-2,AbsoluteValue(dcell[1]));
            h:=L!.elts[dcell[2]];
            g:=cell[2];
            gh:=g*h;
            LO:=Length(orbits);
            z:=0;
            for i in [1..LO] do
                if i>1 then
                    z:=z+Length(grouped[i-1]); #counts number of cell
                fi;
                if AbsoluteValue(dcell[1])=orbits[i] then
                    glist:=(List(grouped[i],y->y[2]));
                    for j in [1..Length(glist)] do
                        g1:=glist[j];
                        if Inverse(g1)*gh in Stab then #getauscht von g1*Inverse(gh) zu Inverse(g1)*gh
                            Add(FL[d][r],z+j);
                            t:=t+1;
                        fi;
                    od;
                fi;
            od;
        od;
        FL[d][r][1]:=t;
    od;
od;

Add(FL,[]);
#Y:=RegularCWComplex(FL);
#cobnd:=Y!.coboundaries;

#return [FL,cobnd];
return FL;
end;
##########################################################################

## Inputs a list U of lists of cells and the dimension n of the lowest dimensional cell and returns the pair [true, Stabilizer group of U] if U is G-balanced, 
## [false, Stabilizer group of U]  otherwise.
PreIsGBalanced:=function(L,e,n,kk) #(L,U,n)->(L,e,n,k)
local orbits, grouped, orb, disj, stable, stableh, orbg, g, h, l, balanced, Stab, k, d, Ud, StabU, U;

U:=kStar(L,e,n,kk);
StabU:=[Identity(L!.elts[1])];
balanced:=true;
for g in L!.elts do
    disj:=true;
    k:=Length(U)-1;
    if g in U[1][1][3] then
    #if g in L!.stabilizer(n,U[1][1]) then
        disj:=false;
    fi;
    for d in [2..k+1] do
        if disj then
            Ud:=U[d];
            orbits:=Set(List(Ud,x->x[1]));
            grouped := List(orbits, y -> Filtered(Ud, x -> x[1] = y));
            for orb in grouped do
                orbg:=List(orb,x->x[2]);
                for h in orbg do
                    for l in orbg do
                        Stab:=orb[1][3];
                        #Stab:=L!.stabilizer(n+d-1,orb[1][1]);
                        if g*h*Inverse(l) in Stab then
                            disj:=false;
                            break;
                        fi;
                    od;
                od;
            od;
        fi;
    od;
    stable:=true;
    if not g in U[1][1][3] then
    #if not g in L!.stabilizer(n,U[1][1]) then
        stable:=false;
    fi;
    for d in [2..k+1] do
        if stable then
            Ud:=U[d];
            orbits:=Set(List(Ud,x->x[1]));
            grouped := List(orbits, y -> Filtered(Ud, x -> x[1] = y));
            for orb in grouped do
                Stab:=orb[1][3];
                orbg:=List(orb,x->x[2]);
                for h in orbg do
                    stableh:=false;
                    for l in orbg do                       
                        #Stab:=L!.stabilizer(n+d-1,orb[1][1]);
                        if g*h*Inverse(l) in Stab then
                            stableh:=true;
                            break;
                        fi;
                    od;
                    if not stableh then
                        stable:=false;
                        break;
                    fi;
                od;
            od;
        fi;
    od;
    if stable then
        Add(StabU,g);
    fi;
    balanced:=stable or disj;
    if not balanced then
        break;
    fi;
od;
StabU:=Group(StabU);
return [balanced,StabU];
#return balanced;
end;

IsGBalanced:=function(L,e,n,k)
return PreIsGBalanced(L,e,n,k)[1];
end;

Stabstar:=function(L,e,n,k)
return PreIsGBalanced(L,e,n,k)[1];
end;
##########################################################################

# d1:=1;;
# n:=1;;
# #FL3:=[[[1,0],[1,0],[1,0]],[[2,1,2],[2,1,3]],[]];;
# #FL3:=[[[1,0]],[[1,1],[1,1]],[]];;
# #FL3:=[[[1,0]],[[1,1],[1,1],[1,1],[1,1]],[[2,1,2],[2,2,3],[2,3,4],[2,4,1]],[]];
# FL3:=[ [ [ 1, 0 ] ], [ [ 1, 1 ], [ 1, 1 ] ], [ [ 2, 1, 2 ], [ 2, 1, 2 ] ], [  ],];;
# Y3:=RegularCWComplex(FL3);;
# bnd:=Y3!.boundaries;;
# cobnd:=Y3!.coboundaries;;


# ShouldRemove:=function(FL,d1,n)
# local d2, d3, cob, V1, V2, V3, bool, bnd, cobnd, Y;

# bnd:=FL;
# Y:=RegularCWComplex(FL);
# cobnd:=Y!.coboundaries;
# bool:=true;
# d2:=d1+1;;
# d3:=d2+1;;
# cob:=1*cobnd[d1][n];

# if not (cobnd[d1][n][1] =2 and bnd[d1][n][1]>0) then 
#     bool:=false;
# fi;
# d2:=d1+1;
# d3:=d2+1;
# cob:=1*cobnd[d1][n];
# if not SortedList(cobnd[d2][cob[2]])= SortedList(cobnd[d2][cob[3]]) then 
#     bool:=false;
# fi;
# V1:=BoundaryOfRegularCWCell(Y,d1-1,n);
# V2:=BoundaryOfRegularCWCell(Y,d1,cobnd[d1][n][2]);
# V3:=BoundaryOfRegularCWCell(Y,d1,cobnd[d1][n][3]);
# if not 1+Size(V1)=Size(Intersection(V2,V3)) then
#     bool:= false; 
# fi;

# return bool;
# end;

EntireBoundary:=function(bnd,d,n)
local V, U, D, tmp, v , cells;

if d=0 then return []; fi;

D:=d+1;
V:=StructuralCopy(bnd[D][n]);
V:=V{[2..Length(V)]};
V:=SSortedList(V);
cells:=List(V,i->[D-2,i]);
D:=D-1;

while D>1 do
tmp:=[];
for v in V do
U:=StructuralCopy(bnd[D][v]);
U:=U{[2..Length(U)]};
Append(tmp,U);
Append(cells,List(U,i->[D-2,i]));   #Changed D-1 to D-2  (May 2022)
od;
tmp:=SSortedList(tmp);
cells:=SSortedList(cells);
V:=tmp;
D:=D-1;
od;

return cells;
end;

Cobound := function(bnd)
local dim, cobnd, b, k, i, j, n;

dim:=PositionProperty(bnd,IsEmpty)-2;
cobnd:=[];; #Cobnd[n+1] contains the info on n-cells.
for n in [0..dim] do
  #k:=2*(n+1)+1;#k:=1+2^(n+1);
  cobnd[n+1]:=List(bnd[n+1],i->[0]);
  for j in [1..Length(bnd[n+2])] do
    b:=bnd[n+2][j];
    k:=Length(b);
    for i in b{[2..k]} do
      cobnd[n+1][i][1]:=cobnd[n+1][i][1]+1;
      Add(cobnd[n+1][i],j);
    od;
  od;
od;
cobnd[dim+1]:=List(bnd[dim+1],a->[0]);
return cobnd;
end;

ShouldRemove := function(FL,d1,n)
local d2, cob, V1, V2, V3, bnd, cobnd, Y;

bnd := FL;

#Y := RegularCWComplex(FL);
cobnd := Cobound(bnd);

if n > Length(bnd[d1]) then
    return false;
fi;

if not (cobnd[d1][n][1] = 2 and bnd[d1][n][1] > 0) then
    return false;
fi;

d2 := d1+1;
cob := 1*cobnd[d1][n];

if Length(cob) < 3 then
    return false;
fi;

if d2 > Length(cobnd) then
    return false;
fi;

if cob[2] > Length(cobnd[d2]) or cob[3] > Length(cobnd[d2]) then
    return false;
fi;

if not SortedList(cobnd[d2][cob[2]]) = SortedList(cobnd[d2][cob[3]]) then
    return false;
fi;


# Print("d1=", d1, " n=", n, "\n");
# Print("d2=", d2, "\n");
# Print("cob=", cob, "\n");
# Print("Length(bnd)=", Length(bnd), "\n");
# Print("Length(bnd[d1])=", Length(bnd[d1]), "\n");
# Print("Length(bnd[d2])=", Length(bnd[d2]), "\n");

V1:=EntireBoundary(bnd,d1-1,n);
V2:=EntireBoundary(bnd,d1,cobnd[d1][n][2]);
V3:=EntireBoundary(bnd,d1,cobnd[d1][n][3]);

if not 1+Size(V1) = Size(Intersection(V2,V3)) then
    return false;
fi;

return true;
end;


JoinCells:=function(FL,d1,n)
local bnd, cobnd, cob, d, d2, d3, m, t, poss, s, pos, u, perm, cnt, F, x, i;

d2:=d1+1;
d3:=d1+2;
bnd:=FL;
cobnd:=Cobound(bnd);
#cobnd:=Y!.coboundaries;
cob:=1*cobnd[d1][n];
##
##REDUCE COBOUNDARIES OF BOUNDARIES OF nTH CELL
if d1>1 then
d:=d1-1;;
for m in bnd[d1][n]{[2..Length(bnd[d1][n])]} do
    t:=1*cobnd[d][m]{[2..Length(cobnd[d][m])]};;
    poss:=Position(t,n);;
    Remove(t,poss);;
    cobnd[d][m]:=Concatenation([Length(t)],t);;
od;;
fi;;
#Print("cbnd ",cobnd,"\n");
##
##COBOUNDARIES OF BOUNDARIES REDUCED

##
##REMOVE nTH CELL, ITS COBOUNDARY, AND ADJUST ITS PRESENCE IN BOUNDARIES
##OF ITS COBOUNDARIES
bnd[d1][n]:=[0];;
cobnd[d1][n]:=[0];;
#Print("cobnd",cobnd,"\n");
#Print("bnd",bnd,"\n");

s:=bnd[d2][cob[2]];;
s:=s{[2..Length(s)]};;
pos:=Position(s,n);;
Remove(s,pos);;


t:=bnd[d2][cob[3]];;
t:=t{[2..Length(t)]};;
pos:=Position(t,n);;
Remove(t,pos);;


u:=SortedList(Concatenation(s,t));;

if IsEmpty(u) then
    bnd[d2][cob[2]]:=[1,0];;
else
    bnd[d2][cob[2]]:=Concatenation([Length(u)],u);;
fi;;

##
##nTH CELL AND ITS PRESENCE REMOVED


##
##FOR SECOND CELL OF DIMENSION n+1  REDUCE THE BOUNDARIES OF ITS 
##COBOUNDARIES

for m in cobnd[d2][cob[3]]{[2..Length(cobnd[d2][cob[3]])]} do
t:=bnd[d3][m]{[2..Length(bnd[d3][m])]};;

pos:=Position(t,cob[3]);;
Remove(t,pos);;
bnd[d3][m]:=Concatenation([Length(t)],t);;
od;;


##
##SECOND CELL  BOUNDARIES OF ITS COBOUNDARIES REDUCED

##
##REMOVE PRESENCE OF SECOND CELL IN COBOUNDARIES OF ITS BOUNDARIES
for m in bnd[d2][cob[3]]{[2..Length(bnd[d2][cob[3]])]} do
if cobnd[d1][m][1]>0 then

t:=cobnd[d1][m]{[2..Length(cobnd[d1][m])]};;

pos:=Position(t,cob[3]);;
t[pos]:=cob[2];;
t:=SSortedList(t);;
cobnd[d1][m]:=Concatenation([Length(t)],t);;

fi;;
od;;


bnd[d2][cob[3]]:=[0];;

cobnd[d2][cob[3]]:=[0];;

perm:=[];  ##perm[d][i] will be the new position of
           ##the old i-th cell of dimension d;

for d in [1..Length(bnd)] do
perm[d]:=[];
cnt:=0;
for n in [1..Length(bnd[d])] do
if bnd[d][n][1]=0 then cnt:=cnt+1;
else
perm[d][n]:=n-cnt;
fi;
od;
od;



for d in [1..Length(bnd)] do
F:=1*Filtered([1..Length(bnd[d])],i->not bnd[d][i][1]=0);
bnd[d]:=1*((bnd[d]){F});

if d>1 then
d1:=d-1;
for x in bnd[d] do
    if x[2]>0 then
        for i in [2..Length(x)] do
            x[i]:=1*perm[d1][x[i]];
        od;
    fi;
od;
fi;
od;

#if all lower-dimensional cells have been deleted, move all the remaining cells downwards
while Length(bnd) > 0 and IsEmpty(bnd[1]) do
    Remove(bnd, 1);
od;

return bnd;
end;
##########################################################################

SimplifyLattice := function(L,kstar,n) #n=dimension of the center of the kstar
local Changed, b, d, j, m, h, g, g1, g2, dim, face, cell, signs, faceorbit, Stab, topcells, involvedcell1, involvedcell2,
ind1, ind2, sign1, sign2, incidence1, incidence2, mergecell1, mergecell2, FLgroups, FLcodim1, permFL, FL, cob;

FL:=ForgetG(L,kstar,n);
permFL:=StructuralCopy(FL);
FL:=StructuralCopy(FL);
FLcodim1:=[1..Length(FL[Length(FL)-2])];
FLgroups:=List([1..Length(FL[Length(FL)-1])], i -> [i]);
dim:=Length(FL)-1;

repeat
    Changed:=false;
    #topcells:=
    for d in [1 .. Length(FL)-2] do
        m:=1;
        while m <= Length(FL[d]) do
            if FL[d][m][1] > 0 and ShouldRemove(FL,d,m) then
                if d=Length(FL)-2 then #for top dimensional cells, we have to keep track of signs
                    mergecell1:=Cobound(FL)[d][m][2];
                    mergecell2:=Cobound(FL)[d][m][3]; #these are the two cells which will be merged
                    #adjust signs
                    b:=FLcodim1[m]; #old index of the cell that will be deleted
                    ind1 := fail;
                    involvedcell1 := fail;
                    ind2 := fail;
                    involvedcell2 := fail; 
                    for j in FLgroups[mergecell1] do
                        if b in permFL[dim][AbsoluteValue(j)]{[2..Length(permFL[dim][AbsoluteValue(j)])]} then
                            #involvedcell1:=permFL[dim][j];
                            ind1:=AbsoluteValue(j);
                            sign1:=SignInt(j);
                            break;
                        fi;
                    od;
                    for j in FLgroups[mergecell2] do
                        if b in permFL[dim][AbsoluteValue(j)]{[2..Length(permFL[dim][AbsoluteValue(j)])]} then
                            #involvedcell2:=permFL[dim][j];
                            ind2:=AbsoluteValue(j);
                            sign2:=SignInt(j);
                            break;
                        fi;
                    od;
                    #Print("Merging the cells ", FLgroups[mergecell1], " and ", FLgroups[mergecell2]," along the face ",b,"\n", " the components involved are ",ind1," and ",ind2,"\n");
                    if ind1 = fail then
                        Error("Could not find face ", b, " in group ", FLgroups[mergecell1]);
                    fi;
                    face:=kstar[dim-1][b];
                    faceorbit:=face[1];
                    g:=face[2];
                    involvedcell1:=kstar[dim][ind1][1];
                    g1:=kstar[dim][ind1][2];
                    involvedcell2:=kstar[dim][ind2][1];
                    g2:=kstar[dim][ind2][2];
                    for cell in Filtered(L!.boundary(n+dim-1,involvedcell1),p->AbsoluteValue(p[1])=faceorbit) do
                        h:=L!.elts[cell[2]];
                        Stab:=kstar[dim-1][b][3];
                        if Inverse(h)*Inverse(g1)*g in Stab then
                            incidence1:=SignInt(cell[1]);
                            break;
                        fi;
                    od;
                    for cell in Filtered(L!.boundary(n+dim-1,involvedcell2),p->AbsoluteValue(p[1])=faceorbit) do
                        h:=L!.elts[cell[2]];
                        Stab:=kstar[dim-1][b][3];
                        if Inverse(h)*Inverse(g2)*g in Stab then
                            incidence2:=SignInt(cell[1]);
                            break;
                        fi;
                    od;
                    if sign1*incidence1*sign2*incidence2 = 1 then
                        FLgroups[mergecell2]:=List(FLgroups[mergecell2], y -> -y);
                    fi;
                    Append(FLgroups[mergecell1],FLgroups[mergecell2]);
                    Remove(FLgroups,mergecell2);
                    Remove(FLcodim1,m);
                fi;

                # NEW: if this JoinCells call will delete a cell from the
                # codim-1 level (FL[Length(FL)-2]) via its d2-side removal,
                # keep FLcodim1 in sync.
                if d = Length(FL)-3 then
                    cob := Cobound(FL)[d][m];
                    FL:=JoinCells(FL,d,m);
                    Remove(FLcodim1, cob[3]);
                else
                    FL:=JoinCells(FL,d,m);
                fi;
                
                # FL:=JoinCells(FL,d,m);
                Changed:=true;               
                m:=1;
            else
                m:=m+1;
            fi;

        od;
    od;
until not Changed;

signs := List([1..Length(FLgroups[1])], x -> 0);
for j in FLgroups[1] do
    signs[AbsoluteValue(j)] := SignInt(j);
od;
return [FL,signs];
end;

HomeomorphicToDisk:=function(L,e,n,k)
local kstar, B, simplified;

B:=false;
kstar:=kStar(L,e,n,k);
simplified:=SimplifyLattice(L,kstar,n)[1];
if Length(simplified)=2 and Length(simplified[1])=1 then
    B:=true;
fi;
return B;
end;
##########################################################################

#Replace the k-star of the e-th n-dimensional cell in L by a single cell
Replace2:=function(L,e,n,k)
local perm, i, j, d, l, topdim, matched, bndnk1, cnt, orbits, bool, StabU, sign, Elts, g, g1, lg, poslg, c, cell, pos,
newbnddup, newbnd, bnd, bndcell, h, instar, candidates, cellstar, incell, Boundary, 
Stab, Dim, action, oldpos, oldindex, RemapBoundary, oldi, signs, found,
kstar, ReduceTermOld, NormalizeBoundaryOld, matches;

Elts := ShallowCopy(L!.elts);
kstar:=kStar(L,e,n,k);

pos := function(g)
local p;
    p := Position(Elts, g);
    if p = fail then
        Add(Elts, g);
        return Length(Elts);
    fi;
    return p;
end;

oldpos := function(g)
local p;
    p := Position(L!.elts, g);
    if p = fail then
        Add(L!.elts, g);
        return Length(L!.elts);
    fi;
    return p;
end;

#ki 09/07
## ReduceTermOld(j,h): j is an OLD (pre-perm) face index at level n+k-1,
## h is an index into Elts. Returns [rep, eps] where rep is the smallest
## index such that Elts[h] and Elts[rep] lie in the same left coset of
## Stab(n+k-1,j).
ReduceTermOld := function(j, h)
local stab, p, t;
    stab := L!.stabilizer(n+k-1, j);
    for p in [1..Length(Elts)] do
        t := Inverse(Elts[h]) * Elts[p];
        if t in stab then
            return [ p, L!.action(n+k-1, j, oldpos(t)) ];
        fi;
    od;
    Error("ReduceTermOld: no representative found (should be unreachable)");
end;

## NormalizeBoundaryOld(bnd): takes a raw list of [±j,h] terms (old
## indices), groups them by (j, canonical right-coset rep), sums signed
## multiplicities (folding in the eps twist), drops classes that cancel
## to zero, and expands surviving classes back into repeated [±j,rep]
## entries. This is what correctly cancels a coset that shows up twice
## with opposite signs -- which the old same-signed-value dedup could
## not do.
NormalizeBoundaryOld := function(bnd)
local out, term, j, sgn, h, red, rep, eps, entry, matched, result, cnt;
    out := [];
    for term in bnd do
        j   := AbsoluteValue(term[1]);
        sgn := SignInt(term[1]);
        h   := term[2];
        red := ReduceTermOld(j, h);
        rep := red[1];
        eps := red[2];
        matched := false;
        for entry in out do
            if entry[1] = j and entry[2] = rep then
                entry[3] := entry[3] + sgn*eps;
                matched := true;
                break;
            fi;
        od;
        if not matched then
            Add(out, [ j, rep, sgn*eps ]);
        fi;
    od;
    out := Filtered(out, x -> x[3] <> 0);
    result := [];
    for entry in out do
        for cnt in [1..AbsoluteValue(entry[3])] do
            Add(result, [ SignInt(entry[3]) * entry[1], entry[2] ]);
        od;
    od;
    return result;
end;
#09/07 ging bis hier

perm:=[];  ##perm[i][j] will be the new position of
           ##the old j-th cell of dimension n+i-1;
for i in [1..k+1] do
    if kstar[i]=[] then break;
    fi;
    perm[i]:=[];
    cnt:=0;
    orbits:=List(kstar[i],x->x[1]);
    for j in [1..L!.dimension(n+i-1)] do
        if j in orbits then cnt:=cnt+1;
        else perm[i][j]:=j-cnt;
        fi;
    od;
od;
#nrremovedtopcells=cnt
topdim:=L!.dimension(n+k)-cnt+1;

#ki #???Fehler: d hat keine Funktion
oldindex := function(p, i, d)
local q;
    q := Position(p, i);
    if q = fail then
        Error("No old index corresponding to new index ",
              i, " in dimension ", d);
    fi;
    return q;
end;

#adjust boundaries of n+k+1-cells
orbits:=List(kstar[k+1],x->x[1]);
signs:=SimplifyLattice(L,kstar,n)[2]; #ki: needed early to correctly sign bndnk1 below
bndnk1:=[];
for j in [1..L!.dimension(n+k+1)] do
    bndnk1[j]:=[];
    for cell in L!.boundary(n+k+1,j) do
        if AbsoluteValue(cell[1]) in orbits then
            for c in [1..Length(kstar[k+1])] do
                if kstar[k+1][c][1]=AbsoluteValue(cell[1]) then
                    g1:=kstar[k+1][c][2];
                    lg:=L!.elts[cell[2]]*Inverse(g1);
                    poslg:=pos(lg);
                    Add(bndnk1[j],[SignInt(cell[1])*signs[c]*topdim,poslg]);
                    break;
                fi;
            od;
        else
            Add(bndnk1[j],[SignInt(cell[1])*perm[k+1][AbsoluteValue(cell[1])],cell[2]]);
        fi;
    od;
od;

newbnddup:=[];
#signs:=SimplifyLattice(L,kstar,n)[2];
for c in [1..Length(kstar[k+1])] do
    cell:=kstar[k+1][c];
    bndcell:=L!.boundary(n+k,cell[1]);
    for bnd in bndcell do
        h:=cell[2]*Elts[bnd[2]];
        instar:=false;
        candidates:=Filtered(kstar[k],y->y[1]=AbsoluteValue(bnd[1]));
        for cellstar in candidates do
            if Inverse(h)*cellstar[2] in cellstar[3] then #changed from h*Inverse(cellstar[2]) to the way it is now
                instar:=true;
                break;
            fi;
        od;
        if not instar then
            #sign:=L!.action(n+k,cell[1],pos(cell[2]));
            #sign:=SignInt(SimplifyLattice(L,kstar,n)[2][c]);
            Add(newbnddup,[signs[c]*bnd[1],pos(h)]);
        fi;
    od;
od;
#Print(newbnddup, "\n");

# newbnd:=[];
# for cell in newbnddup do
#     bool:=true;
#     for incell in Filtered(newbnd,y->y[1]=cell[1]) do
#         if Inverse(Elts[incell[2]])*Elts[cell[2]] in L!.stabilizer(n+k-1,AbsoluteValue(cell[1])) then
#         # changed from if Elts[incell[2]]*Inverse(Elts[cell[2]]) in L!.stabilizer(n+k-1,AbsoluteValue(cell[1])) then
#             bool:=false;
#             break;
#         fi;
#     od;
#     if bool then
#         Add(newbnd,cell);
#     fi;
# od;
newbnd := NormalizeBoundaryOld(newbnddup); #ki


newbnd:=List(newbnd, 
    cell -> [SignInt(cell[1]) * perm[k][AbsoluteValue(cell[1])],cell[2]]);

#ki
RemapBoundary := function(b, p, d, i)
local term, out, newidx;
    out := [];
    for term in b do
        newidx := p[AbsoluteValue(term[1])];
        if newidx = fail then
            Error("Boundary remap failed at d=", d, " i=", i,
                  " old face=", AbsoluteValue(term[1]));
        fi;
        Add(out, [ SignInt(term[1]) * newidx, term[2] ]);
    od;
    return out;
end;


Boundary := function(d,i)
local oldi, b;

    if d<n or d>n+k+1 then
        return L!.boundary(d,i);

    elif d=n then
        # boundaries in degree n go to degree n-1, which was not changed
        oldi := oldindex(perm[d-n+1], i, d);
        return L!.boundary(d, oldi);

    elif d>n and d<n+k then
        oldi := oldindex(perm[d-n+1], i, d);
        b := L!.boundary(d, oldi);
        return RemapBoundary(b, perm[d-n], d, i);

    elif d=n+k then
        if AbsoluteValue(i)<topdim then
            oldi := oldindex(perm[d-n+1], AbsoluteValue(i), d);
            b := L!.boundary(d, oldi);
            b := RemapBoundary(b, perm[d-n], d, AbsoluteValue(i));
            if i<0 then
                b := List(b, t -> [-t[1],t[2]]);
            fi;
            return b;
        elif AbsoluteValue(i)=topdim then
            if i<0 then
                return List(newbnd, t -> [-t[1],t[2]]);
            fi;
            return newbnd;
        fi;

    elif d=n+k+1 then
        return bndnk1[i];
    fi;
end;

StabU:=Stabstar(L,e,n,k);
Stab:=function(d,i)
    if d<n or d>n+k+1 then
        return L!.stabilizer(d,i);
    elif d>=n and d<n+k then
        return L!.stabilizer(d,Position(perm[d-n+1],AbsoluteValue(i)));
    elif d=n+k then
        if AbsoluteValue(i)<topdim then
            return L!.stabilizer(d,Position(perm[d-n+1],AbsoluteValue(i)));
        elif AbsoluteValue(i)=topdim then
            return StabU;
        fi;
    elif d=n+k+1 then
        return L!.stabilizer(d,i);
    fi;
end; 



Dim:=function(d)
local p;
    if d<n or d>=n+k+1 then
        return L!.dimension(d);
    elif d>=n and d<n+k then
        p := Compacted(perm[d-n+1]);
        if Length(p)=0 then
            return 0;
        fi;
        return Maximum(p);
    elif d=n+k then
        p := Compacted(perm[d-n+1]);
        if Length(p)=0 then
            return 1;
        fi;
        return Maximum(p)+1;
    fi;
end;

#changes start here



action := function(d,i,j)
local g, bndcell, bndcell2, Stabbnd, g1, g2, h, eps1, eps2, eps3, found, oldb;
    g := Elts[j];

    if d = n+k and AbsoluteValue(i) = Dim(d) then
        if not (g in StabU) then
            #Print("TOP-CELL: g not in StabU. d=",d," i=",i," g=",g,"\n");
            return 0;
        fi;
        bndcell:=newbnd[1];
        oldb := oldindex(perm[k], AbsoluteValue(bndcell[1]), n+k-1);
        Stabbnd:=L!.stabilizer(n+k-1,oldb); #changed from L!.stabilizer(n+k-1,AbsoluteValue(bndcell[1]));
        g1:=L!.elts[bndcell[2]];       
        found:=false;
        for bndcell2 in Filtered(newbnd,p->AbsoluteValue(p[1])=AbsoluteValue(bndcell[1])) do    
            g2:=L!.elts[bndcell2[2]];
            h:=Inverse(g2)*g*g1;
            if h in Stabbnd then
        #Print(bndcell2,"\n");
        #Print(g2, "\n",h,"\n");
                eps1:=L!.action(n+k-1,oldb,pos(h)); #eps1=1 if g.g1x=g2x, eps1=-1 if g.g1x=-g2x 
                eps2:=SignInt(bndcell[1])*SignInt(bndcell2[1]); #eps2=1 iff g1x and g2x have the same sign in d(star)
                eps3:=eps1*eps2; #eps3=action(star,g)
                found:=true;
                break;
            fi;
        od;
        if not found then
            Print("FAILED: g = ", g, "\n");
            Print("g1 = ", g1, "\n");
            Print("newbnd = ", newbnd, "\n");
            Error("no matching coset found");
        fi;
        return eps3;

    elif d >= n and d <= n+k then
        return L!.action(d, oldindex(perm[d-n+1], AbsoluteValue(i), d), oldpos(g));

    else
        return L!.action(d, AbsoluteValue(i), oldpos(g));
    fi;
end;


return Objectify(HapNonFreeResolution,
        rec(
            dimension := Dim,
            boundary := Boundary,
            stabilizer := Stab,
            action := action,
            elts := Elts,
            group := L!.group,
            homotopy := fail,
            properties := [
                ["characteristic",0],
                ["type","resolution"],
                ["length",1000]
            ]
        )
    );
end;
##########################################################################

#Iteratively simplifies a G-complex L by repeatedly finding a cell e in some
#dimension n whose k-star is
#   (i)  G-balanced      -- IsGBalanced(L,kStar(L,e,n,k),n)[1] = true, and
#   (ii) a disk           -- HomeomorphicToDisk(L,e,n,k) = true,
#and replacing that star by a single new cell via Replace(L,e,n,k).
#
#The search restarts from dimension 0 after every successful replacement,
#since Replace renumbers cells across several dimensions at once and it is
#not safe to keep old indices around. This continues until a full pass over
#the whole complex finds nothing left to replace.
#

ReduceGComplex:=function(L,d)
local CurrentTopDim, Progress, n, e, k, dim, kmax;
 
# L:=arg[1];
# if Length(arg)>1 and arg[2]>0 then
#     MaxK:=arg[2];
# else
#     MaxK:=infinity;
# fi;
 
#Highest dimension in which L currently has any cells.
CurrentTopDim:=function()
local d;
    d:=0;
    while L!.dimension(d+1) > 0 do
        d:=d+1;
    od;
    return d;
end;
 
repeat
    Progress:=false;
    n:=0;
    while n <= CurrentTopDim() and not Progress do
        dim:=L!.dimension(n);
        e:=1;
        while e <= dim and not Progress do
            kmax:=CurrentTopDim()-n;
            k:=1;
            while k <= kmax and not Progress do
                if IsGBalanced(L,e,n,k) and HomeomorphicToDisk(L,e,n,k) and n+k=d then
                    Print("ReduceGComplex: replacing the ",k,
                          "-star of cell ",e," in dimension ",n,"\n");
                    L:=Replace2(L,e,n,k);
                    Progress:=true;
                fi;
                k:=k+1;
            od;
            e:=e+1;
        od;
        n:=n+1;
    od;
until not Progress;
 
return L;
end;
##########################################################################

#Iteratively simplifies a G-complex L by repeatedly finding a cell e in some
#dimension n whose k-star is
#   (i)  G-balanced      -- IsGBalanced(L,kStar(L,e,n,k),n)[1] = true, and
#   (ii) a disk           -- HomeomorphicToDisk(L,e,n,k) = true,
#and replacing that star by a single new cell via Replace(L,e,n,k).
#
#The search restarts from dimension 0 after every successful replacement,
#since Replace renumbers cells across several dimensions at once and it is
#not safe to keep old indices around. This continues until a full pass over
#the whole complex finds nothing left to replace.
#

ReduceGComplex:=function(L,d)
local CurrentTopDim, Progress, n, e, k, dim, kmax;
 
# L:=arg[1];
# if Length(arg)>1 and arg[2]>0 then
#     MaxK:=arg[2];
# else
#     MaxK:=infinity;
# fi;
 
#Highest dimension in which L currently has any cells.
CurrentTopDim:=function()
local d;
    d:=0;
    while L!.dimension(d+1) > 0 do
        d:=d+1;
    od;
    return d;
end;
 
repeat
    Progress:=false;
    n:=0;
    while n <= CurrentTopDim() and not Progress do
        dim:=L!.dimension(n);
        e:=1;
        while e <= dim and not Progress do
            kmax:=CurrentTopDim()-n;
            k:=1;
            while k <= kmax and not Progress do
                if IsGBalanced(L,e,n,k) and HomeomorphicToDisk(L,e,n,k) and n+k=d then
                    Print("ReduceGComplex: replacing the ",k,
                          "-star of cell ",e," in dimension ",n,"\n");
                    L:=Replace2(L,e,n,k);
                    Progress:=true;
                fi;
                k:=k+1;
            od;
            e:=e+1;
        od;
        n:=n+1;
    od;
until not Progress;
 
return L;
end;
