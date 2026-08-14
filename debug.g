k:=2;;
n:=0;;
Elts := ShallowCopy(L!.elts);;
newbnddup:=[];;

pos := function(g)
local p;
    p := Position(Elts, g);
    if p = fail then
        Add(Elts, g);
        return Length(Elts);
    fi;
    return p;
end;

#signs:=SimplifyLattice(L,kstar,n)[2];
# for c in [1..Length(U[k+1])] do
#     cell:=U[k+1][c];
#     bndcell:=L!.boundary(n+k,cell[1]);
#     for bnd in bndcell do
#         h:=cell[2]*Elts[bnd[2]];
#         instar:=false;
#         candidates:=Filtered(U[k],y->y[1]=AbsoluteValue(bnd[1]));
#         for cellstar in candidates do
#             if Inverse(h)*cellstar[2] in cellstar[3] then #changed from h*Inverse(cellstar[2]) to the way it is now
#                 instar:=true;
#                 break;
#             fi;
#         od;
#         if not instar then
#             #sign:=L!.action(n+k,cell[1],pos(cell[2]));
#             #sign:=SignInt(SimplifyLattice(L,2star,0)[2][c]);
#             Add(newbnddup,[signs[c]*bnd[1],pos(h)]);
#         fi;
#     od;
# od;

# newbnd:=[];;
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

#Action(g)
#ausfuehren bevor die Indizes in newbnd aktualisiert werden!
#funktioniert noch nicht weil die Vorzeichen in newbnd noch nicht korrekt sind
# g:=L!.elts[j];
# bndcell:=newbnd[1];
# Stabbnd:=L!.stabilizer(n+k-1,AbsoluteValue(bndcell[1]));
# g1:=L!.elts[bndcell[2]];
# for bndcell2 in Filtered(newbnd,p->AbsoluteValue(p[1])=AbsoluteValue(bndcell[1])) do    
#     g2:=L!.elts[bndcell2[2]];
#     h:=Inverse(g2)*g*g1;
#     if h in Stabbnd then
#         #Print(bndcell2,"\n");
#         #Print(g2, "\n",h,"\n");
#         eps1:=L!.action(n+k-1,AbsoluteValue(bndcell[1]),pos(h)); #eps1=1 if g.g1x=g2x, eps1=-1 if g.g1x=-g2x 
#         eps2:=SignInt(bndcell[1])*SignInt(bndcell2[1]); #eps2=1 iff g1x and g2x have the same sign in d(star)
#         eps3:=eps1*eps2; #eps3=action(star,g)
#         break;
#     fi;
# od;

# EulerChar := function(L)
# local chi, d, i;
#     chi := 0;
#     d := 0;
#     while L!.dimension(d) > 0 do
#         for i in [1..L!.dimension(d)] do
#             chi := chi + (-1)^d / Order(L!.stabilizer(d,i));
#         od;
#         d := d+1;
#     od;
#     return chi;
# end;

# Ka := Replace2(K,  1, 1, 2);;   # step 1: 2-star of cell 1 in dimension 1
# Kb := Replace2(Ka, 2, 1, 1);;   # step 2: 1-star of cell 2 in dimension 1
# Kc := Replace2(Kb, 5, 1, 1);;   # step 3: 1-star of cell 5 in dimension 1
# Print("Chi(K)=",EulerChar(K),"\n");
# Print("Chi(Ka)=",EulerChar(Ka),"\n");
# Print("Chi(Kb)=",EulerChar(Kb),"\n");
# Print("Chi(K1)=",EulerChar(K1),"\n");
# Print("H3(K x Z)=",Homology(TensorWithIntegers(FreeGResolution(K,10)),3),"\n");
# Print("H3(Ka x Z)=",Homology(TensorWithIntegers(FreeGResolution(Ka,10)),3),"\n");
# Print("H3(Kb x Z)=",Homology(TensorWithIntegers(FreeGResolution(Kb,10)),3),"\n");
# Print("H3(K1 x Z)=",Homology(TensorWithIntegers(FreeGResolution(K1,10)),3),"\n");

# #zwischenlager
# OldFaceIndex := function(j)
#     return oldindex(perm[k], j, n+k-1);
# end;

# ## ReduceTerm(j,h) returns [rep, eps]:
# ##   rep is the canonical representative index: the smallest p such
# ##       that Elts[h] * Elts[p]^-1 lies in Stab(n+k-1, oldj) -- i.e.
# ##       Elts[h] and Elts[rep] are in the same RIGHT coset of the
# ##       stabilizer (Stab * Elts[h] = Stab * Elts[rep]). This matches
# ##       the coset convention already used elsewhere in this file
# ##       (e.g. the original dedup loop above, and the "instar" test
# ##       in kstar-style code): terms are identified via h1*h2^-1 in
# ##       Stab, NOT h1^-1*h2. (A left-coset version was tried and
# ##       numerically verified to give the WRONG answer against the
# ##       complex's own boundary data -- e.g. it fails to reproduce
# ##       the known fact that Spiegel*Dreh acts by -1 on the 2-cell's
# ##       own boundary chain in the d8 example. Right cosets + a
# ##       RIGHT group action, as used here and in ActOnNewBoundary,
# ##       reproduce it correctly.)
# ##   eps is the sign such that (Elts[h] tensor 1) = eps*(Elts[rep]
# ##       tensor 1): writing Elts[h] = t*Elts[rep] with
# ##       t := Elts[h]*Elts[rep]^-1 in Stab(oldj), eps is the action
# ##       of t on the face cell.
# ## Because pos()/oldpos() only ever appends new elements at the end
# ## of Elts, the choice of "smallest index" is stable once computed.
# ReduceTerm := function(j, h)
# local oldj, stab, p, t;

#     oldj := OldFaceIndex(j);
#     stab := L!.stabilizer(n+k-1, oldj);

#     for p in [1..Length(Elts)] do
#         t := Inverse(Elts[h]) * Elts[p]; #changed to inverting the first term instead of the second
#         if t in stab then
#             return [ p, L!.action(n+k-1, oldj, oldpos(t)) ];
#         fi;
#     od;

#     Error("ReduceTerm: no representative found (should be unreachable)");
# end;

# ## NormalizeBoundary turns a raw list of terms [±j,h] into a genuine
# ## canonical form for the underlying chain: terms are grouped by
# ## (j, canonical right-coset representative), their signed
# ## multiplicities are summed (including the orientation twist eps
# ## from ReduceTerm), zero-coefficient classes are dropped, and
# ## surviving classes are expanded back into repeated [±j,rep] entries
# ## (matching the convention used elsewhere in this file) before being
# ## sorted. Two chains that are equal as elements of the chain group
# ## ALWAYS produce identical output from this function, and this
# ## equality is preserved under applying any g in G to every term via
# ## RIGHT multiplication first (see ActOnNewBoundary) -- right cosets
# ## Stab*h are automatically preserved under h -> h*g for ANY g,
# ## since (h1*g)*(h2*g)^-1 = h1*h2^-1, with no extra condition on g.
# NormalizeBoundary := function(bnd)
# local out, term, j, sgn, h, red, rep, eps, entry, matched, result, cnt;

#     out := [];  # entries of the form [ j, rep, coeff ]

#     for term in bnd do
#         j   := AbsoluteValue(term[1]);
#         sgn := SignInt(term[1]);
#         h   := term[2];

#         red := ReduceTerm(j, h);
#         rep := red[1];
#         eps := red[2];

#         matched := false;
#         for entry in out do
#             if entry[1] = j and entry[2] = rep then
#                 entry[3] := entry[3] + sgn*eps;
#                 matched := true;
#                 break;
#             fi;
#         od;
#         if not matched then
#             Add(out, [ j, rep, sgn*eps ]);
#         fi;
#     od;

#     out := Filtered(out, x -> x[3] <> 0);

#     result := [];
#     for entry in out do
#         for cnt in [1..AbsoluteValue(entry[3])] do
#             Add(result, [ SignInt(entry[3]) * entry[1], entry[2] ]);
#         od;
#     od;

#     SortBy(result, x -> [ AbsoluteValue(x[1]), x[2], SignInt(x[1]) ]);
#     return result;
# end;

# ## ActOnNewBoundary(g) computes newbnd . g directly: g translates
# ## each basis element (j,h) to (j, h*g) -- RIGHT multiplication,
# ## matching the right-coset convention in ReduceTerm/NormalizeBoundary
# ## above (same face-orbit j, new group-element index). Any orientation
# ## twist or coset-representative change is picked up automatically by
# ## NormalizeBoundary/ReduceTerm, since h*g may land in a different
# ## right coset (with a different, possibly sign-flipped, canonical
# ## representative) than h did. pos() extends Elts if the new group
# ## element hasn't been seen before.
# ActOnNewBoundary := function(g)
# local image, term, j, sgn, newh;

#     image := [];
#     for term in newbnd do
#         j    := AbsoluteValue(term[1]);
#         sgn  := SignInt(term[1]);
#         newh := pos(g*Elts[term[2]]); #changed from right multiplication by g to left multiplication
#         Add(image, [ sgn * j, newh ]);
#     od;

#     return NormalizeBoundary(image);
# end;

# ## action(d,i,j): sign by which Elts[j] acts on cell i of dimension d.
# ## For the new top cell, this is determined by checking whether
# ## g translates the boundary chain newbnd to itself (+1), to its
# ## negative (-1), or to something else (an error - the star was not
# ## actually G-balanced, or an upstream bug).

#Todo: function which takes as input a cell [e,g] with the property that some translate of e is in the star of f
# and decides in which translate of the star the cell is
WhichTranslate:=function(star,cell,n,k)
local bdcell, dim;

currentcell:=cell
g:=L!.elts[cell[2]]
for i in [1..k] do
    dim:=n+k-i+1;
    bdcell:=L!.boundary(dim,AbsoluteValue(currentcell[1]));



end;
