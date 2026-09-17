Read("fred2.g");
Read("tests/helpers.g");

L := Fred2Interval(false);;
U := kStar(L,1,0,1);;
Fred2Check(List(U,Length), [1,2], "star cell counts");
Fred2Check(List(U[2],x->x[1]), [1,2], "star edge orbits");
Fred2Check(List(U[2],x->x[2]), [(),()], "star representatives");
Fred2Check(ForgetG(L,U,0), [[[1,0]],[[1,1],[1,1]],[]],
           "relative incidence lattice");
Fred2Check(IsGBalanced(L,1,0,1), true, "trivial action is balanced");
Fred2Check(Stabstar(L,1,0,1), Group(()), "star stabilizer");
FL := ForgetG(L,U,0);;
Fred2Check(Cobound(FL), [[[2,1,2]],[[0],[0]]], "coboundaries");
Fred2Check(EntireBoundary(FL,0,1), [], "vertex boundary");
Fred2Check(EntireBoundary(FL,1,1), [[0,1]], "edge boundary");
Fred2Check(ShouldRemove(FL,1,1), true, "shared center is removable");
Fred2Check(ShouldRemove(FL,1,2), false, "missing cell is rejected");
Fred2Check(ShouldRemove(FL,2,1), false, "top cell is not removable");
Fred2Check(JoinCells(StructuralCopy(FL),1,1), [[[1,0]],[]],
           "joining edges produces one cell");
Fred2Check(SimplifyLattice(L,U,0), [[[[1,0]],[]],[1,1]],
           "simplification preserves compatible signs");
Fred2Check(ForgetG(L,U,0), FL, "simplification preserves input");
Fred2Check(HomeomorphicToDisk(L,1,0,1), true, "interval star is a disk");
Fred2Check(HomeomorphicToDisk(L,2,0,1), false,
           "endpoint star does not pass the collapse criterion");
L := Fred2Interval(true);;
Fred2Check(SimplifyLattice(L,kStar(L,1,0,1),0)[2], [1,-1],
           "reversed edge changes the merge sign");

# A nontrivial stabilizer produces two distinct translates of one edge.
L := Fred2Interval(false);;
L.elts := [(),(1,2)];;
L.group := Group((1,2));;
L.stabilizer := function(d,i)
    if d = 0 and i = 1 then return Group((1,2)); fi;
    return Group(());
end;;
U := kStar(L,1,0,1);;
Fred2Check(List(U,Length), [1,4], "stabilizer coset enumeration");
Fred2Check(IsGBalanced(L,1,0,1), true, "invariant nontrivial star");
Fred2Check(Stabstar(L,1,0,1), Group((1,2)), "nontrivial star stabilizer");
# Keep the translated edges but give the center a trivial stabilizer:
# translation now overlaps the star without preserving its center.
L.stabilizer := function(d,i) return Group(()); end;;
L.boundary := function(d,i)
    if d = 1 then return [[1,1],[1,2]]; fi;
    return [];
end;;
Fred2Check(IsGBalanced(L,1,0,1), false, "overlapping translate is unbalanced");

# The element table need not generate a cell's full stabilizer. Here it
# contains a subgroup of order 4 while the actual stabilizer has order 12.
L := Fred2Interval(false);;
L.group := Group((1,2,3),(4,5,6,7));;
L.elts := [(),(4,5,6,7)];;
L.stabilizer := function(d,i) return Group((1,2,3),(4,5,6,7)); end;;
Fred2Check(IsGBalanced(L,1,0,1),true,"balance with an incomplete element table");
Fred2Check(Stabstar(L,1,0,1),L.group,"full star stabilizer, including missing 3-factor");
Add(L.elts,(1,2,3));
Fred2Check(Stabstar(L,1,0,1),L.group,"star stabilizer independent of element table growth");

# An edge stabilizer outside the center stabilizer creates an overlapping
# translate even when that element is absent from the element table.
L := Fred2Interval(false);;
L.group := Group((1,2));;
L.stabilizer := function(d,i)
    if d=1 then return Group((1,2)); fi;
    return Group(());
end;;
Fred2Check(IsGBalanced(L,1,0,1),false,"detect overlap outside the element table");

# Noncommuting representatives: the actual translated-edge stabilizer is
# g*Stab(edge)*g^-1, not the untranslated subgroup.
L := Fred2Interval(false);;
L.group := SymmetricGroup(3);;
L.elts := [(),(1,3,2)];;
L.stabilizer := function(d,i)
    if d=0 and i=1 then return Group((1,2)); fi;
    return Group((1,2)^(1,2,3));
end;;
L.boundary := function(d,i)
    if d=1 and i=1 then return [[-1,2],[2,1]]; fi;
    if d=1 and i=2 then return [[1,2],[-3,1]]; fi;
    return [];
end;;
Fred2Check(IsGBalanced(L,1,0,1),true,"balance uses conjugated cell stabilizers");
Fred2Check(Stabstar(L,1,0,1),Group((1,2)),"noncommuting star has the center stabilizer");

# A triangle tests transitive boundary traversal across two dimensions.
FL := [[[1,0],[1,0],[1,0]], [[2,1,2],[2,2,3],[2,1,3]],
       [[3,1,2,3]], []];;
Fred2Check(EntireBoundary(FL,2,1),
           [[0,1],[0,2],[0,3],[1,1],[1,2],[1,3]],
           "entire boundary includes each lower face once");
Print("PASS: fred2 core tests\n");
