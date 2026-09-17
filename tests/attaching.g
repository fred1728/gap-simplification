# Regression tests for stars strictly below the top dimension.
if LoadPackage("hap") = fail then
    Error("HAP and its dependencies must be installed for attaching-boundary tests.");
fi;
Read("fred2.g");
Read("tests/helpers.g");

Fred2TestAttachingBoundaries := function()
local Polygon, BoundarySquared, L, R, reverse, multiplicity, boundary,
      triangle, reflection, g, gi, gri;

# A filled triangle: edges 1 and 2 meet at vertex 1 and are merged.
# Repeated boundary circuits exercise actual attaching multiplicity.
Polygon := function(reverse, multiplicity)
    local edges, circuit, trivial;
    trivial := Group(());
    edges := [[[-1,1],[2,1]], [[1,1],[-3,1]], [[2,1],[-3,1]]];
    circuit := [[1,1],[2,1],[-3,1]];
    if reverse then
        edges[2] := List(edges[2],t->[-t[1],t[2]]);
        circuit[2][1] := -2;
    fi;
    circuit := Concatenation(List([1..multiplicity],i->StructuralCopy(circuit)));
    return rec(
        elts := [()], group := trivial,
        dimension := function(d)
            if d=0 or d=1 then return 3; fi;
            if d=2 then return 1; fi;
            return 0;
        end,
        boundary := function(d,i)
            local b;
            if d=0 then return []; fi;
            if d=1 then b:=edges[AbsoluteValue(i)];
            elif d=2 then b:=circuit;
            else return []; fi;
            return List(b,t->[SignInt(i)*t[1],t[2]]);
        end,
        stabilizer := function(d,i) return trivial; end,
        action := function(d,i,j) return 1; end
    );
end;

# Expand d(d(cell)) over group translates and collect using the lower
# cells' stabilizers and orientation actions. No tensoring with integers:
# distinct translates must remain distinct unless their stabilizer agrees.
BoundarySquared := function(K,d,i)
    local terms, first, second, face, g, p, h, eps, entry;
    terms:=[];
    for first in K!.boundary(d,i) do
        for second in K!.boundary(d-1,AbsoluteValue(first[1])) do
            face:=AbsoluteValue(second[1]);
            g:=K!.elts[first[2]]*K!.elts[second[2]];
            p:=First([1..Length(K!.elts)],
                     p->Inverse(g)*K!.elts[p] in K!.stabilizer(d-2,face));
            if p=fail then
                Add(K!.elts,g);
                p:=Length(K!.elts);
            fi;
            h:=Inverse(g)*K!.elts[p];
            if not h in K!.elts then Add(K!.elts,h); fi;
            eps:=SignInt(first[1])*SignInt(second[1])*
                 K!.action(d-2,face,Position(K!.elts,h));
            entry:=First(terms,t->t[1]=face and t[2]=p);
            if entry=fail then Add(terms,[face,p,eps]);
            else entry[3]:=entry[3]+eps; fi;
        od;
    od;
    return Filtered(terms,t->t[3]<>0);
end;

for reverse in [false,true] do
    for multiplicity in [1,2] do
        L:=Polygon(reverse,multiplicity);
        Fred2Check(BoundarySquared(L,2,1),[],"input triangle boundary squared");
        R:=Replace2(L,1,0,1);
        Fred2Check(List([0..3],R!.dimension),[2,2,1,0],
                   "replacement below top dimension keeps the face");
        boundary:=Concatenation(List([1..multiplicity],i->[[2,1],[-1,1]]));
        Fred2Check(SortedList(R!.boundary(2,1)),SortedList(boundary),
                   "one merged edge per boundary circuit, with multiplicity");
        Fred2Check(SortedList(R!.boundary(2,-1)),
                   SortedList(List(boundary,t->[-t[1],t[2]])),
                   "negative orientation of attaching cell");
        Fred2Check(BoundarySquared(R,2,1),[],"replacement boundary squared");
    od;
od;

# Two copies of the face bound a 3-cell. Replacing a 1-star must update
# both faces while leaving the 3-cell's boundary unchanged.
L:=Polygon(false,1);
triangle:=L.boundary;
L.dimension:=function(d)
    if d=0 or d=1 then return 3; fi;
    if d=2 then return 2; fi;
    if d=3 then return 1; fi;
    return 0;
end;
L.boundary:=function(d,i)
    if d=3 then return [[SignInt(i),1],[-2*SignInt(i),1]]; fi;
    return triangle(d,i);
end;
R:=Replace2(L,1,0,1);
Fred2Check(R!.boundary(3,1),L.boundary(3,1),"higher attaching boundary unchanged");
Fred2Check(BoundarySquared(R,3,1),[],"boundary squared above the attaching cells");

# Reflection exchanges the two constituent edges of the star. Their
# oriented sum is e - r*e, and reflection reverses the new edge.
reflection:=(1,2);
L:=rec(
    elts:=[(),reflection], group:=Group(reflection),
    dimension:=function(d)
        if d=0 or d=1 then return 2; fi;
        if d=2 then return 1; fi;
        return 0;
    end,
    stabilizer:=function(d,i)
        if (d=0 and AbsoluteValue(i)=1) or
           (d=1 and AbsoluteValue(i)=2) or d=2 then
            return Group(reflection);
        fi;
        return Group(());
    end,
    action:=function(d,i,j)
        if (d=2 or (d=1 and AbsoluteValue(i)=2)) and j=2 then return -1; fi;
        return 1;
    end,
    boundary:=function(d,i)
        local b;
        if d=0 then return []; fi;
        if d=1 and AbsoluteValue(i)=1 then b:=[[-1,1],[2,1]];
        elif d=1 then b:=[[2,1],[-2,2]];
        elif d=2 then b:=[[1,1],[-1,2],[-2,1]];
        else return []; fi;
        return List(b,t->[SignInt(i)*t[1],t[2]]);
    end
);
Fred2Check(BoundarySquared(L,2,1),[],"reflection input boundary squared");
R:=Replace2(L,1,0,1);
Fred2Check(Set(R!.boundary(2,1)),Set([[2,1],[-1,1]]),
           "translated constituents give one oriented star attachment");
Fred2Check(R!.action(1,2,2),-1,"reflection reverses the merged edge");
Fred2Check(BoundarySquared(R,2,1),[],"equivariant replacement boundary squared");

# Induce the reflected example to S3, so HAP can request orientation signs
# at elements outside the new edge's stabilizer as well.
L.group:=SymmetricGroup(3);
Append(L.elts,Filtered(Elements(L.group),g->not g in L.elts));
R:=Replace2(L,1,0,1);
for g in Elements(L.group) do
    gi:=Position(R!.elts,g);
    gri:=Position(R!.elts,g*reflection);
    Fred2Check(R!.action(1,2,gi) in [-1,1],true,
               "orientation outside the stabilizer is a sign, never zero");
    Fred2Check(R!.action(1,2,gri),-R!.action(1,2,gi),
               "orientation extension respects the stabilizer character");
od;

# Reducing a triangle again produces a loop with zero cellular boundary.
# Its identity orientation action must still be defined.
R:=ReduceGComplex(Polygon(false,1));
Fred2Check(List([0..3],R!.dimension),[1,1,1,0],"triangle reduction without top restriction");
Fred2Check(R!.boundary(1,1),[],"reduced loop has zero cellular boundary");
Fred2Check(R!.action(1,1,1),1,"orientation action with zero cellular boundary");
Fred2Check(BoundarySquared(R,2,1),[],"full triangle reduction boundary squared");
Print("PASS: attaching-boundary tests\n");
end;

Fred2TestAttachingBoundaries();
