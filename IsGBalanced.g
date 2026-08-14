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
return PreIsGBalanced(L,e,n,k)[2];
end;