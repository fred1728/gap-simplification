CheckResolutionData := function(L, maxd)
local d, i, b, t, s, j, g, ncell, nelt;

    Print("Checking dimensions...\n");
    for d in [0..maxd] do
        if L!.dimension(d) = fail then
            Error("dimension(", d, ") = fail");
        fi;
        if not IsInt(L!.dimension(d)) then
            Error("dimension(", d, ") is not an integer");
        fi;
        if L!.dimension(d) < 0 then
            Error("dimension(", d, ") is negative");
        fi;
    od;

    Print("Checking stabilizers and boundaries...\n");
    for d in [0..maxd] do
        ncell := L!.dimension(d);
        for i in [1..ncell] do

            s := L!.stabilizer(d,i);
            if s = fail then
                Error("stabilizer(", d, ",", i, ") = fail");
            fi;

            if d > 0 then
                b := L!.boundary(d,i);
                if b = fail then
                    Error("boundary(", d, ",", i, ") = fail");
                fi;

                for t in b do
                    if Length(t) <> 2 then
                        Error("bad boundary term at d=", d, " i=", i, ": ", t);
                    fi;

                    if not IsInt(t[1]) then
                        Error("boundary cell index is not integer at d=", d, " i=", i, ": ", t);
                    fi;

                    if AbsoluteValue(t[1]) < 1 or AbsoluteValue(t[1]) > L!.dimension(d-1) then
                        Error("boundary cell index out of range at d=", d, " i=", i, ": ", t);
                    fi;

                    if not IsInt(t[2]) then
                        Error("boundary elt index is not integer at d=", d, " i=", i, ": ", t);
                    fi;

                    if t[2] < 1 or t[2] > Length(L!.elts) then
                        Error("boundary elt index out of range at d=", d, " i=", i, ": ", t);
                    fi;
                od;
            fi;
        od;
    od;

    Print("Checking actions on stabilizers...\n");
    for d in [0..maxd] do
        ncell := L!.dimension(d);
        for i in [1..ncell] do
            s := L!.stabilizer(d,i);
            for g in Elements(s) do
                j := Position(L!.elts, g);
                if j = fail then
                    Error("group element missing from elts at d=", d, " i=", i, ": ", g);
                fi;
                if L!.action(d,i,j) = fail then
                    Error("action(", d, ",", i, ",", j, ") = fail");
                fi;
            od;
        od;
    od;

    Print("All checks passed up to degree ", maxd, ".\n");
end;