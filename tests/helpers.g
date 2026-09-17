Fred2Check := function(actual, expected, label)
    if actual <> expected then
        Print("FAIL: ", label, "\nExpected: ", expected,
              "\nActual: ", actual, "\n");
        Error("Test failed: ", label);
    fi;
end;

# Two oriented edges meeting at vertex 1, with endpoints 2 and 3.
# reverse=true reverses the second edge to exercise orientation tracking.
Fred2Interval := function(reverse)
    local boundaries, trivial;
    trivial := Group(());
    if reverse then
        boundaries := [[[-1,1],[2,1]], [[-1,1],[3,1]]];
    else
        boundaries := [[[-1,1],[2,1]], [[1,1],[-3,1]]];
    fi;
    return rec(
        elts := [()], group := trivial,
        dimension := function(d)
            if d = 0 then return 3; fi;
            if d = 1 then return 2; fi;
            return 0;
        end,
        boundary := function(d,i)
            if d = 0 then return []; fi;
            if d = 1 then return boundaries[i]; fi;
            return [];
        end,
        stabilizer := function(d,i) return trivial; end,
        action := function(d,i,j) return 1; end
    );
end;
