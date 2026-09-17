if LoadPackage("hap") = fail then
    Error("HAP and its dependencies must be installed for integration tests.");
fi;
Read("fred2.g");
Read("tests/helpers.g");
L := Fred2Interval(false);;
R := Replace2(L,1,0,1);;
Fred2Check(List([0..2],R!.dimension), [2,1,0], "replacement dimensions");
Fred2Check(Set(R!.boundary(1,1)), Set([[1,1],[-2,1]]),
           "replacement boundary and renumbering");
Fred2Check(Set(R!.boundary(1,-1)), Set([[-1,1],[2,1]]),
           "negative cell orientation");
Fred2Check(R!.stabilizer(1,1), Group(()), "replacement stabilizer");
Fred2Check(R!.action(1,1,1), 1, "replacement orientation action");
R := ReduceGComplex(Fred2Interval(false));;
Fred2Check(List([0..2],R!.dimension), [2,1,0], "reduction dimensions");
S := ReduceGComplex(R);;
Fred2Check(List([0..2],S!.dimension), [2,1,0], "reduction is idempotent");
Fred2Check(S!.boundary(1,1), R!.boundary(1,1), "stable reduced boundary");
Read("tests/attaching.g");
Print("PASS: fred2 HAP integration tests\n");
