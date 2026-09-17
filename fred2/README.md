# fred2 modules

From the repository root, use the existing entry point:

```gap
LoadPackage("hap");
Read("fred2.g");
```

`fred2.g` loads these modules in dependency order. Function names, signatures,
and implementations are preserved from the original file.

| File | Functions |
| --- | --- |
| `star.g` | `kStar` |
| `forget_action.g` | `ForgetG` |
| `balance.g` | `PreIsGBalanced`, `IsGBalanced`, `Stabstar` |
| `lattice.g` | `EntireBoundary`, `Cobound`, `ShouldRemove`, `JoinCells` |
| `simplify.g` | `SimplifyLattice`, `HomeomorphicToDisk` |
| `replace.g` | `Replace2` |
| `reduce.g` | `ReduceGComplex` |

The similarly named files at the repository root are older implementations.
`setup.g` still loads those files; use `fred2.g` to select this implementation.

Use `ReduceGComplex(L)` to simplify a complex; its top dimension is determined
automatically.

## Sharing a single file

Copy [`fred2-standalone.g`](../fred2-standalone.g) to share all the functions
in one file, without tests or dependencies on the `fred2/` directory.
The recipient needs GAP with HAP installed:

```gap
LoadPackage("hap");
Read("/path/to/fred2-standalone.g");
reduced := ReduceGComplex(L);
```

Here `L` is the input HAP G-complex. The standalone file contains the same
implementations as the modules above, in dependency order. If you modify the
modules, update the standalone copy as well.

## Tests

See the [test guide](../tests/README.md) for the fixtures, individual checks,
homology comparison, and coverage limits.

Inside GAP, with the repository root as the working directory, run:

```gap
Read("tests/core.g");
Read("tests/integration.g"); # Requires HAP
```

Each suite prints a PASS message when successful and leaves GAP running.
Failed checks raise a GAP error; type `quit;` to leave the break loop and
return to your session.

Alternatively, start GAP with the tests from a terminal by supplying the path
to `tests/run.sh`:

```sh
./tests/run.sh             # Standalone core tests; GAP only
./tests/run.sh integration # Replacement and reduction; requires HAP
./tests/run.sh all         # Both suites
```

The terminal runner also leaves GAP open after the tests.
The core suite uses small explicit complexes to check star enumeration,
nontrivial stabilizer cosets, balanced and unbalanced stars, relative incidence
lattices, recursive boundaries, merging, disk detection, and orientation signs.
Without HAP, loading the entry point emits expected unbound-global warnings
for HAP's replacement type and coset helper; the core routines do not need them.

The integration suite checks replacement dimensions, boundary renumbering,
negative orientations, stabilizers, actions, and reduction idempotence on a
subdivided interval. It also runs `tests/attaching.g` to check boundaries above
replaced stars, attaching multiplicities, reflections, and boundary-of-boundary.
It fails explicitly when HAP cannot load. These are small regression fixtures,
not a verification of the algorithms for all G-complexes.

For the slower SL3Zs homology regression, run `Read("tests/sl3zs.g");` inside
GAP. It compares homology before and after unrestricted reduction in degrees
0 through 4 and is separate from the terminal runner's `all` option.

For the Bianchi subdivision regression, run `Read("tests/bianchi.g");` inside
GAP. It reduces before computing homology so that changes to HAP's element
tables cannot hide an incomplete star stabilizer.
