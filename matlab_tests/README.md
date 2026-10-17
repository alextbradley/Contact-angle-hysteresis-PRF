# Figure 11 generator regression

This test checks only the regeneration **script**, not the published physics.

Run from the repository root:

```sh
octave --quiet --no-gui --eval "run('matlab_tests/test_fig11_generator.m')"
```

The regression copies the generator into a temporary directory and supplies a
deterministic `run_droplet` fixture that returns the escape event (5) for all
5 x 11 parameter combinations. It checks that the generator passes the
required eighth argument `logout`, updates **only** the appropriate
`xu0_escape(i,j)` cell, retains the 5 x 11 array, and saves to a new filename
so the deposited reference dataset remains untouched.

A **real scientific regeneration** must use the original
`Auxillary_functions/run_droplet.m`, not the fixture. The original script also
requires the caller to supply `V`, `theta_r`, `N`, `tend` and `logout`.
The script does not document all the historical parameter values, so this
test cannot establish agreement with the paper's reported threshold map.
