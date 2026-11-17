# Execution recipes

These are the 20 GitHub Actions recipes used for the catalogue refits, tracklet simulations, image checks and evidence integration. The run links in each analysis identify its exact commit, executed inputs and logs. The recipes retain the original checkout paths and working directories.

To reproduce one job locally, execute its calculation commands from the indicated working directory with the package and ephemeris versions specified there. Archived observation subsets and checksums are in the relevant result folders. For an Actions replay in a fork, copy the selected recipe to `.github/workflows/` at the repository root and choose the corresponding branch, input data and run parameters. Some follow-ups read a preceding run's artifact by ID. Those transient artifacts can expire. The result folders contain the archived analysis, input selections and input-hash receipts for use with a matching source snapshot.

[Scope and coverage](../../README.md) · [Joined-catalogue verification run](https://github.com/dovvnloading/Open-Space/actions/runs/44039123609)
