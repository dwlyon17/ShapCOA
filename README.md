# ShapCOA

Initial release, version 0.1.0, June 10, 2026.

This R package will compute Shapley Values (SVs) using Component Order of Addition (COA) designs, an approach that works even for huge problems.  It also handles other approaches to SV computation.  Any form of analysis (key drivers, or TURF, for example) can underly the SV analysis, but the package is particularly oriented toward common marketing research applications for SVs.

You can install the latest version of ShapCOA from GitHub using the R function:
```
remotes::install_github("DWLyon17/ShapCOA")
```
This initial version is fully functional for many use cases and includes:
* Generation of Component Order of Addition (COA) designs of an desired size,
* Exact computation of SVs using the combinations approach for small/medium problems (up to 25-30 items),
* Exact computation of SVs using the orderings approach for small problems (up to 12 items or so),
* Accurate approximation of SVs with COAs using the orderings approach, even for huge problems,
* Built-in provisions for key driver regressions,
* Built-in provisions for standard (0/1) TURF problems at any “depth”, 
* Built-in provisions for “TURF on MaxDiff” analyses with several variations,
* Provisions for any user-supplied method of valuing combinations of items or variables,
* Support for “size-limited Shapley Values” in all situations, and
* Support for respondent weights package-wide.

Other features planned for the future (all exist and work, but need thorough testing and some refactoring to be consistent with and interface with the rest of the package) include:
* Built-in provisions for loglinear key driver regressions (with highly optimized C++ code),
* Built-in provisions for many variants of 0/1 TURF, 
* Exact computations via hypergeometric probabilities for both standard TURF and many variants of TURF, and
TURF analysis, and variants of it, proper (searching for best combinations, not just for Shapley Values).

If you are interested in the package, feel free to email Dave Lyon at Dlyon@aurora2000.com to be placed on an e-mail list for notices of updates and changes.  And feel free to suggest priorities on the yet-unincluded features as well.
