# miivfind
This repository contains a Stata program that implements an algorithm for identifying model-implied instrumental variables for structural equation models. The original algorithm was developed for SAS and is presented in Bollen and Bauer (2004). The algorithm for Stata is discussed in Bauldry (2014).

## Revision task list
1. Rewrite syntax command to take specific vectors/matrices as input.
2. Return matrix of IVs
3. Add a selfcheck option (see Stas' code)
4. Reorganized displaying results as a subcommand.
5. Set up tempnames for matrices.

## References
Bauldry, Shawn. 2014. "miivfind: "miivfind: A Program for Identifying Model-Implied Instrumental Variables (MIIVs) for Structural Equation Models in Stata.." *Stata Journal* 14:60-75.

Bollen, Kenneth A. and Daniel J. Bauer. 2004. "Automating the Selection of Model-Implied Instrumental Variables." *Sociological Methods & Research* 32:425-452.