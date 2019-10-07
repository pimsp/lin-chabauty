This code is mainly from Nicolas Mascot, and some code was added by Pim Spelier. See the notes for that at the end.

################################################################################

This package computes Galois representations appearing in the Jacobian of a given curve C, given the char.poly. of the Frobenius at p for p a prime of good reduction of C. See my preprint "Hensel-lifting torsion points on Jacobians and Galois representations" (arXiv 1808.03939) for more information.

It requires the branch "aurel-matmod-ZpXQ" of the Git version of PARI/GP, and the gcc compiler.

At the moment, only the cases of hyperelliptic curves and of smooth plane curves are implemented.

Two examples are provided:
* ExHyper.gp: the representation is afforded by a piece of the 7-torsion of a hyperelliptic curve, namely the modular curve X1(13) (genus 2).
* ExSmooth.gp: the representation is that afforded by the whole 2-torsion of a smooth plane curve, namely the Klein quartic (genus 3).

These examples are intended to demonstrate the functionalities of this package.
To run them, first compile the package by typing "make all", and then start GP and type 
read("ExX.gp")
where X is either "Hyper" or "Smooth".

This package is provided in the hope it will be useful, but comes without any guarantee whatsoever.

Nicolas Mascot, September 1, 2018; updated November 7, 2018.

################################################################################

I added some comments to several .h-files. I also wrote some code for hyperelliptic curves in Hyper2RR.gp, with main accomplishments a function that computes the Abel Jacobimap C -> J, P -> P - infty_+ or P - infty_- for Z/p^eZ points of C, and a function that can compute, with some data given, an upper bound for the number bound affine rational points on C. This last function uses linear Chabauty with Edixhovens method, calculating with divisors in J(Z/p^2 Z) that reduce to 0 modulo p.

For this, see the preprint coming out soon from Edixhoven and Lido.

Pim Spelier
Last updated 21 August, 2019
