# Hydroxyl hydrogen radii modifications

## Background

AMBER-family force field hydroxyl hydrogens are known to have an issue when brought in close proximity where they can "fuse" in some situations.
Hydroxyl hydrogens have zero LJ parameters, and are normally thought to be "protected" by the LJ parameters of the oxygen they are attached to.
However, in our experience, this protection is not adequate; in some cases, when neighboring polar groups (such as an adjacent polar oxygen) come into close proximity, the positive charge on the hydrogen can get sufficiently close to the approaching oxygen that it overcomes the repulsion between the two oxygens (oen of which is attached to and protecting the hydrogen) and the hydrogen then falls into the infinite electrostatic potential well of the approaching oxygen and the simulation crashes.
In our experience this is particular pronounced in the case of multiple polar groups in confined environments (buried in binding sites, or in certain protonation states of Bruce Gibb's "octa acid" supramolecular host, etc.), but also occurs in simpler cases such as association of neutral benzoic acid dimers in water.
We have in the past alleviated this in some cases by artificially creating nonzero LJ parameters for polar hydrogens in specific systems to fix specific crashes, but no change has been made to AMBER-family force fields as a whole because in some sense this would require a refit.
However, parm@Frosst went so far as to introduce a new `HX` atom type which is identical to `HO` except importing nonzero LJ parameters from polar hydrogens attached to nitrogens in order to fix crashes; this type is used in certain cases where multiple hydroxyls are in close proximity.

**Here, we take advantage of the extensive validation we are performing for smirnoff99Frosst to add nonzero LJ parameters for polar hydrogen**.
Our goal is to do so with minimal perturbation to the force field overall, so these are kept quite small, with the thought that we will refit them from scratch at a later date.

## Explanation of the problem

Conceptually, the origin of the problem is simple: An oxygen (A) with an attached hydrogen (H) is approaching another oxygen (B), and the hydrogen is attracted to the target oxygen (B).
If the hydrogen were to overlap perfectly with (B) resulting in an infinite (favorable) electrostatic energy, this would result in a steric clash between A and B.
However, the energy of this steric clash is large but finite, since the separation of A and B is small enough to trigger significant repulsion, but not so small that the potential diverges.

This line of reasoning implies that as a hydrogen (H) approaches its hydrogen bonding partner (B), the energy should fall to a minimum at a preferred "hydrogen bonding" distance, then begin to rise again as (A) and (B) begin to "contact".
However, if (A) and (B) *continue* their approach, at some point the energetics of the favorable (H)-(B) interaction -- which become infinite as the separation goes to zero -- will dominate over the finite repulsive (A)-(B) interaction and (H) will suffer an electrostatic catastrophe and collapse onto (B).

We reproduce this effect for acetic acid dimers in the gas phase in `Dimer energetics.ipynb` and map out the potential for this association event.
In this case, the barrier for this association is very large compared to $k_BT$, suggesting that when this arises in simulations, it may be triggered either by an extremely unusual environment, or by integration error which results in the proton crossing the energetic barrier by "accident" and being unable to recover.

## Approach

Here, we perform three main steps to generate draft polar hydrogen LJ parameters:
1. Generate a set of hydrogen parameters with a `rmin_half` value we choose to be very small and an `epsilon` value selected to ensure that interactions with other atoms outside the radius of the adjoining oxygen are below a specified cutoff threshold; this is done in `hydrogen_radii.ipynb`
2. Run density calculations on TIP3P and modified TIP3P (modified to use the proposed hydrogen parameters, though we do not intend to apply these to water) as a quick and relatively extreme test of whether the proposed new parameters result in a significant perturbation to the density.
3. Check that we can reproduce the electrostatic catastrophe described above for the case of (neutral) acetic acid dimerization in the gas phase, and check that the proposed hydrogen parameters remove this catastrophe. This is done in `Dimer energetics.ipynb`.

After this, we will use the updated hydrogen parameters in validation tests on densities of pure solutions and hydration free energies.

## Manifest
- [`hydrogen_radii.ipynb`][hydrogen_radii.ipynb]: Jupyter notebook estimating what radii should be used for hydrogens to minimize interactions outside oxygen sigma or rmin_half
- [`density.py`](density.py): Script running density calculations for test systems with proposed new hydrogen radii to see how these impact density. Initially testing on TIP3P and modified TIP3P water, where we would NOT be applying these radii but which provides a useful (extreme) test system to see what happens if the parameters are used a great deal.
- [`results`](results): (if present) Output files from `density.py`; generated by running `density.py`
- [`Dimer energetics.ipynb`](Dimer energetics.ipynb): File looking at hydroxyl radii applied to dimers to see what happens to dimer energetics. Using smirnoff99Frosst.ffxml and a modified smirnoff99Frosst.ffxml which adds hydroxyl LJ parameters (which it creates), `hydrogen_radii.ffxml`
- [`dimer_data`](dimer_data): Data relating to `Dimer energetics.ipynb`
- [`smirnoff99Frosst_1.0.6.ffxml`](smirniff99Frosst_1.0.6.ffxml): Version 1.0.6 of SMIRNOFF99Frosst forcefield, which we will add the hydroxyl radii to.
