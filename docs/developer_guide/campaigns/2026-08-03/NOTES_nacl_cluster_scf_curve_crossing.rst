Bug: the "cluster" sub-calculation lands on the wrong SCF solution near the
Na\ :sup:`+`\ /Cl\ :sup:`-` ionic/covalent curve crossing
============================================================================

:Status: Root cause identified, not yet fixed
:Found in: Job 3867 on the ``dev`` dashboard -- Na\ :sup:`+`\ ...Cl\ :sup:`-`
           40-point log-spaced R-scan (1.30-20.00 A), REVDSD-PBEP86-D4/2021
           def2-TZVPPD, ``--charges "1 -1"``, through the M3 ORCA
           N-fragment substep this campaign built. Structures:
           ``~/Sites/training-data/2026-08-03_na-cl-scan-run1/`` (``Na-Cl.sdf``).
:Affects: The **cluster** sub-calculation only (the one full-complex, both
          atoms real, net-neutral SCF -- one of the ``2N + 1`` jobs per point).
          Not the fragment-alone or fragment-in-cluster-basis jobs.
:Reported by: user, from the labeled ``Dimers.extxyz`` -- noticed nonzero
          x/y forces that should vanish by the molecule's linear symmetry,
          plus noisy Chargemol/DDEC6 charges with no clean trend vs. R.

Summary
-------

This is not noise and not a convergence-tolerance problem. A single-
determinant (RHF/KS) description of neutral NaCl has two nearly-degenerate
electronic solutions once R is large enough -- an **ionic** one
(Na\ :sup:`+`\ Cl\ :sup:`-`, the desired diabatic surface for this training
set) and a **covalent** one (Na\ :sup:`.`\ Cl\ :sup:`.`, two neutral atoms).
They cross where the Coulomb stabilization of the ionic state equals the
cost of forming it:

.. math::

   R_c \approx \frac{e^2}{\mathrm{IP}(\mathrm{Na}) - \mathrm{EA}(\mathrm{Cl})}
             = \frac{14.399\ \mathrm{eV \cdot A}}{5.139 - 3.617\ \mathrm{eV}}
             \approx 9.5\ \text{A}

which lands right in the middle of the observed problem window. Near a
crossing like this, the SCF's default (neutral-atom) initial guess can
converge -- self-consistently, no convergence warnings -- to either branch
depending on details of the specific geometry; nothing chains guesses
between R-scan points (checked: no ``%moinp``/``MORead`` in any
``orca.inp``, each point is a fresh default guess), so which branch wins is
essentially decided independently, geometry by geometry. Tighter SCF
thresholds do not fix this: the SCF **is** converging, just to the wrong one
of two locally-stable solutions.

Evidence
--------

Pulled the per-subcalc ORCA logs for job 3867 (``root/5/<config>/2/1/
{cluster,1-alone,2-alone,1-in-cluster,2-in-cluster}/orca.out``) for the
flagged points and several clean neighbors for comparison.

* **Isolated to the cluster job.** The "cluster" SCF energy is non-monotonic
  in R = 6.5-12.2 A (re-deepens instead of relaxing toward the
  E(Na\ :sup:`+`) + E(Cl\ :sup:`-`) asymptote), then snaps back to sane
  behavior at R = 13.13 A. Over the identical R range, the ``1-alone``,
  ``2-alone``, ``1-in-cluster``, and ``2-in-cluster`` energies are all smooth
  to ~1e-6 Hartree -- unsurprising, since each of those is a single,
  unambiguous charge state (+1, -1, or a real fragment plus inert ghost
  functions with no competing electron count), so there is no second
  solution for them to fall into.
* **A virtual orbital is dropping toward becoming an electron-transfer
  acceptor state.** In the cluster job, the LUMO(+1) pair's orbital energy
  drops steadily from +0.0147 Ha at R = 6.52 A through negative values by
  R = 13.13 A -- the textbook signature of approaching this kind of crossing.
* **Genuine SCF symmetry breaking.** The molecule is linear
  (D:math:`_{\infty h}`), so degenerate pi-orbital pairs must be exactly
  degenerate. At R = 12.25 A (config ``1,33``) the normally-exact pair splits
  slightly (-5.4467 vs. -5.4436 eV) -- the same phenomenon that shows up
  downstream as nonzero x/y forces.
* **Total-force check pinpoints the exact broken points.** Sum of forces
  over the two atoms must be exactly zero (translational invariance, no
  external field); this is a stronger, purely computational diagnostic than
  "x/y should vanish by symmetry" and needs no chemistry knowledge to apply.
  Background level (healthy points) is ~1e-4 to 1e-3; three points are
  50-500x above that:

  .. list-table::
     :header-rows: 1

     * - config
       - R (A)
       - \|sum F\|
       - in the crossing window?
     * - ``1,17``
       - 3.99
       - 0.049
       - no -- likely an unrelated, isolated glitch
     * - ``1,25``
       - 6.99
       - 0.066
       - yes (edge)
     * - ``1,27``
       - 8.04
       - 0.054
       - yes

  The energy is smooth through all three points -- only the gradient (and,
  for 1,25/1,27, the Chargemol charge, and the interaction energy itself in
  ``Results.csv``) is affected. The broader charge/energy drift spans
  roughly configs ``1,26``-``1,33`` (R = 7.5-12.25 A) even though only two of
  those points trip the force-sum threshold -- i.e. the crossing can perturb
  a point's density/charge partition measurably before it perturbs the
  gradient enough to trip a simple threshold, so a threshold-based gate
  alone will under-count affected points. ``1,17`` at R = 3.99 A is well
  inside the normal ionic-bonding region, nowhere near :math:`R_c` -- treat
  it as a separate, one-off SCF quirk rather than the same phenomenon until
  shown otherwise.

Why this matters beyond Na/Cl
------------------------------

:math:`R_c` depends only on IP(cation) - EA(anion), so **every** ionic pair
this campaign scans will have some crossing distance, just at a different R.
Worth precomputing :math:`R_c` (or at least flagging "this pair's IP-EA gap
is small enough that :math:`R_c` falls inside the planned scan range") before
each new ion-pair production run, rather than rediscovering it after the
fact per-pair. Ion-water and ion-water-water systems (the angular/trimer test
sets built alongside this scan, see
``~/Sites/mlff-training/2026-08-07_ion-water-angular-and-trimers/``) don't
have this exact failure mode since the "cluster" fragment count/composition
differs, but the same *class* of problem -- a fragment sub-calculation with
more than one accessible electronic state near-degenerate at some geometry --
is worth keeping in mind generally for charged multi-fragment BSSE.

Fix directions (not yet implemented)
-------------------------------------

Because this is "converges to the wrong stable solution," not "fails to
converge," a tighter ``TightSCF -> VeryTightSCF`` will not help. In order of
expected leverage:

#. **Seed the cluster job's guess from the fragment orbitals already being
   computed.** Every BSSE point already runs ``1-alone`` (Na\ :sup:`+`) and
   ``2-alone`` (Cl\ :sup:`-`) as two of the ``2N + 1`` required jobs -- their
   orbitals are unambiguous (single charge state each) and free. Combine
   them into the cluster job's initial guess (ORCA supports building a
   guess from fragment MOs) instead of letting it fall back to the default
   neutral-atom guess. This directly biases every point toward the ionic
   diabatic surface, which is *already* the deliberate target everywhere
   else in this pipeline (the whole point of ``--charges`` fixing
   per-fragment charge), including past :math:`R_c` where the true adiabatic
   ground state would actually be covalent -- that's correct for training-
   set purposes, not a bug to route around.
#. **Add SCF stability analysis** (``!StabPerform`` in current ORCA syntax)
   to the cluster job as a backstop, so a solution that is actually a saddle
   point gets kicked back toward a genuine minimum. Cheaper to add than (1)
   but doesn't by itself guarantee landing on the *ionic* minimum rather
   than the covalent one past :math:`R_c` -- (1) and (2) are complementary,
   not alternatives.
#. **A generic post-hoc QC gate**, independent of the above: check
   :math:`\sum \vec F \approx 0` across the full cluster automatically for
   every point of every charged/multi-fragment production run (cheap, no
   chemistry-specific knowledge needed) and flag/retry anything outside
   tolerance. This is what surfaced the problem here and generalizes to any
   future ionic system, but per the point above about 1,26/1,28-1,33 not all
   tripping the threshold, it should be a first-pass filter, not the only
   check -- pair it with a smoothness check on the CP-corrected energy and/or
   the atomic charges vs. neighboring R points.

Raw data
--------

Dashboard: ``dev`` (``http://localhost:55066``), job 3867, project
``water``. Key files: ``Dimers.extxyz`` (labeled structures/forces/charges),
``Results.csv`` (per-point interaction energies), and per-subcalc ORCA I/O
under ``root/5/<config>/2/1/{cluster,1-alone,2-alone,1-in-cluster,
2-in-cluster}/orca.{inp,out}`` for each ``<config>`` (e.g. ``1,25``).
