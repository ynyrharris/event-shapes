
# Event Shapes

This README documents the details of the Event Shape calculations performed by this package.

Contact [ynyr.harris@physics.ox.ac.uk](mailto:ynyr.harris@physics.ox.ac.uk) with queries.


## Setup

Do this in an environment that has ROOT 6 and Python 3 installed.  On a machine connected to cvmfs, you can use an LCG view to this end
```bash
source /cvmfs/sft.cern.ch/lcg/views/setupView.sh LCG_98python3 x86_64-centos7-gcc8-opt
```

Acquire the repository and build the project.
```bash
git clone <repo url>
```

This repository is intended as a subproject.  To incorporate it into your main project, add a line like `add_subdirectory(event-shapes)` to your top-level CMakeLists.txt, supposing that your project is laid out as
```
- build/
- run/
- src/
  - CMakeLists.txt
  - event-shapes/
  - <other packages>/
```

### Building for Python
```bash
cd build
cmake [-DPYTHON_SETUP=1] ../src # -DPYTHON_SETUP=1 is the default
make -j
make install
```

### Building for C++
```bash
cd build
cmake -DCXX_SETUP=1 ../src
make -j
```


## Introduction

Event Shape measures are intended to give a global view of an event.
In what follows, each event is characterised by a set of four-momentum vectors $`p_i = (\mathbf{p}_i, E_i)`$, where $`i = 1, 2, ..., n`$ index the $`n`$ vectors of the event.
Transverse momentum with respect to a specified axis is denoted $`p_T`$.


## Thrust

Thrust is defined as
```math
T = \max_{\mathbf{t}} \frac{\sum_i |\mathbf{t} \cdot \mathbf{p}_i|}{\sum_i | \mathbf{p}_i|},
```
where the maximising unit vector $`\mathbf{t}`$ is called the *thrust axis*.
The thrust axis corresponds to the direction of greatest energy flow in an event.
Thrust can take values of $`1/2 \leq T \leq 1`$, where 2-jet, pencil-like, events will have larger values, and isotropic events lower values.

We may also define a thrust-like axis in the plane perpendicular to the thrust axis.
This is called the thrust major axis, and the value of thrust calculated with this axis as $`\mathbf{t}`$ is called thrust major, $`T_{\text{major}}`$.
The axis perpendicular to both the thrust and thrust major axes is called the thrust minor axis, and the value of thrust calculated with respect to it is called thrust minor, $`T_{\text{minor}}`$.
The difference between thrust major and thrust minor is called oblateness, $`O = T_{\text{major}} - T_{\text{minor}}`$.

Finally, the broadening of an event is given by
```math
B = \sum_i \frac{|\mathbf{t} \times \mathbf{p}_i|}{\sum_i | \mathbf{p}_i|}.
```



## Sphericity

The sphericity tensor is defined as
```math
S^{\alpha\beta} = \frac{\sum_i |\mathbf{p}_i|^{-1} p_i^{\alpha}p_i^{\beta}}{\sum_i |\mathbf{p}_i|}
```
(Donoghue's definition), where $`\alpha, \beta`$ run over the Cartesian components of the vectors $`p_i`$.
It is geared towards quantifying the shapes of 2 and 3-jet events, as the terminology will make clear.

***Derivation.***
The original notion of a sphericity tensor as first proposed by [Bjorken and Brodsky in 1970, Footnote 7][Bjorken1970],
```math
S^{\alpha\beta} = \frac{\sum_i \left(\frac{3}{2} p_i^{\alpha} p_i^{\beta} - \frac{1}{2} \delta^{\alpha\beta} p_i^2\right)}{\sum_i p_i^2},
```
corresponded essentially to the quadrupole moment of the angular distribution of the $`p_i`$.
The three eigenvalues, $`\lambda_1 \geq \lambda_2 \geq \lambda_3`$, of this momentum ellipsoid are sums of squares of transverse momenta with respect to the three eigenvector directions.
The smallest eigenvalue, $`\lambda_3`$, belongs to the eigenvector along which lies the greatest flow of momentum (the so-called jet axis, particularly in a 2-jet event).

However, as this tensor is not calculable in perturbative QCD, an alternative tensor was proposed by [Donoghue *et al.*][Donoghue1979] which is.
Donoghue's definition is often referred to as the linearised sphericity tensor, since the eigenvalues that result are subsequently linear in the $`p_i`$.
It has a similar geometrical interpretation as the original tensor, as well as the additional desirable properties that: (i) for perfect two-jet events, two of the eigenvalues are zero; (ii) for three-jet events, one of the eigenvalues is zero; and (iii) for events with a spherical momentum distribution, $`\lambda_1 = \lambda_2 = \lambda_3 = 1/3`$.
The normalisation provides that $`\lambda_1 + \lambda_2 + \lambda_3 = 1`$.

Sphericity event shapes are constructed from the eigenvalues.

### 3d

In three dimensions, sphericity itself is defined as
```math
S = \frac{3}{2}(\lambda_2 + \lambda_3).
```
This corresponds to the summed $`p_T^2`$ of an event in the plane that contains the least activity.
The factor of 3/2 arranges that $`0 \leq S \leq 1`$ (as, in a perfectly spherical event, $`\lambda_2 + \lambda_3 = 2/3`$).
A more isotropic event will have an $`S`$ closer to 1, while in a less isotropic event, $`S`$ will take on a lower value.

Aplanarity is defined as
```math
A = \frac{3}{2} \lambda_3
```
to measure the component of momentum outside of the event plane.
It is constrained by the factor of 3/2 to values of $`0 \leq A \leq 1/2`$, where a more planar event has lower values, and a more isotropic event has higher values.

Two more combinations of eigenvalues are drawn from the [Pythia 6.4 manual][Sjostrand2006] to bear on 3 and 4-jet structures.
These are, respectively,
```math
C = 3(\lambda_1 \lambda_2 + \lambda_1 \lambda_3 + \lambda_2\lambda_3),\\
D = 27\lambda_1 \lambda_2 \lambda_3.
```
Both of these quantities take values in $`0 \leq C, D \leq 1`$.

### 2d

Two dimensional sphericity values are defined to have the same interpretation as their three dimensional counterparts.
Since $`\lambda_3 = 0`$ in two dimensions, we have only $`\lambda_1`$ and $`\lambda_2`$ to work with.

The sphericity scalar - quantifying the amount of activity orthogonal to the direction of greatest energy flow - is proportional to the sum of the $`d - 1`$ smallest eigenvalues, thus $`S = 2\lambda_2`$.

As all input vectors here are in a plane, there is no aplanar activity, and thus no notion of aplanarity.

The analog to the three dimensional 3-jet structure can be defined as $`C = 4\lambda_1 \lambda_2`$.
The analog of three-dimensional 4-jet structure $`D`$ is the same as for the 3-jet structure.

The normalisation factor in each case arranges that $`0 < X < 1`$.


## Comparison between thrust and sphericity-based Event Shapes


## References

Here is a list of references to relevant literature on Event Shapes.
They are not necessarily cited in the body above.
- J. D. Bjorken and S. J. Brodsky, *Statistical Model for Electron-Positron Annihilation into Hadrons*, March 1970, [Phys. Rev. D **1**, 1416][Bjorken1970]
- E. Fahri, *Quantum Chromodynamics Test for Jets*, December 1977, [Phys. Rev. Lett. **39**, 1587][Fahri1977] (first definition of a thrust axis)
- A. De Rujula, J. Ellis, E. G. Floratos, M. K. Gaillard, *QCD predictions for hadronic final states in $`e^{+}e^{-}`$ annihilation*, June 1978, [Nucl. Phys. **B** 138, 387][Rujula1978] (arguing for thrust, 'spherocity', and 'acoplanarity' to be used as theoretically rigorous tests of QCD)
- D. P. Barber *et al.*,  *Discovery of Three-Jet Events and a Test of Quantum Chromodynamics at PETRA*, September 1979, [Phys. Rev. Lett. **43**, 830][PETRA1979-1] (thrust axis used to discover the gluon in 3-jet events)

    D. P. Barber *et al.*, *Study of Electron-Positron Collisions at Center-of-Mass Energies of 27.4 and 27.7 GeV at PETRA*, September 1979, [Phys. Rev. Lett. **43**, 901][PETRA1979-2] (thrust and 'spherocity' used to observe 2-jet structure in $`e^{+}e^{-}`$ collision events)

- J. F. Donoghue, F. E. Low and So-Young Pi, *Tensor analysis of hadronic jets in quantum chromodynamics*, December 1979, [Phys. Rev. D **20**, 2759][Donoghue1979]
- A. Banfi, G. P. Salam and G. Zanderighi, *Resummed event shapes at hadron-hadron colliders*, December 2004, [hep-ph/0407287][Banfi2004]
- T. Sjostrand, S. Mrenna and P. Skands, *PYTHIA 6.4 Physics and Manual*, May 2006, [hep-ph/0603175][Sjostrand2006]
- A. Banfi, G. P. Salam and G. Zanderighi, *Phenomenology of event shapes at hadron colliders*, January 2010, [hep-ph/1001.4082][Banfi2010]


[Bjorken1970]: <https://journals.aps.org/prd/abstract/10.1103/PhysRevD.1.1416/>
[Hanson1975]: <https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.35.1609/>
[Fahri1977]: <https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.39.1587/>
[Rujula1978]: <https://doi.org/10.1016/0550-3213(78)90388-7/>
[PETRA1979-1]: <https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.43.830/>
[PETRA1979-2]: <https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.43.901/>
[Donoghue1979]: <https://journals.aps.org/prd/abstract/10.1103/PhysRevD.20.2759/>
[Banfi2004]: <https://arxiv.org/abs/hep-ph/0407287/>
[Sjostrand2006]: <https://arxiv.org/abs/hep-ph/0603175/>
[Banfi2010]: <https://arxiv.org/abs/1001.4082/>