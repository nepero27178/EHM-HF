CURRENT STEPS
1. Finish setup for main SU-Singlet, FakeSU-Singlet run [DONE]
2. Launch main SU-Singlet run on two machines (safety pipeline on Stubborn)
3. Prepare filtered run module
4. Compute AF free energy
5. Finish setup for main AF run
6. Launch main AF run on two machines (safety pipeline on Stubborn)

ALGORITHM
- [x] Use DataFrames.jl to perform data retrieval, selection and plot
- [ ] Time efficient code: write on file every x times, with x reasonably large
- [x] Implement automatic safecopy of data
- [x] Employ symmetry: many calculations are just pure repetitions
- [x] Use normalized K
- [x] Re-use previous value in step for fast convergence
- [x] Restrict to MBZ calculation AND CHECK MULTIPLIERS!
- [x] Insert graphs in report
- [x] Investigate superconducting phase
- [x] Extract and plot Fermi surface
- [x] Fake rigid shift of hopping
- [x] Implement bands renormalization in superconducting phase
- [x] Move RenormalizeHopping => RenormalizeBands (d-wave part of g)
- [ ] Implement triplet simulations
- [x] Correct SU free energy
- [x] Compute bare free energy
- [x] Add bare free energy module
- [x] Compute AF free energy
- [x] New structure for AF simulations
- [x] Discard record-g mode and implement single mode
- [x] Test new general plot function
- [x] Write new filtered run module
- [ ] Add skip in 2D plot
- [ ] Add surface interface for RMP 3D plot
- [x] Clean up setups
- [x] Add coherent neighbors initializer
- [ ] Debug BZ optimization module

HM: AF
- [ ] heatmaps
- [x] scan
- [ ] interesting U scan: vary temperature, vary doping
- [x] write free energy module

HM: SC-SINGLET
- [ ] heatmaps s-wave (negative U too)
- [ ] scan s-wave (negative U too)
- [x] write free energy module

EHM: AF
- [ ] half-filling in depth analysis
- [x] heatmaps: keep doping null
- [x] scan: keep doping null
- [x] write free energy module
- [ ] plot the metallic bands

EHM: SC-SINGLET
- [x] sanity check: compare optimized weights run vs dummy run => use dummy run for now
- [ ] filtered run: re-compute the NaN points => write the filtered run module
- [ ] heatmaps+RMPs s+s*-wave (negative U too) => CHANGE COLUMN NAME, f => fMFT
- [ ] heatmaps+RMPs d-wave => CHANGE COLUMN NAME, f => fMFT
- [ ] heatmaps+RMPs s+s*+d-wave => CHANGE COLUMN NAME, f => fMFT
