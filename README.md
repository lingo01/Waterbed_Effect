
# Waterbed_Effect

Data and code of the article "Stability Degradation Induced by Angle–Voltage Coupling in Power Systems: Bode-Type Fundamental Performance Limitation Analysis"

This source code repository accompanies the following paper:

> [1] X. Peng, Z. Li, C. Fu, P. Yang, Z. Wang, and F. Liu, ‘Stability Degradation Induced by Angle–Voltage Coupling in Power Systems: Bode-Type Fundamental Performance Limitation Analysis’, IEEE Transactions on Automation Science and Engineering, 2026.

The full paper and the source code can be found at: https://github.com/lingo01/Waterbed_Effect

For any questions or uses of the source codes, please feel free to contact the first author, Xiaoyu Peng ([pengxy19@tsinghua.org.cn](mailto:pengxy19@tsinghua.org.cn)), and the corresponding author, Feng Liu ([lfeng@mail.tsinghua.edu.cn](mailto:lfeng@mail.tsinghua.edu.cn)).

**CITATION**: If you use this code in your work, whether directly or indirectly, please cite the above paper.

**LICENSE**: This work is licensed under the MIT License. See the [LICENSE](https://github.com/lingo01/Compositional_Grid_Code/blob/main/LICENSE) file in the repository for details.


# Introduction to the source code

All code is written in `MATLAB` (tested with `MATLAB 2018a` or later).

The workflow is:

1. Load pre-generated network/device data from `data&figure/` (or regenerate network sensitivity matrices for some experiments).
2. Build the closed-loop MIMO model of angle/voltage coupling.
3. Perform frequency-domain sensitivity integral analysis and time-domain simulation verifications.


## Code-Figure Correspondence

The following scripts reproduce the figures in the paper:

| Figure in Paper | MATLAB Script | Description |
|---|---|---|
| Fig.5  | `main_SMIB_GFM.m` | SMIB example (GFM) analysis and plots |
| Fig.6  | `main_SMIB_tuning.m` | SMIB tuning study and plots |
| Fig.7  | `main_fixed_freq.m` | Fixed-frequency case study and plots |
| Fig.8  | `main_multi_inv.m` | Multi-inverter frequency-domain analysis (eigenvalue tracking / directed areas) |
| Fig.9  | `main_multi_inv_simu.m` | Multi-inverter time-domain simulation comparing coupled vs. decoupled models |
| Fig.10 | `main_multi_inv_statistics.m` | Monte-Carlo statistical study of sensitivity integrals under varying operation conditions |


## How to Generate Each Figure

1. Set the MATLAB current folder to the repository root.
3. Run the corresponding script.

Example (reproduce Fig.8):

```matlab
main_multi_inv
```

The scripts generate figures interactively (MATLAB figure windows). Some cases may also load or save intermediate `.mat` data under `data&figure/`.


### Supporting Functions

This repository includes a set of helper functions that are called by the main scripts, e.g.:

- `func_Hnet_generator.m` (network Jacobian / coupling model generation)
- `func_track_eigenvalues.m` (eigenvalue tracking)
- `func_classify_singular.m`, `func_classify_track.m` (mode classification)
- `func_area_bode.m` (directed integral / area computation)
- `func_eig_modi.m` (eigenvalue processing)


## Data Files

Pre-generated data are stored under `data&figure/`, including (depending on case):

- `Hnet_info_*.mat`: network coupling Jacobians and metadata
- `dev_info_*.mat`: device parameters


## Requirements

- `MATLAB 2018a` or newer
- MATLAB Control System Toolbox (for `tf`, `ss`, `lsim`, `evalfr`, etc.)
- MATLAB Symbolic Math Toolbox (for `syms`, `jacobian`, and `vpa` usage)
- MATPOWER (required by `func_Hnet_generator.m` via `loadcase`, `runpf`, `makeYbus`). See [https://matpower.org](https://matpower.org/) for install instruction.


## Notes

- Reproducibility: many scripts use fixed default parameters; Fig.10 uses `rng(2026)` for repeatability.
- Runtime: Fig.8 and Fig.10 can be time-consuming due to wideband frequency sweeps and/or Monte-Carlo sampling.
- Case selection: scripts typically use IEEE test cases (e.g., `case14`, `case39`, `case118`). Ensure the corresponding data files exist under `data&figure/`.


