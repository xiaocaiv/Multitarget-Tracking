# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## How to Run

All `.m` files live in `matlab/`; run from the repository root:

```bash
# Run the full simulation (MATLAB R2024a):
matlab -batch "addpath('matlab'); main"

# Compare single-model CV vs IMM over Monte Carlo runs:
matlab -batch "addpath('matlab'); compare_imm"
```

Alternatively `cd matlab` first and call `matlab -batch "main"`. In the MATLAB GUI add the folder once with `addpath('matlab')`.

### Repository layout

```
matlab/     all simulation code (.m), flat inside this folder
data/       ship-radiated-noise samples for the wav feature mode
docs/       design and implementation notes (Chinese)
archive/    pre-refactor .bak rollbacks, not referenced by any code
```

Data paths are resolved by `repo_path.m`, which derives the repository root from its own file location rather than from the CWD. Never hardcode a relative `.wav` path: `generate_ship_spectrum` falls back to simulated spectra when the file is missing, so a stale relative path silently yields fake feature data (it now raises a warning).

Requires Statistics and Machine Learning Toolbox (`chi2inv`, `poissrnd`) and Signal Processing Toolbox (`spectrogram`, `hilbert`, `butter`, `filtfilt` — only when `use_features=true` with wav import).

## Architecture

Dual-platform cooperative multi-target bearing-only tracking system. Two moving platforms observe bearing angles to maneuvering targets. Each track is filtered by an **IMM-UKF** (interacting multiple model: CV + coordinated turn + constant acceleration), JPDA handles measurement-to-track association, triangulation initializes new tracks, and the two platforms' estimates are fused per model by covariance intersection. Preset truth scenarios (turn / parallel / crossing / speed) live in `scenario_defs.m` and are selected with `params.scenario`.

### Data Flow (per time step)

```
generate_measurements  ──→  z1, z2  (struct arrays with .bearing, .spectrum)
    │
    ├─  imm_predict (all tracks):  model mixing + per-model UKF time update
    │
    ├─  Platform 1:  imm_measurement_predict → JPDA → imm_update → combine
    ├─  Platform 2:  same on a pre-prediction copy (independent association)
    │
    ├─  per-model CI fusion + model-probability fusion (geometric mean)
    ├─  manage_tracks:  delete stale  +  triangulate new + imm_init
    └─  ECM:  move platforms to maximize expected information gain
```

### Module Dependency Order (implementation / modification order)

| Layer | Files | Role |
|-------|-------|------|
| **Params** | `init_parameters.m` | All simulation, IMM, JPDA, GOSPA, and feature parameters; builds `params.mm` |
| **Scenarios** | `scenario_defs.m`, `apply_scenario.m` | Preset target scenarios (initial states, turn/accel segments); `params.scenario` selects one |
| **Ground Truth** | `generate_truth.m`, `motion_cv.m`, `motion_ct.m` | Piecewise CV/CT truth from `params.target_turns` / `params.target_accel` + spectrum templates |
| **Measurements** | `generate_measurements.m` | Struct array: `.bearing` (rad) + `.spectrum` (96×1 normed); includes Poisson clutter |
| **Motion Models** | `motion_cv.m`, `motion_ct.m`, `motion_ca.m`, `build_motion_models.m` | Model set definitions: `F`, `Q`, `f(x,dt)` for CV / CT / CA |
| **IMM Core** | `imm_init.m`, `imm_predict.m`, `imm_measurement_predict.m`, `imm_update.m`, `imm_combine.m`, `ci_fuse.m`, `stabilize_cov.m` | Mixing, per-model predict/update, probability update, CI fusion, covariance repair |
| **UKF Core** | `generate_sigma_points.m`, `measurement_bearing.m`, `ukf_measurement_predict.m`, `ukf_predict.m` | Pure kinematic — **never modify for feature work** |
| **JPDA Core** | `feasible_events.m`, `jpda_associate.m` | Validation matrix → feasible joint events → marginal β probabilities |
| **Orchestration** | `jpda_process_platform.m` | Per-platform: IMM measurement predict → JPDA → per-model update → spectrum update |
| **Management** | `manage_tracks.m`, `triangulate_init.m` | Track birth (triangulation + spectrum consistency check + IMM init) and death |
| **Evaluation** | `gospar.m`, `eval_gospa.m`, `compare_imm.m` | GOSPA metric, per-frame evaluation, model-set ablation (CV / CV+CT / CV+CT+CA) |
| **Entry** | `main.m`, `run_simulation.m`, `run_scenarios.m` | Simulation loop + history storage; single-run wrapper; multi-scenario overview |
| **Output** | `plot_results.m` | 7 figures (6 without features): trajectory (with line-of-sight bearing rays), bearings, position RMSE, velocity RMSE, GOSPA, IMM model probabilities + turn rate, spectrum templates |

### Key Design Decisions

- **Common 7-D state for all models**: `x = [px; py; vx; vy; omega; ax; ay]`. Sharing one state space makes mixing and fusion trivial; the CV model uses neither `omega` nor acceleration, CT uses `omega`, CA uses `(ax, ay)`, and each keeps the unused components unchanged so the IMM can hand them to the owning model. CV is linear (`F`), CT and CA are nonlinear and are propagated by sigma points.
- **IMM-JPDA coupling**: JPDA gates and marginal β are computed from the *probability-weighted combined* prediction (`z_comb`, `S_comb` of `imm_measurement_predict`), then each model is updated with the same β. Model probabilities are updated from the β-weighted innovation likelihood, and skipped when β0 ≥ 0.9 (no usable measurement).
- **Model probability fusion**: the two platform posteriors are combined by normalized geometric mean `mu ∝ sqrt(mu1 .* mu2)`, which corresponds to using each platform's likelihood once.
- **Sequential platform processing**: Platform 1 updated first, Platform 2 operates on a copy taken after prediction. This avoids joint multi-sensor JPDA complexity.
- **Measurements are struct arrays**: Each element has `.bearing` (scalar) and `.spectrum` (96×1 vector). Backward-compatible via `[z1.bearing]` to extract a numeric vector.
- **Feature aid is multiplicative in JPDA**: `L_total = L_kinematic × L_feature`. Features never enter the filter state.
- **Spectrum similarity**: Manual Pearson correlation (no toolbox), sharpened by `exp(gamma × max(0,ρ)²)`. Gamma=5 is conservative; increase for stronger feature reliance.
- **Feature toggle**: Set `params.use_features = false` to revert to pure bearing-only mode. All feature-related code paths are guarded by this flag.
- **Spectrum template caching**: `generate_truth.m` calls `generate_ship_spectrum()` once and caches in `truth.spectrum{tgt}`. `generate_measurements.m` reads from cache — never re-reads wav files per frame.
- **Numerical guards**: `stabilize_cov` clamps negative eigenvalues after UKF/JPDA/CI updates; per-model innovation variances are floored at `R` (theoretical `S ≥ R`). Without these, `alpha=1e-3` with a 7-D state occasionally loses positive definiteness.
- **Single-model baseline**: set `params.imm_models = {'CV'}`. With one model the IMM reduces exactly to the original CV-UKF (verified bit-for-bit on the straight-line scenario), which is the baseline used in `compare_imm.m`.

### Track Struct Fields

```matlab
track.x             % 7×1 combined state [x; y; vx; vy; omega; ax; ay]
track.P             % 7×7 combined covariance
track.imm.nx        % state dimension (7)
track.imm.mu        % [M×1] model probabilities (posterior)
track.imm.mu_prior  % [M×1] predicted model probabilities
track.imm.x         % [7×M] model-conditioned states
track.imm.P         % [7×7×M] model-conditioned covariances
track.id            % unique integer
track.confirmed     % logical (hit_count >= confirm_hits)
track.hit_count     % consecutive associations
track.miss_count    % consecutive misses (deleted at max_misses)
track.spectrum      % 96×1 normalized spectrum (NaN if feature_init==false)
track.feature_init  % logical
```

### Key Tuning Parameters

| Param | Default | Effect |
|-------|---------|--------|
| `Pd` | 0.9 | Detection probability per target per platform |
| `lambda_c` | 2 | Expected clutter per frame per platform |
| `sigma_bearing` | 1° | Bearing measurement noise std |
| `sigma_a_cv` | 0.5 m/s² | CV model acceleration noise (keep near 0.5: smaller values make the gate overconfident and cause track loss) |
| `sigma_a_ct` | 0.4 m/s² | CT model acceleration noise |
| `sigma_omega_cv` | 0.1°/s | CV model turn-rate random walk |
| `sigma_omega_ct` | 0.8°/s | CT model turn-rate random walk (larger = faster turn adaptation, noisier ω) |
| `omega_init_std` | 5°/s | Prior std of a new track's turn rate |
| `accel_init_std` | 0.3 m/s² | Prior std of a new track's acceleration |
| `sigma_jerk` | 0.05 m/s³ | CA model jerk noise |
| `sigma_accel_rw` | 0.1 m/s² | CV/CT acceleration-state random walk |
| `sigma_omega_ca` | 0.1°/s | CA turn-rate random walk (CA does not turn) |
| `imm_models` | `{'CV','CT','CA'}` | Model set; `{'CV'}` / `{'CV','CT'}` for ablations |
| `imm_transition` | `[0.94 0.03 0.03; …]` | 3×3 Markov model-switch matrix (sojourn ≈ 16 s) |
| `imm_mu0` | `[0.4; 0.4; 0.2]` | Initial model probabilities |
| `maneuver_enable` | true | `false` restores the original straight-line scenario |
| `scenario` | 'maneuver' | Preset truth scenario: `maneuver` / `parallel` / `crossing` / `speed` |
| `target_turns` | per scenario | Per-target `[start end deg/s]` turn segments (+ = CCW) |
| `target_accel` | per scenario | Per-target `[start end m/s²]` tangential acceleration segments |
| `feature_gamma` | 8 | Spectrum similarity sharpening; higher = more feature reliance |
| `max_misses` | 2 | Consecutive misses before deletion. With dual-platform OR detection, P(both miss) = (1−Pd)² = 0.01, so 2 frames ≈ 1e-4 false deletion; 3 lets ghost tracks live one frame longer (see `compare_track_mgmt.m`) |
| `confirm_hits` | 3 | Hits before track confirmed. Do **not** raise alone: GOSPA charges unconfirmed tracks as missed, so delaying confirmation hurts |
| `max_new_tracks` | 1 | Max births per frame. 3 lets a fresh ghost be spawned almost every frame and dominates the GOSPA false term; 1 only delays a genuine multi-target birth by one frame |
| `gospa_p` | 2 | GOSPA distance order (1 or 2; 2 = squared-error scale) |
| `gospa_c` | 100 | GOSPA cutoff (m); errors beyond `c` count as `c`. ≈5–10× nominal position error |
| `gospa_alpha` | 2 | GOSPA cardinality penalty; keep 2 so the metric stays symmetric |
| `gospa_use_confirmed` | true | Evaluate confirmed tracks only (`false` = include tentative tracks) |
| `plot_bearing_lines` | false | Overlay platform→target line-of-sight rays on the trajectory figure (set true to show) |
| `plot_bearing_every` | 25 | Snapshot interval (steps) for the bearing rays |

### IMM Motion Model

Model set (`build_motion_models.m`), common state `x = [px; py; vx; vy; ω; ax; ay]`:

```
CV:  px' = px + T·vx,  py' = py + T·vy,  vx' = vx,  vy' = vy,  ω' = ω,  a' = a
CT:  px' = px + sin(ωT)/ω·vx − (1−cos(ωT))/ω·vy
     py' = py + (1−cos(ωT))/ω·vx + sin(ωT)/ω·vy
     vx' = cos(ωT)·vx − sin(ωT)·vy
     vy' = sin(ωT)·vx + cos(ωT)·vy,  ω' = ω,  a' = a
CA:  px' = px + T·vx + T²/2·ax,  py' = py + T·vy + T²/2·ay
     vx' = vx + T·ax,  vy' = vy + T·ay,  a' = a,  ω' = ω
```

CV/CT use the discrete white-noise-acceleration `Q` on (position, velocity) plus a small random walk on the states they do not use; CA uses the white-noise-jerk `Q` on (position, velocity, acceleration). Only CV keeps a linear `F`; CT and CA are propagated through the nonlinear `f(x,dt)` by sigma points.

Each step `imm_predict` performs mixing (`c_j = Σ_i Π(i,j) μ_i`, `μ_{i|j} = Π(i,j)μ_i/c_j`), model-conditioned UKF prediction, and stores `mu_prior = c`. After the two platform updates, `imm_update` recomputes `μ_j ∝ c_j·N(ν_j; 0, S_j)` and the caller fuses the per-model estimates by covariance intersection and the probabilities by `sqrt(μ1·μ2)`.

`generate_truth.m` propagates the same CT kinematics with the commanded turn rate from `params.target_turns` and applies tangential acceleration from `params.target_accel` before each step, so filter and truth share one motion model definition. `truth.speed{tgt}` stores the resulting speed history.

### Scenarios

`scenario_defs.m` holds the preset truth scenarios; `apply_scenario.m` writes their initial states and maneuver segments into `params` (before user overrides, so individual targets can still be overridden). `main.m` exposes `SCENARIO = 'maneuver'` at the top; `compare_imm.m` takes the same `SCENARIO` variable.

| name | description |
|------|-------------|
| `maneuver` | Two targets make 6°/s coordinated turns (156° each) at steps 35–60 / 55–80 |
| `parallel` | Both targets uniform at (3,4) m/s, 400 m apart, constant separation |
| `crossing` | Equal-speed (7.8 m/s) crossing paths, 100 m closest approach at step 51 |
| `speed` | Target 1 uniform 5 m/s; target 2 accelerates 0.06 m/s² then decelerates back, 6→7.8→6 m/s |

`run_scenarios.m` runs `parallel` / `crossing` / `speed` over `N_MC` seeds and plots the three tracking results plus a GOSPA bar chart.

### GOSPA Evaluation

`gospar.m` computes the generalized optimal sub-pattern assignment distance (Rahmathullah, García-Fernández & Svensson, FUSION 2017) between a ground-truth set `X` (`d×n`) and a track-estimate set `Y` (`d×m`):

```
GOSPA_p^p = min_π [ Σ_i d(x_i, y_π(i))^p + (c^p/α)(n + m − 2|π|) ]
```

Each unmatched target or track costs `c^p/α`, so a pair is worth associating only when `d^p < 2c^p/α` (i.e. `d < c` for `α=2, p=2`); larger localization errors saturate at `c`. For `α=2` the metric is symmetric and needs no prior track-to-target correspondence.

`eval_gospa(truth, tracks_history, params)` builds `X_k` from all true target positions and `Y_k` from **confirmed** tracks only (`gospa_use_confirmed`), and returns per-frame `gospa`, the decomposition `loc / missed / false`, the counts `n_matched / n_missed / n_false`, plus averages. The identity `GOSPA^p = loc^p + missed^p + false^p` is the reason to report it alongside RMSE: RMSE only scores the tracks you decided to match, whereas GOSPA also charges for missed targets and false tracks.

```matlab
gos = eval_gospa(truth, tracks_history, params);
fprintf('GOSPA = %.2f m (loc %.2f, missed %.2f, false %.2f)\n', ...
    gos.mean, gos.mean_loc, gos.mean_missed, gos.mean_false);
```

`run_simulation` returns the time-averaged GOSPA as an optional third output and internal data (`tracks_history`, `truth`, `gos`) as an optional fourth output. `sweep_clutter.m` plots GOSPA vs clutter density next to the RMSE curves; `compare_imm.m` uses the same outputs for a paired CV-vs-IMM Monte Carlo. The partial assignment is solved by a self-contained Hungarian routine (no `matchpairs`/toolbox dependency).

### Spectrum Feature System

`generate_ship_spectrum.m` operates in two modes:
- **Simulation** (`wav_file=''`): Synthesizes LOFAR (random line spectra + Gaussian broadening, 10–500Hz) and DEMON (shaft+blade harmonics, 1–50Hz) → 96D vector.
- **Wav import** (valid `.wav` path): STFT → LOFAR extraction; bandpass → Hilbert envelope → FFT → DEMON extraction. Falls back to simulation on read failure.

New tracks created from triangulation require both platform measurements' spectra to have Pearson ρ ≥ 0.4 (consistency check to suppress ghost tracks from clutter pairs).
