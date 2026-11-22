# Camera pose estimation from points (PnP)

MATLAB lab from the **Robot Vision** course (Télécom Physique Strasbourg).

## Objective

Given 6 known 3D points on an object and their projections in a camera image, estimate the object's pose
(rotation **R** and translation **t**) in the camera frame.
Seven methods are implemented and compared on simulated poses, with increasing image noise.

## How it works

`main_canvas.m` runs the whole pipeline, step by step:

1. **Simulation** – each object point is projected into the image with the camera matrix K, then moved by a
   noise of fixed amplitude in a random direction (0 to 2 pixels). All methods are tested on the same noisy points.
2. **Error criteria** – for every estimated pose: 3D error on the object points (mm), translation and rotation
   errors, and reprojection error in the image (pixels).
3. **Pose estimation**, with seven methods:

| # | Method | Principle | File |
|---|---|---|---|
| 1 | Analytical (EPnP) | closed-form solution based on 4 virtual control points | `EPnP/efficient_pnp.m` |
| 2–3 | Linear | the 12 coefficients of [R \| t] solved by SVD, then projection onto a true rotation | `Linear_PnP.m` |
| 4 | Numerical optimization | Levenberg-Marquardt (`lsqnonlin`) on 6 parameters: axis-angle rotation, R = expm([w]×), and translation | `compute_err_persp_proj.m` |
| 5 | Orthogonal Iteration | iterative method minimizing the distance of each point to its line of sight | `OI_algo.m` |
| 6 | Auto-calibration | linear estimation of M = K[R \| t], RQ decomposition, then nonlinear refinement of the intrinsics and pose | `autocalibration_solving.m` |
| 7 | Stereoscopy | second camera, triangulation, 3D-3D registration by SVD, then joint optimization on both images | `pose_3D3D.m`, `compute_err_stereo_persp_proj.m` |

4. **Statistics and plots** – median and maximum errors per noise level, computation time, spatial distribution
   of the errors, and a final comparison of all methods.

## Code organization

```
├── main_canvas.m                      # main script: simulation, the 7 methods, statistics, figures
│
├── compute_err_pts3D.m                # 3D error on the object points
├── compute_err_rep_in_pix.m           # reprojection error in pixels
├── Linear_PnP.m                       # linear method
├── compute_err_persp_proj.m           # cost function for the numerical optimization
├── OI_algo.m                          # Orthogonal Iteration
├── autocalibration_solving.m          # pose + intrinsic parameters from a single image
├── pose_3D3D.m                        # rigid registration between two 3D point sets
├── compute_err_stereo_persp_proj.m    # cost function for the stereo optimization
│
├── skew.m, r2thetau.m, theta2r.m      # rotation utilities
├── plot3Derrors.m                     # 3D display of the errors
├── EPnP/                              # EPnP implementation (third-party, GPLv3)
└── data/data_tutorial_PnP.mat         # simulated poses and camera matrix K
```

## Run

Requirements: MATLAB with the **Optimization Toolbox**.

```matlab
main_canvas   % press a key once the data is loaded
```

## Credits

EPnP: F. Moreno-Noguer, V. Lepetit, P. Fua, *Accurate Non-Iterative O(n) Solution to the PnP Problem*
(code distributed under the GNU GPL v3, see the headers in `EPnP/`).
Orthogonal Iteration: C.-P. Lu, G. D. Hager, E. Mjolsness, *Fast and Globally Convergent Pose Estimation from Video Images*.
