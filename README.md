# DWI preproc and nonlinear registration to T1

<!-- brainlife-app:begin (generated from the App manifest; edit the manifest, not this block) -->

**Preprocesses DWI with MRtrix3 dwipreproc (FSL eddy) and nonlinearly registers it to a T1-weighted image with ANTs, using CSF masks of the mean b0 and T1, rotating the b-vectors with the affine part.**

App for running mrtrix3 dwipreproc and then performing nonlinear registration of dwi to T1 (by running ANTs affine+nonlinear registration of CSF masks of mean b0 and T1). Bvectors are rotated with the affine transform. App is recommended for bringing dwi to T1 space if no phase encoding contrast between b=0 images was collected (so mrtrix3 dwipreproc cannot perform inhomogeneity field estimation).

Run it on brainlife.io: [https://doi.org/10.25663/brainlife.app.111](https://doi.org/10.25663/brainlife.app.111)

| Input | Datatype | Required | Description |
|---|---|---|---|
| `t1` | `neuro/anat/t1w` | yes | T1-weighted image of the same subject; bias-corrected and skull-stripped with fsl_anat and segmented with FSL FAST to get the CSF mask used as registration target. |
| `dwi` | `neuro/dwi` | yes | Raw DWI with bvecs and bvals to preprocess and bring into T1 space. Metadata must include the BIDS fields TotalReadoutTime and PhaseEncodingDirection, required by dwipreproc -rpe_header. |

| Output | Datatype | Description |
|---|---|---|
| `odwi` | `neuro/dwi` | Preprocessed DWI resampled into T1 space with ANTs affine+SyN warps, with the bvals and affine-rotated bvecs. |
| `output` | `raw` | ANTs transforms from b0 to T1: affine_ANTs.mat, warp_ANTs.nii.gz and warp-inverse_ANTs.nii.gz. |

Runtime: container `katealpert/dwi2t1-nonlinear:v0`; 8 CPUs; 24 GB memory; 24 h walltime.

<!-- brainlife-app:end -->

