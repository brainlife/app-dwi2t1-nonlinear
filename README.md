<h1 align="center">DWI preprocessing and nonlinear registration to T1</h1>

<p align="center">
  Preprocessed diffusion data in T1 space, with rotated b-vectors and the ANTs transforms, from a raw DWI and a T1-weighted image.
</p>

<p align="center">
  <a href="https://github.com/brainlife/abcd-spec"><img src="https://img.shields.io/badge/ABCD_Spec-v1.1-green.svg" alt="ABCD Spec v1.1"></a>
  <a href="https://doi.org/10.25663/brainlife.app.111"><img src="https://img.shields.io/badge/Brainlife-brainlife.app.111-blue.svg" alt="Run on brainlife.io"></a>
  <a href="LICENSE"><img src="https://img.shields.io/badge/License-MIT-lightgrey.svg" alt="MIT license"></a>
</p>

## At a glance

<!-- brainlife-app:begin (generated from the App manifest; edit the manifest, not this block) -->

> **Preprocesses DWI with MRtrix3 dwipreproc (FSL eddy) and nonlinearly registers it to a T1-weighted image with ANTs, using CSF masks of the mean b0 and T1, rotating the b-vectors with the affine part.**

App for running mrtrix3 dwipreproc and then performing nonlinear registration of dwi to T1 (by running ANTs affine+nonlinear registration of CSF masks of mean b0 and T1). Bvectors are rotated with the affine transform. App is recommended for bringing dwi to T1 space if no phase encoding contrast between b=0 images was collected (so mrtrix3 dwipreproc cannot perform inhomogeneity field estimation).

| | |
|---|---|
| **Inputs** | `neuro/anat/t1w` · `neuro/dwi` |
| **Outputs** | `neuro/dwi` · `raw` |
| **Container** | `katealpert/dwi2t1-nonlinear:v0` |
| **Resources** | 8 CPUs · 24 GB memory · 24 h walltime |
| **Run it** | [doi.org/10.25663/brainlife.app.111](https://doi.org/10.25663/brainlife.app.111) |

#### Inputs

| Input | Datatype | Required | Description |
|---|---|---|---|
| `t1` | `neuro/anat/t1w` | yes | T1-weighted image of the same subject; bias-corrected and skull-stripped with fsl_anat and segmented with FSL FAST to get the CSF mask used as registration target. |
| `dwi` | `neuro/dwi` | yes | Raw DWI with bvecs and bvals to preprocess and bring into T1 space. Metadata must include the BIDS fields TotalReadoutTime and PhaseEncodingDirection, required by dwipreproc -rpe_header. |

#### Outputs

| Output | Datatype | Description |
|---|---|---|
| `odwi` | `neuro/dwi` | Preprocessed DWI resampled into T1 space with ANTs affine+SyN warps, with the bvals and affine-rotated bvecs. |
| `output` | `raw` | ANTs transforms from b0 to T1: affine_ANTs.mat, warp_ANTs.nii.gz and warp-inverse_ANTs.nii.gz. |

<!-- brainlife-app:end -->

## Running the App

### On brainlife.io

Open the App at [doi.org/10.25663/brainlife.app.111](https://doi.org/10.25663/brainlife.app.111), choose a raw `neuro/dwi` dataset and a T1-weighted image of the same subject, and press **Execute**. The App is recommended for bringing DWI into T1 space when no phase-encoding contrast between b=0 images was collected, so MRtrix3 `dwipreproc` cannot estimate the inhomogeneity field.

The App:

1. bias-corrects and skull-strips the T1 with `fsl_anat`;
2. preprocesses the DWI with MRtrix3 `dwipreproc` (FSL `eddy`, `-rpe_header`) and computes a brain mask and the mean b=0 image;
3. segments the mean b=0 and the T1 with FSL `fast` to obtain CSF masks;
4. registers the b=0 CSF mask to the T1 CSF mask with ANTs (translation, rigid, affine and SyN stages);
5. applies the affine and nonlinear warps to the preprocessed DWI and rotates the b-vectors with the affine part (HCP Pipelines `Rotate_bvecs.sh`).

It outputs a `neuro/dwi` dataset tagged `t1_aligned` and `preprocessed`, and a `raw` dataset tagged `warp` with the transforms.

### Locally

You need [Singularity/Apptainer](https://apptainer.org/) and [`jq`](https://jqlang.github.io/jq/). MRtrix3, ANTs 2.2.0, FSL 5.0.11, Convert3D and HCP Pipelines run inside the `katealpert/dwi2t1-nonlinear:v0` container (built from the [`Dockerfile`](Dockerfile) in this repository).

1. Clone this repository and check out the version brainlife.io runs: `git checkout v0`.
2. Create `config.json` with the paths to your inputs. `main` builds the BIDS JSON sidecars from the `meta` of each entry in `_inputs` (as brainlife.io writes it), and the DWI metadata must include `TotalReadoutTime` and `PhaseEncodingDirection`:

   ```json
   {
       "t1": "input/t1/t1.nii.gz",
       "dwi": "input/dwi/dwi.nii.gz",
       "bvecs": "input/dwi/dwi.bvecs",
       "bvals": "input/dwi/dwi.bvals",
       "_inputs": [
           { "id": "t1", "meta": { "subject": "01" } },
           { "id": "dwi", "meta": { "subject": "01", "TotalReadoutTime": 0.05, "PhaseEncodingDirection": "j-" } }
       ]
   }
   ```

3. Run the App:

   ```bash
   ./main
   ```

## Outputs

```text
dwi.nii.gz       preprocessed DWI resampled into T1 space
dwi.bvals        b-values
dwi.bvecs        b-vectors rotated with the affine transform
output/
├── affine_ANTs.mat            affine transform, b=0 to T1
├── warp_ANTs.nii.gz           nonlinear warp, b=0 to T1
└── warp-inverse_ANTs.nii.gz   inverse warp
```

## Troubleshooting

- **`TotalReadoutTime and PhaseEncodingDirection BIDS fields reqd to run mrtrix3 dwipreproc`**: the DWI metadata lacks the fields `dwipreproc` needs. Add `TotalReadoutTime` and `PhaseEncodingDirection` to the DWI dataset's metadata on brainlife.io (or to its `meta` in `config.json`).
- **`Missing or invalid input images`**: check the `t1`, `dwi`, `bvecs` and `bvals` paths in `config.json`.

## Citing

> [!TIP]
> If you use this App in your work, please cite the App, the software it builds on and brainlife.io (below). GitHub's **Cite this repository** button gives the App citation in APA and BibTeX.

- DWI preproc and nonlinear registration to T1. brainlife.io App. [doi:10.25663/brainlife.app.111](https://doi.org/10.25663/brainlife.app.111)
- Jenkinson et al. FSL. *NeuroImage* 62:782-790 (2012). [doi:10.1016/j.neuroimage.2011.09.015](https://doi.org/10.1016/j.neuroimage.2011.09.015)
- Tournier et al. MRtrix3: A fast, flexible and open software framework for medical image processing and visualisation. *NeuroImage* 202:116137 (2019). [doi:10.1016/j.neuroimage.2019.116137](https://doi.org/10.1016/j.neuroimage.2019.116137)
- Avants et al. A reproducible evaluation of ANTs similarity metric performance in brain image registration. *NeuroImage* 54:2033-2044 (2011). [doi:10.1016/j.neuroimage.2010.09.025](https://doi.org/10.1016/j.neuroimage.2010.09.025)
- Hayashi S, Caron BA, Heinsfeld AS, et al. brainlife.io: a decentralized and open-source cloud platform to support neuroscience research. *Nature Methods* 21:809–813 (2024). [doi:10.1038/s41592-024-02237-2](https://doi.org/10.1038/s41592-024-02237-2)
- Avesani P, McPherson B, Hayashi S, et al. The open diffusion data derivatives, brain data upcycling via integrated publishing of derivatives and reproducible open cloud services. *Scientific Data* 6:69 (2019). [doi:10.1038/s41597-019-0073-y](https://doi.org/10.1038/s41597-019-0073-y)

## Authors

- Kate Alpert

**Contributors:** Franco Pestilli

## Funding

[![Wellcome-226486/Z/22/Z](https://img.shields.io/badge/Wellcome-226486%2FZ%2F22%2FZ-purple.svg)](https://wellcome.org/)
[![NIH-NINDS-U24NS140384](https://img.shields.io/badge/NIH_NINDS-U24NS140384-green.svg)](https://reporter.nih.gov/project-details/11261037)
[![NIH-NINDS-UM1NS132207](https://img.shields.io/badge/NIH_NINDS-UM1NS132207-green.svg)](https://reporter.nih.gov/project-details/11411328)
[![NIH-NIBIB-R01EB029272](https://img.shields.io/badge/NIH_NIBIB-R01EB029272-green.svg)](https://reporter.nih.gov/project-details/10262925)
[![NIH-NIBIB-R01EB030896](https://img.shields.io/badge/NIH_NIBIB-R01EB030896-green.svg)](https://reporter.nih.gov/project-details/10428625)
[![NIH-NIMH-R01MH126699](https://img.shields.io/badge/NIH_NIMH-R01MH126699-green.svg)](https://reporter.nih.gov/project-details/10460628)
[![NIH-NIDCD-R01DC020088](https://img.shields.io/badge/NIH_NIDCD-R01DC020088-green.svg)](https://reporter.nih.gov/project-details/11140284)
[![NSF-BCS-1734853](https://img.shields.io/badge/NSF_BCS-1734853-blue.svg)](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1734853)
[![NSF-BCS-2203524](https://img.shields.io/badge/NSF_BCS-2203524-blue.svg)](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2203524)
[![NSF-IIS-1636893](https://img.shields.io/badge/NSF_IIS-1636893-blue.svg)](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1636893)
[![NSF-OAC-2148729](https://img.shields.io/badge/NSF_OAC-2148729-blue.svg)](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2148729)
[![NSF-IIS-1912270](https://img.shields.io/badge/NSF_IIS-1912270-blue.svg)](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1912270)
[![NSF-IIS-2148700](https://img.shields.io/badge/NSF_IIS-2148700-blue.svg)](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2148700)
[![NSF-OAC-1916518](https://img.shields.io/badge/NSF_OAC-1916518-blue.svg)](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1916518)

<sub>brainlife.io is publicly funded, and acknowledging that support helps keep it running. This App was developed with support from Wellcome, the U.S. National Institutes of Health (including the BRAIN Initiative BRAIN CONNECTS program: The Axonal Projectome EXchange, U24NS140384, and the Center for Mesoscale Connectomics, UM1NS132207) and the U.S. National Science Foundation.</sub>

## License

[MIT](LICENSE). MRtrix3, ANTs and FSL are distributed under their own licenses.
