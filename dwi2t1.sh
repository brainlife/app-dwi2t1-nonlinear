#!/bin/bash

t1=$1
dwi=$2
transformdir=$3
workdir=$4

outdir=$(pwd)

if [[ -z $t1 || ! -e $t1.nii.gz ]]; then
    echo "Missing T1 $t1.nii.gz"
    exit 1
fi
t1=$(pwd)/$t1.nii.gz

if [[ -z $dwi || ! -e $dwi.nii.gz || \
    ! -e $dwi.bval || ! -e $dwi.bvec || ! -e $dwi.json ]]; then
    echo "Missing $dwi image or bval or bvec or json"
    exit 1
fi
dwi=$(pwd)/$dwi

if [[ -z $transformdir ]]; then
    transformdir=output
fi
transformdir=$(pwd)/output

if [[ -z $workdir ]]; then
    workdir=workdir
fi
workdir=$(pwd)/workdir
mkdir -p $workdir
cd $workdir

#prep T1 (from BIDS mrtrix connectome dti preprocessing)
fsl_anat -i $t1 --noseg --nosubcortseg -o fsl
mrconvert fsl.anat/T1_biascorr_brain.nii.gz T1_biascorr_brain.mif

#prep DWI
mrconvert -json_import $dwi.json -fslgrad $dwi.bvec $dwi.bval \
    $dwi.nii.gz  dwi.mif
dwipreproc dwi.mif dwi_corrected.mif -rpe_header -eddy_options " --slm=linear" 
dwi2mask dwi.mif dwi_mask.mif
mrconvert dwi_mask.mif dwi_mask.nii.gz
mrconvert dwi_corrected.mif dwi_corrected.nii.gz \
    -export_grad_fsl dwi_corrected.bvec dwi_corrected.bval

#compute mean bzero
dwiextract dwi.mif -bzero - | \
    mrcalc - 0.0 -max - | \
    mrmath - mean -axis 3 dwi_meanbzero.mif

#instructions from Pierre Beson to use CSF->CSF
#make csf image from b0 and t1
mrconvert dwi_meanbzero.mif dwi_meanbzero.nii.gz
fslmaths dwi_meanbzero.nii.gz -mul dwi_mask.nii.gz dwi_meanbzero_brain.nii.gz
fast -t 2 dwi_meanbzero_brain.nii.gz
fslmaths dwi_meanbzero_brain_seg.nii.gz -uthr 1 dwi_meanbzero_csf.nii.gz

fast -t 1 fsl.anat/T1_biascorr_brain.nii.gz
fslmaths fsl.anat/T1_biascorr_brain_seg.nii.gz -uthr 1 T1_biascorr_csf.nii.gz

#affine+nonlinear reg
#https://github.com/ANTsX/ANTs/wiki/Anatomy-of-an-antsRegistration-call
#from https://www.neuroinf.jp/fmanager/view/737/bah20150723-alex.pdf
f=T1_biascorr_csf.nii.gz
m=dwi_meanbzero_csf.nii.gz
mask=fsl.anat/T1_biascorr_brain_mask.nii.gz
nm=b02t1
imgs=" $f, $m"
its=10000x1111x5
percentage=0.25
syn="20x20x0,0,5"
dim=3

antsRegistration -d $dim  --float 1 \
    -o [${nm},${nm}_diff.nii.gz,${nm}_inv.nii.gz] \
    -r [ $imgs ,1] \
    -m mattes[ $imgs , 1 , 32, regular, 0.05 ] -t translation[ 0.1 ] -c [1000,1.e-8,20] -s 4vox -f 6 -l 1 \
    -m mattes[ $imgs , 1 , 32, regular, 0.1 ] -t rigid[ 0.1 ] -c [1000x1000,1.e-8,20] -s 4x2vox -f 4x2 -l 1 \
    -m mattes[ $imgs , 1 , 32, regular, 0.1 ] -t affine[ 0.1 ] -c [$its,1.e-8,20] -s 4x2x1vox -f 3x2x1 -l 1 \
    -m mattes[ $imgs , 1 , 32 ] -t SyN[ .20, 3, 0 ] -c [ $syn ] -s 1x0.5x0vox -f 4x2x1 -l 1 -u 1 -z 1 -x $mask 

#get FSL-format affine transform matrix
ConvertTransformFile 3 ${nm}0GenericAffine.mat ${nm}0GenericAffine.txt
c3d_affine_tool -itk ${nm}0GenericAffine.txt -ref $f \
    -src $m  -ras2fsl -o ${nm}0GenericAffine_fsl.txt

#apply affine+nonlinear to DWI
antsApplyTransforms -e $dim -i dwi_corrected.nii.gz -r $f \
    -t ${nm}0GenericAffine.mat -t ${nm}1Warp.nii.gz \
    -o $outdir/dwi.nii.gz
cp dwi_corrected.bval $outdir/dwi.bvals

#apply affine reg to bvec
/opt/HCPpipelines/global/scripts/Rotate_bvecs.sh dwi_corrected.bvec \
    ${nm}0GenericAffine_fsl.txt $outdir/dwi.bvecs

if [[ -e dwi.nii.gz && -e dwi.bvals && -e dwi.bvecs ]]; then
    mkdir -p $transformdir
    mv ${nm}0GenericAffine.mat $transformdir/affine_ANTs.mat
    mv ${nm}1Warp.nii.gz $transformdir/warp_ANTs.nii.gz
    mv ${nm}1InverseWarp.nii.gz $transformdir/warp-inverse_ANTs.nii.gz
    exit 0
else
    exit 1
fi
