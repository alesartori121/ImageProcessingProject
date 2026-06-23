**Brain Tumor Edge Detection: Baseline Reproduction and Optimization Strategy**

**Overview**
This repository contains the MATLAB implementation, statistical validation, and algorithmic optimization of the MRI brain tumor edge detection pipeline originally proposed by Zotin et al. The project was developed as an academic study to evaluate the robustness of unsupervised clustering methods on large-scale clinical data and to propose concrete architectural improvements to overcome intrinsic topological limitations.

**Dataset**
The pipeline is validated against the Brain Tumor dataset (Task01: Brain) from the Medical Segmentation Decathlon (MSD) challenge. The dataset comprises 484 multiparametric MRI volumes, including T1, T1Gd, T2, and FLAIR modalities, alongside multi-class expert-annotated ground truth masks.

**Repository Structure and Methodology**

The repository is divided into two main environments to facilitate a direct A/B ablation study:

**1. Baseline Implementation (src/)**
A strict reproduction of the original methodology. It processes T2-weighted MRI slices and relies solely on intensity-based mathematics and topological heuristics.

* **Pre-processing:** Median filtering and robust Otsu-based brain masking.
* **Enhancement:** Balance Contrast Enhancement Technique (BCET) applied strictly within the brain mask to normalize the histogram and isolate hyperintense regions.
* **Segmentation:** Unsupervised Fuzzy C-Means (FCM) clustering ($C=4$) followed by a maximum connected component heuristic to isolate the tumor mass.
* **Edge Extraction & Evaluation:** Canny edge detection validated against the ground truth using rigorous spatial metrics, most notably Pratt's Figure of Merit (FOM).

**2. Optimized Pipeline (opt_src/)**
Developed after batch-processing the baseline and identifying critical failure points caused by the Partial Volume Effect and Cerebrospinal Fluid (CSF) interference. This environment introduces three major innovations:

* **Dynamic Slice Selection:** An automated 3D volume analysis that dynamically extracts the optimal Z-axis slice containing the maximum tumor area, discarding the static central-slice approach.
* **Multimodal Transition (FLAIR):** A radiometric switch from T2 to FLAIR sequence. By leveraging the physical suppression of free water (CSF) inherent to FLAIR imaging, the ambiguity of the FCM clustering is eliminated without the need for fragile topological filters.
* **Dilated Ground Truth Evaluation:** The introduction of a spatial tolerance band (1-pixel radius) in the evaluation metrics to bridge the mathematical rigidity of the Canny operator and the physiological variance of human expert tracing, significantly improving the clinical reliability of the Sensitivity metric.

**References**

* **Methodology:** Zotin, A., et al. (2018). Edge detection in MRI brain tumor images based on fuzzy C-means clustering. *Procedia Computer Science*, 136, 41-50.
* **Dataset:** Antonelli, M., et al. (2022). The Medical Segmentation Decathlon. *Nature Communications*, 13(1), 4128.
