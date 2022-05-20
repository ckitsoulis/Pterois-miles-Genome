# Transposable Elements annotation pipeline

For the transposable elements annotation of the genome assembly, a custom pipeline was used (see *image* below), based on state-of-the-art tools publicly available. Filtering and correction steps between tools were processed with custom python and bash scripts.

These tools are:

1. [RepeatModeler](https://github.com/Dfam-consortium/RepeatModeler)
2. [EDTA](https://github.com/oushujun/EDTA)
3. [DeepTE](https://github.com/LiLabAtVT/DeepTE)
4. [RepeatMasker](http://www.repeatmasker.org)


![Screenshot from 2022-03-24 20-36-45](https://user-images.githubusercontent.com/77507999/159987168-99681bf0-39ca-490e-ac35-f491de8dac45.png)

## Installation
All tools were run in conda environment and the installation instructions for each tool and their dependencies are provided in the corresponding repositories.

## Usage
All commands are provided in a pre-defined order (non-automatic usage for now) in the file *"TE_Annotation_workflow.sh"*.

## References
1.
2. Ou S., Su W., Liao Y., Chougule K., Agda J. R. A., Hellinga A. J., Lugo C. S. B., Elliott T. A., Ware D., Peterson T., Jiang N.envelope, Hirsch C. N.envelope and Hufford M. B.envelope (2019). Benchmarking Transposable Element Annotation Methods for Creation of a Streamlined, Comprehensive Pipeline. Genome Biol. 20(1): 275.
3. Haidong Yan, Aureliano Bombarely, Song Li 2020 DeepTE: a computational method for de novo classification of transposons with convolutional neural network. Bioinformatics, Volume 36, Issue 15, 1 August 2020, Pages 4269–4275.
4. 
