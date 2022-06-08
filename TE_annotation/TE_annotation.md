# Transposable Elements annotation pipeline

For the transposable elements annotation of genome assembly, a custom pipeline was used (see *image* below), based on state-of-the-art tools publicly available. Filtering and correction steps between tools were processed with custom python and bash scripts.

Tools used:

- [RepeatModeler](https://github.com/Dfam-consortium/RepeatModeler)
- [EDTA](https://github.com/oushujun/EDTA)
- [DeepTE](https://github.com/LiLabAtVT/DeepTE)
- [RepeatMasker](http://www.repeatmasker.org)
- [RM_parser](https://github.com/ckitsoulis/Pterois-miles-Genome/blob/main/TE_annotation/RM_parser.py)


![TE_diagram](https://github.com/ckitsoulis/Pterois-miles-Genome/blob/main/TE_annotation/TE_diagram.png)

## Installation
All tools were run in conda environment and the installation instructions for each tool and their dependencies are provided in the corresponding repositories.

## Usage
All commands are provided in a pre-defined order (non-automatic usage for now) in the file *"TE_Annotation_workflow.sh"*.

## References
1. Flynn J.M., Hubley R., Goubert C., Rosen J., Clark A.G., Feschotte C., Smit A.F. (2020). RepeatModeler2 for automated genomic discovery of transposable element families. Proc Natl Acad Sci U S A. 117(17): 9451-9457. doi: 10.1073/pnas.1921046117.
2. Ou S., Su W., Liao Y., Chougule K., Agda J.R.A., Hellinga A.J., Lugo C.S.B., Elliott T.A., Ware D., Peterson T., Jiang N., Hirsch C.N. and Hufford M.B. (2019). Benchmarking Transposable Element Annotation Methods for Creation of a Streamlined, Comprehensive Pipeline. Genome Biol. 20(1): 275.
3. Haidong Yan, Aureliano Bombarely, Song Li 2020 DeepTE: a computational method for de novo classification of transposons with convolutional neural network. Bioinformatics, Volume 36, Issue 15, 1 August 2020, Pages 4269–4275.
4. Tarailo-Graovac M., Chen N. (2009). Using RepeatMasker to identify repetitive elements in genomic sequences. Curr Protoc Bioinformatics.Chapter 4, Unit 4.10. doi: 10.1002/0471250953.bi0410s25
