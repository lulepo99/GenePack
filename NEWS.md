# GenePack 1.0.0

## Major changes 
* Initial release of the GenePack package.
* Each transcript isoform is treated as an individual gene, using Ensembl transcript IDs for flexibility.
* Added S4 classes to represent human gene types, including both protein-coding and non-coding genes.
* Added methods to create gene objects.
* Added validation functions to guarantee robustness.
* Added get and set methods to access and modify gene objects information.
* Added the `lengthProduct` method to compute the length of the gene products.
* Added the `showGeneObject` method to show all the gene objects information.