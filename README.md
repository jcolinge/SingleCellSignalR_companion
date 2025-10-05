# SingleCellSignalR_companion
Companion repository to SingleCellSignalR (Version 2)

&nbsp;

## Contents

This companion repository to `Version 2` of the Bioconductor
[SingleCellSignalR](https://www.bioconductor.org/packages/release/bioc/html/SingleCellSignalR.html)
library, we illustrate its usage with more flexibility than in the vignette or
the on-line documentation. In particular, we provide a basic example that
illustrates how a single-cell data set can be analyzed, which may serve as
starting point to write your own code.

We also provide three additional example scripts addressing particular data
types that are now amenable to SingleCellSignalR analysis:
- Single-cell proteomics data (scProt-MS);
- Bulk RNA-seq data obtained from FACS-separated cell populations;
- PDX bulk RNA-seq data that were aligned against both a human and a murine
  reference genome to separate cancer and stromal components.

These particular data types are discussed in our latest paper about
`version 2` of `SingleCellSignalR`.

Note that `Version 2` of `SingleCellSignalR` is a complete re-design,
which now relies on new S4 classes and on our other library
[BulkSignalR](https://www.bioconductor.org/packages/release/bioc/html/BulkSignalR.html)
that is used as foundation layer. Therefore, the usage of `SingleCellSignalR`
compared to the initial `Version 1` has completely changed. Nonetheless,
the original scoring scheme **LR-score** is still available with the
new implementation.