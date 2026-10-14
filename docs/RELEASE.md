# Release checklist

1. Finalize the manuscript title, author list and software citation in
   `CITATION.cff`. Replace the pending identifiers in the README when available.
2. Upload the experimental files, 151 restart states, four numerical CSVs and
   both derived-data files, the three supplementary workbooks, and the
   supplementary spectrum/summary/recovery files to a Zenodo dataset record, with `SHA256SUMS` and
   the data-layout documentation. Add Figure 1's large file and its explanation.
   Select a data license separately from the MIT code license.
3. Publish the dataset and record its version-specific DOI in the README and
   citation metadata. Preserve the tested directory layout.
4. Upload only the code repository contents to GitHub. When staging with Git,
   inspect the staged files: exclude downloaded data and generated figures.
   This folder is inside another Git checkout locally, so publish it as the
   root of a separate repository rather than publishing the parent research tree.
5. Connect the GitHub repository to Zenodo and create a tagged release, such as
   `v1.0.0`. Archive that release and obtain the software DOI. Link the software
   and dataset records through their related identifiers.
6. Put the software and dataset version DOIs in the manuscript and update the
   software record metadata with the paper DOI once available. The concept DOI
   can additionally identify the evolving project; use version DOIs for the
   precise reproduction artifacts.

The source folder is prepared locally only. GitHub and Zenodo publication remain
manual. Zenodo instructions: [GitHub integration](https://help.zenodo.org/docs/github/)
and [DOI versioning](https://zenodo.org/help/versioning).
