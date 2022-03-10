# AffyCelFiles

This is a low level IO package for Affymetrix .CEL files

**References**

https://media.affymetrix.com/support/developer/powertools/changelog/file-formats.html
https://media.affymetrix.com/support/developer/powertools/changelog/gcos-agcc/cel.html
https://media.affymetrix.com/support/developer/powertools/changelog/gcos-agcc/generic.html
https://media.affymetrix.com/support/developer/powertools/changelog/gcos-agcc/cdf.html
https://media.affymetrix.com/support/developer/powertools/changelog/file-format-pgf.html
https://media.affymetrix.com/support/developer/powertools/changelog/file-format-clf.html
https://media.affymetrix.com/support/developer/powertools/changelog/file-format-mps.html

## Currently supported microarrays

This package is in an early stage and the following microarrays .CEL files are tested:
* Clariom_S_Human
* Mouse430_2
The expression/intensity data is cross checked against R/bioconductor affy and oligo package


# Examples

```
using AffyCelFiles

cel_file = "example_(Mouse430_2).CEL"

cel_data = AffyCelFiles.cel_read(cel_file);

cdf_file = raw"d:\Temp\mouse430_2_libraryfile\CD_Mouse430_2\Full\Mouse430_2\LibFiles\Mouse430_2.cdf"

cdf_data = AffyCelFiles.cdf_read(cdf_file);

intensities = AffyCelFiles.intensities(cel_data,cdf_data);

intensities.pm

intensities.mm

cel_file = "example_(Clariom_S_Human).CEL"

cel_data = AffyCelFiles.cel_read(cel_file);

pgf_file = raw"d:\Temp\clariomShuman\TFS-Assets_LSG_Support-Files_Clariom_S_Human_Analysis-r1\Clariom_S_Human.r1.pgf"

pgf_data = AffyCelFiles.pgf_read(pgf_file);

clf_file  =raw"d:\Temp\clariomShuman\TFS-Assets_LSG_Support-Files_Clariom_S_Human_Analysis-r1\Clariom_S_Human.r1.clf"

clf_data = AffyCelFiles.clf_read(clf_file);

mps_file = raw"d:\Temp\clariomShuman\TFS-Assets_LSG_Support-Files_Clariom_S_Human_Analysis-r1\Clariom_S_Human.r1.mps"

mps_data = AffyCelFiles.mps_read(mps_file);

intensities = AffyCelFiles.intensities(cel_data, pgf_data, clf_data);

intensities = AffyCelFiles.intensities(cel_data, pgf_data, clf_data, mps_data);

intensities.pm

intensities.mm
```



