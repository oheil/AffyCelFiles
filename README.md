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

(more to come)

The expression/intensity data was cross checked against R/bioconductor affy and oligo package.

## Dependencies

#### Julia versions

* Julia 1.0 or above
* CRC32c (Standard library)

## Usage

In general Affymetrix .CEL can always be read in. To read Affymetrix .CEL files with meaningful information you need one of the following additional file/files:
* .cdf
* .pgf + .clf (+ optional .mps)

Those files provide the mapping from a probe location on the chip to a probe_id and the mapping from probeset_ids to probe_ids. For a specific chip type you
can get those files from the ThermoFisher support page, e.g. for the Clariom S human see https://www.thermofisher.com/order/catalog/product/902927?SID=srch-srp-902927
for download of archive TFS-Assets_LSG_Support-Files_Clariom_S_Human_Analysis-r1.zip (you need to have a login) which contains the .pgf,.clf and .mps file for this chip.

For biological meaningful analysis you also need annotation data, which maps probeset_ids to, for example, gene names. This is not part of this package, but
affymetrix annotation data is typically just a file of annotations in a table/csv format which can easily be read in (e.g. with CSV.jl) and mapped using probeset_ids or transcript_ids/meta_probeset_ids (in case of .mps used).

In this early stage this package doesn't provide a very convenient API, but this may change in future.

# Usage examples

```julia
using Pkg
Pkg.add("AffyCelFiles");
#or from github master:
#Pkg.add(url="https://github.com/oheil/AffyCelFiles.jl",rev="master")

using AffyCelFiles
```

Reading .CEL files with corresponding .cdf file:
```julia
cel_file = "example_(Mouse430_2).CEL"
cdf_file = raw"d:\Temp\mouse430_2_libraryfile\CD_Mouse430_2\Full\Mouse430_2\LibFiles\Mouse430_2.cdf"

cel_data = AffyCelFiles.cel_read(cel_file);
cdf_data = AffyCelFiles.cdf_read(cdf_file);

intensities = AffyCelFiles.intensities(cel_data,cdf_data);

#Dict mapping probeset_id to lists of expression values (pm=perfect match, mm=mismatch).
#Depending on the chiptype mm can be empty
intensities.pm
intensities.mm
```
Reading .CEL files with corresponding .pgf and .clf files (you need both):
```julia
cel_file = "example_(Clariom_S_Human).CEL"
pgf_file = raw"d:\Temp\clariomShuman\TFS-Assets_LSG_Support-Files_Clariom_S_Human_Analysis-r1\Clariom_S_Human.r1.pgf"
clf_file  =raw"d:\Temp\clariomShuman\TFS-Assets_LSG_Support-Files_Clariom_S_Human_Analysis-r1\Clariom_S_Human.r1.clf"

cel_data = AffyCelFiles.cel_read(cel_file);
pgf_data = AffyCelFiles.pgf_read(pgf_file);
clf_data = AffyCelFiles.clf_read(clf_file);

intensities = AffyCelFiles.intensities(cel_data, pgf_data, clf_data);

#Dict mapping probeset_id to lists of expression values (pm=perfect match, mm=mismatch).
#Depending on the chiptype mm can be empty
intensities.pm
intensities.mm
```
Providing a .mps file changes the central ids from probeset_ids to meta_probeset_ids/transcript_ids.
Reading .CEL files with corresponding .pgf, .clf files (you need both) and optional a .mps file:
```julia
cel_file = "example_(Clariom_S_Human).CEL"
pgf_file = raw"d:\Temp\clariomShuman\TFS-Assets_LSG_Support-Files_Clariom_S_Human_Analysis-r1\Clariom_S_Human.r1.pgf"
clf_file  =raw"d:\Temp\clariomShuman\TFS-Assets_LSG_Support-Files_Clariom_S_Human_Analysis-r1\Clariom_S_Human.r1.clf"
mps_file = raw"d:\Temp\clariomShuman\TFS-Assets_LSG_Support-Files_Clariom_S_Human_Analysis-r1\Clariom_S_Human.r1.mps"

cel_data = AffyCelFiles.cel_read(cel_file);
pgf_data = AffyCelFiles.pgf_read(pgf_file);
clf_data = AffyCelFiles.clf_read(clf_file);
mps_data = AffyCelFiles.mps_read(mps_file);

intensities = AffyCelFiles.intensities(cel_data, pgf_data, clf_data, mps_data);

#Dict mapping meta_probeset_id/transcript_id to lists of expression values (pm=perfect match, mm=mismatch).
#Depending on the chiptype mm can be empty
intensities.pm
intensities.mm
```










