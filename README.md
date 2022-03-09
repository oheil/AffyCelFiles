# AffyCelFiles
read Affymetrix .CEL files

# Examples
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



