# ATF_RNAsequencing2024
Scripts used for differential expression analysis of processed RNA sequencing count data.

**Format**
The Metadata file should contain columns "name", with each sample grouped by name, and "treatment". The count data file should contain columns "Geneid", and the rest of the columns should match the "name" column in the metadata file; each sample column should contain raw count data.

The script is designed for sample of a protein (ATF3), plasmid control (Aart6), and buffer (Shock) comparisons.
