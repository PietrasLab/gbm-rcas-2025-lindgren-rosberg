#!/bin/bash
cd ../../references/cellranger
tail -n +7 Chromium_Mouse_Transcriptome_Probe_Set_v1.0.1_mm10-2020-A-flex-rcas-v3.csv | cut -d ',' -f 3,2 > Chromium_Mouse_Transcriptome_Probe_Set_v1.0.1_mm10-2020-A-flex-rcas-v3-split.csv
sed -i '' 's/,/|/g' Chromium_Mouse_Transcriptome_Probe_Set_v1.0.1_mm10-2020-A-flex-rcas-v3-split.csv
