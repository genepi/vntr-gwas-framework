#!/bin/bash

# Ensure output folder exists
mkdir -p lpa
unset DX_WORKSPACE_ID
dx cd $DX_PROJECT_CONTEXT_ID:

# Download reference genome files
dx download "Bulk/Exome sequences/Exome OQFE CRAM files/helper_files/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai"
dx download "Bulk/Exome sequences/Exome OQFE CRAM files/helper_files/GRCh38_full_analysis_set_plus_decoy_hla.fa"

while IFS= read -r ID; do
    # Use a safe filename for download
    FNAME=$(basename "$ID" .cram)           # strips .cram
    LOCAL="$FNAME.cram"                      # temporary local CRAM
    OUT="lpa/${FNAME}_lpa.bam"               # final output BAM

    echo "Downloading $ID..."
    dx download "$ID" -o "$LOCAL"            # download the CRAM
    dx download "${ID}.crai" -o "${LOCAL}.crai"

    # Wait 2 seconds to ensure downloads are fully flushed
    sleep 2

    echo "Extracting region chr6:160530483-160665260 to BAM..."
    samtools view -b -T GRCh38_full_analysis_set_plus_decoy_hla.fa \
        -o "$OUT" "$LOCAL" chr6:160530483-160665260

    echo "Deleting temporary CRAM: $LOCAL"
    rm -f "$LOCAL" "$LOCAL.crai"

    dx upload "$OUT" --destination "$DX_PROJECT_CONTEXT_ID:/CRAMS/"
    rm "$OUT"
    echo "Done with $ID"

    # Wait 2 seconds before starting the next iteration
    sleep 2

done < ids_by_ancestry.txt
```