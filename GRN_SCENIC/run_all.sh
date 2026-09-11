#!/bin/bash

SCRIPT="/mnt/nfs/WORKSP/aminata.sow/LAB_DOCS/scripts/loom_analysis.py"

# List of all conditions
CONDITIONS=(
    "SScILD/Upper"
    "SScILD/Low"
    "IPF/Upper"
    "IPF/Low"
    "RA"
    "MS"
    "T1D"
    "SjD"
    "SS"
)

# Run each condition in a separate screen
for COND in "${CONDITIONS[@]}"; do
    SCREEN_NAME=$(echo $COND | tr '/' '_')
    echo "Starting screen: $SCREEN_NAME for condition: $COND"
    screen -dmS $SCREEN_NAME bash -c "
        source ~/miniforge3/etc/profile.d/conda.sh
        conda activate GRN_env
        python $SCRIPT $COND
        echo 'DONE: $COND'
        exec bash
    "
done

echo ""
echo "All conditions started! Check with: screen -ls"
echo "Attach to any screen with: screen -r <name>"
