#!/bin/bash
# This script runs the data preparation pipeline

# Get the project root directory
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# Define paths
SRC_DIR="${PROJECT_ROOT}/src/data_preparation"
PREPARE_SCRIPT="${SRC_DIR}/prepare_all_data.py"
MERGE_SCRIPT="${SRC_DIR}/merge_data.py"

# Print header
echo "=========================================="
echo "Running PrecisionOncology Data Preparation"
echo "=========================================="
echo "Project root: ${PROJECT_ROOT}"
echo "Source directory: ${SRC_DIR}"
echo "=========================================="

# Check if directories and scripts exist
if [ ! -d "${SRC_DIR}" ]; then
    echo "Error: Source directory does not exist: ${SRC_DIR}"
    exit 1
fi

if [ ! -f "${PREPARE_SCRIPT}" ]; then
    echo "Error: Prepare script does not exist: ${PREPARE_SCRIPT}"
    exit 1
fi

if [ ! -f "${MERGE_SCRIPT}" ]; then
    echo "Error: Merge script does not exist: ${MERGE_SCRIPT}"
    exit 1
fi

# Make scripts executable
chmod +x "${PREPARE_SCRIPT}" "${MERGE_SCRIPT}"

# Step 1: Run data preparation
echo -e "\n=========================================="
echo "Step 1: Running data preparation scripts"
echo "=========================================="
python "${PREPARE_SCRIPT}"

# Step 2: Merge data
echo -e "\n=========================================="
echo "Step 2: Merging prepared data"
echo "=========================================="
python "${MERGE_SCRIPT}"

echo -e "\n=========================================="
echo "Data preparation complete!"
echo "=========================================="