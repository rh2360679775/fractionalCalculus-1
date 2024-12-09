#!/bin/bash

# Directory containing the .cpp files (default: current directory)
SOURCE_DIR=${1:-.}

# Output file for the SOURCES list
OUTPUT_FILE="sources.cmake"

# Find all .cpp files, sort them alphabetically, remove the ./ prefix, and save them in CMake-compatible format
echo "# Auto-generated sources list" > "$OUTPUT_FILE"
echo "set(SOURCES" >> "$OUTPUT_FILE"

find "$SOURCE_DIR" -type f -name "*.cpp" | sed 's|^\./||' | sort | while read -r FILE; do
    echo "    ${FILE}" >> "$OUTPUT_FILE"
done

echo ")" >> "$OUTPUT_FILE"

echo "Source list written to $OUTPUT_FILE"
