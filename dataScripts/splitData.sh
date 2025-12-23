#!/usr/bin/env bash

FILES=$(find ../clean/* -type f -size +50M)

for file in $FILES; do
    gzip -S .zip "$file"
    split -b50MB "$file.zip" "$file.zip."
    rm "$file.zip"
done

echo $FILES > "../clean/splitData.txt"