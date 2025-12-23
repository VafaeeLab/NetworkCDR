#!/usr/bin/env bash

FILES=$(cat clean/splitData.txt)

for file in $FILES; do
    if [ -e $file ]; then
        echo "File $file already exists"
        exit 1
    fi

    cat $file.zip.?? > "$file.zip"
    rm $file.zip.??
    gunzip -S .zip "$file.zip"
done


