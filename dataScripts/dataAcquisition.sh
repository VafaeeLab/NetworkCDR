#!/usr/bin/env bash

DP="../rawdata"

# Only keep database- or experiment-backed interactions with high confidence (>0.9)
sed -nE "s/^([^ ]+ [^ ]+) ([0-9]+ ){4}(([0-9]+ [1-9]+)|([1-9]+ [0-9]+)).*$/\1 \3/p" $DP/stringdb-human-raw_ppi.txt > $DP/stringdb-human-known_ppi.txt
grep -E "^[^ ]+ [^ ]+ (9[0-9]{2} [0-9]+)|([0-9]+ 9[0-9]{2})$" $DP/stringdb-human-known_ppi.txt > $DP/stringdb-human-known_ppi-confident.txt

# Get drugs which have target proteins
grep -hE "linkedTargets\":\{\"rows\":\[[A-Za-z0-9\",]+\]" $DP/molecule/*.json > $DP/opentargets-drug.json
sed -iE '$ ! s/$/,/' "$DP/opentargets-drug.json"
sed -i '1i[' $DP/opentargets-drug.json && echo ']' >> $DP/opentargets-drug.json



