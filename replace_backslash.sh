#!/bin/bash

FILE=$1
sed -i 's/\\d/davidd/g' $FILE
sed -i 's/\\n/davidn/g' $FILE
sed -i 's/\x27\\/filesep \x27/g' $FILE
sed -i 's/\\\x27/\x27 filesep/g' $FILE
sed -i 's/\\/\x27 filesep \x27/g' $FILE
sed -i 's/davidd/\\d/g' $FILE
sed -i 's/davidn/\\n/g' $FILE
