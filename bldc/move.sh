#! /bin/bash

i=1
while [ -e "test${i}/" ];
do
    let i=i+1
done
mkdir test$i
mv bldc_* test$i/
