#! /bin/bash

i=1
while [ -e "test${i}/" ];
do
    let i=i+1
done
mkdir test$i
mv bldc_*.ans test$i/
mv bldc_*.bmp test$i/
mv bldc_*.fem test$i/
mv bldc_*.log test$i/
mv bldc_*.mat test$i/
mv bldc_*.pdf test$i/
