.PHONY: all sim move clean

all: clean sim move

sim:
	{ time octave bldc.m; } >> bldc_log.log 2>&1

move:
	./move.sh

clean:
	rm -f bldc_*.ans
	rm -f bldc_*.bmp
	rm -f bldc_*.fem
	rm -f bldc_*.log
	rm -f bldc_*.mat
	rm -f bldc_*.pdf
