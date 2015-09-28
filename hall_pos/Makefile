.PHONY: all sim move clean

all: clean sim

sim:
	{ time octave hall_pos.m; } >> hall_pos_log.log 2>&1

clean:
	rm -f hall_pos.ans
	rm -f hall_pos.fem
	rm -f hall_pos_*.log
	rm -f hall_pos_*.mat
	rm -f hall_pos_*.pdf
	rm -f hall_pos_*.txt
