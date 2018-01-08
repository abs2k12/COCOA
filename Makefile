all: sim2obs/sim2obs fit-king-5 sbpvdp clean_overlap

.PHONY: sim2obs/sim2obs clean
sim2obs/sim2obs:
	cd sim2obs && make

fit-king-5: fit-king5.f
	gfortran -fdefault-real-8 -w fit-king5.f -o fit-king5

sbpvdp: sbpvdp.f
	gfortran sbpvdp.f -o sbp-vdp

clean_overlap: clean_overlap.f
	gfortran clean_overlap.f -o clean_overlap

clean:
	rm -rf fit-king 5 sim2obs/sim2obs
