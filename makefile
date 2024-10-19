# simple makefile for generating figures

all:

figures :=			\
	interp.png		\
	interp_clean.png	\
	interp_noise.png	\
	interp_rot_0.png    	\
	interp_rot_1.png    	\
	interp_rot_2.png    	\

interp.png          : figures/interp.py ; ./$< -o $@ -plotsyms
interp_clean.png    : figures/interp.py ; ./$< -o $@
interp_noise.png    : figures/interp.py ; ./$< -o $@ -nstd 0.1
interp_rot_0.png    : figures/interp.py ; ./$< -o $@ -fc 0.0005
interp_rot_1.png    : figures/interp.py ; ./$< -o $@ -fc 0.00208333
interp_rot_2.png    : figures/interp.py ; ./$< -o $@ -fc 0.00416667


all: ${figures}

clean:
	rm -f ${figures}

