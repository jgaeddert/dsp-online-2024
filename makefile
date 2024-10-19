# simple makefile for generating figures

all:

figures :=			\
	interp.png		\
	interp_noise.png	\
	interp_rot_0.png    \

interp.png          : figures/interp.py ; ./$< -o $@ -plotsyms
interp_noise.png    : figures/interp.py ; ./$< -o $@ -nstd 0.1
interp_rot_0.png    : figures/interp.py ; ./$< -o $@ -fc 0.00208333


all: ${figures}

clean:
	rm -f ${figures}

