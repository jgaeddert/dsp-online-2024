# simple makefile for generating figures

all:

figures :=			\
	interp.png		\
	interp_clean.png	\
	interp_noise.png	\
	interp_rot_0.png    	\
	interp_rot_1.png    	\
	interp_rot_2.png    	\
	interp_rot_3.png    	\
	interp_rot_4.png    	\
	interp_rot_5.png    	\
	partition.png		\
	partition_comp.png	\
	partition_ex_fc.png	\
	partition_ap_fc.png	\
	partition_rot_0.png	\
	partition_rot_1.png	\
	partition_rot_2.png	\
	partition_rot_3.png	\
	partition_rot_3.png	\
	partition_rot_4.png	\
	partition_rot_5.png	\

interp.png          : figures/interp.py    ; ./$< -o $@ -plotsyms
interp_clean.png    : figures/interp.py    ; ./$< -o $@
interp_noise.png    : figures/interp.py    ; ./$< -o $@ -nstd 0.1
interp_rot_0.png    : figures/interp.py    ; ./$< -o $@ -plotcos -plotcor -fc 0
interp_rot_1.png    : figures/interp.py    ; ./$< -o $@ -plotcos -plotcor -fc 0.0005
interp_rot_2.png    : figures/interp.py    ; ./$< -o $@ -plotcos -plotcor -fc 0.00208333
interp_rot_3.png    : figures/interp.py    ; ./$< -o $@ -plotcos -plotcor -fc 0.00416667
interp_rot_4.png    : figures/interp.py    ; ./$< -o $@ -plotcos -plotcor -fc 0.00833333
interp_rot_5.png    : figures/interp.py    ; ./$< -o $@ -plotcos -plotcor -fc 0.02

partition.png       : figures/partition.py ; ./$< -o $@
partition_comp.png  : figures/partition.py ; ./$< -o $@ -plotcomp
partition_ex_fc.png : figures/partition.py ; ./$< -o $@ -plotimag -plotcos -fc 0.0083333
partition_ap_fc.png : figures/partition.py ; ./$< -o $@ -plotimag -plotcos -fc 0.0083333 -fcapprox

partition_rot_0.png : figures/partition.py ; ./$< -o $@ -plotimag -plotcos -plotcor -fc 0
partition_rot_1.png : figures/partition.py ; ./$< -o $@ -plotimag -plotcos -plotcor -fc 0.0005
partition_rot_2.png : figures/partition.py ; ./$< -o $@ -plotimag -plotcos -plotcor -fc 0.00208333
partition_rot_3.png : figures/partition.py ; ./$< -o $@ -plotimag -plotcos -plotcor -fc 0.00416667
partition_rot_4.png : figures/partition.py ; ./$< -o $@ -plotimag -plotcos -plotcor -fc 0.00833333
partition_rot_5.png : figures/partition.py ; ./$< -o $@ -plotimag -plotcos -plotcor -fc 0.02

all: ${figures}

clean:
	rm -f ${figures}

