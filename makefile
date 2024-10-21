# simple makefile for generating figures

all:

figures :=			\
	interp.png		\
	interp_clean.png	\
	interp_noise.png	\
	interp_short_0.png    	\
	interp_short_1.png    	\
	interp_short_2.png    	\
	interp_short_3.png    	\
	interp_long_0.png    	\
	interp_long_1.png    	\
	interp_long_2.png    	\
	interp_long_3.png    	\
	interp_long_4.png    	\
	interp_long_5.png    	\
	interp_long_6.png    	\
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
	qpart_0.png		\
	qpart_1.png		\
	qpart_2.png		\

opts := -N 240 -xticks 30

interp.png          : figures/interp.py    ; ./$< -o $@ ${opts} -plotsyms
interp_clean.png    : figures/interp.py    ; ./$< -o $@ ${opts}
interp_noise.png    : figures/interp.py    ; ./$< -o $@ ${opts} -nstd 0.1

interp_short_0.png  : figures/interp.py    ; ./$< -o $@ ${opts} -plotcor -R 30 -fc 0
interp_short_1.png  : figures/interp.py    ; ./$< -o $@ ${opts} -plotcor -R 30 -fc 0
interp_short_2.png  : figures/interp.py    ; ./$< -o $@ ${opts} -plotcor -R 30 -fc 0 -nstd 0.1 -scale=0.95
interp_short_3.png  : figures/interp.py    ; ./$< -o $@ ${opts} -plotcor -R 30 -fc 0 -nstd 0.3 -scale=0.7

interp_long_0.png   : figures/interp.py    ; ./$< -o $@ ${opts} -plotcor -fc 0 -nstd 0.3 -scale=0.7
interp_long_1.png   : figures/interp.py    ; ./$< -o $@ ${opts} -plotcos -plotcor -fc 0
interp_long_2.png   : figures/interp.py    ; ./$< -o $@ ${opts} -plotcos -plotcor -fc 0.0005
interp_long_3.png   : figures/interp.py    ; ./$< -o $@ ${opts} -plotcos -plotcor -fc 0.00208333
interp_long_4.png   : figures/interp.py    ; ./$< -o $@ ${opts} -plotcos -plotcor -fc 0.00416667
interp_long_5.png   : figures/interp.py    ; ./$< -o $@ ${opts} -plotcos -plotcor -fc 0.00833333
interp_long_6.png   : figures/interp.py    ; ./$< -o $@ ${opts} -plotcos -plotcor -fc 0.02

partition.png       : figures/partition.py ; ./$< -o $@ ${opts} -plotsyms
partition_comp.png  : figures/partition.py ; ./$< -o $@ ${opts} -plotcomp
partition_ex_fc.png : figures/partition.py ; ./$< -o $@ ${opts} -plotimag -plotcos -fc 0.0083333
partition_ap_fc.png : figures/partition.py ; ./$< -o $@ ${opts} -plotimag -plotcos -fc 0.0083333 -fcapprox
partition_rot_0.png : figures/partition.py ; ./$< -o $@ ${opts} -plotimag -plotcos -plotcor -fc 0
partition_rot_1.png : figures/partition.py ; ./$< -o $@ ${opts} -plotimag -plotcos -plotcor -fc 0.0005
partition_rot_2.png : figures/partition.py ; ./$< -o $@ ${opts} -plotimag -plotcos -plotcor -fc 0.00208333
partition_rot_3.png : figures/partition.py ; ./$< -o $@ ${opts} -plotimag -plotcos -plotcor -fc 0.00416667
partition_rot_4.png : figures/partition.py ; ./$< -o $@ ${opts} -plotimag -plotcos -plotcor -fc 0.00833333
partition_rot_5.png : figures/partition.py ; ./$< -o $@ ${opts} -plotimag -plotcos -plotcor -fc 0.02

qpart_0.png : %.png : plot_qpartition.py qpartition.py ; ./$< -o $*
qpart_1.png : %.png : plot_qpartition.py qpartition.py ; ./$< -o $* -fc 0.0005
qpart_2.png : %.png : plot_qpartition.py qpartition.py ; ./$< -o $* -fc 0.0005 -dt 17

all: ${figures}

clean:
	rm -f ${figures}

