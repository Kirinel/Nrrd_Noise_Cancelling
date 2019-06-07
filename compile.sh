#!/bin/bash

SCIVIS=/Users/yanchuanqi/UChicago/Winter_2019/Scientific_Visualization/scivis-work/scivis-2019

# cmd="cc -O2 -W -Wall -DTEEM_32BIT=1  -Wl,-prebind \
#   -I$SCIVIS/tpz/include -o conegen conegen.c \
#   -L$/Users/yanchuanqi/UChicago/Spring_2019/vpdata/teem-src/arch/darwi.64 -lteem \
#   -lz -lpthread -lm"
cmd="cc -g -Wall -O2 -o mean_stdv mean_stdv.c -I$SCIVIS/tpz/include -L$SCIVIS/tpz/lib-osx -ltpz -lm"
echo $cmd
eval $cmd
