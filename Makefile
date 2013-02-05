.PHONY: ldlsolver.js
ldlsolver.js: generated.js post.js
	copy /b generated.js+post.js $@
generated.js: ldlsolver.c .\SuiteSparse\lib\libldl.bc .\SuiteSparse\lib\libcamd.bc
	emcc -O2 -Wall -I.\SuiteSparse\include
	-s TOTAL_MEMORY=200000000  -s EXPORTED_FUNCTIONS="['_malloc','_free','_getMatrices','_freeMatrices']" \
	$^ -o $@
test: ldlsolver.c
	gcc -O2 -Wall -std=c99 -I.\SuiteSparse\include -DMAIN $^ -L.\SuiteSparse\lib -lcamd -lldl -o $@