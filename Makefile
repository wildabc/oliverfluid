.PHONY: ldlsolver.js
ldlsolver.js: ldlsolver.c .\SuiteSparse\lib\lib.bc
	emcc -O2 -Wall -I.\SuiteSparse\include \
	-s TOTAL_MEMORY=200000000  -s EXPORTED_FUNCTIONS="['_malloc','_free','_initSolver','_solve']" \
	$^ -o $@
	copy /b $@+post.js $@
