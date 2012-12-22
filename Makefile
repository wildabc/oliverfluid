.PHONY: ldlsolver.js
ldlsolver.js: generated.js post.js
	copy /b generated.js+post.js $@
generated.js: ldlsolver.c .\SuiteSparse\lib\lib.bc
	emcc -O2 -Wall -I.\SuiteSparse\include \
	-s TOTAL_MEMORY=200000000  -s EXPORTED_FUNCTIONS="['_malloc','_free','_getMatrices','_freeMatrices']" \
	$^ -o $@