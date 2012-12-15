ldlsolver.js: ldlsolver.c
	emcc -O2 -s TOTAL_MEMORY=20000000 $^ -o $@
	copy /b $@+post.js $@

