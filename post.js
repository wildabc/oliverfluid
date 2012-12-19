Module.newFloat64Array = function(size){
	var numBytes = size*Float64Array.BYTES_PER_ELEMENT;
	var ptr = Module._malloc(numBytes);
    return new Float64Array(Module.HEAPU8.buffer, ptr, size);
}
Module.freeFloat64Array = function(array){
	if(array !== undefined && 'byteOffset' in array){
		Module._free(array.byteOffset);
	}
}
Module.initSolver = Module.cwrap('initSolver','',['number','number']);
Module.csolve = Module.cwrap('solve','',['number','number']);
Module.solve = function(p, div){
	Module.csolve(p.byteOffset, div.byteOffset);
}
Module.getMatrices = function(){
	var f = Module.cwrap('getMatrices', '', ['number']);
	var ptr = Module._malloc(4*7);
	var a = new Int32Array(Module.HEAPU8.buffer, ptr, 7);
	f(ptr);
	var N = a[6];
	Lp = new Int32Array(Module.HEAPU8.buffer, a[0], N+1);
	Li = new Int32Array(Module.HEAPU8.buffer, a[1], Lp[N]);
	Lx = new Float64Array(Module.HEAPU8.buffer, a[2], Lp[N]);
	D = new Float64Array(Module.HEAPU8.buffer, a[3], N);
	P = new Int32Array(Module.HEAPU8.buffer, a[4], N);
	Pinv = new Int32Array(Module.HEAPU8.buffer, a[5], N);
}
Module.jsolve = function(p, div, width, height){
	var k = 0;
	var l = width+3;
	var N = width*height-1;
	var x = new Float64Array(N);
	for(var j = 1 ; j <= height; j++){
		for(var i = 1 ; i <= width; i++){
			x[Pinv[k++]] = div[l++];
		}
		l += 2;
	}

	for(var j = 0 ; j < N ; j++){
		var imax = Lp[j+1];
		for(var i = Lp[j] ; i < imax; i++){
			x[Li[i]] -= Lx[i]*x[j];
		}
	}

	for(var j = 0 ; j < N ; j++){
		x[j] /= D[j];
	}

	for(var j = N-1 ; j >= 0 ; j--){
		var imax = Lp[j+1];
		for(var i = Lp[j] ; i < imax; i++){
			x[j] -= Lx[i]*x[Li[i]];
		}
	}

	k = 0;
	l = width+3;
	for(var j = 1 ; j <= height; j++){
		for(var i = 1 ; i <= width; i++){
			p[l++] = x[Pinv[k++]];
		}
		l += 2;
	}
	p[width+1+(height+1)*(width+2)] = 0;
}