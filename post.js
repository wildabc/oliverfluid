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

Module.getMatrices = function(width,height){
	var f = Module.cwrap('getMatrices', '', ['number']);
	var ptr = Module._malloc(4*7);
	var a = new Int32Array(Module.HEAPU8.buffer, ptr, 7);
	f(ptr);
	var N = a[6];
	Module.matrices = new Object();
	Module.matrices.Lp = new Int32Array(Module.HEAPU8.buffer, a[0], N+1);
	Module.matrices.Li = new Int32Array(Module.HEAPU8.buffer, a[1], Module.matrices.Lp[N]);
	Module.matrices.Lx = new Float64Array(Module.HEAPU8.buffer, a[2], Module.matrices.Lp[N]);
	Module.matrices.D = new Float64Array(Module.HEAPU8.buffer, a[3], N);
	Module.matrices.P = new Int32Array(Module.HEAPU8.buffer, a[4], N);
	Module.matrices.Pinv = new Int32Array(Module.HEAPU8.buffer, a[5], N);

	var x2p = new Int32Array(N);
	var k = 0;
	var l = width+3;
	for(var j = 1 ; j <= height; j++){
		for(var i = 1 ; i <= width; i++){
			x2p[Module.matrices.Pinv[k++]] = l++;
		}
		l += 2;
	}
	Module.matrices.x2p = x2p;
}

Module.jsolve = function(p, div, width, height){
	var Lp = Module.matrices.Lp;
	var Li = Module.matrices.Li;
	var Lx = Module.matrices.Lx;
	var D = Module.matrices.D;
	var x2p = Module.matrices.x2p;

	var N = width*height-1;
	var x = new Float64Array(N);
	for(var k = 0 ; k < N; k++){
			x[k] = div[x2p[k]];
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

	for(var k = 0 ; k < N; k++){
			p[x2p[k]] = x[k];
	}
	p[width+height*(width+2)] = 0;
}