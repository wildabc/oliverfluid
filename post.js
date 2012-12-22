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

Module.getMatrices = function(width,height){
	var getMatrices = Module.cwrap('getMatrices', '', ['number','number','number']);
	var freeMatrices = Module.cwrap('freeMatrices', '', []);
	var ptr = Module._malloc(4*7);
	var a = new Int32Array(Module.HEAPU8.buffer, ptr, 7);
	getMatrices(width,height,ptr);
	var N = a[6];
	var matrices = new Object();
	var Lp = new Int32Array(Module.HEAPU8.buffer, a[0], N+1);
	var Li = new Int32Array(Module.HEAPU8.buffer, a[1], Lp[N]);
	var Lx = new Float64Array(Module.HEAPU8.buffer, a[2], Lp[N]);
	var D = new Float64Array(Module.HEAPU8.buffer, a[3], N);
	var P = new Int32Array(Module.HEAPU8.buffer, a[4], N);
	var Pinv = new Int32Array(Module.HEAPU8.buffer, a[5], N);
	//copy arrays out
	matrices.Lp = new Int32Array(Lp);
	matrices.Li = new Int32Array(Li);
	matrices.Lx = new Float64Array(Lx);
	matrices.D = new Float64Array(D);
	matrices.P = new Int32Array(P);
	matrices.Pinv = new Int32Array(Pinv);

	var x2p = new Int32Array(N);
	var k = 0;
	var l = width+3;
	for(var j = 1 ; j <= height; j++){
		for(var i = 1 ; i <= width; i++){
			x2p[matrices.Pinv[k++]] = l++;
		}
		l += 2;
	}
	matrices.x2p = x2p;

	// clean up
	Module._free(ptr);
	freeMatrices();

	return matrices;
}

Module.solve = function(p, div, width, height){
	var Lp = selected.matrices.Lp;
	var Li = selected.matrices.Li;
	var Lx = selected.matrices.Lx;
	var D = selected.matrices.D;
	var x2p = selected.matrices.x2p;

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