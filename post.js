Module.newFloat64Array = function(size){
	var numBytes = size*Float64Array.BYTES_PER_ELEMENT;
	var ptr = Module._malloc(numBytes);
    return new Float64Array(Module.HEAPF64.buffer, ptr, size);
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