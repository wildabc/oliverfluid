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
