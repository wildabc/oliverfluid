#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "ldl.h"
#include "amd.h"

/* ---------Macros borrowed from ldl demo------------------------------------ */
/* -------------------------------------------------------------------------- */
/* ALLOC_MEMORY: allocate a block of memory */
/* -------------------------------------------------------------------------- */
#define ALLOC_MEMORY(p,type,size) \
p = (type *) calloc ((((size) <= 0) ? 1 : (size)) , sizeof (type)) ; \
if (p == (type *) NULL) \
{ \
    printf ("out of memory\n") ; \
    exit(1) ; \
}
/* -------------------------------------------------------------------------- */
/* FREE_MEMORY: free a block of memory */
/* -------------------------------------------------------------------------- */
#define FREE_MEMORY(p,type) \
if (p != (type *) NULL) \
{ \
    free (p) ; \
    p = (type *) NULL ; \
}

#define IX(i,j) ((i)+(WIDTH)*(j))
#define VALID(i,j) ((i)>=0 && (i)<WIDTH && (j)>=0 && (j)<HEIGHT)

static int WIDTH, HEIGHT, N, SIZE;
static int *Ap, *Ai, *Lp, *Li, *Parent, *Lnz, *Flag, *Pattern, *P, *Pinv;
static double *Ax, *Lx, *D, *Y;

void constructMatrices(int width, int height){
	WIDTH = width;
	HEIGHT = height;
	N = width*height-1;
	SIZE = (width+2)*(height+2);

	clock_t begin = clock();

	int Anz = 5*(width-2)*(height-2)+4*(width-2+height-2)*2+3*3-2;
	ALLOC_MEMORY(Ap, int, N+1);
	ALLOC_MEMORY(Ai, int, Anz);
	ALLOC_MEMORY(Ax, double, Anz);

	// construct matix A
	int p = 0, c, r, d;
	for(int j=0; j<height; ++j){
		for(int i=0; i<width; ++i){
			c = IX(i,j);
			if(c==N){
				break;
			}
			Ap[c] = p;
			d = 0;

			if(VALID(i,j-1)){
				r = IX(i,j-1);
				if(r<N){
					Ai[p] = r;
					Ax[p++] = -1;
				}
				++d;
			}
			if(VALID(i-1,j)){
				r = IX(i-1,j);
				if(r<N){
					Ai[p] = r;
					Ax[p++] = -1;
				}
				++d;
			}
			if(VALID(i+1,j)){
				r = IX(i+1,j);
				if(r<N){
					Ai[p] = r;
					Ax[p++] = -1;
				}
				++d;
			}
			if(VALID(i,j+1)){
				r = IX(i,j+1);
				if(r<N){
					Ai[p] = r;
					Ax[p++] = -1;
				}
				++d;
			}

			Ai[p] = c;
			Ax[p++] = d;
		}
	}
	Ap[N] = Anz;

	ALLOC_MEMORY(Parent, int, N);
	ALLOC_MEMORY(Lnz, int, N);
	ALLOC_MEMORY(Flag, int, N);
	ALLOC_MEMORY(Lp, int, N+1);
	ALLOC_MEMORY(Pattern, int, N);
	ALLOC_MEMORY(P, int, N);
	ALLOC_MEMORY(Pinv, int, N);
	ALLOC_MEMORY(D, double, N);
	ALLOC_MEMORY(Y, double, N);

	// LDL factorization
	amd_order(N, Ap, Ai, P, NULL, NULL);
	ldl_symbolic(N, Ap, Ai, Lp, Parent, Lnz, Flag, P, Pinv);
	ALLOC_MEMORY(Li, int, Lp[N]);
	ALLOC_MEMORY(Lx, double, Lp[N]);
	ldl_numeric(N, Ap, Ai, Ax, Lp, Parent, Lnz, Li, Lx, D, Y, Flag, Pattern, P, Pinv);

	clock_t end = clock();
	printf("%d*%d: NNZ of A:%d; NNZ of L:%d; time of construction:%gs.\n", width, height, Anz, Lp[N]+N, (float)(end-begin)/CLOCKS_PER_SEC);
}

void getMatrices(int width, int height, void **p){
	constructMatrices(width, height);

	p[0] = Lp;
	p[1] = Li;
	p[2] = Lx;
	p[3] = D;
	p[4] = P;
	p[5] = Pinv;
	p[6] = (void *)N;
}

void freeMatrices(){
	// free all
	FREE_MEMORY(Ap, int);
	FREE_MEMORY(Ai, int);
	FREE_MEMORY(Ax, double);
	FREE_MEMORY(Parent, int);
	FREE_MEMORY(Lnz, int);
	FREE_MEMORY(Flag, int);
	FREE_MEMORY(Lp, int);
	FREE_MEMORY(Pattern, int);
	FREE_MEMORY(P, int);
	FREE_MEMORY(Pinv, int);
	FREE_MEMORY(D, double);
	FREE_MEMORY(Y, double);
	FREE_MEMORY(Li, int);
	FREE_MEMORY(Lx, double);
}