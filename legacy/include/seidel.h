//
// Patched Seidel header used by our benchmark CLI.
//
// Upstream palmerc/Seidel ships `seidel.h` with a small static SEGSIZE (200),
// which makes the implementation reject polygons with >200 vertices.
//
// We keep the upstream submodule unmodified and provide this header earlier on
// the include path so that `#include "seidel.h"` resolves here.
//
#ifndef Triangulate_seidel_h
#define Triangulate_seidel_h

#ifndef SEGSIZE
// Must be >= max polygon vertex count used in benchmarks.
// Paper tables use n up to 10,000, so 20,000 gives comfortable slack.
#define SEGSIZE 20000
#endif

/* Functions */
extern int triangulate_polygon(int ncontours, int cntr[], double (*vertices)[2], int (*triangles)[3]);
extern int is_point_inside_polygon(double *);

extern int read_segments(const char *, int *);
extern int initialise(int n);
extern int construct_trapezoids(int);
extern int monotonate_trapezoids(int);
extern int triangulate_monotone_polygons(int, int, int (*)[3]);

#endif


