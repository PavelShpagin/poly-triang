// Seidel triangulation baseline CLI (palmerc/Seidel, C implementation).
//
// Input: .poly (n then n lines x y), simple polygon without holes, CCW or CW.
// Output: .tri and a single stdout line:
//   seidel,vertices=n,triangles=t,time_ms=...
//
// NOTE: The upstream implementation uses statically allocated arrays sized by SEGSIZE
// in seidel.h. We compile with a larger SEGSIZE via -DSEGSIZE=... in build.sh.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>

#include "seidel.h"

static double signed_area(const double (*v)[2], int n) {
  double area = 0.0;
  for (int i = 0; i < n; ++i) {
    const int j = (i + 1) % n;
    area += v[i][0] * v[j][1] - v[j][0] * v[i][1];
  }
  return area * 0.5;
}

static void print_usage(const char* prog) {
  fprintf(stderr, "Usage: %s --input <polygon.poly> --output <output.tri>\n", prog);
}

static double now_ms(void) {
  struct timespec ts;
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts);
  return (double)ts.tv_sec * 1000.0 + (double)ts.tv_nsec / 1e6;
}

int main(int argc, char** argv) {
  const char* input = NULL;
  const char* output = NULL;
  for (int i = 1; i < argc; ++i) {
    if ((strcmp(argv[i], "--input") == 0 || strcmp(argv[i], "-i") == 0) && i + 1 < argc) {
      input = argv[++i];
    } else if ((strcmp(argv[i], "--output") == 0 || strcmp(argv[i], "-o") == 0) && i + 1 < argc) {
      output = argv[++i];
    } else {
      print_usage(argv[0]);
      return 1;
    }
  }
  if (!input || !output) {
    print_usage(argv[0]);
    return 1;
  }

  FILE* f = fopen(input, "r");
  if (!f) {
    fprintf(stderr, "Error: cannot open input file: %s\n", input);
    return 1;
  }

  int n = 0;
  if (fscanf(f, "%d", &n) != 1 || n < 3) {
    fclose(f);
    fprintf(stderr, "Error: malformed polygon file (n)\n");
    return 1;
  }

  // Upstream expects vertices[0] unused, indices start at 1.
  double (*vertices)[2] = (double (*)[2])calloc((size_t)(n + 1), sizeof(double[2]));
  if (!vertices) {
    fclose(f);
    fprintf(stderr, "Error: out of memory\n");
    return 1;
  }

  for (int i = 1; i <= n; ++i) {
    if (fscanf(f, "%lf %lf", &vertices[i][0], &vertices[i][1]) != 2) {
      fclose(f);
      free(vertices);
      fprintf(stderr, "Error: malformed polygon file (coords)\n");
      return 1;
    }
  }
  fclose(f);

  // Ensure CCW for outer contour (as required by Seidel code).
  // Compute signed area using 1-based array.
  double (*tmp0)[2] = (double (*)[2])calloc((size_t)n, sizeof(double[2]));
  if (!tmp0) {
    free(vertices);
    fprintf(stderr, "Error: out of memory\n");
    return 1;
  }
  for (int i = 0; i < n; ++i) {
    tmp0[i][0] = vertices[i + 1][0];
    tmp0[i][1] = vertices[i + 1][1];
  }
  if (signed_area(tmp0, n) < 0.0) {
    // Reverse vertices 1..n
    for (int i = 1; i <= n / 2; ++i) {
      const int j = n - i + 1;
      const double x = vertices[i][0];
      const double y = vertices[i][1];
      vertices[i][0] = vertices[j][0];
      vertices[i][1] = vertices[j][1];
      vertices[j][0] = x;
      vertices[j][1] = y;
    }
  }
  free(tmp0);

  int cntr[1];
  cntr[0] = n;

  // Upper bound on triangles for simple polygon without holes: n-2.
  int (*triangles)[3] = (int (*)[3])calloc((size_t)(n * 2), sizeof(int[3]));
  if (!triangles) {
    free(vertices);
    fprintf(stderr, "Error: out of memory\n");
    return 1;
  }

  const double t0 = now_ms();
  const int rc = triangulate_polygon(1, cntr, vertices, triangles);
  const double t1 = now_ms();
  if (rc != 0) {
    free(vertices);
    free(triangles);
    fprintf(stderr, "Error: seidel triangulation failed (rc=%d)\n", rc);
    return 1;
  }

  const int tri_count = n - 2;

  FILE* out = fopen(output, "w");
  if (!out) {
    free(vertices);
    free(triangles);
    fprintf(stderr, "Error: cannot open output file: %s\n", output);
    return 1;
  }

  fprintf(out, "# vertices\n");
  fprintf(out, "%d\n", n);
  for (int i = 1; i <= n; ++i) {
    fprintf(out, "%.10g %.10g\n", vertices[i][0], vertices[i][1]);
  }
  fprintf(out, "# triangles\n");
  fprintf(out, "%d\n", tri_count);
  for (int i = 0; i < tri_count; ++i) {
    // Seidel indices are 1-based.
    fprintf(out, "%d %d %d\n", triangles[i][0] - 1, triangles[i][1] - 1, triangles[i][2] - 1);
  }
  fclose(out);

  free(vertices);
  free(triangles);

  printf("seidel,vertices=%d,triangles=%d,time_ms=%.6f\n", n, tri_count, (t1 - t0));
  return 0;
}

