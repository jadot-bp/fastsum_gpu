void one_x_one(double complex U[], int pos[], int mu, int nu, int Nt, int Ns,
               double complex plaq[3][3]);

int* allocatePaths(int N, int M);
void one_x_one_Paths(int* fPath, int* bPath);
double tracePathWilsonLoop(double complex U[], int pos[], int Nt, int Ns,
                           int NPath, int pathLength, int* fPath, int* bPath,
                           double complex plaq[3][3]);

int calculatePathSize(double maxR, double maxT);
void wRT(double maxR, double maxT, int* fPaths, int* bPaths, int pathSize);
