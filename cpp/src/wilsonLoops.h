// double one_x_one(dc U[], int pos[], int mu, int nu, int Nt, int Ns);

int* allocatePaths(int N, int M);
void one_x_one_Paths(int* fPath, int* bPath);
float tracePathWilsonLoop(wc U[], int pos[], int Nt, int Ns, int NPath,
                          int pathLength, int* fPath, int* bPath);

int calculatePathSize(double maxR, double maxT);
void wRT(double maxR, double maxT, int* fPaths, int* bPaths, int pathSize);
