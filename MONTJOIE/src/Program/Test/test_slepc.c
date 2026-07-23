#include <slepcnep.h>
#include <stdio.h>

static char help[] = "Usage : ./exe num_matrices coef_fi.dat\n";

void ReadMatrix(char* name, Mat *A)
{
  FILE* fp; char* ligne; size_t len;
  fp = fopen(name, "r");
  while (fgetc(fp) != '\n');
  
  PetscInt m, n, nz;
  fscanf(fp, "%d %d %d", &m, &n, &nz);
  PetscInt* nnz = malloc(m*sizeof(PetscInt));
  for (int i = 0; i < m; i++)
    nnz[i] = 0;

  double val_r, val_i;
  for (int k = 0; k < nz; k++)
    {
      int i, j;
      fscanf(fp, "%d %d %lf %lf\n", &i, &j, &val_r, &val_i);
      nnz[i-1]++;
    }

  fclose(fp);

  fp = fopen(name, "r");
  while (fgetc(fp) != '\n');

  fscanf(fp, "%d %d %d", &m, &n, &nz);
  /* printf("coucou %d %d %d\n", m, n, nz); */
  
  MatCreate(PETSC_COMM_WORLD, A);
  MatSetSizes(*A, m, n, m, n);
  MatSetType(*A, MATSEQAIJ);

  MatSeqAIJSetPreallocation(*A, PETSC_DEFAULT, nnz);
  MatSetUp(*A);
  
  for (int i = 0; i < m; i++)
    {
      int i0, j;
      PetscInt* ind; PetscScalar* data;
      ind = malloc(nnz[i]*sizeof(PetscInt));
      data = malloc(nnz[i]*sizeof(PetscScalar));
      for (int k = 0; k < nnz[i]; k++)
        {
          fscanf(fp, "%d %d %lf %lf\n", &i0, &j, &val_r, &val_i);
          ind[k] = j-1; data[k] = val_r + I * val_i;
        }
      
      MatSetValues(*A, 1, &i, nnz[i], ind, data, INSERT_VALUES);
      free(ind); free(data);
    }
  
  MatAssemblyBegin(*A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*A, MAT_FINAL_ASSEMBLY);

  free(nnz);
}

int main(int argc, char **argv)
{
  if (argc < 3)
    {
      printf("Provide at least three arguments\n");
      printf("Usage ./test.x number_of_matrices file_data.dat\n");
      printf("number_of_matrices : number of matrices \n");
      printf("file_coef : coefficients for the function f_i and matrix\n");

      exit(0);
    }
  
  int nb_mat = atoi(argv[1]);  
  PetscErrorCode ierr;
  
  ierr = SlepcInitialize(&argc, &argv, (char*)0, help);
  if (ierr)
    return ierr;
  
  NEP            nep;             /* nonlinear eigensolver context */
  Mat            F, J, *A;
  FN* f; RG rg; PetscInt nev;
  NEPType        type;
  
  A = malloc(nb_mat*sizeof(Mat));
  f = malloc(nb_mat*sizeof(FN));
  
  ierr = NEPCreate(PETSC_COMM_WORLD, &nep); CHKERRQ(ierr);

  FILE* fp;
  fp = fopen(argv[2], "r");

  char nom[256];
  PetscScalar coef_num[2], coef_denom[3];  PetscInt deg_n;
  coef_num[1] = 0.0;
  coef_denom[1] = 0.0; double nr, ni, omega, dr, di;
  for (int k = 0; k < nb_mat; k++)
    {
      fscanf(fp, "%d %lf %lf %lf %lf %lf %s", &deg_n, &nr, &ni, &omega, &dr, &di, nom);
      ReadMatrix(nom, &A[k]);
      coef_num[0] = nr + I*ni; coef_denom[0] = omega*omega;
      coef_denom[2] = dr + I*di;
      
      ierr = FNCreate(PETSC_COMM_WORLD, &f[k]); CHKERRQ(ierr);
      ierr = FNSetType(f[k], FNRATIONAL); CHKERRQ(ierr);
      ierr = FNRationalSetNumerator(f[k], deg_n+1, coef_num); CHKERRQ(ierr);
      if (omega == 0.0)
        ierr = FNRationalSetDenominator(f[k], 0, coef_denom);
      else
        ierr = FNRationalSetDenominator(f[k], 3, coef_denom);
    }
  
  fclose(fp);

  ierr = NEPSetType(nep, NEPNLEIGS); CHKERRQ(ierr);
  ierr = NEPSetProblemType(nep, NEP_RATIONAL); CHKERRQ(ierr);
  NEPSetWhichEigenpairs(nep, NEP_TARGET_MAGNITUDE);
  NEPRIISetLagPreconditioner(nep, 0);
  ierr = NEPGetRG(nep, &rg);CHKERRQ(ierr);
  ierr = RGSetType(rg, RGINTERVAL);CHKERRQ(ierr);
  ierr = RGIntervalSetEndpoints(rg, 0.1, 10.0, 0.01, 10.0); CHKERRQ(ierr);
  ierr = NEPSetTarget(nep, 0.58+I*0.02); CHKERRQ(ierr);
  /* PetscScalar shifts[1]; shifts[0] = 2.0+I*0.1;
     NEPNLEIGSSetRKShifts(nep, 1, shifts); */

  ierr = NEPSetSplitOperator(nep, nb_mat, A, f, DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);

  ierr = NEPSetFromOptions(nep); CHKERRQ(ierr);
  ierr = NEPSolve(nep); CHKERRQ(ierr);
  ierr = NEPGetType(nep, &type);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, " Solution method: %s\n", type); CHKERRQ(ierr);
  ierr = NEPGetDimensions(nep, &nev, NULL, NULL); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, " Number of requested eigenvalues: %D\n", nev); CHKERRQ(ierr);
  
  ierr = PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_INFO_DETAIL); CHKERRQ(ierr);
  ierr = NEPConvergedReasonView(nep, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  ierr = NEPErrorView(nep, NEP_ERROR_BACKWARD, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  ierr = PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  ierr = NEPDestroy(&nep); CHKERRQ(ierr);
    
  for (int k = 0; k < nb_mat; k++)
    {
      MatDestroy(&A[k]);
      FNDestroy(&f[k]);
    }
  
  free(A);

  ierr = SlepcFinalize();
  return ierr;
}
