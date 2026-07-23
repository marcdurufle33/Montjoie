#define MONTJOIE_WITH_TWO_DIM
#define MONTJOIE_WITH_THREE_DIM

#define MONTJOIE_WITH_NODAL_H1
#define MONTJOIE_WITH_NODAL_DG

#include "Elliptic/Helmholtz/MontjoieHelmholtz.hxx"

using namespace Montjoie;

Real_wp threshold;

template<class T>
bool isnan(const complex<T>& x)
{
  if (isnan(real(x)) || isnan(imag(x)))
    return true;
  
  return false;
}

template<class T>
void GetRand(T & x)
{
  x = T(rand())/RAND_MAX;
}

template<class T>
void GetRand(complex<T> & x)
{
  int type = rand()%3;
  if (type == 0)
    x = complex<T>(0, rand())/Real_wp(RAND_MAX);
  else if (type == 1)
    x = complex<T>(rand(), 0)/Real_wp(RAND_MAX);
  else
    x = complex<T>(rand(), rand())/Real_wp(RAND_MAX);
}

template<class T>
void GetRandCplx(T & x)
{
  x = T(rand())/RAND_MAX;
}

template<class T>
void GetRandCplx(complex<T> & x)
{
  x = complex<T>(rand(), rand())/Real_wp(RAND_MAX);
}

template<class T>
void GenerateRandomVector(Vector<T>& x, int m)
{
  x.Reallocate(m);
  for (int i = 0; i < m; i++)
    GetRand(x(i));  
}

void ChangeOrder(Vector<string>& lines_data_file, int r)
{
  for (int i = 0; i < lines_data_file.GetM(); i++)
    if (lines_data_file(i).substr(0, 19)  == "OrderDiscretization")
      lines_data_file(i) = "OrderDiscretization = " + to_str(r);
}

void ChangePath(Vector<string>& lines_data_file, const string& path_mesh)
{
  for (int i = 0; i < lines_data_file.GetM(); i++)
    if (lines_data_file(i).substr(0, 8)  == "MeshPath")
      lines_data_file(i) = "MeshPath = " + path_mesh;
}

void ChangeMesh(Vector<string>& lines_data_file, const string& file_mesh)
{
  for (int i = 0; i < lines_data_file.GetM(); i++)
    if (lines_data_file(i).substr(0, 8)  == "FileMesh")
      lines_data_file(i) = "FileMesh = " + file_mesh;
}

template<class TypeEquation>         
void ConstructProblem(EllipticProblem<TypeEquation>& var,
                      const string& input_file, const string& name_element,
		      int r, string file_mesh = string())
{
  var.InitIndices(20);
  Vector<string> lines_data_file;
  if (input_file.compare("NONE"))
    ReadLinesFile(input_file, lines_data_file, var.comm_group_mode);

  if (r >= 0)
    ChangeOrder(lines_data_file, r);
  
  if (file_mesh.size() > 2)
    {
      ChangePath(lines_data_file, "./");
      ChangeMesh(lines_data_file, file_mesh);
    }
  
  ReadInputFile(lines_data_file, var);
  
  bool split_mesh = true;
  if (var.GetNbProcPerMode() == 1)
    split_mesh = false;
  
  var.ComputeMeshAndFiniteElement(name_element, split_mesh);

  var.PerformOtherInitializations();
  
  var.ComputeMassMatrix();
  
  var.SetCurrentModeNumber(0);
  var.ComputeQuasiPeriodicPhase();
}

template<class TypeEquation,
         class MatrixSparse, class MatrixSparseRef>
void CheckGlobalProduct(EllipticProblem<TypeEquation>& var,
                        const MatrixSparse& A, const MatrixSparseRef& Aref)
{
  typedef typename MatrixSparseRef::entry_type T;
  T one, zero;
  SetComplexZero(zero); SetComplexOne(one);
  T alpha, beta;
  GetRandCplx(alpha); GetRandCplx(beta);
  
  if (var.GetNbProcPerMode() == 1)
    {
      int nodl = var.GetNbDof();
      Vector<T> Ones(nodl), Ah_Ones(nodl), B;
      GenerateRandomVector(B, nodl);
      for (int i = 0; i < nodl; i++)
        {
          Ones.Fill(0); Ones(i) = one;
          A.MltVector(Ones, Ah_Ones);
          for (int j = 0; j < nodl; j++)
            if ( (abs(Ah_Ones(j) - Aref(i, j)) > threshold) || (isnan(Ah_Ones(j))))
              {
                cout << "Matrix vector product incorrect (sequential) " << endl;
                DISP(i); DISP(j);
                DISP(Ah_Ones(j)); DISP(Aref(i, j));
                abort();
              }

          // product with transpose
          A.MltAddVector(one, SeldonTrans, Ones, zero, Ah_Ones);
          for (int j = 0; j < nodl; j++)
            if ( (abs(Ah_Ones(j) - Aref(j, i)) > threshold) || (isnan(Ah_Ones(j))))
              {
                cout << "Matrix vector product incorrect (transpose sequential) " << endl;
                DISP(i); DISP(j);
                DISP(Ah_Ones(j)); DISP(Aref(j, i));
                abort();
              }
          
          // product with conjugate transpose
          A.MltAddVector(one, SeldonConjTrans, Ones, zero, Ah_Ones);
          for (int j = 0; j < nodl; j++)
            if ( (abs(Ah_Ones(j) - conjugate(Aref(j, i))) > threshold) || (isnan(Ah_Ones(j))))
              {
                cout << "Matrix vector product incorrect (transpose conjugate sequential) " << endl;
                DISP(i); DISP(j);
                DISP(Ah_Ones(j)); DISP(Aref(j, i));
                abort();
              }
        }

      // testing with scalars and random values
      Copy(B, Ah_Ones);
      GenerateRandomVector(Ones, nodl);
      A.MltAddVector(alpha, Ones, beta, B);
      MltAdd(alpha, Aref, Ones, beta, Ah_Ones);
      for (int i = 0; i < nodl; i++)
        if (abs(Ah_Ones(i) - B(i)) > threshold || (isnan(Ah_Ones(i)))) 
          {
            cout << "Matrix vector product incorrect (random sequential) " << endl;
            DISP(i); DISP(Ah_Ones(i)); DISP(B(i));
            abort();
          }      

      A.MltAddVector(alpha, SeldonTrans, Ones, beta, B);
      MltAdd(alpha, SeldonTrans, Aref, Ones, beta, Ah_Ones);
      for (int i = 0; i < nodl; i++)
        if (abs(Ah_Ones(i) - B(i)) > threshold || (isnan(Ah_Ones(i)))) 
          {
            cout << "Matrix vector product incorrect (transpose random sequential) " << endl;
            DISP(i); DISP(Ah_Ones(i)); DISP(B(i));
            abort();
          }      

      A.MltAddVector(alpha, SeldonConjTrans, Ones, beta, B);
      MltAdd(alpha, SeldonConjTrans, Aref, Ones, beta, Ah_Ones);
      for (int i = 0; i < nodl; i++)
        if (abs(Ah_Ones(i) - B(i)) > threshold || (isnan(Ah_Ones(i)))) 
          {
            cout << "Matrix vector product incorrect (transpose conjugate random sequential) " << endl;
            DISP(i); DISP(Ah_Ones(i)); DISP(B(i));
            abort();
          }      
    }
  else
    {
      int nodl = var.GetNbDof();
      Vector<T> Ones(nodl), Ah_Ones(nodl);
      Ones.Fill(zero); Ah_Ones.Fill(zero);
      for (int i = 0; i < var.GetNbGlobalDof(); i++)
        {
          for (int j = 0; j < nodl; j++)
            {
              if (var.mesh_num.GlobDofNumber_Subdomain(j) == i)
                Ones(j) = one;
              else
                Ones(j) = zero;
            }
          
          A.MltVector(Ones, Ah_Ones);
          for (int j = 0; j < nodl; j++)
            {
              int irow = var.mesh_num.GlobDofNumber_Subdomain(j);
              int icol = i;
              if ((abs(Ah_Ones(j) - Aref(irow, i)) > threshold) || (isnan(Ah_Ones(j))))
                {
                  cout << "Matrix vector product incorrect (parallel) " << endl;
                  DISP(irow); DISP(icol);
                  DISP(Ah_Ones(j)); DISP(Aref(irow, i));
                  abort();
                }
            }
          
          // testing product with transpose
          A.MltAddVector(one, SeldonTrans, Ones, zero, Ah_Ones);
          for (int j = 0; j < nodl; j++)
            {
              int irow = var.mesh_num.GlobDofNumber_Subdomain(j);
              int icol = i;
              if ((abs(Ah_Ones(j) - Aref(i, irow)) > threshold) || (isnan(Ah_Ones(j))))
                {
                  cout << "Matrix vector product incorrect (transpose parallel) " << endl;
                  DISP(irow); DISP(icol);
                  DISP(Ah_Ones(j)); DISP(Aref(i, irow));
                  abort();
                }
            }

          // testing product with conjugate transpose
          A.MltAddVector(one, SeldonConjTrans, Ones, zero, Ah_Ones);
          for (int j = 0; j < nodl; j++)
            {
              int irow = var.mesh_num.GlobDofNumber_Subdomain(j);
              int icol = i;
              if ((abs(Ah_Ones(j) - conjugate(Aref(i, irow))) > threshold) || (isnan(Ah_Ones(j))))
               {
                 cout << "Matrix vector product incorrect (conjugate transpose parallel) " << endl;
                 DISP(irow); DISP(icol);
                 DISP(Ah_Ones(j)); DISP(Aref(i, irow));
                 abort();
               }
            }
        }

      // testing with scalars and random values
      srand(1024); // needed to synchronize random generators in all processors
      GetRandCplx(alpha); GetRandCplx(beta);
      Vector<T> Bref, Xref, B(nodl);
      GenerateRandomVector(Bref, var.GetNbGlobalDof());
      GenerateRandomVector(Xref, var.GetNbGlobalDof());
      for (int i = 0; i < nodl; i++)
        {
          Ones(i) = Xref(var.mesh_num.GlobDofNumber_Subdomain(i));
          B(i) = Bref(var.mesh_num.GlobDofNumber_Subdomain(i));
        }
      
      A.MltAddVector(alpha, Ones, beta, B);
      MltAdd(alpha, Aref, Xref, beta, Bref);
      for (int i = 0; i < nodl; i++)
        if (abs(B(i) - Bref(var.mesh_num.GlobDofNumber_Subdomain(i))) > threshold || (isnan(B(i)))) 
          {
            cout << "Matrix vector product incorrect (random parallel) " << endl;
            DISP(i); DISP(Bref(var.mesh_num.GlobDofNumber_Subdomain(i))); DISP(B(i));
            abort();
          }      

      A.MltAddVector(alpha, SeldonTrans, Ones, beta, B);
      MltAdd(alpha, SeldonTrans, Aref, Xref, beta, Bref);
      for (int i = 0; i < nodl; i++)
        if (abs(B(i) - Bref(var.mesh_num.GlobDofNumber_Subdomain(i))) > threshold || (isnan(B(i)))) 
          {
            cout << "Matrix vector product incorrect (transpose random parallel) " << endl;
            DISP(i); DISP(Bref(var.mesh_num.GlobDofNumber_Subdomain(i))); DISP(B(i));
            abort();
          }      

      A.MltAddVector(alpha, SeldonConjTrans, Ones, beta, B);
      MltAdd(alpha, SeldonConjTrans, Aref, Xref, beta, Bref);
      for (int i = 0; i < nodl; i++)
        if (abs(B(i) - Bref(var.mesh_num.GlobDofNumber_Subdomain(i))) > threshold || (isnan(B(i)))) 
          {
            cout << "Matrix vector product incorrect (transpose conjugate random parallel) " << endl;
            DISP(i); DISP(Bref(var.mesh_num.GlobDofNumber_Subdomain(i))); DISP(B(i));
            abort();
          }      

    }  
}
                        
template<class TypeEquation>
void CheckProductMatrixVector(EllipticProblem<TypeEquation>& var_test,
                              const string& input_file, const string& name_element, int r1, int r2)
{
  typedef typename TypeEquation::Complexe Complexe;
  GlobalGenericMatrix<Complexe> nat_mat;
  nat_mat.SetCoefMass(0.92); nat_mat.SetCoefStiffness(2.4); nat_mat.SetCoefDamping(1.3);
  
  for (int r = r1; r <= r2; r++)
    {
      DISP(r);
      // parallel problem
      EllipticProblem<TypeEquation> var;
      ConstructProblem(var, input_file, name_element, r);
      
      FemMatrixFreeClass<Complexe, TypeEquation> A(var);
      var.AddMatrixWithBC(A, nat_mat);
      
      // sequential problem
      EllipticProblem<TypeEquation> var_only;
      DistributedMatrix<Complexe, General, ArrayRowSparse> Aref;
      var_only.comm_group_mode = MPI_COMM_SELF;
      if (var.GetNbProcPerMode() == 1)
        ConstructProblem(var_only, input_file, name_element, r);
      else
        ConstructProblem(var_only, input_file, name_element, r, "initial.mesh");
      
      var_only.AddMatrixWithBC(Aref, nat_mat);
      // Aref.WriteText("Aref.dat");
      
      CheckGlobalProduct(var, A, Aref);
    }
}

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  
  threshold = 1e4*epsilon_machine;

  {
    // testing continuous finite element and Helmholtz equation
    //EllipticProblem<HelmholtzEquation<Dimension2> > vars;
    
    //CheckProductMatrixVector(vars, "src/Program/Unit/Hyperbolic/carre.ini",
    //		     string ("TRIANGLE_LOBATTO"), 1, 10);
  }

  {
    // testing continuous finite element and Laplace equation
    EllipticProblem<LaplaceEquation<Dimension3> > vars;
    
    CheckProductMatrixVector(vars, "src/Program/Unit/Hyperbolic/cube.ini",
    			     string ("TETRAHEDRON_LOBATTO"), 1, 5);
  }

  {
    // testing continuous finite element and Laplace equation
    EllipticProblem<LaplaceEquation<Dimension2> > vars;
    
    CheckProductMatrixVector(vars, "src/Program/Unit/Hyperbolic/carre.ini",
			     string ("TRIANGLE_LOBATTO"), 1, 10);
  }
  
  {
    // testing SIPG and Helmholtz equation
    //EllipticProblem<HelmholtzEquationSipg<Dimension2> > vars;
    
    //CheckProductMatrixVector(vars, "src/Program/Unit/Hyperbolic/carre.ini", 1, 10);
  }
  
  MPI_Barrier(MPI_COMM_WORLD);
  int rank_proc; MPI_Comm_rank(MPI_COMM_WORLD, &rank_proc);
  if (rank_proc == 0)
    {
      cout << endl << endl;
      cout << "All tests successfully passed" << endl;
    }
  
  return FinalizeMontjoie();
}
