#include "Hyperbolic/Acoustic/MontjoieAcoustic.hxx"

#include "FiniteElement/Hexahedron/OptHexahedronLobatto.hxx"
#include "FiniteElement/Hexahedron/OptHexahedronLobatto.cxx"

#ifdef MONTJOIE_WITH_RAJA
#include <RAJA/RAJA.hpp>
#endif

using namespace Montjoie;

#include "mlt_test3D.cc"

typedef Matrix3_3sym MatrixN_Nsym;
typedef Dimension3 DimN;


Vector<Vector<MatrixN_Nsym> > Glob_matrice;

void GenerateCode(const Matrix<double, General, RowSparse>& A, const string& Uh, const string& Vh, const string& file_c)
{
  ofstream file_out(file_c.data());
  file_out.precision(15);
  long* ptr = A.GetPtr();
  int* ind = A.GetInd();
  double* data = A.GetData();
  for (int i = 0; i < A.GetM(); i++)
    {
      file_out<<Vh<<"("<<i<<") = ";
      for (int j = ptr[i]; j < ptr[i+1]; j++)
	{
	  if (data[j] < 0)
	    file_out<<data[j];
	  else
	    file_out<<"+"<<data[j];
	  
	  file_out<<"*"<<Uh<<"("<<ind[j]<<")";
	}
      file_out<<";\n  ";
    }
  file_out.close();
}


template<int p, int i, int j>
class AcousticTripleLoop
{
public :
  template<int m, class Prop,
	   class T0, class T1, class T2>
  static inline void
  Mlt(const TinyMatrix<T0, Prop, m, m>& A,
      const Vector<T1>& x, Vector<T2>& y)
  {
    y((p-1)*m+i-1) += A(i-1, j)*x((p-1)*m+j);
    AcousticTripleLoop<p, i, j-1>::Mlt(A, x, y);
  }
};

template<int p, int i>
class AcousticTripleLoop<p, i, 0>
{
public :
  template<int m, class Prop,
	   class T0, class T1, class T2>
  static inline void
  Mlt(const TinyMatrix<T0, Prop, m, m>& A,
      const Vector<T1>& x, Vector<T2>& y)
  {}
};


template<int p, int i>
class AcousticDoubleLoop
{
public :
  template<int m, class Prop,
	   class T0, class T1, class T2>
  static inline void
  Mlt(const TinyMatrix<T0, Prop, m, m>& A,
      const Vector<T1>& x, Vector<T2>& y)
  {
    y((p-1)*m + i-1) = A(i-1, 0)*x((p-1)*m);
    AcousticTripleLoop<p, i, m-1>::Mlt(A, x, y);
    AcousticDoubleLoop<p, i-1>::Mlt(A, x, y);
  }
  
};

template<int p>
class AcousticDoubleLoop<p, 0>
{
public :
  template<int m, class Prop,
	   class T0, class T1, class T2>
  static inline void
  Mlt(const TinyMatrix<T0, Prop, m, m>& A,
      const Vector<T1>& x, Vector<T2>& y)
  {}
};

template<int p>
class AcousticLoop
{
public :
  template<class T0, class Vector1, int dim_N>
  static inline void MltGeom(Vector<TinyMatrix<T0, Symmetric, dim_N, dim_N> >& mat_dfj,
			     const Vector1& Vh, Vector1& BhVh)
  {
    AcousticDoubleLoop<p, dim_N>::Mlt(mat_dfj(p-1), Vh, BhVh);
    AcousticLoop<p-1>::MltGeom(mat_dfj, Vh, BhVh);
  }
  

  
};

template<>
class AcousticLoop<0>
{
public :
  template<class T0, class Vector1, int dim_N>
  static inline void MltGeom(Vector<TinyMatrix<T0, Symmetric, dim_N, dim_N> >& mat_dfj, const Vector1& Vh, Vector1& BhVh)
  {}
  
};

// B = beta*B + alpha*Kh*x
template<class TypeEquation, int m>
void MltAddTest(Real_wp alpha, int level, HyperbolicProblem<TypeEquation>& var,
		VectReal_wp& X, Real_wp beta, VectReal_wp& B, OptHexahedronLobatto<m>& FaceBasis)
{
  typedef typename TypeEquation::TypeEquationStationary::Dimension Dimension;
  Mesh<Dimension>& mesh = var.var_harmonic.mesh;
  MeshNumbering<Dimension>& mesh_num = var.var_harmonic.GetMeshNumbering(0);
  
  int nb_elt = mesh.GetNbElt();
  //int order = mesh_num.GetOrder();
  int dim_N = Dimension::dim_N;
  const ElementReference<Dimension, 1>& Fb = var.var_harmonic.GetReferenceElementH1(0);
  
  int nb_dof_elt = Fb.GetNbDof();  
  
  Mlt(beta, B);
  
  VectReal_wp Uh(nb_dof_elt);
  VectReal_wp Vh(dim_N*nb_dof_elt);
  VectReal_wp BhVh(dim_N*nb_dof_elt);
  VectReal_wp Prod_Uh(nb_dof_elt);    
  
  // GenerateCode(Fb.rht_loc, "Uh", "Vh", "rht.c");
  //GenerateCode(Fb.rh_loc, "Vh", "Uh", "rh.c");
  //exit(0);
  // DISP(Fb.rht_weight); DISP(Fb.rh_weight); exit(0);

  Real_wp* BhVh_ptr = BhVh.GetData();
  TinyVector<Real_wp, Dimension::dim_N> tmp, vh_loc;
  int num_dof;
  for (int i = 0; i < nb_elt; i++)
   {
     for (int j = 0; j < nb_dof_elt; j++)
       {
         num_dof = mesh_num.Element(i).GetNumberDof(j);
         Uh(j) = X(num_dof);
       }
     
     // produit avec R^T
     ApplyRhTransposeExp(FaceBasis, Uh, Vh); // code en dur (deroulage manuel)
     //FaceBasis.ApplyRhTranspose(Uh, Vh); // C++ deroule les boucles
     //Fb.ApplyRhTranspose(Uh, Vh); // ici on met des boucles (produit matrice vecteur creux)
     
     const Vector<MatrixN_Nsym>& mat_dfj = Glob_matrice(i);
     // multiplication by J_i / omega_k DF_i^{-1} \mu DF_i^{*-1}
     // AcousticLoop<(m+1)*(m+1)*(m+1)>::MltGeom(Glob_matrice(i), Vh, BhVh);
     for (int j = 0; j < nb_dof_elt; j++)
       { 
	 TinyVectorLoop<Dimension::dim_N>::CopyVector(Vh, j, tmp); 
	 
	 //Mlt(Glob_matrice(i)(j), tmp, vh_loc);
	 Mlt(mat_dfj(j), tmp, vh_loc);
	 
	 int offset = 3*j;
	 BhVh_ptr[offset] = vh_loc(0);
	 BhVh_ptr[offset+1] = vh_loc(1);
	 BhVh_ptr[offset+2] = vh_loc(2);
	 //TinyVectorLoop<Dimension::dim_N>::CopyVector(vh_loc, j, BhVh);
       }

     // integration against \nabla \varphi
     ApplyRhExp(FaceBasis, BhVh, Prod_Uh);
     //FaceBasis.ApplyRh(BhVh, Prod_Uh);
     //Fb.ApplyRh(BhVh, Prod_Uh);

     // adding interactions
     for (int j = 0; j < nb_dof_elt; j++)
       {
	 num_dof = mesh_num.Element(i).GetNumberDof(j);
	 B(num_dof) += alpha*Prod_Uh(j);
       }
   }

 
 var.var_harmonic.ImposeNullDirichletCondition(B);
 //int test_input; cout << "waiting" << endl; cin >> test_input;
}


template<class Dimension>
void GenerateColorElement(Mesh<Dimension>& mesh,
                          Vector<Vector<int> >& color_element)
{
  // NeighborVertex(i) : liste des elements autour du sommet i
  Vector<Vector<int> > NeighborVertex(mesh.GetNbVertices());
  for (int i = 0; i < mesh.GetNbElt(); i++)
    {
      int nv = mesh.Element(i).GetNbVertices();
      for (int j = 0; j < nv; j++)
        NeighborVertex(mesh.Element(i).numVertex(j)).PushBack(i);
    }

  // NeighborElt(i) : liste des elements voisins de l'element i
  Vector<Vector<int> > NeighborElt(mesh.GetNbElt());
  for (int i = 0; i < mesh.GetNbElt(); i++)
    {
      int nv = mesh.Element(i).GetNbVertices();
      Vector<int> num;
      for (int j = 0; j < nv; j++)
        {
          int p = mesh.Element(i).numVertex(j);
          for (int k = 0; k < NeighborVertex(p).GetM(); k++)
            if (NeighborVertex(p)(k) != i)
              num.PushBack(NeighborVertex(p)(k));
        }
      
      RemoveDuplicate(num);
      NeighborElt(i) = num;
    }
  
  NeighborVertex.Clear();
  int nb_colors = 8;
  bool valid_coloring = false;
  Vector<int> ColorElt(mesh.GetNbElt());
  Vector<int> NbEltPerColor(nb_colors);
  while (!valid_coloring)
    {
      ColorElt.Fill(-1);      
      NbEltPerColor.Zero();

      // boucle sur les elements
      valid_coloring = true;
      for (int i = 0; i < mesh.GetNbElt(); i++)
        {
          int icc = -1;
          // on cherche une couleur admissible pour l'element i
          int nb_elt_min = 0;
          for (int j = 0; j < nb_colors; j++)
            {
              bool admissible_color = true;
              for (int k = 0; k < NeighborElt(i).GetM(); k++)
                {
                  if (ColorElt(NeighborElt(i)(k)) == j)
                    admissible_color = false;
                }
              
              if (admissible_color)
                {
                  if ((icc == -1) || (NbEltPerColor(j) < nb_elt_min))
                    {
                      icc = j;
                      nb_elt_min = NbEltPerColor(j);
                    }
                }
            }
          
          if (icc == -1)
            {
              cout << "No coloring possible with " << nb_colors << " colors " << endl;
              nb_colors++;
              NbEltPerColor.Reallocate(nb_colors);
              valid_coloring = false;
              break;
            }          
          else
            {
              ColorElt(i) = icc;
              NbEltPerColor(icc)++;
            }
        }      
    }

  cout << "The mesh is colored with " << nb_colors << " colors " << endl;
  
  DISP(NbEltPerColor);
  mesh.WriteColor("color.mesh", ColorElt);

  color_element.Reallocate(nb_colors);
  for (int i = 0; i < nb_colors; i++)
    color_element(i).Reallocate(NbEltPerColor(i));
  
  NbEltPerColor.Zero();
  for (int i = 0; i < mesh.GetNbElt(); i++)
    {
      int icc = ColorElt(i);
      color_element(icc)(NbEltPerColor(icc)) = i;
      NbEltPerColor(icc)++;
    }
  
}

// B = beta*B + alpha*Kh*x
template<class TypeEquation, int m>
void MltAddOpenMP(Real_wp alpha, int level, const Vector<Vector<int> >& color_element,
                  HyperbolicProblem<TypeEquation>& var,
                  VectReal_wp& X, Real_wp beta, VectReal_wp& B, OptHexahedronLobatto<m>& FaceBasis)
{
  typedef typename TypeEquation::TypeEquationStationary::Dimension Dimension;
  //Mesh<Dimension>& mesh = var.var_harmonic.mesh;
  MeshNumbering<Dimension>& mesh_num = var.var_harmonic.GetMeshNumbering(0);
  
  //int nb_elt = mesh.GetNbElt();
  //int order = mesh_num.GetOrder();
  int dim_N = Dimension::dim_N;
  const ElementReference<Dimension, 1>& Fb = var.var_harmonic.GetReferenceElementH1(0);
  
  //  Mlt(beta, B);
#pragma omp parallel for
  for (int i = 0; i < B.GetM(); i++)
    B(i) *= beta;
  
  // GenerateCode(Fb.rht_loc, "Uh", "Vh", "rht.c");
  //GenerateCode(Fb.rh_loc, "Vh", "Uh", "rh.c");
  //exit(0);
  // DISP(Fb.rht_weight); DISP(Fb.rh_weight); exit(0);

  for (int nc = 0; nc < color_element.GetM(); nc++)
    {
      int Ne = color_element(nc).GetM();
      //Vector<int> DofCount(X.GetM());
      //DofCount.Zero();
      
#pragma omp parallel
      {
        int nb_dof_elt = Fb.GetNbDof();  
        VectReal_wp Uh(nb_dof_elt);
        VectReal_wp Vh(dim_N*nb_dof_elt);
        VectReal_wp BhVh(dim_N*nb_dof_elt);
        VectReal_wp Prod_Uh(nb_dof_elt);    
        Real_wp* BhVh_ptr = BhVh.GetData();
        TinyVector<Real_wp, Dimension::dim_N> tmp, vh_loc;
        int num_dof;
        
#pragma omp for
        //For elements with same color
        for (int i0 = 0; i0 < Ne; i0++)
          {
            
            int i = color_element(nc)(i0);
            for (int j = 0; j < nb_dof_elt; j++)
              {
                num_dof = mesh_num.Element(i).GetNumberDof(j);
                Uh(j) = X(num_dof);
                //DofCount(num_dof)++;
              }
            
            // produit avec R^T
            //ApplyRhTransposeOpenMP(FaceBasis, Uh, Vh); // code en dur (deroulage manuel)
            //FaceBasis.ApplyRhTranspose(Uh, Vh); // C++ deroule les boucles
            Fb.ApplyRhTranspose(Uh, Vh); // ici on met des boucles (produit matrice vecteur creux)
            
            const Vector<MatrixN_Nsym>& mat_dfj = Glob_matrice(i);
            // multiplication by J_i / omega_k DF_i^{-1} \mu DF_i^{*-1}
            // AcousticLoop<(m+1)*(m+1)*(m+1)>::MltGeom(Glob_matrice(i), Vh, BhVh);
            for (int j = 0; j < nb_dof_elt; j++)
              { 
                TinyVectorLoop<Dimension::dim_N>::CopyVector(Vh, j, tmp); 
                
                //Mlt(Glob_matrice(i)(j), tmp, vh_loc);
                Mlt(mat_dfj(j), tmp, vh_loc);
                
                int offset = 3*j;
                BhVh_ptr[offset] = vh_loc(0);
                BhVh_ptr[offset+1] = vh_loc(1);
                BhVh_ptr[offset+2] = vh_loc(2);
                //TinyVectorLoop<Dimension::dim_N>::CopyVector(vh_loc, j, BhVh);
              }
            
            // integration against \nabla \varphi
            //ApplyRhOpenMP(FaceBasis, BhVh, Prod_Uh);
            //FaceBasis.ApplyRh(BhVh, Prod_Uh);
            Fb.ApplyRh(BhVh, Prod_Uh);
            
            // adding interactions
            for (int j = 0; j < nb_dof_elt; j++)
              {
                num_dof = mesh_num.Element(i).GetNumberDof(j);
                B(num_dof) += alpha*Prod_Uh(j);
              }
          }
      }
    } 
  
  var.var_harmonic.ImposeNullDirichletCondition(B);
}

template<class TypeEquation, int m>
void MltAddRAJA(Real_wp alpha, int level, const Vector<Vector<int> >& color_element,
                  HyperbolicProblem<TypeEquation>& var,
                  VectReal_wp& X, Real_wp beta, VectReal_wp& B, OptHexahedronLobatto<m>& FaceBasis)
{
  typedef typename TypeEquation::TypeEquationStationary::Dimension Dimension;
  //Mesh<Dimension>& mesh = var.var_harmonic.mesh;
  MeshNumbering<Dimension>& mesh_num = var.var_harmonic.GetMeshNumbering(0);
  
  //int nb_elt = mesh.GetNbElt();
  //int order = mesh_num.GetOrder();
  int dim_N = Dimension::dim_N;
  
  
  //  Mlt(beta, B);
#pragma omp parallel for
  for (int i = 0; i < B.GetM(); i++)
    B(i) *= beta;
  
  //GenerateCode(Fb.rht_loc, "Uh", "Vh", "rht.c");
  //GenerateCode(Fb.rh_loc, "Vh", "Uh", "rh.c");
  //exit(0);
  //DISP(Fb.rht_weight); DISP(Fb.rh_weight); exit(0);
  Real_wp* X_ptr = X.GetData();
  Real_wp* B_ptr = B.GetData();
  for (int nc = 0; nc < color_element.GetM(); nc++)
    {
      int Ne = color_element(nc).GetM();
      //Vector<int> DofCount(X.GetM());
      //DofCount.Zero();
       
      using exec3 = RAJA::omp_parallel_for_exec;                                                                                        RAJA::forall< exec3 >(RAJA::RangeSegment(0, Ne), [&] (int i0) {
          const ElementReference<Dimension, 1>& Fb = var.var_harmonic.GetReferenceElementH1(0);
          int nb_dof_elt = Fb.GetNbDof();  
          VectReal_wp Uh(nb_dof_elt);
          VectReal_wp Vh(dim_N*nb_dof_elt);
          VectReal_wp BhVh(dim_N*nb_dof_elt);
          VectReal_wp Prod_Uh(nb_dof_elt);    
           Real_wp* BhVh_ptr = BhVh.GetData();
          //Real_wp* Vh_ptr = Vh.GetData();
          Real_wp* Uh_ptr = Uh.GetData();
          Real_wp* Prod_Uh_ptr = Prod_Uh.GetData();
          TinyVector<Real_wp, Dimension::dim_N> tmp, vh_loc;
          int num_dof;
          
          int i = color_element(nc)(i0);
          for (int j = 0; j < nb_dof_elt; j++)
            {
              num_dof = mesh_num.Element(i).GetNumberDof(j);
              Uh_ptr[j] = X_ptr[num_dof];
              //DofCount(num_dof)++;
            }
          
          // produit avec R^T
          //ApplyRhTransposeOpenMP(FaceBasis, Uh, Vh); // code en dur (deroulage manuel)
          //FaceBasis.ApplyRhTranspose(Uh, Vh); // C++ deroule les boucles
          Fb.ApplyRhTranspose(Uh, Vh); // ici on met des boucles (produit matrice vecteur creux)
          
          const Vector<MatrixN_Nsym>& mat_dfj = Glob_matrice(i);
          // multiplication by J_i / omega_k DF_i^{-1} \mu DF_i^{*-1}
          // AcousticLoop<(m+1)*(m+1)*(m+1)>::MltGeom(Glob_matrice(i), Vh, BhVh);
          
          for (int j = 0; j < nb_dof_elt; j++)
            { 
              TinyVectorLoop<Dimension::dim_N>::CopyVector(Vh, j, tmp); 
              
              //Mlt(Glob_matrice(i)(j), tmp, vh_loc);
              Mlt(mat_dfj(j), tmp, vh_loc);
              
              int offset = 3*j;
              BhVh_ptr[offset] = vh_loc(0);
              BhVh_ptr[offset+1] = vh_loc(1);
              BhVh_ptr[offset+2] = vh_loc(2);
              //TinyVectorLoop<Dimension::dim_N>::CopyVector(vh_loc, j, BhVh);
            }
          
          // integration against \nabla \varphi
          //ApplyRhOpenMP(FaceBasis, BhVh, Prod_Uh);
          //FaceBasis.ApplyRh(BhVh, Prod_Uh);
          Fb.ApplyRh(BhVh, Prod_Uh);
            
          // adding interactions
          for (int j = 0; j < nb_dof_elt; j++)
            {
              num_dof = mesh_num.Element(i).GetNumberDof(j);
              B_ptr[num_dof] += alpha*Prod_Uh_ptr[j];
            }
        });
    }
  
  var.var_harmonic.ImposeNullDirichletCondition(B);
}

#if defined(RAJA_ENABLE_CUDA)                                                                            
const int CUDA_BLOCK_SIZE = 256;                                                                         
#endif 

template<class TypeEquation, int m>
void MltAddCUDA(  Real_wp alpha, int level, const Vector<Vector<int> >& color_element, 
                  HyperbolicProblem<TypeEquation>& var, VectReal_wp& X, 
                  Real_wp beta, VectReal_wp& B, OptHexahedronLobatto<m>& FaceBasis  )
{
  typedef typename TypeEquation::TypeEquationStationary::Dimension Dimension;
  //Mesh<Dimension>& mesh = var.var_harmonic.mesh;
  MeshNumbering<Dimension>& mesh_num = var.var_harmonic.GetMeshNumbering(0);
  
  //int nb_elt = mesh.GetNbElt();
  //int order = mesh_num.GetOrder();
  int dim_N = Dimension::dim_N;
  
  
  //  Mlt(beta, B);
#pragma omp parallel for
  for (int i = 0; i < B.GetM(); i++)
    B(i) *= beta;
  
  //GenerateCode(Fb.rht_loc, "Uh", "Vh", "rht.c");
  //GenerateCode(Fb.rh_loc, "Vh", "Uh", "rh.c");
  //exit(0);
  //DISP(Fb.rht_weight); DISP(Fb.rh_weight); exit(0);
  Real_wp* X_ptr = X.GetData();
  Real_wp* B_ptr = B.GetData();
  for (int nc = 0; nc < color_element.GetM(); nc++)
    {
      int Ne = color_element(nc).GetM();
      //Vector<int> DofCount(X.GetM());
      //DofCount.Zero();
       
      RAJA::forall<RAJA::cuda_exec<CUDA_BLOCK_SIZE> >(RAJA::RangeSegment(0, Ne), [&] (int i0) {
          const ElementReference<Dimension, 1>& Fb = var.var_harmonic.GetReferenceElementH1(0);
          int nb_dof_elt = Fb.GetNbDof();  
          VectReal_wp Uh(nb_dof_elt);
          VectReal_wp Vh(dim_N*nb_dof_elt);
          VectReal_wp BhVh(dim_N*nb_dof_elt);
          VectReal_wp Prod_Uh(nb_dof_elt);    
          Real_wp* BhVh_ptr = BhVh.GetData();
          Real_wp* Vh_ptr = Vh.GetData();
          Real_wp* Uh_ptr = Uh.GetData();
          Real_wp* Prod_Uh_ptr = Prod_Uh.GetData();
          TinyVector<Real_wp, Dimension::dim_N> tmp, vh_loc;
          int num_dof;
          
          int i = color_element(nc)(i0);
          for (int j = 0; j < nb_dof_elt; j++)
            {
              num_dof = mesh_num.Element(i).GetNumberDof(j);
              Uh_ptr[j] = X_ptr[num_dof];
              //DofCount(num_dof)++;
            }
          
          // produit avec R^T
          //ApplyRhTransposeOpenMP(FaceBasis, Uh, Vh); // code en dur (deroulage manuel)
          //FaceBasis.ApplyRhTranspose(Uh, Vh); // C++ deroule les boucles
          Fb.ApplyRhTranspose(Uh, Vh); // ici on met des boucles (produit matrice vecteur creux)
          
          const Vector<MatrixN_Nsym>& mat_dfj = Glob_matrice(i);
          // multiplication by J_i / omega_k DF_i^{-1} \mu DF_i^{*-1}
          // AcousticLoop<(m+1)*(m+1)*(m+1)>::MltGeom(Glob_matrice(i), Vh, BhVh);
          
          for (int j = 0; j < nb_dof_elt; j++)
            { 
              TinyVectorLoop<Dimension::dim_N>::CopyVector(Vh, j, tmp); 
              
              //Mlt(Glob_matrice(i)(j), tmp, vh_loc);
              Mlt(mat_dfj(j), tmp, vh_loc);
              
              int offset = 3*j;
              BhVh_ptr[offset] = vh_loc(0);
              BhVh_ptr[offset+1] = vh_loc(1);
              BhVh_ptr[offset+2] = vh_loc(2);
              //TinyVectorLoop<Dimension::dim_N>::CopyVector(vh_loc, j, BhVh);
            }
          
          // integration against \nabla \varphi
          //ApplyRhOpenMP(FaceBasis, BhVh, Prod_Uh);
          //FaceBasis.ApplyRh(BhVh, Prod_Uh);
          Fb.ApplyRh(BhVh, Prod_Uh);
            
          // adding interactions
          for (int j = 0; j < nb_dof_elt; j++)
            {
              num_dof = mesh_num.Element(i).GetNumberDof(j);
              B_ptr[num_dof] += alpha*Prod_Uh_ptr[j];
            }
        });
    }
  
  var.var_harmonic.ImposeNullDirichletCondition(B);
}





template<class TypeEquation, class RefElt>
void RunAll(HyperbolicProblem<TypeEquation>& var, string name_element, RefElt& quad)
{  
  typedef typename TypeEquation::TypeEquationStationary::Dimension Dimension;
  Mesh<Dimension>& mesh = var.var_harmonic.mesh;
  MeshNumbering<Dimension>& mesh_num = var.var_harmonic.GetMeshNumbering(0);

  var.var_harmonic.ComputeMeshAndFiniteElement(name_element);
  quad.ConstructFiniteElement(dynamic_cast<const HexahedronLobatto&>(var.var_harmonic.GetReferenceElementH1(0)));
  var.var_harmonic.PerformOtherInitializations();
  var.ComputeRightHandSide();
  var.var_harmonic.ComputeMassMatrix();
  
  var.var_harmonic.ComputeQuasiPeriodicPhase();
  
  //var.var_harmonic.TreatDirichletCondition();
  
  var.ComputeMassMatrix();
  var.InitTimeIterations();
  
  Real_wp t; 
  VectReal_wp U0, U1, Ah_u1;
  int nb_iter = 0;
  
  int N = mesh_num.GetNbDof();
  srand(0);
  U0.Reallocate(N); U0.FillRand(); U0 *= 1e-100;
  U1.Reallocate(N); U1.FillRand(); U1 *= 1e-100;
  //U0.Zero(); U1.Zero();
  Ah_u1.Reallocate(mesh_num.GetNbDof()); Ah_u1.Zero();
  
  int nb_dof = mesh_num.Element(0).GetNbDof();
  // calcul des petites matrices DF_e^{-1} DF_e^{-T}
  //on stocke
  Glob_matrice.Reallocate(mesh.GetNbElt());
  Vector<VectReal_wp> Glob_jacobian(mesh.GetNbElt());
  for (int i = 0; i < mesh.GetNbElt(); i++)
    {
      Glob_matrice(i).Reallocate(nb_dof);
      Glob_jacobian(i).Reallocate(nb_dof);
      const VectReal_wp& omega = var.var_harmonic.WeightsND(i);
      bool affine = mesh.IsElementAffine(i);
      if (affine)
        {
          MatrixN_Nsym A;
          MltTrans(var.var_harmonic.Glob_DFjm1(i)(0),
                   var.var_harmonic.Glob_DFjm1(i)(0), A);

          Mlt(1.0/var.var_harmonic.Glob_jacobian(i)(0), A);
          for (int j = 0; j < nb_dof; j++)
            {
              Glob_matrice(i)(j) = A;
              Glob_jacobian(i)(j) = var.var_harmonic.Glob_jacobian(i)(0)*omega(j);
              Mlt(omega(j), Glob_matrice(i)(j));
            }          
        }
      else
        for (int j = 0; j < nb_dof; j++)
          {
            Glob_jacobian(i)(j) = var.var_harmonic.Glob_jacobian(i)(j);
            MltTrans(var.var_harmonic.Glob_DFjm1(i)(j),
                     var.var_harmonic.Glob_DFjm1(i)(j), Glob_matrice(i)(j));
            
            Mlt(1.0/var.var_harmonic.Glob_jacobian(i)(j)*pow(omega(j), 2), Glob_matrice(i)(j));
          }

      var.var_harmonic.Glob_DFjm1(i).Clear();
    }
  
  t = var.GetInitialTime();
  Real_wp dt = var.GetTimeStep();
  DISP(var.GetTimeStep()); DISP(var.GetFinalTime());
  if (var.GetTimeStep() <= Real_wp(0))
    abort();

  int test_input;
  cout << "Voulez vous ecrire les donnes pour test_montjoie ? (1- Oui)" << endl;
  cin >> test_input; DISP(test_input);
  if (test_input == 1)
    {
      ofstream file_out("matrice.dat"); int nb_elt = mesh.GetNbElt();
      file_out.write(reinterpret_cast<char*>(&N), sizeof(int));
      file_out.write(reinterpret_cast<char*>(&nb_elt), sizeof(int));
      file_out.write(reinterpret_cast<char*>(&nb_dof), sizeof(int));
      for (int i = 0; i < nb_elt; i++)
        for (int j = 0; j < nb_dof; j++)
          {
            int num_dof = mesh_num.Element(i).GetNumberDof(j);
            file_out.write(reinterpret_cast<char*>(&Glob_matrice(i)(j)(0, 0)), sizeof(Real_wp));
            file_out.write(reinterpret_cast<char*>(&Glob_matrice(i)(j)(0, 1)), sizeof(Real_wp));
            file_out.write(reinterpret_cast<char*>(&Glob_matrice(i)(j)(0, 2)), sizeof(Real_wp));
            file_out.write(reinterpret_cast<char*>(&Glob_matrice(i)(j)(1, 1)), sizeof(Real_wp));
            file_out.write(reinterpret_cast<char*>(&Glob_matrice(i)(j)(1, 2)), sizeof(Real_wp));
            file_out.write(reinterpret_cast<char*>(&Glob_matrice(i)(j)(2, 2)), sizeof(Real_wp));
            file_out.write(reinterpret_cast<char*>(&Glob_jacobian(i)(j)), sizeof(Real_wp));
            file_out.write(reinterpret_cast<char*>(&num_dof), sizeof(int));
          }
      
      file_out.close();
    }

  Glob_jacobian.Clear();
  

  const VectReal_wp& invDh = var.Glob_mat_Dh->GetInverseDiagonalDhPlusdtSh();
  const VectReal_wp& Sh = var.Glob_mat_Dh->GetDiagonalSh();
  if (test_input == 1)
    invDh.Write("invDh.dat");
  if (test_input == 1)
    Sh.Write("Sh.dat");
  if ((Sh.GetM() == 0)  || (invDh.GetM() == 0))
    {
      cout << "Mass matrices not computed" <<endl;
      abort();
    }

  if (test_input == 1)
    var.sparse_vector_source.Write("Fn.dat");

  Vector<Vector<int> > color_element;
  GenerateColorElement(mesh, color_element);
  
  cout<<"type 1 for raja and 0 for no raja, 2 for openmp"<<endl;
  int RAJA; cin>>RAJA; DISP(RAJA); DISP(var.GetFinalTime()); DISP(t);

  glob_chrono.Reset(VirtualTimer::ALL);
  glob_chrono.Start(VirtualTimer::ALL);
  glob_chrono.Reset(VirtualTimer::SCHEME);
  glob_chrono.Reset(VirtualTimer::STIFFNESS);
  glob_chrono.Reset(VirtualTimer::FLUX);
  glob_chrono.Reset(VirtualTimer::COMM);
  glob_chrono.Reset(VirtualTimer::OUTPUT);
  
  while (t <= (var.GetFinalTime() + var.epsilon_time))
    {
      glob_chrono.Start(VirtualTimer::OUTPUT);
      if (var.var_harmonic.GetPrintLevel() >= 2)
        var.WriteSnapshot(nb_iter, t, U1);

      glob_chrono.Stop(VirtualTimer::OUTPUT);
      
      glob_chrono.Start(VirtualTimer::FLUX);
      
      Ah_u1.Zero();
      // we compute the source term Fh(n)
      if (t < var.GetFinalTimeSource())
	var.AddPrimitiveSourceAtTime(1.0, t, 1, Ah_u1);

      glob_chrono.Stop(VirtualTimer::FLUX);
      glob_chrono.Start(VirtualTimer::STIFFNESS);
      
      // now we compute Ah_u1 = Fh_n + Ah u_n
      if (RAJA == 2)
        MltAddOpenMP(-dt, -1, color_element, var, U1, dt, Ah_u1, quad);
      else if (RAJA == 1) 
        {
#ifdef RAJA_ENABLE_CUDA
          MltAddCUDA(-dt, -1, color_element, var, U1, dt, Ah_u1, quad);
#else
          MltAddRAJA(-dt, -1, color_element, var, U1, dt, Ah_u1, quad);
#endif
        }
      else
        MltAddTest(-dt, -1, var, U1, dt, Ah_u1, quad);
      
      glob_chrono.Stop(VirtualTimer::STIFFNESS);
      glob_chrono.Start(VirtualTimer::SCHEME);
      

      Real_wp* Ah_u1_ptr = Ah_u1.GetData();
      Real_wp* U0_ptr = U0.GetData();
      Real_wp* U1_ptr = U1.GetData();

  if (RAJA == 1)
    {
      
      // dans Ah_u1, tu peux rajouter Sh (Un + Un-1)
#ifdef MONTJOIE_WITH_RAJA
      RAJA::RangeSegment L2(0, mesh_num.GetNbDof());
      using exec2 = RAJA::omp_parallel_for_exec;
      RAJA::forall<exec2>(L2 , [&] (int i)
                          {
                            Ah_u1_ptr[i] += Sh(i)*(U0(i) - U1(i));
                          });
#else
      for (int i = 0; i < N; i++)
        Ah_u1(i) += Sh(i)*(U0(i) - U1(i));
#endif
    }
  else if (RAJA == 2)
    {
#pragma omp parallel
      {
#pragma omp for
        for (int i = 0; i < N; i++)
          Ah_u1_ptr[i] +=Sh(i)*(U0_ptr[i] - U1_ptr[i]);
      }
    }
  else
    {
      for (int i = 0; i < N; i++)
        Ah_u1(i) += Sh(i)*(U0(i) - U1(i));
    }
  
  glob_chrono.Stop(VirtualTimer::SCHEME);				      
  glob_chrono.Start(VirtualTimer::COMM);
  
  // une fois qu'on a tout rajoute, on assemble
  var.Assemble(Ah_u1);
  
  glob_chrono.Stop(VirtualTimer::COMM);
  glob_chrono.Start(VirtualTimer::SCHEME);
  
  // avancee du schema
  if (RAJA == 1)
    {      
#ifdef MONTJOIE_WITH_RAJA
      RAJA::RangeSegment L3(0, mesh_num.GetNbDof());
      using exec3 = RAJA::omp_parallel_for_exec;
      RAJA::forall<exec3>(L3 , [&] (int i)
        {
      Real_wp unp1_l = 2.0*U1(i) - U0(i) + dt*invDh(i)*Ah_u1(i);
      U0_ptr[i] = U1_ptr[i]; U1_ptr[i] = unp1_l;
    });
#else
  Real_wp mass, unp1;
  for (int i = 0; i < N; i++)
    {
      mass = invDh(i);
      unp1 = 2.0*U1(i) - U0(i) + dt*mass*Ah_u1(i);
      U0(i) = U1(i); U1(i) = unp1;
    }
#endif
    }
  else if (RAJA == 2)
    {
#pragma omp parallel
      {
#pragma omp for
      for (int i = 0; i < N; i++)
        {
      Real_wp unp1_l = 2.0*U1_ptr[i] - U0_ptr[i] + dt*invDh(i)*Ah_u1(i);
      U0_ptr[i] = U1_ptr[i]; U1_ptr[i] = unp1_l;
    }
    }   
      
    }
  else
    {
      Real_wp mass, unp1;
      for (int i = 0; i < N; i++)
        {
      mass = invDh(i);
      unp1 = 2.0*U1(i) - U0(i) + dt*mass*Ah_u1(i);
      U0(i) = U1(i); U1(i) = unp1;
    }
    }
  
      glob_chrono.Stop(VirtualTimer::SCHEME);
      
      nb_iter++;
      t = var.GetInitialTime() + var.GetTimeStep()*nb_iter;
    }

  glob_chrono.Stop(VirtualTimer::ALL);
  glob_chrono.DisplayTime(VirtualTimer::ALL, " to complete time iterations ");
  glob_chrono.DisplayTime(VirtualTimer::OUTPUT, " to write the solution ");
  glob_chrono.DisplayTime(VirtualTimer::FLUX, " to compute source terms ");
  glob_chrono.DisplayTime(VirtualTimer::STIFFNESS, " to compute the product with Kh ");
  glob_chrono.DisplayTime(VirtualTimer::SCHEME, " to advance the time scheme ");
}


template<class TypeEquation>
void RunAll(HyperbolicProblem<TypeEquation>& var, string file_name, string name_element)
{
  var.var_harmonic.InitIndices(PhysicalConstant::nb_max_indices);
  var.var_harmonic.SetTypeEquation("ACOUSTIC");

  ReadInputFile(file_name, var);

  int r = var.var_harmonic.GetMeshNumbering(0).GetOrder();

  switch(r)
    {
    case 1 :
      {
	OptHexahedronLobatto<1> quad;
	RunAll(var, name_element, quad);
	break;
      }
    case 2 :
      {
        OptHexahedronLobatto<2> quad;
	RunAll(var, name_element, quad);
	break;
      }
    case 3 :
      {
	OptHexahedronLobatto<3> quad;
	RunAll(var, name_element, quad);
	break;
        }
    case 4:
      {
	OptHexahedronLobatto<4> quad;
	RunAll(var, name_element, quad);
	break;
      }
      /*
    case 5:
      {
	OptHexahedronLobatto<5> quad;
	RunAll(var, name_element, quad);
	break;
      }
      */
    default:
      cout << "Not implemented" << endl;
      abort();
    }
  
}



int main(int argc, char **argv) 
{  
  InitMontjoie(argc, argv);
  
  if (argc>1)
    {
      // we add the default path to the file name given
      string file_name_data;
      file_name_data = string(argv[1]);
      // DISP(file_name_data);
      
      // we get the type of element selected by the user, and type of equation
      string type_element, type_equation;
      getElement_Equation(file_name_data, type_element, type_equation);
      // cout<<" Equation "<<type_equation<<" Type Element "<<type_element<<endl;

      //HyperbolicProblem<AcousticEquation<Dimension2> > Vars;
      HyperbolicProblem<AcousticEquation<Dimension3> > Vars;

      cout<<" Acoustic Solver with Lobatto quadrangles Qk"<<endl; 
      RunAll(Vars, file_name_data, type_element);
      cout<<" we destroy the variables "<<endl;
      
    }
  else
    {
      cout<<"This code needs a data file in argument"<<endl;
      cout<<"time2D.x nom_fichier"<<endl;
      cout<<"is a good syntax"<<endl;
    }

  return FinalizeMontjoie();
}
