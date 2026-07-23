#define MONTJOIE_WITH_TWO_DIM
#define MONTJOIE_WITH_THREE_DIM

#define MONTJOIE_WITH_NODAL_H1
#define MONTJOIE_WITH_NODAL_DG
#define MONTJOIE_WITH_ORTHO_DG
#define MONTJOIE_WITH_LEGENDRE_DG

#define MONTJOIE_WITH_NODAL_HCURL
#define MONTJOIE_WITH_HP_HCURL

#define MONTJOIE_WITH_NODAL_HDIV
#define MONTJOIE_WITH_HP_HDIV

#define SELDON_WITH_PRECONDITIONING

#include "Elliptic/Helmholtz/MontjoieLaplace.hxx"
#include "Elliptic/Maxwell/MontjoieMaxwell3D.hxx"

#include "FiniteElement/Pyramid/PyramidHdivOther.cxx"

#include "Elliptic/Helmholtz/HelmholtzHdiv.hxx"
#include "Elliptic/Helmholtz/HelmholtzHdiv.cxx"
#include "Hyperbolic/Aeroacoustic/MontjoieTimeGalbrun.hxx"

using namespace Montjoie;

int rank_proc(0), nb_proc(1);

template<class Dimension>
class InputDispersionVariables : public InputDataProblem_Base
{
  public :
    typedef typename Dimension::R_N R_N;
    //! kmax : nombre d'onde utilisee pour la sortie sur la sphere
    Real_wp kmax, threshold_rational, kh1_order, kh2_order;
    //! order de discretisation
    int order;
    //! ordre d'integration
    int order_integration;
    int angle_deg, print_level;
    string name_file;
    //! maillage regulier ou pas ?
    bool regular_mesh;
    bool  mesh_created;
    int  irregular_mesh;  

    //! base vectorielle du motif periodique 
    Vector<R_N> vdirector_per;
    R_N k1_plot, k2_plot; // wave vectors for the two ends of the plot
    Real_wp alpha_penal, delta_penal;

    //! solveur pour les valeurs propres
    int type_solver_eigenval;
    enum {LAPACK, ITERATIVE_POWER, ARNOLDI};

    int nb_examples, order_dispersion;
    VectReal_wp Xsi1Points, Xsi2Points, Xsi3Points;
    bool search_taylor_dvt, find_order;
    string file_name_plot; // file where plot is stored
    int nb_points_plot; // number of points for the plot
    int type_hybrid_mesh;
    int nb_div_for_moyenne;
    bool display_constant_sphere;
    string file_sphere;
    string file_sphere_output;
    bool courbe_log;
    Real_wp log_min; //!< valeur minimal du log(h) pour la courbe de dispersion
    bool use_warburton_trick;
    int nb_eigenvalues_wanted;

    InputDispersionVariables()
    {
      kmax = 0;
      threshold_rational = 1e-10;
      kh1_order = 0;
      kh2_order = 0;
      order = 1;
      order_integration = 3;
      angle_deg = 90;
      print_level = 1;
      name_file = "toto.mesh";
      regular_mesh = true;
      mesh_created = false;

      nb_examples = 0;
      order_dispersion = 0;
      search_taylor_dvt = false;
      find_order = false;
      file_name_plot = "plot.dat"; // file where plot is stored
      nb_points_plot = 0; // number of points for the plot

      type_hybrid_mesh = 0;
      nb_div_for_moyenne = 3;
      display_constant_sphere = false;
      file_sphere = string("sphere.mesh");
      file_sphere_output = string("sphere_out.mesh");
      courbe_log = false;
      log_min = 0;

      alpha_penal = 0; delta_penal = 0;
      type_solver_eigenval = LAPACK;
      nb_eigenvalues_wanted = 4;
    }


    void SetInputData(const string& description_field, const VectString& parameters)
    {
      if (!description_field.compare("OrderDiscretization"))
      {
        to_num(parameters(0), order);
      }
      else if (!description_field.compare("TypeMesh"))
      {

        if (Dimension::dim_N == 2)
        {
          type_hybrid_mesh = Mesh<Dimension2>::TRIANGULAR_MESH;
          if (!parameters(0).compare("TRI"))
            type_hybrid_mesh = Mesh<Dimension2>::TRIANGULAR_MESH;
          else if (!parameters(0).compare("QUAD"))
            type_hybrid_mesh = Mesh<Dimension2>::QUADRILATERAL_MESH;
          else if (!parameters(0).compare("RADAU"))
            type_hybrid_mesh = Mesh<Dimension2>::RADAU_MESH;
        }
        else
        {
          type_hybrid_mesh = Mesh<Dimension3>::TETRAHEDRAL_MESH;
          if (!parameters(0).compare("TETRA"))
            type_hybrid_mesh = Mesh<Dimension3>::TETRAHEDRAL_MESH;
          else if (!parameters(0).compare("HEXA"))
            type_hybrid_mesh = Mesh<Dimension3>::HEXAHEDRAL_MESH;
          else if (!parameters(0).compare("PYRAMID"))
            type_hybrid_mesh = Mesh<Dimension3>::PYRAMID_MESH;
          else if (!parameters(0).compare("HYBRID"))
            type_hybrid_mesh = Mesh<Dimension3>::HYBRID_MESH;
          else if (!parameters(0).compare("WEDGE"))
            type_hybrid_mesh = Mesh<Dimension3>::WEDGE_MESH;
          else if (!parameters(0).compare("CRUSH"))
            type_hybrid_mesh = Mesh<Dimension3>::CRUSH_MESH;
        }        
      }
      else if (!description_field.compare("NbDivEvalMoyenne"))
      {
        nb_div_for_moyenne = to_num<int>(parameters(0));
      } 
      else if (!description_field.compare("DisplayMesh"))
      {
        if (!parameters(0).compare("YES"))
        {
          display_constant_sphere = true;
          file_sphere = parameters(1);
          file_sphere_output = parameters(2);
        }
      } 
      else if (!description_field.compare("DisplayPlot"))
      {
        int nb = 0;
        for (int k = 0; k < Dimension::dim_N; k++)
          to_num(parameters(nb++), k1_plot(k));

        for (int k = 0; k < Dimension::dim_N; k++)
          to_num(parameters(nb++), k2_plot(k));

        if (!parameters(nb).compare("LOG"))
        {
          courbe_log = true;
          nb++;
          log_min = to_num<Real_wp>(parameters(nb++));
        }

        file_name_plot = parameters(nb++);
        to_num(parameters(nb++), nb_points_plot);
      }
      else if (!description_field.compare("OrderIntegration"))
      {
        to_num(parameters(0), order_integration);
      } 
      else if (!description_field.compare("FileMesh"))
      {
        regular_mesh = false;
        if (!parameters(0).compare("REGULAR"))
        {
          regular_mesh = true;
          if (parameters.GetM() > 1)
          {
            if (!parameters(1).compare("EQUILATERAL"))
              angle_deg = 60;
            else
              to_num(parameters(1), angle_deg);
          }
        }
        else
        {
          name_file = parameters(0);
#ifdef MONTJOIE_WITH_TWO_DIM
          vdirector_per.Reallocate(2);
          to_num(parameters(1), vdirector_per(0)(0));
          to_num(parameters(2), vdirector_per(0)(1));
          to_num(parameters(3), vdirector_per(1)(0));
          to_num(parameters(4), vdirector_per(1)(1));
#else
          vdirector_per.Reallocate(3);
          to_num(parameters(1), vdirector_per(0)(0));
          to_num(parameters(2), vdirector_per(0)(1));
          to_num(parameters(3), vdirector_per(0)(2));
          to_num(parameters(4), vdirector_per(1)(0));
          to_num(parameters(5), vdirector_per(1)(1));
          to_num(parameters(6), vdirector_per(1)(2));
          to_num(parameters(7), vdirector_per(2)(0));
          to_num(parameters(8), vdirector_per(2)(1));
          to_num(parameters(9), vdirector_per(2)(2));
#endif
        }
      }
      else if (!description_field.compare("IrregularMesh"))
      {
        if (!parameters(0).compare("YES"))
          irregular_mesh = 1;
        else if (!parameters(0).compare("LIGHT"))
          irregular_mesh = 2;
        else
          irregular_mesh = 0;
      }
      else if (!description_field.compare("MaxWaveNumber"))
      {
        to_num(parameters(0), kmax);
      }
      else if (!description_field.compare("ThresholdRational"))
      {
        to_num(parameters(0), threshold_rational);
        // DISP(threshold_rational);
      }
      else if (!description_field.compare("NbExamples"))
      {
        to_num(parameters(0), nb_examples);
      }
      else if (!description_field.compare("Xsi1Examples"))
      {
        Xsi1Points.Reallocate(nb_examples);
        for (int i = 0; i < nb_examples; i++)
          to_num(parameters(i), Xsi1Points(i));
      }
      else if (!description_field.compare("Xsi2Examples"))
      {
        Xsi2Points.Reallocate(nb_examples);
        for (int i = 0; i < nb_examples; i++)
          to_num(parameters(i), Xsi2Points(i));
      }
      else if (!description_field.compare("Xsi3Examples"))
      {
        Xsi3Points.Reallocate(nb_examples);
        for (int i = 0; i < nb_examples; i++)
          to_num(parameters(i), Xsi3Points(i));
      }
      else if (!description_field.compare("TaylorDevelopment"))
      {
        if (!parameters(0).compare("YES"))
          search_taylor_dvt = true;
        else
          search_taylor_dvt = false;
      }
      else if (!description_field.compare("PrintLevel"))
      {
        to_num(parameters(0), print_level);
      }
      else if (!description_field.compare("FindOrder"))
      {
        if (!parameters(0).compare("YES"))
        {
          find_order = true;
          to_num(parameters(1), kh1_order);
          to_num(parameters(2), kh2_order);
        }
        else
          find_order = false;
      }
      else if (!description_field.compare("UseWarburtonTrick"))
      {
        use_warburton_trick = false;
        if (!parameters(0).compare("YES"))
          use_warburton_trick = true;
      }
      else if (!description_field.compare("PenalizationDG"))
      {
        alpha_penal = to_num<Real_wp>(parameters(0));
        delta_penal = to_num<Real_wp>(parameters(1));
      }
      else if (!description_field.compare("TypeSolver"))
      {
        if (!parameters(0).compare("POWER"))
          type_solver_eigenval = ITERATIVE_POWER;
        else if (!parameters(0).compare("ARNOLDI"))
        {
          type_solver_eigenval = ARNOLDI;
          if (parameters.GetM() > 1)
            nb_eigenvalues_wanted = to_num<int>(parameters(1));
          else
            nb_eigenvalues_wanted = 5;
        }
      }
    }

};


// returns true if ptA is the translated of ptB
// INPUT
// ptA , ptB : the two vertices
// OUTPUT 
// n1 : the components of the translation
//    : ptA is translated from ptB by the vector n1(0)*(1, 0, 0) + n1(1)*(0, 1, 0) + ...
  template<class VectorInt, class Rn>
bool VertexTranslated(const Rn& ptA, const Rn& ptB, VectorInt& n1)
{
  Rn w;

  // w  = ptA - ptB
  w = ptA - ptB;

  // now, the question is :
  // components of w are integers ?
  bool translated = true;
  for (int k = 0; k < w.GetM(); k++)
  {
    Real_wp w0 = w(k) - round(w(k));
    if (abs(w0) <= Rn::threshold)
      n1(k) = toInteger(round(w(k)));
    else
      translated = false;
  }

  return translated;
}


  template<class T, int m>
bool PointInsideUnitCube(const TinyVector<T, m>& point)
{
  bool inside = true;
  for (int k = 0; k < m; k++)
    if ((point(k) < -TinyVector<T, m>::threshold)
        ||(point(k) > T(1)-TinyVector<T, m>::threshold))
      inside = false;

  return inside;
}

#ifdef MONTJOIE_WITH_TWO_DIM
// motif de base en 2-D (quadrilatere)
  template<class Vector1>
void ConstructBaseMotif(const Vector1& vdirector_per, Mesh<Dimension2>& mesh_elt)
{
  mesh_elt.ReallocateVertices(4);
  Real_wp zero(0); R2 V0(zero,zero);
  mesh_elt.Vertex(0) = V0;
  mesh_elt.Vertex(1) = vdirector_per(0);
  mesh_elt.Vertex(2) = vdirector_per(0);
  Add(Real_wp(1), vdirector_per(1), mesh_elt.Vertex(2));
  mesh_elt.Vertex(3) = vdirector_per(1);

  mesh_elt.ReallocateElements(1);
  mesh_elt.Element(0).InitQuadrangular(0, 1, 2, 3, 1);

  mesh_elt.ReorientElements();
  mesh_elt.FindConnectivity();
}
#endif

#ifdef MONTJOIE_WITH_THREE_DIM
// motif de base en 3-D (hexahedre)
  template<class Vector1>
void ConstructBaseMotif(const Vector1& vdirector_per, Mesh<Dimension3>& mesh_elt)
{
  mesh_elt.ReallocateVertices(8);
  Real_wp zero(0); R3 V0(zero,zero,zero);
  mesh_elt.Vertex(0) = V0;
  mesh_elt.Vertex(1) = vdirector_per(2);
  mesh_elt.Vertex(3) = vdirector_per(1);
  mesh_elt.Vertex(2) = vdirector_per(1);
  Add(Real_wp(1), vdirector_per(2), mesh_elt.Vertex(2));
  mesh_elt.Vertex(4) = vdirector_per(0);
  mesh_elt.Vertex(5) = vdirector_per(0);
  Add(Real_wp(1), vdirector_per(2), mesh_elt.Vertex(5));
  mesh_elt.Vertex(7) = vdirector_per(0);
  Add(Real_wp(1), vdirector_per(1), mesh_elt.Vertex(7));
  mesh_elt.Vertex(6) = vdirector_per(0);
  Add(Real_wp(1), vdirector_per(1), mesh_elt.Vertex(6));
  Add(Real_wp(1), vdirector_per(2), mesh_elt.Vertex(6));

  mesh_elt.ReallocateElements(1);
  mesh_elt.Element(0).InitHexahedral(0, 1, 2, 3, 4, 5, 6, 7, 1);

  mesh_elt.ReorientElements();
  mesh_elt.FindConnectivity();
}
#endif


#ifdef MONTJOIE_WITH_TWO_DIM
// construction of mesh
  template<class Vector1>
void ConstructMesh(Mesh<Dimension2>& mesh, int type_mesh, bool regular_mesh,
    const Real_wp& angle_deg, const string& name_file, 
    int irregular_mesh, Vector1& vdirector_per)
{
  if (regular_mesh)
  {
    // regular mesh
    Real_wp teta = Real_wp(angle_deg)/Real_wp(180)*pi_wp;
    Real_wp sin_teta = sin(teta), cos_teta = cos(teta);

    // two main directions of periodicity in 2-D
    vdirector_per.Reallocate(2); 
    vdirector_per(0)(0) = Real_wp(1);  vdirector_per(0)(1) = Real_wp(0);
    vdirector_per(1)(0) = cos_teta; vdirector_per(1)(1) = sin_teta;

    TinyVector<int, 4> ref_boundary; ref_boundary.Fill(1);
    mesh.CreateRegularMesh(R2(-1, -1), R2(2, 2), TinyVector<int, 2>(4, 4), 1,
        ref_boundary, type_mesh);

    // bending the mesh 
    for (int i = 0; i < mesh.GetNbVertices(); i++)
    {
      mesh.Vertex(i)(0) += cos_teta*mesh.Vertex(i)(1);
      mesh.Vertex(i)(1) *= sin_teta;
    }
  }
  else
  {
    // mesh is read in a file
    mesh.Read(name_file);
  }

  if (irregular_mesh == 1)
  {
    mesh.SplitIntoTriangles();
    mesh.SplitIntoQuadrilaterals();
  }

  if (rank_proc == 0)
    mesh.Write("toto.mesh");
}
#endif

#ifdef MONTJOIE_WITH_THREE_DIM
// construction of mesh
// construction of mesh
  template<class Vector1>
void ConstructMesh(Mesh<Dimension3>& mesh, int type_mesh, bool regular_mesh,
    const Real_wp& angle_deg, const string& name_file, 
    int irregular_mesh, Vector1& vdirector_per)
{
  if (regular_mesh)
  {
    // regular mesh

    // three main directions of periodicity in 3-D
    vdirector_per.Reallocate(3); 
    vdirector_per(0)(0) = Real_wp(1); vdirector_per(0)(1) = Real_wp(0); vdirector_per(0)(2) = Real_wp(0);
    vdirector_per(1)(0) = Real_wp(0); vdirector_per(1)(1) = Real_wp(1); vdirector_per(1)(2) = Real_wp(0);
    vdirector_per(2)(0) = Real_wp(0); vdirector_per(2)(1) = Real_wp(0); vdirector_per(2)(2) = Real_wp(1);

    TinyVector<int, 6> ref_boundary; ref_boundary.Fill(1);
    if (irregular_mesh == 1)
      type_mesh = Mesh<Dimension3>::TETRAHEDRAL_MESH;

    mesh.CreateRegularMesh(R3(-1, -1, -1), R3(2, 2, 2), TinyVector<int, 3>(4, 4, 4), 1,
        ref_boundary, type_mesh);
  }
  else
  {
    // mesh is read in a file
    mesh.Read(name_file);
  }

  if (irregular_mesh == 1)
    mesh.SplitIntoHexahedra();

  // irregular mesh moving the central point
  if (irregular_mesh == 2)
  {
    Vector<string> param(1); param(0) = string("YES");
    mesh.SetInputData(string("SlightModificationOnRegularMesh"), param);

    // three main directions of periodicity in 3-D
    vdirector_per.Reallocate(3); 
    vdirector_per(0)(0) = Real_wp(1); vdirector_per(0)(1) = Real_wp(0); vdirector_per(0)(2) = Real_wp(0);
    vdirector_per(1)(0) = Real_wp(0); vdirector_per(1)(1) = Real_wp(1); vdirector_per(1)(2) = Real_wp(0);
    vdirector_per(2)(0) = Real_wp(0); vdirector_per(2)(1) = Real_wp(0); vdirector_per(2)(2) = Real_wp(1);

    TinyVector<int, 6> ref_boundary; ref_boundary.Fill(1);
    //mesh.CreateRegularMesh(R3(0, 0, 0), R3(1, 1, 1), 3, 1,
    mesh.CreateRegularMesh(R3(-1, -1, -1), R3(2, 2, 2), TinyVector<int, 3>(7, 7, 7), 1,				 
        ref_boundary, type_mesh);

  }

  // sorting faces
  mesh.ReallocateBoundariesRef(mesh.GetNbBoundary());
  for (int i = 0; i < mesh.GetNbBoundary(); i++)
  {
    int nv = mesh.Boundary(i).GetNbVertices();
    IVect num(nv);
    for (int k = 0; k < nv; k++)
      num(k) = mesh.Boundary(i).numVertex(k);

    if (num(nv-1) < num(1))
    {
      // swapping numbers
      int itmp = num(nv-1);
      num(nv-1) = num(1);
      num(1) = itmp;
    }

    mesh.BoundaryRef(i).Init(num, mesh.Boundary(i).GetReference());
  }

  mesh.FindConnectivity();

  if (rank_proc == 0)
    mesh.Write("toto.mesh");
}
#endif


//! class to find dispersion relations
/*!
  Here we denote "original" dofs/elements/edges, the dofs/elements/edges 
  which belong to the original pattern (e.g. the square [0, 1]^2)
  All the other dofs are called periodic dofs.
  When computing the dispersion relation, we will obtain an eigenvalue 
  problem concerning only original dofs
  */
template<class Dimension>
class DispersionSolver
{
  public :
    typedef typename Dimension::R_N R_N;

    //! considered problem
    DistributedProblem<Dimension>& var;

    //! Input data
    InputDispersionVariables<Dimension> data;

    //! Relation between periodic dofs and originals dofs
    /*!
      PeriodicityDof(i, 0) = number of the original dof obtained by translation of periodic dof i
      PeriodicityDof(i, 1:dim_N) = translation vector between the two dofs
      */
    Matrix<int> PeriodicityDof;

    //! list of original dof numbers
    IVect Num_OriginalDofs;

    //! true is the dof is periodic
    VectBool IsDofPeriodic;

    //! original pattern
    Mesh<Dimension> mesh_elt;

    //! phase de quasi-periodicite
    VectComplex_wp phase;

    //! matrices apres periodisation
#ifdef SELDON_WITH_MPI
    DistributedMatrix<Complex_wp, General, ArrayRowSparse> Kh_per, Dh_per, Sh_per;
#else
    Matrix<Complex_wp, General, ArrayRowSparse> Kh_per, Dh_per, Sh_per;
#endif
    //EigenProblemSparse<Complex_wp, Matrix<Complex_wp, General, ArrayRowSparse>,
    //MatrixMumps<Complex_wp> > var_eig;

    typename Dimension::MatrixN_Nsym Aperiodicity, InvAper;

    Real_wp mean_edge, volume_elt;
    R_N flow;
    bool is_galbrun_formulation;

    template<class TypeEquation>
      DispersionSolver(EllipticProblem<TypeEquation>& var_) : var(var_)
    {
      mean_edge = 1;
      volume_elt = 1;
      is_galbrun_formulation = false;
    }

    void TreatPeriodicityVertices(const Mesh<Dimension>& mesh, int nb_dof_vertex)
    {
      Vector<R_N> CoorVertex(mesh.GetNbVertices());
      R_N rhs, coor;
      TinyVector<int, 3> n1;

      IVect NumOriginal_Vertices(mesh.GetNbVertices());
      int nb_original_vertices = 0;
      Vector<bool> IsVertexPeriodic(mesh.GetNbVertices());
      IsVertexPeriodic.Fill(false);
      // finding all original vertices
      for (int i = 0; i < mesh.GetNbVertices(); i++)
      {
        // expressing vertex in the base of vdirector_per
        for (int j = 0; j < Dimension::dim_N; j++)
          rhs(j) = DotProd(data.vdirector_per(j), mesh.Vertex(i));

        Mlt(InvAper, rhs, coor);

        // storing this result
        CoorVertex(i) = coor;

        if (PointInsideUnitCube(coor))
        {
          NumOriginal_Vertices(nb_original_vertices++) = i;
          for (int k = 0; k < nb_dof_vertex; k++)
            IsDofPeriodic(i*nb_dof_vertex + k) = false;
        }
        else
          IsVertexPeriodic(i) = true;
      }

      // then finding translation of periodic vertices 
      for (int i = 0; i < mesh.GetNbVertices(); i++)
        if (IsVertexPeriodic(i))
        {
          for (int j = 0; j < nb_original_vertices; j++)
          {
            if (VertexTranslated(CoorVertex(i), CoorVertex(NumOriginal_Vertices(j)), n1))
            {
              for (int k = 0; k < nb_dof_vertex; k++)
              {
                int ni = i*nb_dof_vertex + k;
                int nj = NumOriginal_Vertices(j)*nb_dof_vertex + k;
                PeriodicityDof(ni, 0) = nj;
                for (int k = 0; k < Dimension::dim_N; k++)
                  PeriodicityDof(ni, k+1) = n1(k);
              }
            }
          }
        }
    }

    void TreatPeriodicityEdges(const Mesh<Dimension>& mesh, int nb_dof_vertex, int nb_dof_edge)
    {
      Vector<R_N> CoorEdge(mesh.GetNbEdges());
      R_N rhs, coor, center;
      TinyVector<int, 3> n1;

      IVect NumOriginal_Edges(mesh.GetNbEdges());
      int nb_original_edges = 0;
      Vector<bool> IsEdgePeriodic(mesh.GetNbEdges());
      IsEdgePeriodic.Fill(false);
      int offset = mesh.GetNbVertices()*nb_dof_vertex;
      // finding all original edges
      for (int i = 0; i < mesh.GetNbEdges(); i++)
      {
        int ne1 = mesh.GetEdge(i).numVertex(0);
        int ne2 = mesh.GetEdge(i).numVertex(1);
        center = Real_wp(0.5)*(mesh.Vertex(ne1) + mesh.Vertex(ne2));

        // expressing center in the base of vdirector_per
        for (int j = 0; j < Dimension::dim_N; j++)
          rhs(j) = DotProd(data.vdirector_per(j), center);

        Mlt(InvAper, rhs, coor);

        // storing this result
        CoorEdge(i) = coor;

        if (PointInsideUnitCube(coor))
        {
          NumOriginal_Edges(nb_original_edges++) = i;
          for (int k = 0; k < nb_dof_edge; k++)
            IsDofPeriodic(offset + i*nb_dof_edge + k) = false;
        }
        else
          IsEdgePeriodic(i) = true;
      }

      // then finding translation of periodic edges
      for (int i = 0; i < mesh.GetNbEdges(); i++)
        if (IsEdgePeriodic(i))
        {
          for (int j = 0; j < nb_original_edges; j++)
          {
            if (VertexTranslated(CoorEdge(i), CoorEdge(NumOriginal_Edges(j)), n1))
            {
              for (int k = 0; k < nb_dof_edge; k++)
              {
                int ni = offset + i*nb_dof_edge + k;
                int nj = offset + NumOriginal_Edges(j)*nb_dof_edge + k;
                PeriodicityDof(ni, 0) = nj;
                for (int k = 0; k < Dimension::dim_N; k++)
                  PeriodicityDof(ni, k+1) = n1(k);
              }
            }
          }
        }
    }


    void TreatPeriodicityFaces(Mesh<Dimension2>& mesh, int nb_dof_vertex, int nb_dof_edge,
        int nb_dof_tri, int nb_dof_quad, int& offset)
    {
      offset = nb_dof_vertex*mesh.GetNbVertices() + nb_dof_edge*mesh.GetNbEdges();
    }


    void TreatPeriodicityFaces(Mesh<Dimension3>& mesh, int nb_dof_vertex, int nb_dof_edge,
        int nb_dof_tri, int nb_dof_quad, int& offset)
    {    
      R_N rhs, coor, center;
      TinyVector<int, 3> n1;

      // then finding original and periodic dofs associated with faces    
      int nb_tri = mesh.GetNbTriangles();
      int nb_quad = mesh.GetNbQuadrangles();
      offset = nb_dof_vertex*mesh.GetNbVertices() + nb_dof_edge*mesh.GetNbEdges();
      if (((nb_tri > 0)&&(nb_dof_tri > 0)) || ((nb_quad > 0)&&(nb_dof_quad > 0)) )
      {
        Vector<R_N> CoorFace(mesh.GetNbFaces());

        IVect NumOriginal_Faces(mesh.GetNbFaces());
        int nb_original_faces = 0;
        Vector<bool> IsFacePeriodic(mesh.GetNbFaces());
        IsFacePeriodic.Fill(false);

        IVect offset_face(mesh.GetNbBoundary()+1);
        offset_face(0) = offset;
        // finding all original faces
        for (int i = 0; i < mesh.GetNbBoundary(); i++)
        {
          int ne1 = mesh.Boundary(i).numVertex(0);
          int ne2 = mesh.Boundary(i).numVertex(1);
          int ne3 = mesh.Boundary(i).numVertex(2);
          center = mesh.Vertex(ne1) + mesh.Vertex(ne2) + mesh.Vertex(ne3);

          int nv = mesh.Boundary(i).GetNbVertices();
          int nb_dof = 0;
          if (nv == 3)
          {
            // triangle
            Mlt(Real_wp(1)/Real_wp(3), center);
            nb_dof = nb_dof_tri;
          }
          else
          {
            // quadrangle
            int ne4 = mesh.Boundary(i).numVertex(3);
            center = Real_wp(0.25)*(center + mesh.Vertex(ne4));
            nb_dof = nb_dof_quad;
          }

          // expressing center in the base of vdirector_per
          for (int j = 0; j < Dimension::dim_N; j++)
            rhs(j) = DotProd(data.vdirector_per(j), center);

          Mlt(InvAper, rhs, coor);

          // storing this result
          CoorFace(i) = coor;

          if (PointInsideUnitCube(coor))
          {
            NumOriginal_Faces(nb_original_faces++) = i;
            for (int k = 0; k < nb_dof; k++)
              IsDofPeriodic(offset + k) = false;
          }
          else
            IsFacePeriodic(i) = true;

          offset += nb_dof;
          offset_face(i+1) = offset;
        }

        // then finding translation of periodic faces
        for (int i = 0; i < mesh.GetNbFaces(); i++)
          if (IsFacePeriodic(i))
          {
            int nb_dof = offset_face(i+1) - offset_face(i);
            for (int j = 0; j < nb_original_faces; j++)
            {
              int nf = NumOriginal_Faces(j);
              if (VertexTranslated(CoorFace(i), CoorFace(nf), n1))
              {
                for (int k = 0; k < nb_dof; k++)
                {
                  int ni = offset_face(i) + k;
                  int nj = offset_face(nf) + k;
                  PeriodicityDof(ni, 0) = nj;
                  for (int k = 0; k < Dimension::dim_N; k++)
                    PeriodicityDof(ni, k+1) = n1(k);
                }

                // checking that vertices are ok (no rotation)
                for (int k = 0; k < mesh.Boundary(i).GetNbVertices(); k++)
                {
                  int ne1 = mesh.Boundary(i).numVertex(k);
                  int ne2 = mesh.Boundary(nf).numVertex(k);
                  TinyVector<int, 3> n2;
                  if (!VertexTranslated(mesh.Vertex(ne1), mesh.Vertex(ne2), n2))
                  {
                    //cout << "case not handled " << endl;
                    //abort();
                  }
                  else
                  {
                    if (n2 != n1)
                    {
                      //cout << "case not handled " << endl;
                      //abort();
                    }                                  
                  }
                }
              }
            }
          }
      }    
    }


    void TreatPeriodicityElements(Mesh<Dimension>& mesh,
        MeshNumbering<Dimension>&mesh_num, int offset)
    {
      R_N rhs, coor, center;
      TinyVector<int, 3> n1;
      int order = mesh_num.GetOrder();
      //DISP(order);
      //DISP(mesh_num.number_map.GetNbDofVertex(order));
      //DISP(mesh_num.number_map.GetNbDofEdge(order));
      //DISP(mesh_num.number_map.GetNbDofTriangle(order));
      //DISP(mesh_num.number_map.GetNbDofQuadrangle(order));

      // we find all the elements inside this pattern
      IVect NumOriginal_Elements(mesh.GetNbElt());
      int nb_original_elt = 0;

      Vector<R_N> CenterElt(mesh.GetNbElt());    
      Vector<bool> IsElementPeriodic(mesh.GetNbElt());
      IsElementPeriodic.Fill(true);
      IVect OffsetElt(mesh.GetNbElt()+1); OffsetElt.Fill(0);
      OffsetElt(0) = offset;
      for (int i = 0; i < mesh.GetNbElt(); i++)
      {
        // computing barycenter of the element
        center.Fill(Real_wp(0));
        for (int j = 0; j < mesh.Element(i).GetNbVertices(); j++)
          Add(Real_wp(1), mesh.Vertex(mesh.Element(i).numVertex(j)), center);

        Real_wp coef = Real_wp(1)/Real_wp(mesh.Element(i).GetNbVertices());
        Mlt(coef, center);

        // expressing barycenter in the base of vdirector_per
        for (int j = 0; j < Dimension::dim_N; j++)
          rhs(j) = DotProd(data.vdirector_per(j), center);

        Mlt(InvAper, rhs, coor);

        // storing this result
        CenterElt(i) = coor;

        int nb_dof = mesh_num.number_map.GetNbDofElement(order, mesh.Element(i));
        if (PointInsideUnitCube(coor))
        {
          // the element is original
          NumOriginal_Elements(nb_original_elt) = i;
          for (int k = 0; k < nb_dof; k++)
            IsDofPeriodic(offset + k) = false;

          nb_original_elt++;
        }
        else
          IsElementPeriodic(i) = true;

        offset += nb_dof;
        OffsetElt(i+1) = offset;
      }

      // then finding translation of periodic elements
      for (int i = 0; i < mesh.GetNbElt(); i++)
        if (IsElementPeriodic(i))
        {
          int nb_dof = OffsetElt(i+1) - OffsetElt(i);
          for (int j = 0; j < nb_original_elt; j++)
          {
            if (VertexTranslated(CenterElt(i), CenterElt(NumOriginal_Elements(j)), n1))
            {
              for (int k = 0; k < nb_dof; k++)
              {
                int ni = OffsetElt(i) + k;
                int nj = OffsetElt(NumOriginal_Elements(j)) + k;
                PeriodicityDof(ni, 0) = nj;
                for (int k = 0; k < Dimension::dim_N; k++)
                  PeriodicityDof(ni, k+1) = n1(k);
              }
            }
          }
        }
    }

    // for nodal H1 only
    void CheckPeriodicity()
    {
      Mesh<Dimension>& mesh = var.mesh;
      MeshNumbering<Dimension>& mesh_num = var.mesh_num;
      Vector<R_N> s, CoordinateDofs(mesh_num.GetNbDof());  
      SetPoints<Dimension> PointsElem; R_N rhs, coor;
      for (int i = 0; i < mesh.GetNbElt(); i++)
      {
        mesh.GetVerticesElement(i, s);
        var.FjElemNodal(s, PointsElem, mesh, i);
        int nb_dof = mesh_num.Element(i).GetNbDof();
        for (int j = 0; j < nb_dof; j++)
        {
          for (int k = 0; k < Dimension::dim_N; k++)
            rhs(k) = DotProd(PointsElem.GetPointNodal(j), data.vdirector_per(k));

          Mlt(InvAper, rhs, coor);            
          CoordinateDofs(mesh_num.Element(i).GetNumberDof(j)) = coor;
        }         
      }

      TinyVector<int, 3> n1, n2;
      for (int i = 0; i < mesh_num.GetNbDof(); i++)
        if (IsDofPeriodic(i))
        {
          int j = PeriodicityDof(i, 0);
          int node = Num_OriginalDofs(j);
          if (VertexTranslated(CoordinateDofs(i), CoordinateDofs(node), n1))
          {
            for (int k = 0; k < Dimension::dim_N; k++)
              n2(k) = PeriodicityDof(i, k+1);

            if (n1 != n2)
            {
              DISP(i); DISP(node);
              DISP(CoordinateDofs(i)); DISP(CoordinateDofs(node));
              DISP(n1); DISP(n2);
              abort();
            }
          }
          else
          {
            abort();
          }
        }
    }

    // constructs informations on periodicity of dofs of the mesh
    void GetPeriodicityDofs()
    {    
      // mass matrix associated with periodic directions
      for (int i = 0; i < Dimension::dim_N; i++)
        for (int j = i; j < Dimension::dim_N; j++)
          Aperiodicity(i,j) = DotProd(data.vdirector_per(i), data.vdirector_per(j));

      GetInverse(Aperiodicity, InvAper);

      // construction of basic pattern
      ConstructBaseMotif(data.vdirector_per, mesh_elt);

      // allocating arrays;
      MeshNumbering<Dimension>& mesh_num = var.GetMeshNumbering(0);
      int order = mesh_num.GetOrder();
      IsDofPeriodic.Reallocate(mesh_num.GetNbDof());
      IsDofPeriodic.Fill(true);
      int nb_original_dof = 0;
      Num_OriginalDofs.Reallocate(mesh_num.GetNbDof());
      PeriodicityDof.Reallocate(mesh_num.GetNbDof(), Dimension::dim_N+1);
      PeriodicityDof.Fill(-1);

      // then finding original and periodic dofs associated with vertices
      int nb_dof_vertex = mesh_num.number_map.GetNbDofVertex(order);    
      if (nb_dof_vertex > 0)
        TreatPeriodicityVertices(var.mesh, nb_dof_vertex);

      // then finding original and periodic dofs associated with vertices
      int nb_dof_edge = mesh_num.number_map.GetNbDofEdge(order);    
      if (nb_dof_edge > 0)
        TreatPeriodicityEdges(var.mesh, nb_dof_vertex, nb_dof_edge);

      int offset;
      int nb_dof_tri = mesh_num.number_map.GetNbDofTriangle(order);    
      int nb_dof_quad = mesh_num.number_map.GetNbDofQuadrangle(order);
      TreatPeriodicityFaces(var.mesh, nb_dof_vertex, nb_dof_edge,
          nb_dof_tri, nb_dof_quad, offset);

      TreatPeriodicityElements(var.mesh, var.GetMeshNumbering(0), offset);

      // counting original dofs
      nb_original_dof = 0;
      IVect InverseDof(mesh_num.GetNbDof()); InverseDof.Fill(-1);
      for (int i = 0; i < mesh_num.GetNbDof(); i++)
      {
        if (IsDofPeriodic(i))
        {
          if (PeriodicityDof(i, 0) == -1)
          {
            DISP(i);
            abort();
          }
        }
        else
        {
          InverseDof(i) = nb_original_dof;
          nb_original_dof++;
        }
      }

      // updating PeriodicityDof
      Num_OriginalDofs.Reallocate(nb_original_dof);
      nb_original_dof = 0;
      for (int i = 0; i < mesh_num.GetNbDof(); i++)
      {
        if (IsDofPeriodic(i))
        {
          int nv = PeriodicityDof(i, 0);
          int nt = InverseDof(nv);
          PeriodicityDof(i, 0) = nt;
        }
        else
        {
          Num_OriginalDofs(nb_original_dof++) = i;
          PeriodicityDof(i, 0) = InverseDof(i);
          for (int k = 0; k < Dimension::dim_N; k++)
            PeriodicityDof(i, k+1) = 0;
        }
      }

      // checking dof coordinates (for H1)
      // CheckPeriodicity();

      // volume of elementary element
      typename Dimension::MatrixN_N dfj;
      for (int i = 0; i < Dimension::dim_N; i++)
        for (int j = 0; j < Dimension::dim_N; j++)
          dfj(i, j) = data.vdirector_per(j)(i);

      mean_edge = Real_wp(1)/pow(nb_original_dof, Real_wp(1)/Real_wp(3));
      volume_elt = Det(dfj); 
    }

#ifdef SELDON_WITH_MPI
    void GetDispersionError(const R_N& kwave, DistributedMatrix<Complex_wp, General, ArrayRowSparse>& Dh,
        DistributedMatrix<Complex_wp, General, ArrayRowSparse>& Kh,
        DistributedMatrix<Complex_wp, General, ArrayRowSparse>& Sh,
        Real_wp& error_dispersion, bool first_formulation)
#else
      void GetDispersionError(const R_N& kwave, Matrix<Complex_wp, General, ArrayRowSparse>& Dh,
          Matrix<Complex_wp, General, ArrayRowSparse>& Kh,
          Matrix<Complex_wp, General, ArrayRowSparse>& Sh,
          Real_wp& error_dispersion, bool first_formulation)
#endif
      {
        Real_wp omega, arg, zero(0), omega2, one(1);
        Complex_wp czero(zero, zero), cone(one, zero);
        R_N ur;
        omega = Norm2(kwave);
        omega2 = omega*omega;
        // DISP(kwave); DISP(omega);


        bool presence_Sh = false;
        if(Sh.GetNonZeros()>0)
        {
          cout << "There is a flow" << endl;
          presence_Sh = true;
        }


        // phase for periodic dofs = exp(i k \cdot x)
        int Nvol = IsDofPeriodic.GetM(), Nper = Num_OriginalDofs.GetM();
        phase.Reallocate(Nvol);
        phase.Fill(cone);
        for (int i = 0; i < IsDofPeriodic.GetM(); i++)
          if (IsDofPeriodic(i))
          {
            ur.Fill(zero);
            for (int k = 0; k < Dimension::dim_N; k++)
              Add(Real_wp(PeriodicityDof(i, k+1)), data.vdirector_per(k), ur); 

            arg = DotProd(kwave, ur);
            phase(i) = Complex_wp(cos(arg), sin(arg));
          }

        // periodic matrices
        int nb_unknowns = var.nb_unknowns;
        Dh_per.Reallocate(Nper*nb_unknowns, Nper*nb_unknowns);
        Kh_per.Reallocate(Nper*nb_unknowns, Nper*nb_unknowns);
        Dh_per.Fill(czero);
        Kh_per.Fill(czero);
        Complex_wp vloc;
        if(presence_Sh)
        {
          Sh_per.Reallocate(Nper*nb_unknowns, Nper*nb_unknowns);
          Sh_per.Fill(czero);
        }

        // on assemble les matrices periodiques
        for (int i1 = 0; i1 < Nper; i1++)
        {
          // on recupere le ddl original numero i
          int i = Num_OriginalDofs(i1);
          for (int ik = 0; ik < nb_unknowns; ik++)
            for (int j1 = 0; j1 < Dh.GetRowSize(i+ik*Nvol); j1++)
            {
              int jn = Dh.Index(i+ik*Nvol, j1);
              // j : numero du ddl, jk : numero d'inconnue pour la colonne consideree
              int j = jn%Nvol, jk = jn/Nvol;

              if (IsDofPeriodic(j))
              {
                // cas d'un ddl periodique, on rajoute phase * valeur
                int num_dof = PeriodicityDof(j, 0);
                vloc = phase(j)*Dh.Value(i+ik*Nvol,j1);
                Dh_per.AddInteraction(i1+ik*Nper, num_dof+jk*Nper, vloc);
              }
              else
              {
                // cas d'un ddl original, on rajoute valeur
                int num_dof = PeriodicityDof(j,0);
                vloc = Complex_wp(Dh.Value(i+ik*Nvol,j1));
                //vloc = Complex_wp(Dh.Value(i+ik*Nvol,j1), zero);
                Dh_per.AddInteraction(i1+ik*Nper, num_dof+jk*Nper, vloc);
              }
            }

          // meme chose pour la matrice Kh
          for (int ik = 0; ik < nb_unknowns; ik++)
            for (int j1 = 0; j1 < Kh.GetRowSize(i+ik*Nvol); j1++)
            {
              int jn = Kh.Index(i+ik*Nvol, j1);
              int j = jn%Nvol, jk = jn/Nvol;

              if (IsDofPeriodic(j))
              {
                int num_dof = PeriodicityDof(j,0);
                vloc = phase(j)*Kh.Value(i+ik*Nvol,j1);
                Kh_per.AddInteraction(i1+ik*Nper, num_dof+jk*Nper, vloc);
              }
              else
              {
                int num_dof = PeriodicityDof(j,0);
                vloc = Complex_wp(Kh.Value(i+ik*Nvol,j1));
                //vloc = Complex_wp(Kh.Value(i+ik*Nvol,j1), zero);
                Kh_per.AddInteraction(i1+ik*Nper, num_dof+jk*Nper, vloc);
              }
            }

          // meme chose pour la matrice Sh
          if(presence_Sh)
          {
            for (int ik = 0; ik < nb_unknowns; ik++)
              for (int j1 = 0; j1 < Sh.GetRowSize(i+ik*Nvol); j1++)
              {
                int jn = Sh.Index(i+ik*Nvol, j1);
                int j = jn%Nvol, jk = jn/Nvol;

                if (IsDofPeriodic(j))
                {
                  int num_dof = PeriodicityDof(j,0);
                  vloc = phase(j)*Sh.Value(i+ik*Nvol,j1);
                  Sh_per.AddInteraction(i1+ik*Nper, num_dof+jk*Nper, vloc);
                }
                else
                {
                  int num_dof = PeriodicityDof(j,0);
                  vloc = Complex_wp(Sh.Value(i+ik*Nvol,j1));
                  //vloc = Complex_wp(Sh.Value(i+ik*Nvol,j1), zero);
                  Sh_per.AddInteraction(i1+ik*Nper, num_dof+jk*Nper, vloc);
                }
              }
          }
        }


        // now we are looking for eigenvalues near omega
        // pour l'ordre 2, on cherche pres de omega^2
        Complex_wp sigma;
        if (first_formulation)
          sigma = Iwp*omega;
        else
          sigma = omega2;

        if (data.type_solver_eigenval == data.LAPACK)
        {
          Matrix<Complex_wp> Dh_dense, Kh_dense, Ah_dense;
          if(!presence_Sh)
          { 
            // conversion sur des matrices denses
            Copy(Dh_per, Dh_dense);
            Copy(Kh_per, Kh_dense);
            //Dh_dense.Write("DhPer.dat");
            //Kh_dense.Write("KhPer.dat");

            // je calcule Ah = D_h^{-1} K_h
            Ah_dense.Reallocate(Nper*nb_unknowns, Nper*nb_unknowns);

            // Lapack solver with dense matrices
            // now we compute Dh_per^{-1} Kh_per
            IVect ipivot(Nper*nb_unknowns);
            GetLU(Dh_dense, ipivot);
            VectComplex_wp Xh(Nper*nb_unknowns), Yh(Xh.GetM());
            for (int i = 0; i < Nper*nb_unknowns; i++)
            {
              for (int j = 0; j < Nper*nb_unknowns; j++)
                Xh(j) = Kh_dense(j, i);

              SolveLU(Dh_dense, ipivot, Xh);

              for (int j = 0; j < Nper*nb_unknowns; j++)
                Ah_dense(j, i) = Xh(j);
            }

            // we get eigenvalues in Yh
            GetEigenvalues(Ah_dense, Yh);

            Real_wp err_minimal = Real_wp(100), coef(0);
            error_dispersion = Real_wp(1);

            // boucle sur toutes les valeurs propres 
            // pour calculer l'erreur de dispersion
            // dans le cas Helmholtz/Maxwell, on cherche la valeur propre la plus proche
            // de omega (formulation d'ordre 1), et de omega^2 (formulation d'ordre 2)
            for (int i = 0; i < Nper*nb_unknowns; i++)
            {
              arg = abs(Yh(i));
              //DISP(i); DISP(arg); DISP(arg-omega); 
              if (first_formulation)
                coef = abs(omega-arg);
              else
                coef = abs(omega2-arg);

              if (coef < err_minimal)
              {
                if (first_formulation)
                  error_dispersion = arg-omega;
                else
                  error_dispersion = arg-omega2;

                err_minimal = coef;
              }
            }
          }
          else // presence d'un flow
          {
            //Matrix<Complex_wp> Dh_dense, Kh_dense, Ah_dense, Sh_dense; 
            Matrix<Complex_wp> Sh_dense;
            Copy(Dh_per, Dh_dense);
            Copy(Kh_per, Kh_dense);
            Copy(Sh_per, Sh_dense);

            GetInverse(Dh_dense);
            Ah_dense.Reallocate(2*Nper*nb_unknowns, 2*Nper*nb_unknowns);
            Ah_dense.Fill(czero);
            Matrix<Complex_wp> prod_Dh_Kh;
            Copy(Dh_dense, prod_Dh_Kh);
            Mlt(Dh_dense, Kh_dense, prod_Dh_Kh);
            Matrix<Complex_wp> prod_Dh_Sh; 
            Copy(Dh_dense, prod_Dh_Sh);
            Mlt(Dh_dense, Sh_dense, prod_Dh_Sh);

            for(int i=0; i<Nper*nb_unknowns; i++)
            {
              Ah_dense(i,Nper*nb_unknowns+i) = -1.0; // Ou plus ???
              for(int j=0; j<Nper*nb_unknowns; j++)
              {
                Ah_dense(Nper*nb_unknowns+i,Nper*nb_unknowns+j) = prod_Dh_Sh(i,j);
                Ah_dense(Nper*nb_unknowns+i,j) = -1.0*prod_Dh_Kh(i,j);
              }
            }
            VectComplex_wp Yh(2*Nper*nb_unknowns);
            GetEigenvalues(Ah_dense, Yh);
            DISP(kwave);Yh.Write("vp.dat");abort();
            Real_wp err_minimal(100);
            error_dispersion = Real_wp(1);
            for(int j=0; j<2*Nper*nb_unknowns; j++)
            {
              //arg = abs(Yh(j));
              Complex_wp vp = Yh(j);
              //On evalue la dispersion theorique : omega = M\cdot k +/- |k|
              R_N M; M(0) = flow(0); M(1) = flow(1);
              Real_wp convection = M(0)*kwave(0)+M(1)*kwave(1);
              // Multiplier par Iwp ? Oui ou non ?
              // formulation : ordre 2, mais relation de dispersion ordre 1 (obtenue en factorisant un polynome de degre 2) 
              // Rappel : omega = norme(kwave), omega2 = omega^2
              Complex_wp err_eig1 = vp-convection-omega;
              Complex_wp err_eig2 = vp-convection+omega;
              Complex_wp err_eig_normal = vp-convection;
              //DISP(Yh(j));// DISP(err_eig1); DISP(err_eig2);


              if(abs(err_eig1) < err_minimal)
              {
                error_dispersion = abs(err_eig1);
                err_minimal = abs(err_eig1);
              }
              if(abs(err_eig2) < err_minimal)
              {
                error_dispersion = abs(err_eig2);
                err_minimal = abs(err_eig2);
              }
              if(is_galbrun_formulation && abs(err_eig_normal) < err_minimal)
              {
                error_dispersion = abs(err_eig_normal);
                err_minimal = abs(err_eig_normal);
              }

            }
          }

        }
        else if (data.type_solver_eigenval == data.ITERATIVE_POWER)
        {
          // iterative power
#ifdef SELDON_WITH_MPI
          DistributedMatrix<Complex_wp, General, ArrayRowSparse> Ah(Kh_per);
#else
          Matrix<Complex_wp, General, ArrayRowSparse> Ah(Kh_per);
#endif
          Add(-sigma, Dh_per, Ah);

          IVect perm(Ah.GetM()); perm.Fill();

#ifdef SELDON_WITH_PASTIX
          MatrixPastix<double> mat_lu;
          mat_lu.HideMessages();
          mat_lu.FindOrdering(Ah, perm, true);
#endif

          IlutPreconditioning<Complex_wp> mat_ilut;
          mat_ilut.SetPrintLevel(1);
          mat_ilut.SetFactorisationType(mat_ilut.ILUT);
          mat_ilut.SetPivotThreshold(0.1);
          mat_ilut.SetDroppingThreshold(0.0);

          mat_ilut.FactorizeMatrix(perm, Ah);

          VectComplex_wp xn(Kh_per.GetM()), yn(Kh_per.GetM());
          for (int i = 0; i < yn.GetM(); i++)
            yn(i) = Complex_wp(rand(), rand());

          Mlt(Real_wp(1e-8), yn); xn.Zero();

          bool test_loop = true; Real_wp lambda = 0.1, lambda_prec = 0;
          int nb_max_iter = 1000, nb_iter = 0;
          while (test_loop)
          {
            lambda_prec = lambda;

            xn = yn;
            MltAdd(cone, Kh_per, xn, czero, yn);
            MltAdd(sigma, Dh_per, xn, cone, yn);

            mat_ilut.Solve(yn);
            lambda = Norm2(yn)/Norm2(xn);
            if (abs(lambda - lambda_prec)/abs(lambda) < 10*sqrt(lambda)*epsilon_machine)
              test_loop = false;

            Mlt(one/Norm2(yn), yn);

            nb_iter++;
            if (nb_iter > nb_max_iter)
              test_loop = false;
          }

          if (nb_iter > nb_max_iter)
          {
            cout << "Research of closest eigenvalue to omega failed " << endl;
            DISP(lambda); DISP(lambda_prec);
            DISP(abs(lambda - lambda_prec)/abs(lambda));
          }

          if (real(yn(0))*real(xn(0)) < 0)
            lambda  = -lambda;

          arg = abs(sigma)*abs( (one + lambda)/(one - lambda));
          error_dispersion = arg-abs(sigma);
        }
        else
        {
          // Dh_per.WriteText("Dh.dat");
          // Kh_per.WriteText("Kh.dat");
          int N = Kh_per.GetM();

          // Arpack-like algorithm
#ifdef SELDON_WITH_MPI
          DistributedMatrix<Complex_wp, General, ArrayRowSparse> Ah(Kh_per);
#else
          Matrix<Complex_wp, General, ArrayRowSparse> Ah(Kh_per);
#endif
          Add(-sigma, Dh_per, Ah);

          IVect perm(Ah.GetM()); perm.Fill();
          VectComplex_wp Ones(Ah.GetM()), AhOnes(Ah.GetM());
          for (int i = 0; i < N; i++)
            Ones(i) = Complex_wp(Real_wp(i), zero);

          AhOnes.Fill(czero);
          MltAdd(cone, Ah, Ones, czero, AhOnes);

          // GetLU(Ah);

#ifdef SELDON_WITH_PASTIX
          MatrixPastix<double> mat_lu;
          mat_lu.HideMessages();
          mat_lu.FindOrdering(Ah, perm, true);
#endif

          IlutPreconditioning<Complex_wp> mat_ilut;
          mat_ilut.SetPrintLevel(1);
          mat_ilut.SetFactorisationType(mat_ilut.ILUT);
          mat_ilut.SetPivotThreshold(0.1);
          mat_ilut.SetDroppingThreshold(0.0);

          mat_ilut.FactorizeMatrix(perm, Ah, true);

          mat_ilut.Solve(AhOnes);

          for (int i = 0; i < N; i++)
            if (abs(AhOnes(i) - Complex_wp(Real_wp(i), zero)) > 1e-4)
            {
              DISP(i); DISP(AhOnes(i));
              cout << "Problem during inversion " << endl;
              //abort();
            }

          VectComplex_wp xn(N), yn(N);
          for (int i = 0; i < yn.GetM(); i++)
            yn(i) = Complex_wp(rand(), rand());

          Mlt(cone/Norm2(yn), yn); xn.Fill(0);

          int nb_eigenval = data.nb_eigenvalues_wanted;
          Vector<VectComplex_wp> Vn(nb_eigenval);
          Vn(0) = yn;
          Matrix<Complex_wp> Hn(nb_eigenval, nb_eigenval);
          Hn.Fill(czero);
          for (int k = 1; k <= Hn.GetM(); k++)
          {
            xn = Vn(k-1);
            MltAdd(Complex_wp(one, zero), Kh_per, xn, Complex_wp(zero, zero), yn);
            MltAdd(sigma, Dh_per, xn, Complex_wp(one, zero), yn);

            //mat_ilut.Solve(yn);
            SolveLU(Ah, yn);

            for (int j = 0; j < k; j++)
            {
              Hn(j, k-1) = DotProdConj(Vn(j), yn);  
              Add(-Hn(j, k-1), Vn(j), yn);
            }

            Real_wp norm_yn = Norm2(yn);
            Mlt(Complex_wp(one/norm_yn, zero), yn);

            if (k < Hn.GetM())
            {
              Hn(k, k-1) = Complex_wp(norm_yn, zero);
              Vn(k) = yn;
            }
          }

          //DISP(Hn);
          VectComplex_wp lambda_all(Hn.GetM());
          GetEigenvalues(Hn, lambda_all);
          //DISP(lambda_all);

          Real_wp lambda = abs(lambda_all(0));
          for (int i = 1; i < lambda_all.GetM(); i++)
            lambda = max(lambda, abs(lambda_all(i)));

          arg = abs(sigma)*abs( (one + lambda)/(one - lambda));
          error_dispersion = arg-abs(sigma);
        }
      }

};


  template<class Dimension, class Matrix1>
void GetDispersionRelation(DispersionSolver<Dimension>& var_disp,
    MultivariatePolynomial<Real_wp>& dispersion_relation, Matrix1& Dh, Matrix1& Kh, Matrix1& Sh)
{
  typedef typename Dimension::R_N R_N;

  InputDispersionVariables<Dimension>& data = var_disp.data;

  /****************************
   * Declaration de variables *
   ****************************/

  bool first_formulation = var_disp.var.FirstOrderFormulation();
  if(Sh.GetNonZeros()>0 && first_formulation)
  {
    Add(Real_wp(1), Sh, Kh);
    Sh.Clear();
  }

  Real_wp arg, coef; 
  R_N kwave, ur;
  int order = data.order;

  Vector<R_N> Points; 
  Real_wp xsi1, xsi2, xsi3, theta, phi;
  Matrix<Real_wp> VanDerMonde;
  Real_wp zero(0), error_dispersion;
  cout<<" Number of original dof " << var_disp.Num_OriginalDofs.GetM() << endl;


  /************************************************
   * On cherche l'ordre de l'erreur de dispersion *
   ************************************************/


  if (data.find_order)
  {
    cout<<"We compute the order of dispersion " << endl;

    Real_wp err1, err2;
#ifdef SELDON_WITH_MPI
    double err1d; MPI_Status status;
    if (rank_proc == 0)
    {
#endif
      // calcul de abs(omega^2 - k^2) pour k = kh1_order
      kwave.Fill(zero); kwave(0) = data.kh1_order;
      var_disp.GetDispersionError(kwave, Dh, Kh, Sh, err1, first_formulation);

      err1 = abs(err1); 
      cout << "k1 = " << kwave << endl;
      cout << "|omega^2 - ||k1||^2 | = " << err1 << endl;

#ifdef SELDON_WITH_MPI
      err1d = toDouble(err1);
      if (nb_proc > 1)
      {
        // we send the error found to first processor
        MPI_Ssend(&err1d, 1, MPI_DOUBLE, 1, 16, MPI_COMM_WORLD);
      }
    }
#endif

#ifdef SELDON_WITH_MPI
    if (rank_proc == 1)
    {
#endif
      // calcul de abs(omega - k) pour k = kh2_order
      kwave(0) = data.kh2_order;
      var_disp.GetDispersionError(kwave, Dh, Kh, Sh, err2, first_formulation);

      err2 = abs(err2); 
      cout << "k2 = " << kwave << endl;
      cout << "|omega^2 - ||k2||^2 | = " << err2 << endl;

#ifdef SELDON_WITH_MPI
      if (nb_proc > 1)
      {
        // we retrive the error found by first processor
        MPI_Recv(&err1d, 1, MPI_DOUBLE, 0, 16, MPI_COMM_WORLD, &status);
        err1 = err1d;
      }
#endif

      // on calcule l'ordre avec un ratio de logs
      error_dispersion = (log(err2)-log(err1))/(log(data.kh2_order)-log(data.kh1_order));
      data.order_dispersion = toInteger(round(error_dispersion));
      if (first_formulation)
        cout<<"We found a dispersion order of "<<(error_dispersion-1)<<endl;
      else
        cout<<"We found a dispersion order of "<<(error_dispersion-2)<<endl;

#ifdef SELDON_WITH_MPI
    }
#endif

    // calcul de la constante dans l'erreur de dispersion
    // we have the relation
    // omega/k = 1 + C1 k^p   (with p = order_dispersion - 1/2)
    // we want to write it as 
    // omega/k = 1 + C2 K^p   with K = (6k)/(2 pi r)
    // so that K = 1 will correspond to a mesh with 6 points per wavelength
    // therefore C2 = C1 (k/K)^p = C1 [(2 pi r)/6]^p
    // err2 = omega - k2, hence err2/k2 = omega/k2 - 1 = C1 k2^p
    Real_wp k2 = Norm2(kwave);
    Real_wp C1 = 0.5*err2/pow(k2, error_dispersion);
    Real_wp C2 = 0;
    if (first_formulation)
      C2 = C1*pow(pi_wp*order/Real_wp(3), error_dispersion - 1);      
    else
      C2 = C1*pow(pi_wp*order/Real_wp(3), error_dispersion - 2);      

    cout << "If we /-write the dispersion relation as omega/k = 1 + C K^p" << endl;
    cout << "where K = (6 k)/(2 pi r) " << endl;
    cout << "The constant C is equal to " << C2*var_disp.mean_edge << endl;
  }
  else
  {
    data.order_dispersion = 2*(order+1);
  }


  /*************************************
   * Calcul de la courbe de dispersion *
   *************************************/


  if (data.nb_points_plot > 0)
  {
    int nb_points_effective = data.nb_points_plot;
    Matrix<double> result_plot(nb_points_effective, 3);
    result_plot.Fill(0);

#ifdef SELDON_WITH_MPI
    nb_points_effective = (data.nb_points_plot-rank_proc+nb_proc-1)/nb_proc;
    DISP(nb_points_effective);
    if (rank_proc > 0)
      result_plot.Reallocate(nb_points_effective, 3);

    int nplot = 0;
    for (int i = rank_proc; i < data.nb_points_plot; i += nb_proc)
#else
      for (int i = 0; i < data.nb_points_plot; i++)
#endif
      {
        // kwave va decrire des valeurs regulierement espacees
        // entre k1_plot et k2_plot (en echelle log)
        Real_wp lambda = Real_wp(i)/Real_wp(data.nb_points_plot-1);
        if (data.courbe_log)
        {
          Real_wp ls = lambda;
          lambda = pow(10, data.log_min*(1.0-ls));
          cout<<"Lambda "<<lambda<<endl;
        }
        kwave.Fill(zero);
        Add(Real_wp(1)-lambda, data.k1_plot, kwave);
        Add(lambda, data.k2_plot, kwave);
        Real_wp omega = Norm2(kwave);

        Real_wp err1(0);
        // on calcule err1 = |omega^2 - k1^2|
        var_disp.GetDispersionError(kwave, Dh, Kh, Sh, err1, first_formulation);

        cout << "Wave vector " << kwave << endl << "Error : " << err1 << endl;

        // ratio_omega = omega/k
        Real_wp ratio_omega = 0;
        if (first_formulation)
          ratio_omega = (err1+omega)/omega;
        else
          ratio_omega = sqrt((err1+omega*omega)/(omega*omega));

        Real_wp err_relative = abs(Real_wp(1) - ratio_omega);

        if (rank_proc > 0)
        {
#ifdef SELDON_WITH_MPI
          result_plot(nplot, 0) = toDouble(omega);
          result_plot(nplot, 1) = toDouble(ratio_omega);
          //result_plot(nplot, 2) = err_relative;
          result_plot(nplot, 2) = toDouble(err1);
#endif
        }
        else
        {
          result_plot(i, 0) = toDouble(omega);
          result_plot(i, 1) = toDouble(ratio_omega);
          //result_plot(i, 2) = err_relative;
          result_plot(i, 2) = toDouble(err1);
        }
#ifdef SELDON_WITH_MPI
        nplot++;
#endif
      }

    if (rank_proc == 0)
    {
      cout << "On calcule la courbe de dispersion "<<endl;
      cout << "Ecriture sur le fichier " << data.file_name_plot << endl;
      ofstream file_out(data.file_name_plot.data());
      file_out.precision(15);

#ifdef SELDON_WITH_MPI
      for (int i = 1; i < nb_proc; i++)
      {
        nb_points_effective = (data.nb_points_plot-i+nb_proc-1)/nb_proc;
        Matrix<double> result_other(nb_points_effective, 3);
        nplot = 0; MPI_Status status;
        MPI_Recv(result_other.GetData(), 3*nb_points_effective, MPI_DOUBLE, i, 3, MPI_COMM_WORLD, &status);
        for (int j = i; j < data.nb_points_plot; j+= nb_proc)
        {
          for (int k = 0; k < 3; k++)
            result_plot(j, k) = result_other(nplot, k);

          nplot++;
        }
      }
#endif
      result_plot.WriteText(file_out);
      file_out.close();
    }
    else
    {
#ifdef SELDON_WITH_MPI
      MPI_Ssend(result_plot.GetData(), 3*nb_points_effective, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD);
#endif
    }
  }


  /*****************************************
   * Calcul de la constante sur une sphere *
   *****************************************/


#ifdef MONTJOIE_WITH_THREE_DIM
  if (data.display_constant_sphere)
  {
    cout << "On calcule les constantes sur un maillage de sphere "<< endl;
    Mesh<Dimension> mesh_sphere;
    mesh_sphere.Read(data.file_sphere);

    Vector<Real_wp> const_err(mesh_sphere.GetNbVertices());
#ifdef SELDON_WITH_MPI
    for (int i = 0; i < mesh_sphere.GetNbVertices(); i+= nb_proc)
#else
      for (int i = 0; i < mesh_sphere.GetNbVertices(); i++)
#endif
      {
        kwave = mesh_sphere.Vertex(i);
        Mlt(data.kmax/Norm2(kwave), kwave);

        Real_wp err1(0);
        var_disp.GetDispersionError(kwave, Dh, Kh, Sh, err1, first_formulation);

        const_err(i) = abs(err1);
        Mlt(const_err(i), mesh_sphere.Vertex(i));
      }

    mesh_sphere.Write(data.file_sphere_output);
  }
#endif


  /**************************************************************************
   * Calcul de l'erreur de dispersion pour des k demandes par l'utilisateur *
   **************************************************************************/


  int Np = data.order_dispersion, nb = 0;
  for (int n = 0; n < data.nb_examples; n++)
  {
    kwave(0) = data.Xsi1Points(n);
    kwave(1) = data.Xsi2Points(n);
    if (Dimension::dim_N == 3)
      kwave(2) = data.Xsi3Points(n);

    var_disp.GetDispersionError(kwave, Dh, Kh, Sh, error_dispersion, first_formulation);

    error_dispersion = error_dispersion/(Norm2(kwave)*Norm2(kwave));
    int nb_original_dof = var_disp.Num_OriginalDofs.GetM();
    arg = pi_wp*pow(Real_wp(nb_original_dof/Dimension::dim_N)/var_disp.volume_elt,
        Real_wp(Real_wp(1)/Real_wp(Dimension::dim_N)))/3;

    arg /= Norm2(kwave);
    coef = pow(arg, Real_wp(Np-1));
    DISP(kwave); DISP(error_dispersion); DISP(error_dispersion*coef);
  }

  return;


  /********************************************************
   * Calcul du terme preponderant du developpement limite *
   ********************************************************/


  // Le vandermonde c'est pour identifier les coefficient aij 
  // du terme preponderant du developpement limite de l'erreur de dispersion
  // omega = k + \sum_{i,j} aij kx^i ky^j kz^{p-i-j}

  int nb_polynoms = Np + 1;
  if (Dimension::dim_N == 3)
    nb_polynoms = (Np+1)*(Np+2)/2;

  Points.Reallocate(nb_polynoms);
  VanDerMonde.Reallocate(Points.GetM(), Points.GetM());
  if (Dimension::dim_N == 2)
  {
    for (int i = 0; i <= Np; i++)
    {
      xsi1 = data.kmax*Real_wp(i)/Real_wp(Np);
      xsi2 = sqrt(data.kmax*data.kmax-xsi1*xsi1);
      for (int j = 0; j <= Np; j++)
        VanDerMonde(i,j) = pow(xsi1, Real_wp(j))*pow(xsi2, Real_wp(Np-j));

      Points(i)(0) = xsi1;
      Points(i)(1) = xsi2;
    }
  }
  else
  {
    srand(time(NULL));
    for (int i = 0; i < nb_polynoms; i++)
    {
      theta = 1e-9*rand()*pi_wp*2;
      phi = 1e-9*rand()*pi_wp;
      xsi1 = data.kmax*cos(theta)*cos(phi);
      xsi2 = data.kmax*cos(theta)*sin(phi);
      xsi3 = data.kmax*sin(theta);
      int icol = 0;
      for (int k = 0; k <= Np; k++)
        for (int l = 0; l <= (Np-k); l++)
          VanDerMonde(nb,icol++) = pow(xsi1, Real_wp(k))*pow(xsi2, Real_wp(l))*pow(xsi3, Real_wp(Np-k-l));

      Points(nb)(0) = xsi1;
      Points(nb)(1) = xsi2;
      Points(nb)(2) = xsi3;
      nb++;
    }
  }

  // DISP(VanDerMonde);
  IVect ipivot_vdm(VanDerMonde.GetM()); GetLU(VanDerMonde, ipivot_vdm);

  // DISP(Points);
  // DISP(vec_phi);

  int nb_points = Points.GetM();
  VectReal_wp coef_relation(nb_points); coef_relation.Fill(zero);
  // dispersion_relation.SetOrder(dim_N, 0);
  // dispersion_relation(0) = Real_wp(0);

  if (data.search_taylor_dvt)
  {
    cout<<"We compute dispersion coefficients "<<endl;
    cout << "On calcule le developpement de Taylor de l'erreur de dispersion " << endl;
    for (int n = 0; n < nb_points; n++)
    {
      cout<<"Point "<<n<<endl;

      var_disp.GetDispersionError(Points(n), Dh, Kh, Sh, error_dispersion, first_formulation);

      // DISP(Points(n)); DISP(error_dispersion);
      coef_relation(n) = error_dispersion;
    }
  }


  /**************************************************************
   * Calcul d'une constante moyenne dans l'erreur de dispersion *
   **************************************************************/


  // DISP(coef_relation);
  // L2 norm of the dispersion error
  VectReal_wp norm_L2(1); norm_L2(0) = zero;
  Globatto<Real_wp> gauss;
  gauss.ConstructQuadrature(data.order_integration, gauss.QUADRATURE_GAUSS);

  if (data.search_taylor_dvt)
  {
    SolveLU(VanDerMonde, ipivot_vdm, coef_relation);
    // DISP(coef_relation);
    cout<<"we compute rational approximation of each coefficient"<<endl;
    // coef_relation.PrintRational(threshold_rational);

    /* if (dim_N == 2)
       {
       for (int i = 0; i < gauss.nb_points_quadrature; i++)
       {      
       theta = pi_wp*2*gauss.Points(i);
       xsi1 = cos(theta); xsi2 = sin(theta);
       arg = zero;
       for (int j = 0; j <= Np; j++)
       arg += coef_relation(j)*pow(xsi1, Real_wp(j))*pow(xsi2, Real_wp(Np-j));

       coef = arg*arg;
       norm_L2(0) += gauss.Weights(i)*coef;
       }
       }
       else
       {

       for (int i = 0; i < gauss.nb_points_quadrature; i++)
       for (int j = 0; j < gauss.nb_points_quadrature; j++)
       {      
       theta = pi_wp*gauss.Points(i);
       phi = pi_wp*2*gauss.Points(j);
       xsi1 = cos(theta)*cos(phi); xsi2 = cos(theta)*sin(phi); xsi3 = sin(theta);
       arg = zero;
       nb = 0;
       for (int k = 0; k <= Np; k++)
       for (int l = 0; l <= (Np-k); l++)
       arg += coef_relation(nb++)*pow(xsi1, Real_wp(k))*pow(xsi2, Real_wp(l))*pow(xsi3, Real_wp(Np-k-l));

       coef = arg*arg;
       norm_L2(0) += gauss.Weights(i)*gauss.Weights(j)*coef*abs(xsi3);
       }
       norm_L2(0) *= pi_wp/2;
       }*/
  }
  else
  {
    if (Dimension::dim_N == 2)
    {
      for (int i = 0; i < gauss.GetNbPointsQuad(); i++)
      {      
        theta = pi_wp*2*gauss.Points(i);
        xsi1 = data.kmax*cos(theta);
        xsi2 = data.kmax*sin(theta);
        kwave(0) = xsi1; kwave(1) = xsi2;
        var_disp.GetDispersionError(kwave, Dh, Kh, Sh, error_dispersion, first_formulation);
        arg = Real_wp(1)/data.kmax;
        coef = pow(arg, Real_wp(2*Np));
        coef *= error_dispersion*error_dispersion;
        // DISP(i); DISP(kwave); DISP(error_dispersion/(kmax*kmax));
        // DISP(arg); DISP(coef); DISP(error_dispersion); DISP(norm_L2(0));
        norm_L2(0) += gauss.Weights(i)*coef; // DISP(norm_L2(0));
      }
    }
    else
    {
      if (Dimension::dim_N <= 2)
        return;

      // Quelle est l'erreur moyenne de dispersion pour || k || = kmax ?
      // moyenne sur toutes les directions possibles
      cout<<"Calcul de l'erreur moyenne de dispersion "<<endl;
      cout << "Balayage des angles teta, phi sur [0, 2pi] et [0,pi] " 
        << "avec " << data.nb_div_for_moyenne + 1 << " points sur chaque angle " << endl;

      int nb_div = data.nb_div_for_moyenne;
      Real_wp err1 = 0, err2  = 0;
      for (int i = 0; i <= nb_div; i++)
        for (int j = 0; j <= nb_div; j++)
        {
          theta = 2.0*pi_wp*Real_wp(i)/Real_wp(nb_div);
          phi = pi_wp*Real_wp(j)/Real_wp(nb_div);

          kwave(0) = data.kmax*cos(theta)*cos(phi);
          kwave(1) = data.kmax*cos(theta)*sin(phi);
          kwave(2) = data.kmax*sin(theta);

          // err2 est ici la valeur propre de D^{-1} K la plus proche de omega
          var_disp.GetDispersionError(kwave, Dh, Kh, Sh, err2, first_formulation);

          // DISP(kwave); DISP(err2);
          err1 += abs(err2);
        }

      err1 /= (nb_div+1)*(nb_div+1);

      cout << "The mean value of |omega^2 - ||k||^2 |  for all directions and ||k|| = " << data.kmax << endl;
      cout << "is equal to " << err1 << endl;

      // DISP(Np);
      Real_wp C1 = 0;
      if (first_formulation)
        C1 = 0.5*err1/pow(data.kmax,Np);
      else
        C1 = err1/pow(data.kmax,Np);

      Real_wp C2 = C1*pow(pi_wp*order/Real_wp(3), Np - 1);      

      cout << "If we write the dispersion relation as omega/k = 1 + C K^p" << endl;
      cout << "where K = (6 k)/(2 pi r) " << endl;
      cout << "The mean value of constant C is equal to " << C2 << endl;

    }
  }

  // for (int i = 0; i < coef_relation.GetM(); i++)
  //  coef_relation(i) *= coef_relation(i);
  // PrintRational(coef_relation, threshold);
  // coef = Real_wp(1)/(kmax*kmax);
  // ScalePolynomVariable(dispersion_relation, coef);

  // DISP(dispersion_relation);
  // Real_wp epsilon_machine = GetPrecisionMachine(Real_wp(1));
  // Real_wp threshold = sqrt(epsilon_machine); DISP(threshold);
  // CompressPolynom(dispersion_relation, threshold);
  // dispersion_relation.PrintRational(threshold);

  norm_L2(0) = sqrt(norm_L2(0));
  // DISP(norm_L2); 
  // PrintRational(norm_L2, threshold_rational); 
  cout<<"Norm after renormalization of h "<<endl;
  int nb_original_dof = var_disp.Num_OriginalDofs.GetM();
  arg = pi_wp*pow(Real_wp(nb_original_dof)/var_disp.volume_elt, Real_wp(Real_wp(1)/Real_wp(Dimension::dim_N)))/3;
  DISP(arg);
  coef = pow(arg, Real_wp(Np-2));
  norm_L2(0) *= coef; 
  DISP(norm_L2);
}


  template<class Dimension>
void ConstructReferenceElement_And_AssembleMatrix(DispersionSolver<Dimension>& var_disp,
    DistributedMatrix<Complex_wp, General, ArrayRowSparse>& Dh,
    DistributedMatrix<Complex_wp, General, ArrayRowSparse>& Kh,
    DistributedMatrix<Complex_wp, General, ArrayRowSparse>& Sh,
    const string& name_element)
{
  DistributedProblem<Dimension>& var = var_disp.var;
  InputDispersionVariables<Dimension>& data = var_disp.data;
  MeshNumbering<Dimension>& mesh_num = var.GetMeshNumbering(0);
 
  //var.InitIndices(10);
  var.mesh_data.Clear();

  // construction du maillage
  ConstructMesh(var.mesh, data.type_hybrid_mesh, data.regular_mesh,
      data.angle_deg, data.name_file, data.irregular_mesh, data.vdirector_per);

  for (int i = 1; i < var.mesh.GetNbReferences(); i++)
    var.mesh.SetBoundaryCondition(i, BoundaryConditionEnum::LINE_NEUMANN);

  int order = var_disp.data.order;
  var.mesh.ResizeNbReferences(6);
  mesh_num.SetOrder(order);
  DISP(mesh_num.GetOrder());

  var.GetOutputProblem().grid_to_be_computed = false;
  var.alpha_penalization = data.alpha_penal;
  var.delta_penalization = data.delta_penal;  
  var.ComputeMeshAndFiniteElement(name_element, false);


  // calcul des matrices
  int nb_unknowns = var.nb_unknowns;
  Dh.Reallocate(nb_unknowns*mesh_num.GetNbDof(), nb_unknowns*mesh_num.GetNbDof());

  GlobalGenericMatrix<Complex_wp> nat_mat;
  nat_mat.SetCoefMass(Real_wp(1)); nat_mat.SetCoefStiffness(Real_wp(0)); nat_mat.SetCoefDamping(Real_wp(0));
  var.ComputeMassMatrix();
  var.GetComputationProblem().AddMatrixFEM(Dh, nat_mat);

  Kh.Reallocate(Dh.GetM(), Dh.GetM());
  nat_mat.SetCoefMass(Real_wp(0)); nat_mat.SetCoefStiffness(Real_wp(1)); nat_mat.SetCoefDamping(Real_wp(0));
  var.GetComputationProblem().AddMatrixFEM(Kh, nat_mat);

  Sh.Reallocate(Dh.GetM(), Dh.GetM());
  nat_mat.SetCoefMass(Real_wp(0)); nat_mat.SetCoefStiffness(Real_wp(0)); nat_mat.SetCoefDamping(Real_wp(1));
  var.GetComputationProblem().AddMatrixFEM(Sh, nat_mat);
  if (data.print_level > 3)
  {
    var.CheckContinuity();
  }


  // For Galbrun, we need to evaluate the flow
  try
  {
    //VarGalbrun_Eq<TypeEquation>& g = dynamic_cast<VarGalbrun_Eq<TypeEquation>& >(var_);
    VarGalbrunIndex_Base<Dimension>& g = dynamic_cast<VarGalbrunIndex_Base<Dimension>& >(var);
    cout << "Galbrun : " << g.eval_flow.GetM() << endl;
    var_disp.flow = g.eval_flow(0)(0);
    var_disp.is_galbrun_formulation = true;
  }
  catch (const bad_cast& e)
  {
  }

  Dh.WriteText("mat_Dh.dat");
  Kh.WriteText("mat_Kh.dat");
  Sh.WriteText("mat_Sh.dat");
}


template<class Dimension>
void RunAll(DispersionSolver<Dimension>& var_disp, const string& file_input,
            const string& name_element)
{
  DistributedProblem<Dimension>& var = var_disp.var;

  var.InitIndices(50);

  // lecture du fichier de donnees
  ReadInputFile(file_input, var_disp.data);
  ReadInputFile(file_input, var.GetOutputProblem());

  Vector<string> param(1);
  if (var_disp.data.use_warburton_trick)
  {
    param(0) = "YES";
    var.SetInputData("UseWarburtonTrick", param);
  }

  // matrice de masse et rigidite
  // en ordre 2 : chercher omega et X tel que (- omega^2 Dh + Kh) X = 0 
  //              avec conditions de quasi-periodicite
  DistributedMatrix<Complex_wp, General, ArrayRowSparse>  Dh, Kh, Sh;
  MultivariatePolynomial<Real_wp> dispersion_relation;

  // calcul des matrices
  var.exit_if_no_boundary_condition = false;  
  ConstructReferenceElement_And_AssembleMatrix(var_disp, Dh, Kh,Sh, name_element);
  cout<<"global matrices computed"<<endl;

  // ddls periodiques
  var_disp.GetPeriodicityDofs();

  cout<<"Periodicity treated"<<endl;
  // calcul de la relation de dispersion
  GetDispersionRelation(var_disp, dispersion_relation, Dh, Kh, Sh);

}


int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);

#ifdef SELDON_WITH_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank_proc);
  MPI_Comm_size(MPI_COMM_WORLD, &nb_proc);
#endif

  if (argc != 2)
  {
    cout<<"entrez le nom du fichier de donnees"<<endl;
    return -1;
  }

  string input_file(argv[1]);

  cout.precision(15);
  srand(time(NULL));
  DISP(epsilon_machine);
  R2::threshold = 100*epsilon_machine;
  R3::threshold = 100*epsilon_machine;

  // we get the type of element selected by the user, and type of equation
  string type_element, type_equation;
  getElement_Equation(input_file, type_element, type_equation);

  int dim_N = 3;
  if ((type_element.find("TRIANGLE") == 0) || (type_element.find("QUADRANGLE") == 0))
    dim_N = 2;

  cout << "TypeEquation = " << type_equation << endl;
  if (dim_N == 2)
  {
    InputDispersionVariables<Dimension2> data; 
    ReadInputFile(input_file, data);

    if (type_equation == "LAPLACE_DIV")
    {
      // formulation H(div) de Helmholtz
      EllipticProblem<LaplaceEquationHdiv<Dimension2> > var_harmonic;
      DispersionSolver<Dimension2> var(var_harmonic);
      RunAll(var, input_file, type_element);
    }
    else if (type_equation == "LAPLACE_DG")
    {
      // formulation LDG de Helmholtz
      EllipticProblem<LaplaceEquationDG<Dimension2> > var_harmonic;
      DispersionSolver<Dimension2> var(var_harmonic);
      RunAll(var, input_file, type_element);
    }
    else if (type_equation == "HARMONIC_GALBRUN")
    {
      // Formulation pq galbrun
      EllipticProblem<HarmonicGalbrunEquationDG<Dimension2> > var_harmonic;
      DispersionSolver<Dimension2> var(var_harmonic);
      RunAll(var, input_file, type_element);
    }
    else if (type_equation == "HARMONIC_GALBRUN_H1")
    {
      // Formulation agumentee galbrun
      EllipticProblem<HarmonicGalbrunEquationH1> var_harmonic;
      DispersionSolver<Dimension2> var(var_harmonic);
      RunAll(var, input_file, type_element);
    }
    else if (type_equation == "HARMONIC_GALBRUN_MODEL_DIV")
    {
      // Formulation pv galbrun
      EllipticProblem<HarmonicGalbrunEquation<Dimension2> > var_harmonic;
      DispersionSolver<Dimension2> var(var_harmonic);
      RunAll(var, input_file, type_element);
    }
    else if (type_equation == "HARMONIC_GALBRUN_MODEL_SIPG")
    {
      // Fromulation sipg galbrun
      EllipticProblem<HarmonicGalbrunEquationSipg<Dimension2> > var_harmonic;
      DispersionSolver<Dimension2> var(var_harmonic);
      RunAll(var, input_file, type_element);
    }
    else
    {
      // formulation H1 de Helmholtz
      EllipticProblem<LaplaceEquation<Dimension2> > var_harmonic;
      DispersionSolver<Dimension2> var(var_harmonic);
      RunAll(var, input_file, type_element);
    }
  }
  else
  {
    InputDispersionVariables<Dimension3> data; 
    ReadInputFile(input_file, data);

    if (type_equation == "MAXWELL_DG")
    {
      // formulation LDG pour Maxwell
      EllipticProblem<StaticMaxwellEquation_3D_DG> var_harmonic;
      DispersionSolver<Dimension3> var(var_harmonic);
      RunAll(var, input_file, type_element);
    }
    else if (type_equation == "MAXWELL_3D")
    {
      // formulation H(curl) pour Maxwell
      EllipticProblem<StaticMaxwellEquation_3D> var_harmonic;
      DispersionSolver<Dimension3> var(var_harmonic);
      RunAll(var, input_file, type_element);
    }
    else if (type_equation == "LAPLACE_DIV")
    {
      // formulation H(div) pour Helmholtz
      EllipticProblem<LaplaceEquationHdiv<Dimension3> > var_harmonic;
      DispersionSolver<Dimension3> var(var_harmonic);
      RunAll(var, input_file, type_element);
    }
    else if (type_equation == "LAPLACE_DG")
    {
      // formulation LDG pour Helmholtz
      EllipticProblem<LaplaceEquationDG<Dimension3> > var_harmonic;
      DispersionSolver<Dimension3> var(var_harmonic);
      RunAll(var, input_file, type_element);
    }
    else
    {
      // formulation H1 pour Helmholtz
      EllipticProblem<LaplaceEquation<Dimension3> > var_harmonic;
      DispersionSolver<Dimension3> var(var_harmonic);
      RunAll(var, input_file, type_element);
    }
  }

  return FinalizeMontjoie();
}
