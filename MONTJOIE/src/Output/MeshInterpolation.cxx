#ifndef MONTJOIE_FILE_MESH_INTERPOLATION_CXX

namespace Montjoie
{

  //! default constructor
  template<class Dimension>
  MeshInterpolation_Base<Dimension>::MeshInterpolation_Base()
    : proj(), proj_surf()
  {
    nb_subdivisions = 0;
    type_mesh = VOLUME_MESH;
    nb_points_all_quadrature = 0;
    rank_processor = 0;
  }


  //! returns the size of memory used by the object in bytes
  template<class Dimension>
  size_t MeshInterpolation_Base<Dimension>::GetMemorySize() const
  {
    size_t taille = proj.GetMemorySize() + sizeof(*this);
    taille += proj_surf.GetMemorySize();
    taille += ref_surf.GetMemorySize() + Seldon::GetMemorySize(Nodle) +
      Seldon::GetMemorySize(points_div)+ Seldon::GetMemorySize(points_div_surf)
      + Seldon::GetMemorySize(weights_surf) + step_subdiv.GetMemorySize() + 
      Seldon::GetMemorySize(file_name) +  NumElement.GetMemorySize() + 
      LocalPositionBoundary.GetMemorySize() + TypeBoundary.GetMemorySize();

    taille += Seldon::GetMemorySize(NormaleBoundary) + Seldon::GetMemorySize(PointsBoundary)
      + Seldon::GetMemorySize(PointsNodalBoundary) + Seldon::GetMemorySize(WeightsBoundary)
      + Seldon::GetMemorySize(DsBoundary) + Seldon::GetMemorySize(DFiBoundary)
      + Seldon::GetMemorySize(invDFiBoundary) + Seldon::GetMemorySize(DFiNodalBoundary)
      + Seldon::GetMemorySize(invDFiNodalBoundary) + AllPoints.GetMemorySize()
      + AllNormale.GetMemorySize() + AllWeights.GetMemorySize();
    
    taille += nb_points_quadrature_per_proc.GetMemorySize() + offset_quadrature_per_proc.GetMemorySize();
    return taille;
  }
  

  //! the interpolation mesh is obtained by a regular subdivision
  template<class Dimension>
  void MeshInterpolation_Base<Dimension>::SetRegularSubdivisions(int r)
  {
    nb_subdivisions = r;
    step_subdiv.Reallocate(r+1);
    step_subdiv.Fill();
    Mlt(Real_wp(1)/r, step_subdiv);
    SetRegularPoints(r);
  }
  
  
  //! the interpolation mesh is obtained with Gauss-Lobatto points as subdivisions
  template<>
  void MeshInterpolation_Base<Dimension2>::SetLobattoSubdivisions(int r)
  {
    nb_subdivisions = r;
    ComputeGaussLobatto(step_subdiv, weights_surf(0), r);
    points_div_surf.Reallocate(1);
    points_div_surf(0) = step_subdiv;
  }
  
    
  //! the interpolation mesh is obtained with Gauss-Lobatto points as subdivisions
  template<>
  void MeshInterpolation_Base<Dimension3>::SetLobattoSubdivisions(int r)
  {
    nb_subdivisions = r;
    ComputeGaussLobatto(step_subdiv, weights_surf(0), r);
    points_div_surf.Reallocate(2);
    points_div_surf(1).Reallocate((r+1)*(r+1));
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r; j++)
	points_div_surf(1)(i*r + j).Init(step_subdiv(i), step_subdiv(j));    
  }
  
    
  //! the interpolation mesh is obtained by a regular subdivision  
  template<>
  void MeshInterpolation_Base<Dimension2>::SetRegularPoints(int r)
  {
    points_div_surf.Reallocate(1);
    points_div_surf(0) = step_subdiv;
  }

  
  //! the interpolation mesh is obtained by a regular subdivision
  template<>
  void MeshInterpolation_Base<Dimension3>::SetRegularPoints(int r)
  {
    points_div_surf.Reallocate(2);
    points_div_surf(0).Reallocate((r+1)*(r+2)/2);
    points_div_surf(1).Reallocate((r+1)*(r+1));
    int nbt = 0;
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r; j++)
        {
          points_div_surf(1)(i*r + j).Init(step_subdiv(i), step_subdiv(j));
          if (i+j <= r)
            {
              points_div_surf(0)(nbt).Init(step_subdiv(i), step_subdiv(j));
              nbt++;
            }
        }
  }
  

  //! interpolation on Gauss points is required
  template<>
  void MeshInterpolation_Base<Dimension2>::SetGaussQuadrature(int r, int type)
  {
    nb_subdivisions = r;
    if (type == Globatto<Real_wp>::QUADRATURE_GAUSS)
      ComputeGaussLegendre(step_subdiv, weights_surf(0), r);
    else
      ComputeGaussLobatto(step_subdiv, weights_surf(0), r);
    
    points_div_surf.Reallocate(1);
    points_div_surf(0) = step_subdiv;
  }
  

  //! interpolation on Gauss points is required
  template<>
  void MeshInterpolation_Base<Dimension3>::SetGaussQuadrature(int r, int type)
  {
    points_div_surf.Reallocate(2);
    
    nb_subdivisions = r;
    // quadrature points for the square
    VectReal_wp weights1d; Matrix<int> NumNodes2D, coor;
    if (type == Globatto<Real_wp>::QUADRATURE_GAUSS)
      ComputeGaussLegendre(step_subdiv, weights1d, r);
    else
      ComputeGaussLobatto(step_subdiv, weights1d, r);
        
    MeshNumbering<Dimension2>::ConstructQuadrilateralNumbering(r, NumNodes2D, coor);
    points_div_surf(1).Reallocate((r+1)*(r+1));
    weights_surf(1).Reallocate((r+1)*(r+1));
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r; j++)
        {
          points_div_surf(1)(NumNodes2D(i, j)).Init(step_subdiv(i), step_subdiv(j));
          weights_surf(1)(NumNodes2D(i, j)) = weights1d(i)*weights1d(j);
        }
    
    // quadrature points for the unit triangle
    TriangleQuadrature::ConstructQuadrature(2*r, points_div_surf(0), weights_surf(0));
  }
  

  //! extracts a surface mesh and computes normales, jacobian matrices, etc
  void MeshInterpolation<Dimension2>::
  ComputeSurfaceMesh(const IVect& ref_cond, const Mesh<Dimension2>& mesh,
                     Mesh<Dimension2>& mesh_subdiv, int ref_domain)
  {
    SurfacicMesh<Dimension2> mesh_surf;
    
    // extraction of the surfacic mesh
    mesh.GetBoundaryMesh(1, mesh_surf, ref_cond, ref_domain);
    
    // some local variables
    int num_elem, num_loc;
    
    // subdivision of the initial mesh
    Vector<IVect> NodleVol, NodleSurf;
    mesh_subdiv = mesh_surf;
    mesh_subdiv.SubdivideMesh(step_subdiv, points_div_surf, points_div, NodleVol, NodleSurf);
    int nb_points = step_subdiv.GetM();
    
    // copmputation of the normale, projection of points
    int N = mesh_surf.ListeBoundaries.GetM();
    NumElement.Reallocate(N);
    LocalPositionBoundary.Reallocate(N);        
    Nodle.Reallocate(N);
    TypeBoundary.Reallocate(N); TypeBoundary.Fill(0);    
    NormaleBoundary.Reallocate(N);
    PointsBoundary.Reallocate(N);
        
    // loop on each face
    for (int i1 = 0; i1 < N; i1++)
      {
        //int i = mesh_surf.ListeBoundaries(i1);
        //num_boundary = i;
        num_elem = mesh_surf.NumElement(i1);
        
        NumElement(i1) = num_elem;
        num_loc = mesh_surf.NumLocalBoundary(i1);
        LocalPositionBoundary(i1) = num_loc;
                
	Nodle(i1).Reallocate(nb_points);
	if (mesh.Element(num_elem).GetOrientationEdge(num_loc))
	  {
	    for (int j = 0; j < nb_points; j++)
	      Nodle(i1)(j) = NodleSurf(i1)(j);
	  }
	else
	  {
	    for (int j = 0; j < nb_points; j++)
	      Nodle(i1)(j) = NodleSurf(i1)(nb_points-1-j);
	  }
	
        this->StoreGeometricDataSurface(mesh, i1, num_elem, num_loc);
      }
  }
  

  //! internal function
  void MeshInterpolation<Dimension2>::
  StoreGeometricDataSurface(const Mesh<Dimension2>& mesh, int i1, int num_elem, int num_loc)
  {
    // points and jacobian matrices
    SetPoints<Dimension2> pts;
    SetMatrices<Dimension2> mat;
    VectR2 s;    

    // computation of transformations Fi and DFi
    const ElementGeomReference<Dimension2>& Fb = mesh.GetReferenceElement(num_elem);
    mesh.GetVerticesElement(num_elem, s);
    Fb.FjElemNodal(s, pts, mesh, num_elem);

    // loop on subdivided points
    int nb_points = step_subdiv.GetM();
    PointsBoundary(i1).Reallocate(nb_points);
    NormaleBoundary(i1).Reallocate(nb_points);
       
    R2 point_loc, point_glob, normale_fj;
    Matrix2_2 mat_dfj, dfjm1; Real_wp dsj;
    for (int j = 0; j < nb_points; j++)
      {
        Fb.GetLocalCoordOnBoundary(num_loc, step_subdiv(j), point_loc);
        
        Fb.Fj(s, pts, point_loc, point_glob, mesh, num_elem);
        Fb.DFj(s, pts, point_loc, mat_dfj, mesh, num_elem);
        GetInverse(mat_dfj, dfjm1);
        Fb.GetNormale(dfjm1, normale_fj, dsj, num_loc);
            	    
        PointsBoundary(i1)(j) = point_glob;
        NormaleBoundary(i1)(j) = normale_fj;
      }
  }
  
  
#ifdef MONTJOIE_WITH_THREE_DIM
  //! extracts a surface mesh and computes normales, jacobian matrices, etc
  void MeshInterpolation<Dimension3>::
  ComputeSurfaceMesh(const IVect& ref_cond, const Mesh<Dimension3>& mesh,
                     Mesh<Dimension3>& mesh_subdiv, int ref_domain)
  {
    SurfacicMesh<Dimension3> mesh_surf;
    
    // extraction of the surfacic mesh
    mesh.GetBoundaryMesh(1, mesh_surf, ref_cond, ref_domain);
    
    // some local variables
    int num_elem, num_loc;
    
    // subdivision of the initial mesh
    Vector<IVect> NodleVol, NodleSurf;
    Matrix<int> FacesPoints_Rotation_Tri, FacesPoints_Rotation_Quad;
    mesh_subdiv = mesh_surf;
    mesh_subdiv.SubdivideMesh(step_subdiv, points_div_surf, points_div,
			      NodleVol, NodleSurf);
    
    // rotation of points        
    MeshNumbering<Dimension3>::
      GetRotationTriangularFace(points_div_surf(0), FacesPoints_Rotation_Tri);    
    
    MeshNumbering<Dimension3>::
      GetRotationQuadrilateralFace(points_div_surf(1), FacesPoints_Rotation_Quad);
    
    // computation of the normale, projection of points
    int N = mesh_surf.ListeBoundaries.GetM();
    NumElement.Reallocate(N);
    LocalPositionBoundary.Reallocate(N);    
    PointsBoundary.Reallocate(N);
    NormaleBoundary.Reallocate(N);
    Nodle.Reallocate(N);
    TypeBoundary.Reallocate(N); TypeBoundary.Fill(0);
    
    // loop on each face
    for (int i1 = 0; i1 < N; i1++)
      {
	int i = mesh_surf.ListeBoundaries(i1);
        num_elem = mesh_surf.NumElement(i1);
        TypeBoundary(i1) = mesh.BoundaryRef(i).GetHybridType();
        bool tri = (mesh.BoundaryRef(i).GetNbVertices() == 3);
        NumElement(i1) = num_elem;
        num_loc = mesh_surf.NumLocalBoundary(i1);
        LocalPositionBoundary(i1) = num_loc;
        
	int nb_points = points_div_surf(TypeBoundary(i1)).GetM();
	int rot = mesh.Element(num_elem).GetOrientationFace(num_loc);
	
	// we take opposite orientation
	rot = MeshNumbering<Dimension3>::GetOppositeOrientationFace(rot, tri);
        
	Nodle(i1).Reallocate(nb_points);
	for (int j = 0; j < nb_points; j++)
	  {
	    int j1 = -1;
	    if (tri)
	      j1 = FacesPoints_Rotation_Tri(rot, j);
	    else
	      j1 = FacesPoints_Rotation_Quad(rot, j);
	    
	    Nodle(i1)(j) = NodleSurf(i1)(j1);
	  }
        
        this->StoreGeometricDataSurface(mesh, i1, num_elem, num_loc);
      }
  }
  

    
  //! internal function
  void MeshInterpolation<Dimension3>::
  StoreGeometricDataSurface(const Mesh<Dimension3>& mesh, int i1, int num_elem, int num_loc)
  {
    // points and jacobian matrices
    SetPoints<Dimension3> pts;
    SetMatrices<Dimension3> mat;
    VectR3 s;

    // computation of transformations Fi and DFi    
    mesh.GetVerticesElement(num_elem, s);
    const ElementGeomReference<Dimension3>& Fb = mesh.GetReferenceElement(num_elem);
    Fb.FjElemNodal(s, pts, mesh, num_elem);
        
    // loop on subdivided points
    int nb_points = points_div_surf(TypeBoundary(i1)).GetM();
    PointsBoundary(i1).Reallocate(nb_points);
    NormaleBoundary(i1).Reallocate(nb_points);
    Nodle(i1).Reallocate(nb_points);

    R3 point_loc, point_glob, normale_fj;
    Matrix3_3 mat_dfj, dfjm1; Real_wp dsj;
    for (int j = 0; j < nb_points; j++)
      {
        Fb.GetLocalCoordOnBoundary(num_loc,
                                   points_div_surf(TypeBoundary(i1))(j), point_loc);
        
        Fb.Fj(s, pts, point_loc, point_glob, mesh, num_elem);
        Fb.DFj(s, pts, point_loc, mat_dfj, mesh, num_elem);
        GetInverse(mat_dfj, dfjm1);
        Fb.GetNormale(dfjm1, normale_fj, dsj, num_loc);
            	    
        PointsBoundary(i1)(j) = point_glob;
        NormaleBoundary(i1)(j) = normale_fj;
      }
  }  
#endif

#ifdef SELDON_WITH_MPI
  //! quadrature points for each processor are gathered,
  //! so that each processor knows all the points
  template<class Dimension>
  void MeshInterpolation_Base<Dimension>::GatherQuadraturePoints(const MPI_Comm& comm, bool assemble)
#else
  template<class Dimension>
  void MeshInterpolation_Base<Dimension>::GatherQuadraturePoints(bool assemble)
#endif
  {
#ifdef SELDON_WITH_MPI
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    int rank_proc; MPI_Comm_rank(comm, &rank_proc);
    rank_processor = rank_proc;
#else
    int nb_proc(1), rank_proc(0);
#endif
    
    int N = 0;
    for (int i = 0; i < this->GetNbBoundary(); i++)
      N += this->GetNbPointsQuadrature(i);
    
    nb_points_all_quadrature = N;
    Vector<R_N> PointsProc(N), NormaleProc(N);
    VectReal_wp WeightsProc(N);
    
    N = 0;
    for (int i = 0; i < this->GetNbBoundary(); i++)
      {
        int nb_points = this->GetNbPointsQuadrature(i);
        for (int j = 0; j < nb_points; j++)
          {
            PointsProc(N+j) = PointsBoundary(i)(j);
            NormaleProc(N+j) = NormaleBoundary(i)(j);
            WeightsProc(N+j) = WeightsBoundary(i)(j);
          }
        
        N += nb_points;
      }
    
    //DISP(nb_proc); DISP(N);
    if ((nb_proc > 1) && (assemble))
      {
#ifdef SELDON_WITH_MPI
        // we retrieve the total number of quadrature points for all processors
        int nb_all_proc = 0;
        MPI_Allreduce(&nb_points_all_quadrature, &nb_all_proc, 1,
                      MPI_INTEGER, MPI_SUM, comm);
        
        // allocating arrays containing all quadrature points
        AllPoints.Reallocate(nb_all_proc);
        AllNormale.Reallocate(nb_all_proc);
        AllWeights.Reallocate(nb_all_proc);
        
        nb_points_quadrature_per_proc.Reallocate(nb_proc);
        offset_quadrature_per_proc.Reallocate(nb_proc+1);
        MPI_Allgather(&N, 1, MPI_INTEGER, nb_points_quadrature_per_proc.GetData(),
                      1, MPI_INTEGER, comm);
        
        offset_quadrature_per_proc(0) = 0;
        for (int i = 0; i < nb_proc; i++)
          offset_quadrature_per_proc(i+1) = offset_quadrature_per_proc(i)
            + nb_points_quadrature_per_proc(i);
        
        int offset = offset_quadrature_per_proc(rank_proc);
        for (int j = 0; j < N; j++)
          {
            AllPoints(offset+j) = PointsProc(j);
            AllNormale(offset+j) = NormaleProc(j);
            AllWeights(offset+j) = WeightsProc(j);
          }
        
        Vector<int64_t> xtmp;
        for (int i = 0; i < nb_proc; i++)
          if (nb_points_quadrature_per_proc(i) > 0)
            {
              MpiBcast(comm, reinterpret_cast<Real_wp*>(&AllPoints(offset_quadrature_per_proc(i))),
                       xtmp, Dimension::dim_N*nb_points_quadrature_per_proc(i), i);
              
              MpiBcast(comm, reinterpret_cast<Real_wp*>(&AllNormale(offset_quadrature_per_proc(i))),
                       xtmp, Dimension::dim_N*nb_points_quadrature_per_proc(i), i);
              
              MpiBcast(comm, &AllWeights(offset_quadrature_per_proc(i)), xtmp,
                       nb_points_quadrature_per_proc(i), i);
            }        
#endif
      }
    else
      {
        nb_points_quadrature_per_proc.Reallocate(nb_proc);
        offset_quadrature_per_proc.Reallocate(nb_proc+1);
        nb_points_quadrature_per_proc(rank_proc) = PointsProc.GetM();
        offset_quadrature_per_proc.Fill(0);
        for (int i = rank_proc+1; i <= nb_proc; i++)
          offset_quadrature_per_proc(i) = PointsProc.GetM();
        
        AllPoints = PointsProc;
        AllNormale = NormaleProc;
        AllWeights = WeightsProc;
      }
                  
  }
    
} // namespace Montjoie

#define MONTJOIE_FILE_MESH_INTERPOLATION_CXX
#endif
