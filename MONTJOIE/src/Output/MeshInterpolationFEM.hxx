#ifndef MONTJOIE_FILE_MESH_INTERPOLATION_FEM_HXX

namespace Montjoie
{

  class MeshInterpolationFEM_Base
  {
  protected:
    template<class T>
    static void SetContribU(const Vector<Vector<TinyVector<T, 1> > >& u,
                           int nb_points_quad, const IVect& offset, Vector<Vector<T> >& u_boundary);

    template<class T>
    static void SetContribU(const Vector<Vector<TinyVector<T, 2> > >& u,
                           int nb_points_quad, const IVect& offset, Vector<Vector<T> >& u_boundary);

    template<class T>
    static void SetContribU(const Vector<Vector<TinyVector<T, 3> > >& u,
                           int nb_points_quad, const IVect& offset, Vector<Vector<T> >& u_boundary);
    
  };
  
  //! base class to do interpolation on subdivided meshes
  template<class Dimension>
  class MeshInterpolationFEM : public MeshInterpolation<Dimension>,
    public MeshInterpolationFEM_Base
  {
  protected :
    typedef typename Dimension::DimensionBoundary DimensionBoundary;
    typedef typename Dimension::MatrixN_N MatrixN_N;
    typedef typename Dimension::R_N R_N;
    
  public :

    void StoreGeometricDataSurface(const DistributedProblem<Dimension>& var,
				   const Mesh<Dimension>& mesh,
				   int i1, int num_elem, int num_loc);
    
    void ComputeSurfaceMesh(const IVect& ref_cond, const Mesh<Dimension>& mesh,
                            Mesh<Dimension>& mesh_subdiv,
			    const DistributedProblem<Dimension>& var, int ref_domain = 0);
    
    template<class Complexe>
    void ComputeEnHnOnBoundary(const VarProblem<Dimension>& var, const Vector<Complexe>& U0,
                               Vector<Complexe>& trace_En, Vector<Complexe>& trace_Hn,
                               bool assemble = true, bool compute_H = true) const;

    template<class T>
    void AddSourceBoundary(const VarProblem<Dimension>& var,
                           const Vector<T>& g_source, Vector<T>& rhs) const;
    
  };
  
} // namespace Montjoie

#define MONTJOIE_FILE_MESH_INTERPOLATION_FEM_HXX
#endif
