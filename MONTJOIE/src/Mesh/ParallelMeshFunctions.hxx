#ifndef MONTJOIE_FILE_PARALLEL_MESH_FUNCTIONS_HXX

namespace Montjoie
{

  template<class Dimension>
  struct ParamParallelMesh
  {
    IVect num;
    VectReal_wp all_param;
    typename Dimension::VectR_N PointsEdgeRef, TranslatPeriodicEdge;
    IVect NumBoundary, NumDofNeighbor, NumElem, NumDof;
    Vector<IVect> MinimalProc; IVect OffsetDofV, NumDofPML;
    int dg_form; IVect nodl_mesh;
    Vector<IVect> offset_matching;
    
  };

  template<class Dimension>
  void SendMeshToProcessor(const Mesh<Dimension>& glob_mesh,
                           const Vector<MeshNumbering<Dimension>* >& glob_mesh_num,
			   const IVect& NumElement, const IVect& Epart, const IVect& NumLoc,
                           Mesh<Dimension>& mesh, int proc, const MPI_Comm& comm, ParamParallelMesh<Dimension>& param, int tag = 48);                           

  template<class Dimension>
  void RecvMeshFromProcessor(Mesh<Dimension>& mesh, Vector<MeshNumbering<Dimension>* >& mesh_num,
			     Vector<IVect>& ConnecEdge, IVect& MatchingProc, Vector<IVect>& MatchingDofs,
                             ParamParallelMesh<Dimension>& param, int proc, MPI_Comm& comm, int tag = 48);
  
}

#define MONTJOIE_FILE_PARALLEL_MESH_FUNCTIONS_HXX
#endif
