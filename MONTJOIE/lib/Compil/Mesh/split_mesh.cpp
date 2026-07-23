
#include "MontjoieFlag.hxx"

#include "Mesh/MontjoieMeshHeader.hxx"
#include "Mesh/MontjoieMeshInline.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#ifdef SELDON_WITH_MPI
#include "Mesh/ParallelMeshFunctions.cxx"
#endif
#endif

namespace Montjoie
{

#ifdef SELDON_WITH_MPI
  
  SELDON_EXTERN template void SendMeshToProcessor(const Mesh<Dimension2>& glob_mesh,
                                                  const Vector<MeshNumbering<Dimension2>* >& glob_mesh_num,
                                                  const IVect& NumElement, const IVect& Epart, const IVect& NumLoc,
                                                  Mesh<Dimension2>& mesh, int proc, const MPI_Comm& comm, ParamParallelMesh<Dimension2>& Param, int tag);

  SELDON_EXTERN template void SendMeshToProcessor(const Mesh<Dimension3>& glob_mesh,
                                                  const Vector<MeshNumbering<Dimension3>* >& glob_mesh_num,
                                                  const IVect& NumElement, const IVect& Epart, const IVect& NumLoc,
                                                  Mesh<Dimension3>& mesh, int proc, const MPI_Comm& comm,
                                                  ParamParallelMesh<Dimension3>& Param, int tag);
  
  SELDON_EXTERN template void RecvMeshFromProcessor(Mesh<Dimension2>& mesh, Vector<MeshNumbering<Dimension2>* >& mesh_num,
                                                    Vector<IVect>& ConnecEdge, IVect& MatchingProc, Vector<IVect>& MatchingDofs,
                                                    ParamParallelMesh<Dimension2>& Param, int proc, MPI_Comm& comm, int tag);
  
  SELDON_EXTERN template void RecvMeshFromProcessor(Mesh<Dimension3>& mesh, Vector<MeshNumbering<Dimension3>* >& mesh_num,
                                                    Vector<IVect>& ConnecFace, IVect& MatchingProc, Vector<IVect>& MatchingDofs,
                                                    ParamParallelMesh<Dimension3>& Param, int proc, MPI_Comm& comm, int tag);
#endif
  
}
