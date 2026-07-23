#ifndef MONTJOIE_FILE_VAR_PROBLEM_1D_INLINE_CXX

namespace Montjoie
{
  
  /*****************
   * VarProblem_1D *
   *****************/

  
  inline void VarProblem_1D::GetMiomega(Real_wp& m_iomega) 
  {
    m_iomega = 1.0;
  }
  
  
  inline void VarProblem_1D::GetMiomega(Complex_wp& m_iomega)
  {
    m_iomega = Complex_wp(0, -omega); 
  }
  
  
  inline void VarProblem_1D::SetOmega(const Real_wp& om)
  {
    omega = om;
    omega2 = om*om;
  }
    
  
  inline bool VarProblem_1D::FirstOrderFormulation() const
  {
    return mixed_formulation; 
  }
  
  
  inline void VarProblem_1D::SetFirstOrderFormulation(bool first)
  {
    mixed_formulation = first;
  }

  
  inline Real_wp VarProblem_1D::GetXmin() const
  {
    return this->xmin_;
  }
  
  
  inline Real_wp VarProblem_1D::GetXmax() const
  {
    return this->xmax_; 
  }
    
  
  inline int VarProblem_1D::GetNbDofBoundaries(int i)
  {
    return 2;
  }
  
  
  inline int VarProblem_1D::GetNbDof() const
  {
    return this->nodl;
  }
  
  
  inline int VarProblem_1D::GetNbRows() const
  {
    return this->nodl;
  }
  
  
  inline int VarProblem_1D::GetNbElt() const
  {
    return this->mesh.GetNbElt();
  }
  
  
  inline int VarProblem_1D::GetPrintLevel() const
  {
    return this->print_level;
  }
  

  inline VectString& VarProblem_1D::GetMeshData(int n)
  {
    return mesh_data(n);
  }


  inline void VarProblem_1D::SetMeshData(const Vector<VectString>& param)
  {
    mesh_data = param; 
  }
  
  
  inline GridInterpolation<Dimension1>& 
  VarProblem_1D::GetSectionGrid()
  {
    return var_section; 
  }
  
  
  inline const GridInterpolation<Dimension1>& 
  VarProblem_1D::GetSectionGrid() const 
  {
    return var_section; 
  }
  
  
  inline ParamOutputClass& VarProblem_1D::GetOutputParameters()
  {
    return output_section_param; 
  }
  
  
  inline const ParamOutputClass& VarProblem_1D::GetOutputParameters() const
  {
    return output_section_param; 
  }

  
  inline const VectReal_wp& VarProblem_1D::GetCoordinateDof() const
  {
    return CoorDofs; 
  }
  
  
  inline const Real_wp& VarProblem_1D::GetCoordinateDof(int i) const
  {
    return CoorDofs(i); 
  }
  
  
  inline const IVect& VarProblem_1D::GetOverlappedDofNumber() const
  {
    return OverlapDofNumber_Subdomain; 
  }


  template<class MatrixSparse>
  inline void VarProblem_1D::SetPeriodicCondition(MatrixSparse& A) 
  {
  }
  
  
  template<class MatrixSparse>
  inline void VarProblem_1D::SetDirichletCondition(MatrixSparse& A)
  {
  }


  inline int VarProblem_1D::GetLeftBoundaryCondition() const
  {
    return condition_left_side; 
  }
  
  
  inline int VarProblem_1D::GetRightBoundaryCondition() const
  {
    return condition_right_side; 
  }

  
  inline bool VarProblem_1D::EnvelopeToCompute() const
  {
    return calcul_enveloppe; 
  }
  
  
  inline void VarProblem_1D::SetBoundaryCondition(int left_cond, int right_cond)
  {
    condition_left_side = left_cond;
    condition_right_side = right_cond; 
  }

  
  inline void VarProblem_1D::SetPulsation(const Real_wp& om)
  { 
    omega = om; 
    omega2 = om*om;
  }


  inline void VarProblem_1D
  ::GetInternalNodesElement(int i, int nb_dof_loc, int& nb_dof_edges, int& nb_dof_int,
                            Vector<int>& intern_node) const
  {
    abort();
  }

}

#define MONTJOIE_FILE_VAR_PROBLEM_1D_INLINE_CXX
#endif



