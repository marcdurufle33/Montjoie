#ifndef MONTJOIE_FILE_GENERALIZED_IMPEDANCE_MODEL_CXX

namespace Montjoie
{
  
  //! Sets values of attributes to default values
  void VarGeneralizedImpedance_Base::InitDefaultValues()
  {
    nodl_mesh = 0;
    order_gibc = 1;
    nb_ddl_volume_on_surface = 0;
    pade_approx_sqrt = false;
    epsilon_pade = 0.0;
    use_unsymmetric_algo = true;
  }
  

  size_t VarGeneralizedImpedance_Base::GetMemorySize() const
  {
    size_t taille = sizeof(*this) + gibc_references.GetMemorySize()
      + ref_cond.GetMemorySize() + DdlVol.GetMemorySize()
      + Seldon::GetMemorySize(NodleVol) + MatchingNumber_Subdomain.GetMemorySize()
      + Seldon::GetMemorySize(MatchingDofOrig_Subdomain) + CoefAl.GetMemorySize();
    
    return taille;
  }
  
    
  //! modifying gibc conditions with a line of the data file
  void VarGeneralizedImpedance_Base
  ::SetInputData(const string& description_field, const VectString& parameters)
  {
    if (description_field == "GibcABC")
      {
        if (parameters(0) == "PADE")
          {
            if (parameters.GetM() <= 3)
              {
                cout << "Not enough parameters in GibcABC" << endl;
                cout << "Write for instance GibcAbc = PADE alpha n epsilon" << endl;
                cout << "Current parameters are " << endl << parameters << endl;
                abort();
              }
            
            pade_approx_sqrt = true;
            Real_wp alpha = to_num<Real_wp>(parameters(1))*pi_wp/180;
            int n = to_num<int>(parameters(2));
            epsilon_pade = to_num<Real_wp>(parameters(3));
            ComputePadeCoefficientsSqrt(alpha, n, CoefC0, CoefAl, CoefBl);
          }
      }
    else if (description_field == "GibcSymmetric")
      {
        if (parameters(0) == "YES")
          use_unsymmetric_algo = false;
        else
          use_unsymmetric_algo = true;
      }
  }
  

  //! releases memory used by the object
  void VarGeneralizedImpedance_Base::Clear()
  {
    DdlVol.Clear();
    MatchingNumber_Subdomain.Clear();
    MatchingDofOrig_Subdomain.Clear();
  }

  
  //! sets column numbers for impedance boundary condition
  void VarGeneralizedImpedance_Base::SetModifiedColNumbers(int num)
  {
    int Nvol = var_problem.offset_dof_unknown(1);
    int Nsurf = this->DdlVol.GetM();    
    int Nvol_mesh = nodl_mesh;
    
    for (int m = 0; m < var_problem.nb_unknowns; m++)
      for (int i = 0; i < DdlVol.GetM(); i++)
        var_boundary.NewColumnNumbers_Impedance(this->DdlVol(i) + m*Nvol)
          = m*Nvol + Nvol_mesh + num*Nsurf + i;
  }
  

  //! sets row numbers for impedance boundary condition  
  void VarGeneralizedImpedance_Base::SetModifiedRowNumbers(int num)
  {
    int Nvol = var_problem.offset_dof_unknown(1);
    int Nsurf = this->DdlVol.GetM();    
    int Nvol_mesh = nodl_mesh;
    
    for (int m = 0; m < var_problem.nb_unknowns; m++)
      for (int i = 0; i < DdlVol.GetM(); i++)
        var_boundary.NewRowNumbers_Impedance(this->DdlVol(i) + m*Nvol)
          = m*Nvol + Nvol_mesh + num*Nsurf + i;
  }
  
  
  //! internal function used to retrieve dofs on a boundary
  template<class Dimension>
  void VarGeneralizedImpedance_Base::
  FindLocalDofs(int i, int ne, int num_loc, VarProblem<Dimension>& var,
                IVect& IndexDof, const ElementReference_Dim<Dimension>& Fb)
  {
    IVect Nodle = var.GetDofNumberOnElement(ne, 0);
    int nb_dof_face = Fb.GetNbDofBoundary(num_loc);
    NodleVol(i).Reallocate(nb_dof_face);
    for (int j = 0; j < nb_dof_face; j++)
      {
        int num_dof_loc = Fb.GetLocalNumber(num_loc, j);
        int num_dof = Nodle(num_dof_loc);
        NodleVol(i)(j) = num_dof;
        if (num_dof >= 0)
          IndexDof(num_dof) = 1;
      }
  }

  
  //! initialisation of impedance condition
  template<class Dimension>
  void VarGeneralizedImpedance_Base
  ::FindDofsGibc(VarProblem<Dimension>& var, const IVect& Epart)
  {
    ref_cond.Reallocate(var.mesh.GetNbReferences()+1);
    ref_cond.Fill(0);
    
    if (gibc_references.GetM() <= 0)
      return;
    
    // counting the number of edges/faces in Gamma
    for (int i = 0; i < gibc_references.GetM(); i++)
      ref_cond(gibc_references(i)) = 1;
    
    int nb_faces = 0;
    for (int i = 0; i < var.mesh.GetNbBoundaryRef(); i++)
      if (ref_cond(var.mesh.BoundaryRef(i).GetReference()) == 1)
        nb_faces++;
    
    // fills NodleVol and IndexDof
    nodl_mesh = var.GetMeshNumbering(0).GetNbDof();
    NodleVol.Reallocate(nb_faces);
    IVect IndexDof(nodl_mesh);
    IndexDof.Fill(0);
    nb_faces = 0;
    for (int i = 0; i < var.mesh.GetNbBoundaryRef(); i++)
      if (ref_cond(var.mesh.BoundaryRef(i).GetReference()) == 1)
        {
          int num_face = i;
          int ne = var. mesh.Boundary(num_face).numElement(0);	  
          if (ne >= 0)
            {
              int num_loc = var.mesh.Element(ne).GetPositionBoundary(num_face);
	      const ElementReference_Dim<Dimension>& Fb = var.GetReferenceElement(ne);
	      FindLocalDofs(nb_faces, ne, num_loc, var, IndexDof, Fb);
	    }
	  
	  nb_faces++;
        }
    
    // counting the number of dofs on Gamma
    nb_ddl_volume_on_surface = 0;
    for (int i = 0; i < IndexDof.GetM(); i++)
      if (IndexDof(i) > 0)
        nb_ddl_volume_on_surface++;
    
    // filling DdlVol
    DdlVol.Reallocate(nb_ddl_volume_on_surface);
    nb_ddl_volume_on_surface = 0;
    for (int i = 0; i < IndexDof.GetM(); i++)
      if (IndexDof(i) > 0)
        DdlVol(nb_ddl_volume_on_surface++) = i;
    
    /* if (nb_proc > 1)
      {
        int nb_edges_ref = NumFaceMinus.GetM();
        NodleVolMinus.Reallocate(nb_edges_ref);
        NodleVolPlus.Reallocate(nb_edges_ref);
        ProcVolMinus.Reallocate(nb_edges_ref);
        ProcVolPlus.Reallocate(nb_edges_ref);
        for (int i = 0; i < nb_edges_ref; i++)
          {
            int nb_dof_edge = NodleSurf(i).GetM();
            NodleVolMinus(i).Reallocate(nb_dof_edge);
            NodleVolPlus(i).Reallocate(nb_dof_edge);
            int ne_minus = mesh.BoundaryRef(NumFaceMinus(i)).numElement(0);
            int ne_plus = mesh.BoundaryRef(NumFacePlus(i)).numElement(0);
            ProcVolMinus(i) = Epart(ne_minus);
            ProcVolPlus(i) = Epart(ne_plus);
            for (int j = 0; j < nb_dof_edge; j++)
              {
                int num_dof = NodleSurf(i)(j);
                NodleVolMinus(i)(j) = DdlVolMinus(num_dof);
                NodleVolPlus(i)(j) = DdlVolPlus(num_dof);
	      }
	  }
          }     */
  }  
  
} // end namespace

#define MONTJOIE_FILE_GENERALIZED_IMPEDANCE_MODEL_CXX
#endif

  
