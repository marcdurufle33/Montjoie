

// version simplifiee de ConstructFiniteElement (dans VarProblem.cxx)
void ModeEs_Solver::ConstructFiniteElement(const string& name_elt)
{
  // previous finite elements are cleared if present
  this->ClearFiniteElement();

  // orders present in the mesh
  TinyVector<IVect, 4> order;
  this->mesh_num.GetOrder(order);

  Vector<bool> change_elt(this->mesh.GetNbElt());
  change_elt.Fill(true);

  this->AddFiniteElement(name_elt, order, change_elt, this->mesh, this->mesh_num,
                         0, false);

  this->mesh_num.GetOrderQuadrature(order);
  RemoveDuplicate(order(0));

  int type_integration_edge, type_integration_tri, type_integration_quad;
  Vector<MeshNumbering<Dimension2>* > other_mesh_num;
  this->GetTypeIntegrationBoundary(this->mesh, this->mesh_num, other_mesh_num, type_integration_edge,
           type_integration_tri, type_integration_quad);

  this->mesh_num.number_map.ConstructQuadrature2D(order, type_integration_edge);
}

// traitement de Dirichlet, quasi-periodique et ddl conductors
// À COMMENTER POUR BIEN COMPRENDRE
// Certaines choses sont enlevées (bords, conducteurs notamment),
 // les rechercher des Cable_Droit/Mesh.cxx si besoin
void ModeEs_Solver::RenumberDofs()
{
  Vector<bool> is_dirichlet(this->mesh.GetNbReferences()+1);
  is_dirichlet.Fill(false);
  for (int i = 0; i < this->ref_bord.GetM(); i++)
    is_dirichlet(ref_bord(i)) = true;

  // on met le ddl 0 en Dirichlet (dans le cas quasi-periodique avec m <> 0)
  Vector<int> DofDirichlet(this->mesh_num.GetNbDof());
  DofDirichlet.Fill(-1);

  int nb_dof_dirichlet = 0;


  // on renumerote
  if (nb_dof_dirichlet > 0)
    {
      Vector<int> IndexDof(this->mesh_num.GetNbDof());
      IndexDof.Fill(-1);
      int nodl = 0;
      for (int i = 0; i < this->mesh_num.GetNbDof(); i++)
        {
          if (DofDirichlet(i) == -1)
            IndexDof(i) = nodl++;
        }

      //DISP(IndexDof);

      for (int i = 0; i < this->mesh.GetNbElt(); i++)
        {
          //DISP(i); DISP(this->mesh_num.Element(i).GetNodle());
          for (int j = 0; j < this->mesh_num.Element(i).GetNbDof(); j++)
            {
              int num_dof = this->mesh_num.Element(i).GetNumberDof(j);
              if (IndexDof(num_dof) <= -2)
                {
                  this->mesh_num.Element(i).SetNumberDof(j, -IndexDof(num_dof)-2);
                }
              else
                this->mesh_num.Element(i).SetNumberDof(j, IndexDof(num_dof));
            }

          // DISP(this->is_dof_with_phase(i));
          //DISP(this->mesh_num.Element(i).GetNodle());
        }

      this->mesh_num.SetNbDof(nodl);
    }

  // boucle sur les aretes de bord pour lister les ddls Dirichlets
  nb_dof_dirichlet = 0;
  DofDirichlet.Fill(-1);
  for (int i = 0; i < this->mesh.GetNbBoundaryRef(); i++)
    {
      int ref = this->mesh.BoundaryRef(i).GetReference();
      if (is_dirichlet(ref))
        {
          int num_elem = this->mesh.BoundaryRef(i).numElement(0);
          int num_loc = this->mesh.Element(num_elem).GetPositionBoundary(i);
          const ElementReference<Dimension2, 1>& Fb = this->GetReferenceElementH1(num_elem);
          const IVect& Nodle = this->mesh_num.Element(num_elem).GetNodle();
          int nb_dof = Fb.GetNbDofBoundary(num_loc);
          for (int j = 0; j < nb_dof; j++)
            {
              int num_dof_loc = Fb.GetLocalNumber(num_loc, j);
              int num_dof = Nodle(num_dof_loc);
              if ((num_dof >= 0) && (DofDirichlet(num_dof) == -1))
                {
                  DofDirichlet(num_dof) = 0;
                  nb_dof_dirichlet++;
                }
            }
        }
    }

  cout << "Number of Dirichlet dofs = " << nb_dof_dirichlet << endl;

  int nb_dof_border = 0;
  IndexDirichlet.Reallocate(this->mesh_num.GetNbDof());
  IndexDirichlet.Fill(-1);
  IndexLocalBorder.Reallocate(this->mesh_num.GetNbDof());
  IndexLocalBorder.Fill(-1);
  rhoInfini = Complex_wp(0.0, 0.0);
  DofKeptDir.Reallocate(this->mesh_num.GetNbDof()-nb_dof_dirichlet);
  nb_dof_dirichlet = 0;
  for (int i = 0; i < this->mesh_num.GetNbDof(); i++)
    if (DofDirichlet(i) == -1)
      {
        IndexDirichlet(i) = nb_dof_dirichlet;
        DofKeptDir(nb_dof_dirichlet) = i;
        nb_dof_dirichlet++;
      }

  for (int i = 0; i < this->mesh.GetNbBoundaryRef(); i++)
    {
      int ref = this->mesh.BoundaryRef(i).GetReference();
      if (ref == ref_bord(0))
        {
          int num_elem = this->mesh.BoundaryRef(i).numElement(0);
          int num_loc = this->mesh.Element(num_elem).GetPositionBoundary(i);
          const ElementReference<Dimension2, 1>& Fb = this->GetReferenceElementH1(num_elem);
          const IVect& Nodle = this->mesh_num.Element(num_elem).GetNodle();
          if (realpart(rhoInfini) == 0.0)
            {
              ref_water = this->mesh.Element(num_elem).GetReference();
              rhoInfini = this->ref_epsilon(ref_water)*this->ref_mu(ref_water)
                          + Iwp * this->ref_sigma(ref_water) * this->ref_mu(ref_water) / this->omega;
            }
          int nb_dof = Fb.GetNbDofBoundary(num_loc);
          for (int j = 0; j < nb_dof; j++)
            {
              int num_dof_loc = Fb.GetLocalNumber(num_loc, j);
              int num_dof = Nodle(num_dof_loc);
              // cout << "Arête " << i << ", ddl local " << j << endl;
              // DISP(num_dof_loc);
              // DISP(num_dof);
              if (IndexLocalBorder(num_dof) == -1)
                {
                  IndexLocalBorder(num_dof) = nb_dof_border;
                  NumGlobBorder.PushBack(num_dof);
                  nb_dof_border++;
                }
            }
        }
    }
}
