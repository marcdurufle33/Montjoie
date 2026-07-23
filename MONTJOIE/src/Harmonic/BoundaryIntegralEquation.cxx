#ifndef MONTJOIE_FILE_BOUNDARY_INTEGRAL_EQUATION_CXX

#include "BoundaryIntegralEquation.hxx"

namespace Montjoie
{

  /*********************************
   * BoundaryIntegralEquation_Base *
   *********************************/


  //! default constructor
  BoundaryIntegralEquation_Base::BoundaryIntegralEquation_Base()
  {
    omega = 2.0*pi_wp; omega2 = omega*omega;
    frequency = 1.0;
    default_order = 1;
    nb_unknowns_u = 1; nodl = 0;
    print_level = -1;

    order_over_quadrature = 0;
    type_quadrature_element = 0;

    order_integration_singular = 1;
  }
  
  
  //! a line of the data file is interpreted
  void BoundaryIntegralEquation_Base::SetInputData(const string& description_field, const VectString& parameters)
  {
    if (!description_field.compare("FileMesh")) 
      {	
	mesh_data.Clear();
	mesh_data.PushBack(parameters);
      }
    else if (!description_field.compare("AdditionalMesh")) 
      {
	mesh_data.PushBack(parameters);
      }
    else if (!description_field.compare("MateriauDielec"))
      {
	if (parameters.GetM() <= 2)
	  {
	    cout << "In SetInputData of VarProblem" << endl;
	    cout << "MateriauDielec needs more parameters, for instance :" << endl;
	    cout << "MateriauDielec = ref ISOTROPE value" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	int ref = to_num<int>(parameters(0));
	Vector<string> parameters_(parameters.GetM()-1);
	for (int i = 1; i < parameters.GetM(); i++)
	  parameters_(i-1) = parameters(i);
	
	this->SetIndices(ref, parameters_);
      }
    else if (!description_field.compare("PhysicalMedia"))
      {
	if (parameters.GetM() <= 3)
	  {
	    cout << "In SetInputData of VarProblem" << endl;
	    cout << "PhysicalMedia needs more parameters, for instance :" << endl;
	    cout << "PhysicalMedia = nom_media ref ISOTROPE value" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }
        
        string name_media = parameters(0);
        int ref = to_num<int>(parameters(1));
	Vector<string> parameters_(parameters.GetM()-2);
	for (int i = 2; i < parameters.GetM(); i++)
	  parameters_(i-2) = parameters(i);
	
	this->SetPhysicalIndex(name_media, ref, parameters_);        
      }
    else if (!description_field.compare("OrderGeometry"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of MeshNumbering_Base" << endl;
	    cout << "OrderGeometry needs more parameters, for instance :" << endl;
	    cout << "OrderGeometry = r" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	int r = to_num<int>(parameters(0));
	default_order = r;
      }
    else if (!description_field.compare("OrderDiscretization"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of MeshNumbering_Base" << endl;
	    cout << "OrderDiscretization needs more parameters, for instance :" << endl;
	    cout << "OrderDiscretization = r" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	if ((!parameters(0).compare("MAX_EDGE")) || (!parameters(0).compare("MEAN_EDGE")))
	  {
	  }
	else
	  {
	    int r = to_num<int>(parameters(0));
	    default_order = r;
	  }
      }
    else if (description_field == "IntegrationFarPanel")
      {
	order_over_quadrature = to_num<int>(parameters(0));
	type_quadrature_element = to_num<int>(parameters(1));	
      }
    else if (!description_field.compare("PrintLevel"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of VarProblem" << endl;
	    cout << "PrintLevel needs more parameters, for instance :" << endl;
	    cout << "PrintLevel = level" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	print_level = to_num<int>(parameters(0));
      }
    else if (!description_field.compare("Frequency"))
      {
	if (parameters.GetM() <= 1)
	  {
	    cout << "In SetInputData of VarHarmonic" << endl;
	    cout << "Frequency needs more parameters, for instance :" << endl;
	    cout << "Frequency = a b" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }
	
	this->omega = to_num<Real_wp>(parameters(0))*pi_wp*2.0
	  + to_num<Real_wp>(parameters(1));

	// omega is the pulsation
	// omega2 is the square of omega
	this->omega2 = this->omega*this->omega;
	this->frequency = this->omega/(2.0*pi_wp);
      }
    else if (!description_field.compare("PhysicalFrequency"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of VarHarmonic" << endl;
	    cout << "PhysicalFrequency needs more parameters, for instance :" << endl;
	    cout << "PhysicalFrequency = f0" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	to_num(parameters(0), this->frequency);
	this->frequency *= 1.0/PhysicalConstant::speed_light; // DISP(frequency);
	this->omega = 2.0*pi_wp*this->frequency; 
	
	// omega is the pulsation
	// omega2 is the square of omega
	this->omega2 = this->omega*this->omega;
      }
    else if (description_field == "OrderIntegrationSingular")
      this->order_integration_singular = to_num<int>(parameters(0));
  }
  

  void BoundaryIntegralEquation_Base
  ::ConstructAll(const string& input_file, const string& name_elt)
  {
    Vector<string> lines_data_file;
    if (input_file.compare("NONE"))
      ReadLinesFile(input_file, lines_data_file);
    
    ConstructAll(lines_data_file, name_elt);
  }


  void BoundaryIntegralEquation_Base
  ::ConstructAll(const Vector<string>& lines_data_file, const string& name_elt)
  {
    this->InitIndices(PhysicalConstant::nb_max_indices);

    ReadInputFile(lines_data_file, *this);

    ComputeMeshAndFiniteElement(name_elt);
    
    ComputeGeometryQuantities();
  }


  void BoundaryIntegralEquation_Base
  ::RunAll(const string& input_file, const string& name_elt)
  {
    // on inialise les calculs
    ConstructAll(input_file, name_elt);

    // on calcule la matrice elements finis
    Matrix<Complex_wp> A;
    ComputeIntegralMatrix(A);
    A.Write("mat.dat");
    
    // on calcule le second membre
    VectComplex_wp b;
    ComputeRightHandSide(b);
    b.Write("rhs.dat");

    // on resout le systeme A x = b
    Vector<int> pivot;
    GetLU(A, pivot);
    SolveLU(A, pivot, b);

    b.Write("sol.dat");
    DISP(Norm2(b));
    // on ecrit la solution
    WriteSolution(b);
  }

  void BoundaryIntegralEquation_Base
  ::ComputeIntegralMatrix(Matrix<Complex_wp>& A) const
  {
    A.Reallocate(nodl, nodl);
    Matrix<Complex_wp> mat_elem;
    Vector<int> num_row, num_col;
    cout << "Computing the dense matrix..." << endl;
    int old_percent = 0, new_percent = 0;
    for (int i = 0; i < this->GetNbBoundary(); i++)
      {
        // displays a progress bar
	new_percent = toInteger(round(Real_wp(i)/(this->GetNbBoundary()-1)*80));
        if (this->print_level >= 2)
          for (int percent = old_percent; percent < new_percent; percent++)
            { cout<<"#"; cout.flush(); }
	
	old_percent = new_percent;
        
        for (int j = 0; j < this->GetNbBoundary(); j++)
          {
            this->ComputeElementaryMatrix(i, j, num_row, num_col, mat_elem);
            
            for (int i2 = 0; i2 < num_row.GetM(); i2++)
              for (int j2 = 0; j2 < num_col.GetM(); j2++)
                A(num_row(i2), num_col(j2)) += mat_elem(i2, j2);
          }
      }

    if (this->print_level >= 2)
      cout << endl;
    
    cout << "Computation done" << endl;
  }


  /***********************
   * GeometryQuantityBEM *
   ***********************/


  template<class Dimension>
  void GeometryQuantityBEM<Dimension>
  ::ComputeGeometry(const Mesh<Dimension>& mesh,
		    const Vector<const ElementReference_Dim<DimensionB>* >& fem)
  {
    VectR3 s;
    SetPoints<Dimension> Points;
    SetMatrices<Dimension> Mat;
    TinyMatrix<Real_wp, General, 3, 2> mat_DF;
    int offset = 0, offsetQ = 0;
    for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
      {
	const ElementReference_Dim<DimensionB>& Fb = *fem(i);
	offset += Fb.GetNbPointsNodalElt();
	offsetQ += Fb.GetNbPointsQuadratureInside();
      }

    PointsNodal.Reallocate(offset);
    JacobNodal.Reallocate(offset);
    NormaleNodal.Reallocate(offset);
    MatDfNodal.Reallocate(offset);

    PointsQuad.Reallocate(offsetQ);
    JacobQuad.Reallocate(offsetQ);
    NormaleQuad.Reallocate(offsetQ);
    MatDfQuad.Reallocate(offsetQ);

    offset = 0; offsetQ = 0;
    VectR3 dF_dx, dF_dy, dF_dx_quad, dF_dy_quad;
    R3 vec_u;
    OffsetNodal.Reallocate(mesh.GetNbBoundaryRef()+1);
    OffsetQuad.Reallocate(mesh.GetNbBoundaryRef()+1);
    OffsetNodal(0) = 0; OffsetQuad(0) = 0;
    for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
      {
	int nb_vert = mesh.BoundaryRef(i).GetNbVertices();
	s.Reallocate(nb_vert);
	for (int j = 0; j < nb_vert; j++)
	  s(j) = mesh.Vertex(mesh.BoundaryRef(i).numVertex(j));

	const ElementReference_Dim<DimensionB>& Fb = *fem(i);
	const ElementGeomReference<DimensionB>& Fb_geom = fem(i)->GetGeometricElement();
	Fb_geom.FjElemQuadrature(s, Points, mesh, i, mesh.BoundaryRef(i));

	bool is_curved = mesh.IsBoundaryCurved(i);
	Fb_geom.DFjElemNodal(Points.GetPointNodal(), is_curved, dF_dx, dF_dy);
	Fb_geom.DFjElemQuadrature(Points.GetPointNodal(), is_curved, dF_dx, dF_dy,
				  dF_dx_quad, dF_dy_quad);

	for (int j = 0; j < Fb.GetNbPointsNodalElt(); j++)
	  {
	    PointsNodal(offset + j) = Points.GetPointNodal(j);
	    
	    TimesProd(dF_dx(j), dF_dy(j), vec_u);
	    JacobNodal(offset + j) = Norm2(vec_u);

	    Mlt(1.0/JacobNodal(offset + j), vec_u);
            
	    NormaleNodal(offset + j) = vec_u;
            	    
	    SetCol(dF_dx(j), 0, mat_DF); SetCol(dF_dy(j), 1, mat_DF);
	    MatDfNodal(offset+j) = mat_DF;
	  }

	for (int j = 0; j < Fb.GetNbPointsQuadratureInside(); j++)
	  {
	    PointsQuad(offsetQ + j) = Points.GetPointQuadrature(j);
	      
	    TimesProd(dF_dx_quad(j), dF_dy_quad(j), vec_u);
	    JacobQuad(offsetQ + j) = Norm2(vec_u);

	    Mlt(1.0/JacobQuad(offsetQ + j), vec_u);
            
	    NormaleQuad(offsetQ + j) = vec_u;

	    SetCol(dF_dx_quad(j), 0, mat_DF); SetCol(dF_dy_quad(j), 1, mat_DF);
            //DISP(i); DISP(j); DISP(mat_DF);
	    MatDfQuad(offsetQ + j) = mat_DF;
	  }

	offset += Fb.GetNbPointsNodalElt();
	offsetQ += Fb.GetNbPointsQuadratureInside();
	OffsetNodal(i+1) = offset;
	OffsetQuad(i+1) = offsetQ;
      }

    //DISP(PointsQuad); DISP(NormaleQuad); DISP(JacobQuad); 
  }


  template<class Dimension>
  void GeometryQuantityBEM<Dimension>
  ::FillJacobianMatQuad(int i, Vector<R_N>& PointsQuadI,
			Vector<TinyMatrix<Real_wp, General, 3, 2> >& MatDfI) const
  {
    int N = OffsetQuad(i+1) - OffsetQuad(i);
    PointsQuadI.SetData(N, const_cast<R_N*>(&PointsQuad(OffsetQuad(i))));
    MatDfI.SetData(N, const_cast<TinyMatrix<Real_wp, General, 3, 2>* >(&MatDfQuad(OffsetQuad(i))));
  }


  template<class Dimension>
  void GeometryQuantityBEM<Dimension>
  ::FillNormaleQuad(int i, Vector<R_N>& PointsQuadI, Vector<R_N>& NormaleQuadI) const
  {
    int N = OffsetQuad(i+1) - OffsetQuad(i);
    PointsQuadI.SetData(N, const_cast<R_N*>(&PointsQuad(OffsetQuad(i))));
    NormaleQuadI.SetData(N, const_cast<R_N*>(&NormaleQuad(OffsetQuad(i))));
  }


  template<class Dimension>
  void GeometryQuantityBEM<Dimension>
  ::FillNormaleDfQuad(int i, Vector<R_N>& PointsQuadI, Vector<R_N>& NormaleQuadI,
		      Vector<TinyMatrix<Real_wp, General, 3, 2> >& MatDfI) const
  {
    int N = OffsetQuad(i+1) - OffsetQuad(i);
    PointsQuadI.SetData(N, const_cast<R_N*>(&PointsQuad(OffsetQuad(i))));
    NormaleQuadI.SetData(N, const_cast<R_N*>(&NormaleQuad(OffsetQuad(i))));
    MatDfI.SetData(N, const_cast<TinyMatrix<Real_wp, General, 3, 2>* >(&MatDfQuad(OffsetQuad(i))));
  }
  
  
  /********************************
   * BoundaryIntegralEquation_Dim *
   ********************************/


  //! default constructor
  template<class Dimension>
  BoundaryIntegralEquation_Dim<Dimension>::BoundaryIntegralEquation_Dim() : mesh_num(mesh)
  {
  }
  
  
  //! modifies the object with a line of the data file
  template<class Dimension>
  void BoundaryIntegralEquation_Dim<Dimension>
  ::SetInputData(const string& description_field, const VectString& parameters)
  {
    BoundaryIntegralEquation_Base::SetInputData(description_field, parameters);
    mesh.SetInputData(description_field, parameters);
    mesh_num.SetInputData(description_field, parameters);

    if (!description_field.compare("IncidentAngle"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of VarHarmonic" << endl;
	    cout << "IncidentAngle needs more parameters, for instance :" << endl;
	    cout << "IncidentAngle = teta" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }
	
	// in 2-D teta is the incident angle of the plane wave
	// in 3-D, the wave vector is oriented along
	// (sin(teta)*cos(phi) ; sin(teta)*sin(phi) ; cos(teta))
	Real_wp teta = 0.0, phi = 0.0;
	if (Dimension::dim_N >= 2)
	  teta = to_num<Real_wp>(parameters(0))*pi_wp/180;
	    
	if (Dimension::dim_N == 3)
	  if (parameters.GetM() > 1)
	    phi = to_num<Real_wp>(parameters(1))*pi_wp/180;
	    
	SetIncidentAngle(this->omega, this->wave_vector, teta, phi);
      }

  }
  

  template<class Dimension>
  void BoundaryIntegralEquation_Dim<Dimension>
  ::ConstructMesh(const Vector<string>& param, Mesh<Dimension>& mesh)
  {
    int nb = 0;
    if (param(nb) == "PLATE")
      {
        nb++;
        int nx = to_num<int>(param(nb++));
        int ny = to_num<int>(param(nb++));
        Mesh<Dimension2> carre;
        carre.CreateRegularMesh(R2(0, 0), R2(1, 1), TinyVector<int, 2>(nx, ny),
                                1, TinyVector<int, 4>(1, 1, 1, 1), carre.QUADRILATERAL_MESH);
        
        R3 pt0(to_num<Real_wp>(param(3)), to_num<Real_wp>(param(4)), to_num<Real_wp>(param(5)));
        R3 ptA(to_num<Real_wp>(param(6)), to_num<Real_wp>(param(7)), to_num<Real_wp>(param(8)));
        R3 ptB(to_num<Real_wp>(param(9)), to_num<Real_wp>(param(10)), to_num<Real_wp>(param(11)));
        nb += 9;
        
        mesh.ReallocateVertices(carre.GetNbVertices());
        mesh.ReallocateBoundariesRef(carre.GetNbElt());
        for (int i = 0; i < carre.GetNbVertices(); i++)
          {
            Real_wp x = carre.Vertex(i)(0);
            Real_wp y = carre.Vertex(i)(1);

            R3 pt = pt0 + x*(ptA - pt0) + y*(ptB-pt0);
            mesh.Vertex(i) = pt;
          }

        for (int i = 0; i < carre.GetNbElt(); i++)
          {
            int n1 = carre.Element(i).numVertex(0);
            int n2 = carre.Element(i).numVertex(1);
            int n3 = carre.Element(i).numVertex(2);
            int n4 = carre.Element(i).numVertex(3);
            int ref = carre.Element(i).GetReference();
            mesh.BoundaryRef(i).InitQuadrangular(n1, n2, n3, n4, ref);            
          }

        mesh.FindConnectivity();
      }
    else if (param(nb) == "PAVE")
      {
        nb++;
        int nx = to_num<int>(param(nb++));
        int ny = to_num<int>(param(nb++));
        int nz = to_num<int>(param(nb++));
        Mesh<Dimension2> carre_xy, carre_xz, carre_yz;
        mesh.Clear();
        
        carre_xy.CreateRegularMesh(R2(0, 0), R2(1, 1), TinyVector<int, 2>(nx, ny),
                                   1, TinyVector<int, 4>(1, 1, 1, 1), carre_xy.QUADRILATERAL_MESH);
        
        carre_xz.CreateRegularMesh(R2(0, 0), R2(1, 1), TinyVector<int, 2>(nx, nz),
                                   1, TinyVector<int, 4>(1, 1, 1, 1), carre_xy.QUADRILATERAL_MESH);

        carre_yz.CreateRegularMesh(R2(0, 0), R2(1, 1), TinyVector<int, 2>(ny, nz),
                                   1, TinyVector<int, 4>(1, 1, 1, 1), carre_xy.QUADRILATERAL_MESH);
        
        Real_wp xmin = to_num<Real_wp>(param(nb++));
        Real_wp xmax = to_num<Real_wp>(param(nb++));
        Real_wp ymin = to_num<Real_wp>(param(nb++));
        Real_wp ymax = to_num<Real_wp>(param(nb++));
        Real_wp zmin = to_num<Real_wp>(param(nb++));
        Real_wp zmax = to_num<Real_wp>(param(nb++));
        
        mesh.ReallocateVertices(2*(nx*ny + nx*nz + ny*nz));
        mesh.ReallocateBoundariesRef(2*((nx-1)*(ny-1) + (nx-1)*(nz-1) + (ny-1)*(nz-1)));

        int off = (nx-1)*(ny-1), offn = nx*ny;
        for (int i = 0; i < carre_xy.GetNbVertices(); i++)
          {
            Real_wp x = carre_xy.Vertex(i)(0);
            Real_wp y = carre_xy.Vertex(i)(1);

            R3 pt(xmin + (xmax-xmin)*x, ymin + (ymax-ymin)*y, zmin);
            mesh.Vertex(i) = pt;
            pt(2) = zmax;
            mesh.Vertex(offn + i) = pt;
          }        

        for (int i = 0; i < carre_xy.GetNbElt(); i++)
          {
            int n1 = carre_xy.Element(i).numVertex(0);
            int n2 = carre_xy.Element(i).numVertex(1);
            int n3 = carre_xy.Element(i).numVertex(2);
            int n4 = carre_xy.Element(i).numVertex(3);
            int ref = carre_xy.Element(i).GetReference();
            mesh.BoundaryRef(i).InitQuadrangular(n1, n2, n3, n4, ref);
            mesh.BoundaryRef(off+i).InitQuadrangular(offn + n1, offn + n2, offn + n3, offn + n4, ref);
          }

        int off0 = 2*off; off = off0 + (nx-1)*(nz-1);
        int offn0 = 2*offn; offn = offn0 + nx*nz;
        for (int i = 0; i < carre_xz.GetNbVertices(); i++)
          {
            Real_wp x = carre_xz.Vertex(i)(0);
            Real_wp y = carre_xz.Vertex(i)(1);

            R3 pt(xmin + (xmax-xmin)*x, ymin, zmin + (zmax-zmin)*y);
            mesh.Vertex(offn0 + i) = pt;
            pt(1) = ymax;
            mesh.Vertex(offn + i) = pt;
          }        

        for (int i = 0; i < carre_xz.GetNbElt(); i++)
          {
            int n1 = carre_xz.Element(i).numVertex(0);
            int n2 = carre_xz.Element(i).numVertex(1);
            int n3 = carre_xz.Element(i).numVertex(2);
            int n4 = carre_xz.Element(i).numVertex(3);
            int ref = carre_xz.Element(i).GetReference();
            mesh.BoundaryRef(off0+i).InitQuadrangular(offn0 + n1, offn0 + n2, offn0 + n3, offn0 + n4, ref);
            mesh.BoundaryRef(off+i).InitQuadrangular(offn + n1, offn + n2, offn + n3, offn + n4, ref);
          }
        
        off0 = off + (nx-1)*(nz-1); off = off0 + (ny-1)*(nz-1);
        offn0 = offn + nx*nz; offn = offn0 + ny*nz;
        for (int i = 0; i < carre_yz.GetNbVertices(); i++)
          {
            Real_wp x = carre_yz.Vertex(i)(0);
            Real_wp y = carre_yz.Vertex(i)(1);

            R3 pt(xmin, ymin + (ymax-ymin)*x, zmin + (zmax-zmin)*y);
            mesh.Vertex(offn0 + i) = pt;
            pt(0) = xmax;
            mesh.Vertex(offn + i) = pt;
          }        

        for (int i = 0; i < carre_yz.GetNbElt(); i++)
          {
            int n1 = carre_yz.Element(i).numVertex(0);
            int n2 = carre_yz.Element(i).numVertex(1);
            int n3 = carre_yz.Element(i).numVertex(2);
            int n4 = carre_yz.Element(i).numVertex(3);
            int ref = carre_yz.Element(i).GetReference();
            mesh.BoundaryRef(off0+i).InitQuadrangular(offn0 + n1, offn0 + n2, offn0 + n3, offn0 + n4, ref);
            mesh.BoundaryRef(off+i).InitQuadrangular(offn + n1, offn + n2, offn + n3, offn + n4, ref);
          }
        
        mesh.FindConnectivity();
        mesh.RemoveDuplicateVertices();
      }
    else
      mesh.Read(param(nb++));

    R3 center;
    if (param.GetM() > nb)
      {
        if (param(nb) == "Center")
          {
            nb++;
            Real_wp x = to_num<Real_wp>(param(nb++));
            Real_wp y = to_num<Real_wp>(param(nb++));
            Real_wp z = to_num<Real_wp>(param(nb++));
            center.Init(x, y, z);
          }
      }
    
    CheckOrientationNormales(mesh, center);
  }
  

  template<class Dimension>
  void BoundaryIntegralEquation_Dim<Dimension>
  ::CheckOrientationNormales(Mesh<Dimension>& mesh, const R3& center)
  {
    R3 normale; bool change_elt = false;
    R3 center_elt;
    for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
      {
        int nb_vert = mesh.BoundaryRef(i).GetNbVertices();
        int n1 = mesh.BoundaryRef(i).numVertex(0);
        int n2 = mesh.BoundaryRef(i).numVertex(1);
        int n3 = mesh.BoundaryRef(i).numVertex(2);
        int n4 = mesh.BoundaryRef(i).numVertex(nb_vert-1);
        int ref = mesh.BoundaryRef(i).GetReference();
        R3 vec_u = mesh.Vertex(n2) - mesh.Vertex(n1);
        R3 vec_v = mesh.Vertex(n4) - mesh.Vertex(n1);
        center_elt = mesh.Vertex(n1) + mesh.Vertex(n2) + mesh.Vertex(n3);
        if (nb_vert == 4)
          center_elt += mesh.Vertex(n4);
        
        center_elt *= Real_wp(1)/nb_vert;
        TimesProd(vec_u, vec_v, normale);
        R3 vec_w = center_elt - center;
        if (DotProd(normale, vec_w) < 0)
          {
            change_elt = true;
            if (nb_vert == 3)
              mesh.BoundaryRef(i).InitTriangular(n1, n3, n2, ref);
            else
              mesh.BoundaryRef(i).InitQuadrangular(n1, n4, n3, n2, ref);
          }
      }
    
    mesh.FindConnectivity();
  }

  
  template<class Dimension>
  void BoundaryIntegralEquation_Dim<Dimension>
  ::ComputeMeshAndFiniteElement(const string& name_elt)
  {
    DISP(this->mesh_data(0));
    if (this->mesh_data.GetM() > 0)
      ConstructMesh(this->mesh_data(0), this->mesh);

    this->mesh.Write("test.mesh");
    DISP(this->mesh.GetNbEdges());
    
    this->ConstructFiniteElement(name_elt);
    this->mesh_num.NumberMesh(true);
    
    nodl = nb_unknowns_u * this->mesh_num.GetNbDof();
    DISP(nodl); DISP(this->mesh.GetNbEdges());
  }


  /********************************
   * BoundaryIntegralEquation_Fem *
   ********************************/


  template<class Dimension, int type>
  BoundaryIntegralEquation_Fem<Dimension, type>::BoundaryIntegralEquation_Fem()
  {
    tri_element = NULL;
    quad_element = NULL;
  }


  template<class Dimension, int type>
  BoundaryIntegralEquation_Fem<Dimension, type>::~BoundaryIntegralEquation_Fem()
  {
    if (tri_element != NULL)
      delete tri_element;

    if (quad_element != NULL)
      delete quad_element;
  }

  
  //! construction of finite elements
  template<class Dimension, int type>
  void BoundaryIntegralEquation_Fem<Dimension, type>::ConstructFiniteElement(const string& name_elt)
  {
    // ordre constant
    int r = this->mesh_num.GetOrder();
    int r1 = max(r, 1);
    int rgeom = this->mesh.GetGeometryOrder();
    DISP(r); DISP(rgeom);
    int rquad = r + this->order_over_quadrature;
    int type_quad = this->type_quadrature_element;
    DISP(rquad); DISP(type_quad);
    if (this->mesh.GetNbTrianglesRef() > 0)
      {
	int id_elt = VarFiniteElementEnum<DimensionB, type>::GetIdentityNumber(name_elt, 0);
	tri_element  = VarFiniteElementEnum<DimensionB, type>::GetNewReferenceElement(r, id_elt);
	tri_element->ConstructFiniteElement(r, min(r1, rgeom), rquad, type_quad);
        tri_element->ConstructNumberMap(this->mesh_num.number_map, ElementReference_Base::CONTINUOUS);
      }

    if (this->mesh.GetNbQuadranglesRef() > 0)
      {
	int id_elt = VarFiniteElementEnum<DimensionB, type>::GetIdentityNumber(name_elt, 1);
	quad_element  = VarFiniteElementEnum<DimensionB, type>::GetNewReferenceElement(r, id_elt);
	quad_element->ConstructFiniteElement(r, min(r1, rgeom), rquad, type_quad);
        quad_element->ConstructNumberMap(this->mesh_num.number_map, ElementReference_Base::CONTINUOUS);
      }

    reference_element.Reallocate(this->mesh.GetNbBoundaryRef());
    for (int i = 0; i < this->mesh.GetNbBoundaryRef(); i++)
      {
	if (this->mesh.BoundaryRef(i).IsTriangular())
	  reference_element(i) = tri_element;
	else
	  reference_element(i) = quad_element;
      }

    int type_integration_edge(0), type_integration_tri(0), type_integration_quad(0);
    //this->GetTypeIntegrationBoundary(this->mesh, this->mesh_num, *this, type_integration_edge,
    //				     type_integration_tri, type_integration_quad);

    TinyVector<IVect, 4> order;
    this->mesh_num.GetOrder(order);
    
    if (Dimension::dim_N == 2)
      {
	this->mesh_num.number_map.ConstructQuadrature2D(order, type_integration_edge);
      }
    else
      {
#ifdef MONTJOIE_WITH_THREE_DIM        
        this->mesh_num.number_map.ConstructQuadrature3D(order, type_integration_tri, type_integration_quad);
#endif
      }    
  }


  template<class Dimension, int type>
  void BoundaryIntegralEquation_Fem<Dimension, type>
  ::ComputeGeometryQuantities()
  {
    Vector<const ElementReference_Dim<DimensionB>* > obj_fem(this->mesh.GetNbBoundaryRef());
    for (int i = 0; i < this->mesh.GetNbBoundaryRef(); i++)
      obj_fem(i) = &this->GetReferenceElement(i);
    
    this->var_geom.ComputeGeometry(this->mesh, obj_fem);    
  }

  /**************************************
   * BoundaryIntegralEquation_Maxwell3D *
   **************************************/

  BoundaryIntegralEquation_Maxwell3D::BoundaryIntegralEquation_Maxwell3D()
  {
    polarization.Init(Real_wp(1), Real_wp(0), Real_wp(0));
  }

  
  void BoundaryIntegralEquation_Maxwell3D::InitIndices(int N)
  {
    ref_epsilon.Reallocate(N);
    ref_mu.Reallocate(N);
    ref_sigma.Reallocate(N);
  }

  
  void BoundaryIntegralEquation_Maxwell3D::SetIndices(int ref, const Vector<string>& param)
  {
    ref_epsilon(ref) = to_num<Complex_wp>(param(1));
    ref_mu(ref) = to_num<Complex_wp>(param(2));
    ref_sigma(ref) = to_num<Complex_wp>(param(3));
  }

  
  void BoundaryIntegralEquation_Maxwell3D
  ::SetPhysicalIndex(const string& name, int ref, const Vector<string>& param)
  {
    if (name == "epsilon")
      ref_epsilon(ref) = to_num<Complex_wp>(param(1));
    else if (name == "mu")
      ref_mu(ref) = to_num<Complex_wp>(param(1));
    else if (name == "sigma")
      ref_sigma(ref) = to_num<Complex_wp>(param(1));
  }

  
  void BoundaryIntegralEquation_Maxwell3D
  ::SetInputData(const string& description_field, const VectString& parameters)
  {
    BoundaryIntegralEquation_Fem<Dimension3, 3>::SetInputData(description_field, parameters);
    if (description_field == "Polarization")
      {
	for (int k = 0; k < 3; k++)
	  polarization(k) = to_num<Real_wp>(parameters(k));
      }
  }

  
  void BoundaryIntegralEquation_Maxwell3D
  ::ComputeElementaryMatrix(int i, int j, Vector<int>& row, Vector<int>& col,
			    Matrix<Complex_wp>& mat_elem) const
  {
    const ElementReference<Dimension2, 3>& Fb = this->GetReferenceElement(i);
    
    int N = Fb.GetNbDof();
    int Nquad = Fb.GetNbPointsQuadratureInside();
    Matrix<R2> ValuePhi(N, Nquad);
    Matrix<Real_wp> DivPhi(N, Nquad);
    VectR2 phi(N); VectReal_wp div_phi(N);
    for (int k = 0; k < Nquad; k++)
      {
	Fb.GetValuePhiOnQuadraturePoint(k, phi);
	Fb.GetDivPhiOnQuadraturePoint(k, div_phi);
	for (int i = 0; i < N; i++)
	  {
	    ValuePhi(i, k) = phi(i);
	    DivPhi(i, k) = div_phi(i);
	  }
      }

    row.Reallocate(N); col.Reallocate(N);
    for (int i1 = 0; i1 < N; i1++)
      {
	row(i1) = this->mesh_num.Element(i).GetNumberDof(i1);
	col(i1) = this->mesh_num.Element(j).GetNumberDof(i1);
      }

    /* if (j == 0)
      {
        DISP(i); DISP(this->mesh.BoundaryRef(i).numVertex(0));
        DISP(this->mesh.BoundaryRef(i).numVertex(1));
        DISP(this->mesh.BoundaryRef(i).numVertex(2));
        DISP(this->mesh.BoundaryRef(i).numVertex(3));
        int n1 = this->mesh.BoundaryRef(i).numVertex(0);
        DISP(this->mesh.Vertex(n1));
        DISP(row);
      }
    */

    mat_elem.Reallocate(N, N);
    mat_elem.Zero();
    if (i == j)
      {
	// elements confondus
	VectR3 PointsQuadI;
	Vector<TinyMatrix<Real_wp, General, 3, 2> > MatDfI;
	this->var_geom.FillJacobianMatQuad(i, PointsQuadI, MatDfI);

	R3 ptA = this->mesh.Vertex(this->mesh.BoundaryRef(j).numVertex(0));
	R3 ptB = this->mesh.Vertex(this->mesh.BoundaryRef(j).numVertex(1));
	R3 ptC = this->mesh.Vertex(this->mesh.BoundaryRef(j).numVertex(2));
        R3 ptD;
        if (Fb.GetNbVertices() == 4)
          ptD = this->mesh.Vertex(this->mesh.BoundaryRef(j).numVertex(3));
	
	for (int k = 0; k < Nquad; k++)
	  {
	    // Duffy transformation
	    // we create a small mesh around the singular point
	    Mesh<Dimension2> mesh_sing;
	    int nb_vert = this->mesh.BoundaryRef(i).GetNbVertices();
	    mesh_sing.ReallocateVertices(nb_vert + 1);
	    if (nb_vert == 3)
	      {
		mesh_sing.Vertex(0).Init(Real_wp(0), Real_wp(0));
		mesh_sing.Vertex(1).Init(Real_wp(1), Real_wp(0));
		mesh_sing.Vertex(2).Init(Real_wp(0), Real_wp(1));
	      }
	    else
	      {
		mesh_sing.Vertex(0).Init(Real_wp(0), Real_wp(0));
		mesh_sing.Vertex(1).Init(Real_wp(1), Real_wp(0));
		mesh_sing.Vertex(2).Init(Real_wp(1), Real_wp(1));
		mesh_sing.Vertex(3).Init(Real_wp(0), Real_wp(1));
	      }

	    mesh_sing.Vertex(nb_vert) = Fb.PointsND(k);
	    mesh_sing.ReallocateElements(nb_vert);
	    for (int i1 = 0; i1 < nb_vert; i1++)
	      {
		int n1 = i1, n2 = (i1+1)%nb_vert, n3 = nb_vert;
		mesh_sing.Element(i1).InitTriangular(n3, n1, n2, 1);
	      }

            if ((i == 12) && (k == 0))
              mesh_sing.Write("sing.mesh");
            
	    VectR2 s(3);
	    Globatto<Real_wp> gauss_Y; R2 vec_u, vec_v;
	    int Ny = this->order_integration_singular+1;
	    gauss_Y.ConstructQuadrature(Ny-1);
	    // loop on triangles
	    VectR2 PointsQuadY(nb_vert*Ny*Ny);
	    VectReal_wp WeightsQuadY(nb_vert*Ny*Ny);
	    int nb = 0;
            ofstream file_out; file_out.precision(15);
            if ((i == 12) && (k == 0))
              file_out.open("points_sing.dat");
            
	    for (int i1 = 0; i1 < mesh_sing.GetNbElt(); i1++)
	      {
		mesh_sing.GetVerticesElement(i1, s);
                //if ((i == 12) && (k == 0))
                //{ DISP(i1); DISP(s); }
		for (int l1 = 0; l1 < Ny; l1++)
		  for (int l2 = 0; l2 < Ny; l2++)
		    {
		      Real_wp x = gauss_Y.Points(l1);
		      Real_wp y = gauss_Y.Points(l2);
                      //if ((i == 12) && (k == 0))
                      //{ DISP(l1); DISP(l2); DISP(x); DISP(y); }
		      PointsQuadY(nb) = (1.0-x)*s(0) + x*(1.0-y)*s(1) + x*y*s(2);
                      //if ((i == 12) && (k == 0))
                      //{ DISP(PointsQuadY(nb)); }
		      vec_u = s(1) - s(0) + y*(s(2)-s(1));
		      vec_v = x*(s(2) - s(1));
		      Real_wp detY = abs(vec_u(0)*vec_v(1) - vec_u(1)*vec_v(0));
		      WeightsQuadY(nb) = gauss_Y.Weights(l1)*gauss_Y.Weights(l2)*detY;
                      if ((i == 12) && (k == 0))
                        file_out << PointsQuadY(nb)(0) << " " << PointsQuadY(nb)(1) << " " << WeightsQuadY(nb) << '\n';
		      nb++;
		    }
	      }

            if ((i == 12) && (k == 0))
              file_out.close();

            /*if (i == 12)
              {
                DISP(k); DISP(Fb.PointsND(k));
                DISP(PointsQuadY); DISP(WeightsQuadY);
                } */
            //Complex_wp val(0, 0);
	    TinyMatrix<Real_wp, General, 2, 3> mat_df_t;
	    TinyMatrix<Real_wp, General, 2, 2> mat22;
	    for (int l = 0; l < PointsQuadY.GetM(); l++)
	      {
		Fb.ComputeValuesPhiRef(PointsQuadY(l), phi);
		Fb.ComputeDivPhiRef(PointsQuadY(l), div_phi);

		Real_wp x = PointsQuadY(l)(0), y = PointsQuadY(l)(1);
		R3 pt_reel;
                if (Fb.GetNbVertices() == 3)
                  pt_reel = (1.0-x-y)*ptA + x*ptB + y*ptC;
                else
                  pt_reel = (1.0-x)*(1.0-y)*ptA + x*(1.0-y)*ptB + x*y*ptC + (1.0-x)*y*ptD;
                
		Real_wp dist = PointsQuadI(k).Distance(pt_reel);
		Complex_wp phase(cos(this->omega*dist), sin(this->omega*dist));
		Complex_wp green_kernel = phase / (4.0*pi_wp*dist);
		Complex_wp poids = Fb.WeightsND(k)*WeightsQuadY(l)*green_kernel;
		Complex_wp poids2 = -poids / this->omega2;

		Transpose(MatDfI(0), mat_df_t);
		Mlt(mat_df_t, MatDfI(0), mat22);
		for (int i1 = 0; i1 < N; i1++)
		  for (int j1 = 0; j1 < N; j1++)
		    {
		      Mlt(mat22, ValuePhi(i1, k), vec_u);
		      Complex_wp vloc = poids*DotProd(vec_u, phi(j1));
                      //if ((i1 == 3) && (j1 == 3))
                      //val += vloc;
                      
		      vloc += poids2*DivPhi(i1, k)*div_phi(j1);
		      mat_elem(i1, j1) += vloc;
		    }
	      }

            // if (i == 12) { DISP(val / (0.25*Fb.WeightsND(k))); exit(0); }
            
	  }

	PointsQuadI.Nullify(); MatDfI.Nullify();
      }
    else
      {
	VectR3 PointsQuadI, PointsQuadJ;
	Vector<TinyMatrix<Real_wp, General, 3, 2> > MatDfI, MatDfJ;
	this->var_geom.FillJacobianMatQuad(i, PointsQuadI, MatDfI);
	this->var_geom.FillJacobianMatQuad(j, PointsQuadJ, MatDfJ);

	// cas de deux elements loins
	TinyMatrix<Real_wp, General, 2, 3> mat_df_t;
	TinyMatrix<Real_wp, General, 2, 2> mat22; R2 vec_u;
	for (int k = 0; k < Nquad; k++)
          for (int l = 0; l < Nquad; l++)
	    {
	      Real_wp dist = PointsQuadI(k).Distance(PointsQuadJ(l));
	      Complex_wp phase(cos(this->omega*dist), sin(this->omega*dist));
	      Complex_wp green_kernel = phase / (4.0*pi_wp*dist);
	      Complex_wp poids = Fb.WeightsND(k)*Fb.WeightsND(l)*green_kernel;
	      Complex_wp poids2 = -poids / this->omega2;

	      Transpose(MatDfJ(l), mat_df_t);
	      Mlt(mat_df_t, MatDfI(k), mat22);
	      for (int i1 = 0; i1 < N; i1++)
		for (int j1 = 0; j1 < N; j1++)
		  {
		    Mlt(mat22, ValuePhi(i1, k), vec_u);
		    Complex_wp vloc = poids*DotProd(vec_u, ValuePhi(j1, l));
		    vloc += poids2*DivPhi(i1, k)*DivPhi(j1, l);
		    mat_elem(i1, j1) += vloc;
		  }
	    }
        
	PointsQuadI.Nullify(); PointsQuadJ.Nullify();
	MatDfI.Nullify(); MatDfJ.Nullify();
      }

    {
      const IVect& neg_dof = this->mesh_num.Element(i).GetNegativeDofNumber();
      for (int i1 = 0; i1 < neg_dof.GetM(); i1++)
	{
	  int i2 = neg_dof(i1);
	  for (int j1 = 0; j1 < mat_elem.GetN(); j1++)
	    mat_elem(i2, j1) = -mat_elem(i2, j1);
	}
    }

    {
      const IVect& neg_dof = this->mesh_num.Element(j).GetNegativeDofNumber();
      for (int j1 = 0; j1 < neg_dof.GetM(); j1++)
	{
	  int j2 = neg_dof(j1);
	  for (int i1 = 0; i1 < mat_elem.GetN(); i1++)
	    mat_elem(i1, j2) = -mat_elem(i1, j2);
	}
    }
  }
  

  void BoundaryIntegralEquation_Maxwell3D
  ::ComputeRightHandSide(Vector<Complex_wp>& b) const
  {
    b.Reallocate(mesh_num.GetNbDof());
    b.Zero();
    VectComplex_wp feval, contrib;
    DISP(wave_vector); 
    for (int i = 0; i < this->mesh.GetNbBoundaryRef(); i++)
      {
	const ElementReference<Dimension2, 3>& Fb = this->GetReferenceElement(i);
	int N = Fb.GetNbDof();
	int Nquad = Fb.GetNbPointsQuadratureInside();

	VectR3 PointsQuadI, NormaleQuadI;
	Vector<TinyMatrix<Real_wp, General, 3, 2> > MatDfI;
	this->var_geom.FillNormaleDfQuad(i, PointsQuadI, NormaleQuadI, MatDfI);
	
	TinyVector<Complex_wp, 3> Einc, vec_u;
	TinyVector<Complex_wp, 2> vec_v;
	feval.Reallocate(2*Nquad); contrib.Reallocate(N);
	for (int j = 0; j < Nquad; j++)
	  {
	    Real_wp arg = DotProd(wave_vector, PointsQuadI(j));
	    Complex_wp u_inc(cos(arg), sin(arg));
	    Einc = u_inc*polarization;
	    TimesProd(Einc, NormaleQuadI(j), vec_u);
	    MltTrans(MatDfI(j), vec_u, vec_v);
	    Real_wp poids = Fb.WeightsND(j);
	    feval(2*j) = poids*vec_v(0);
	    feval(2*j+1) = poids*vec_v(1);
	  }

        //DISP(i); DISP(NormaleQuadI); DISP(MatDfI); DISP(feval);
        
	Fb.ApplyCh(feval, contrib);
	const IVect& neg_dof = this->mesh_num.Element(i).GetNegativeDofNumber();
	for (int j = 0; j < neg_dof.GetM(); j++)
	  contrib(neg_dof(j)) = -contrib(neg_dof(j));

        //DISP(contrib);
	for (int j = 0; j < N; j++)
	  b(mesh_num.Element(i).GetNumberDof(j)) += contrib(j);

	PointsQuadI.Nullify(); NormaleQuadI.Nullify(); MatDfI.Nullify();
      }

    DISP(Norm2(b));
  }

  
  void BoundaryIntegralEquation_Maxwell3D
  ::WriteSolution(const Vector<Complex_wp>& b) const
  {
    Vector<int> NbFct(mesh.GetNbVertices()); NbFct.Zero();
    Vector<TinyVector<Complex_wp, 3> > ValueJ(mesh.GetNbVertices());
    for (int i = 0; i < this->mesh.GetNbBoundaryRef(); i++)
      {
	const ElementReference<Dimension2, 3>& Fb = this->GetReferenceElement(i);
	int N = Fb.GetNbDof();

	const IVect& neg_dof = this->mesh_num.Element(i).GetNegativeDofNumber();
	VectComplex_wp contrib(N);
	for (int j = 0; j < N; j++)
	  contrib(j) = b(mesh_num.Element(i).GetNumberDof(j));

	for (int j = 0; j < neg_dof.GetM(); j++)
	  contrib(neg_dof(j)) = -contrib(neg_dof(j));
	
	int nb_vert = mesh.BoundaryRef(i).GetNbVertices();
	VectR2 s(nb_vert);
	if (nb_vert == 3)
	  {
	    s(0).Init(Real_wp(0), Real_wp(0));
	    s(1).Init(Real_wp(1), Real_wp(0));
	    s(2).Init(Real_wp(0), Real_wp(1));
	  }
	else
	  {
	    s(0).Init(Real_wp(0), Real_wp(0));
	    s(1).Init(Real_wp(1), Real_wp(0));
	    s(2).Init(Real_wp(1), Real_wp(1));
	    s(3).Init(Real_wp(0), Real_wp(1));
	  }

	TinyVector<Complex_wp, 2> vecJ;
	TinyVector<Complex_wp, 3> vecJreel;
	VectR2 phi;
	for (int j = 0; j < nb_vert; j++)
	  {
	    Fb.ComputeValuesPhiRef(s(j), phi);
	    vecJ.Zero();
	    for (int k = 0; k < N; k++)
	      vecJ += contrib(k)*phi(k);

	    Mlt(var_geom.GetDfNodal(i, j), vecJ, vecJreel);
	    vecJreel *= Real_wp(1) / var_geom.GetJacobNodal(i, j);
	    
	    int nv = mesh.BoundaryRef(i).numVertex(j);
            if (nv == 0)
              {
                DISP(s(j)); DISP(i); DISP(j);
                DISP(vecJ); DISP(vecJreel);
              }
            
	    NbFct(nv)++;
	    ValueJ(nv) += vecJreel;
	  }
      }

    for (int i = 0; i < mesh.GetNbVertices(); i++)
      Mlt(Real_wp(1)/NbFct(i), ValueJ(i));

    Vector<Real_wp> ValueJmod(mesh.GetNbVertices());
    DISP(ValueJ(0));
    for (int i = 0; i < mesh.GetNbVertices(); i++)
      ValueJmod(i) = Norm2(ValueJ(i));

    DISP(ValueJmod(0));
    this->mesh.Write("test.mesh");
    //ofstream file_out("test.vtk", ios::app);
    //file_out.precision(15);
    Dimension3 dim;
    
    WriteMedit(ValueJmod, "test.bb", dim, true); 
    //file_out.close();
  }
  
}

#define MONTJOIE_FILE_BOUNDARY_INTEGRAL_EQUATION_CXX
#endif

