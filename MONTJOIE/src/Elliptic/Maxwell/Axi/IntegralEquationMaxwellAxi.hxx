#ifndef MONTJOIE_FILE_INTEGRAL_EQUATION_MAXWELL_AXI_HXX

namespace Montjoie
{  
  //! base class for integral equations for Maxwell's equations in axisymmetric domain
  template<class TypeElement, class TypeEquation>
  class VarIntegralEquationAxi : public VarIntegralEquation_Base<TypeElement, TypeEquation>
  {
  public :    
    typedef typename TypeEquation::Complexe Complexe;
    
    // singular integration
    enum {QUADRATURE_GAUSS , QUADRATURE_HEXA_POLAR, QUADRATURE_DUFFY,
	  QUADRATURE_POLAR, QUADRATURE_GAUSS_SQUARE};
    
    //! file name to write currents
    string name_file_Jb, name_file_Jt, name_file_Kb, name_file_Kt,
      name_file_Jnorme, name_file_Knorme;
    
    //! output in 2-D or 3-D ?
    int type_output_file_ie;
    Mesh<Dimension2> mesh_refined_ie; //!< 2-D surfacic mesh for output
    Mesh<Dimension3> mesh_boundary_inside;
    IVect Index_VertexSurf_to_Vertex; VectReal_wp AngleVertex;
    Matrix<Real_wp, General, ArrayRowSparse> matmass_Jt_IE; //!< mass matrix for integral equation
    Matrix<Real_wp, General, ArrayRowSparse> matmass_Jb_IE; //!< mass matrix for integral equation
    Matrix<Real_wp, General, ArrayRowSparse> Proj_Ch2_IE;
    //!< coupling matrix for integral equation
    Matrix<Real_wp, General, ArrayRowSparse> Proj_Ch1_IE;
    //!< coupling matrix for integral equation
    Matrix<Real_wp, General, ArrayRowSparse> matcrossed_Ch1; 
    //!< coupling matrix for integral equation
    
    // variables used to compute matrix 
    // quadrature formulas for IE (IE = Integral Equation)
    Globatto<Real_wp> lob_IE; //!< Gauss-Lobatto interpolation of integral equation
    Globatto<Real_wp> gauss_reg_IE; 
    //!< quadrature formula for regular parts of integral equation
    Globatto<Real_wp> gauss_sy_IE; 
    //!< quadrature formula for singular parts of integral equation (y = inner integral)
    Globatto<Real_wp> gauss_sx_IE; 
    //!< quadrature formula for singular parts of integral equation (x = outer integral)
    Globatto<Real_wp> gauss_phi_IE; 
    //!< quadrature formula to integrate over phi
    Globatto<Real_wp> gauss_rx_IE; 
    //!< quadrature formula to axis-singular parts of integral equation (x = outer integral)
    Globatto<Real_wp> gauss_ry_IE; 
    //!< quadrature formula to axis-singular parts of integral equation (y = inner integral)
    Globatto<Real_wp> gauss_syPhi_IE; //!< quadrature formula
    // ri (radius) or zi (coordinate z), or t = (tx,tz) (tangential vector)
    VectReal_wp Ri_gauss_IE; //!< radius for regular quadrature points
    VectReal_wp Zi_gauss_IE; //!< z-coordinate for regular quadrature points
    VectReal_wp Tx_gauss_IE; //!< x-coordinate of tangential vector for regular quadrature points
    VectReal_wp Tz_gauss_IE; //!< z-coordinate of tangential vector for regular quadrature points
    // integration weights
    VectReal_wp Ds_gauss_IE; //!< lineic element for regular quadrature points
    VectReal_wp Weight_gauss_IE; //!< integration weights (regular parts)
    
    // orders of integration
    int order_integration_phi; //!< order of integration over phi
    int order_integration_radiusX; //!< order of integration for axis-singular parts
    int order_integration_radiusY; //!< order of integration for axis-singular parts
    // orders of integration automatically chosen ?
    bool automatic_integration_phi; //!< order of integration over phi automatically chosen ?
    bool automatic_integration_radius; 
    //!< order of integration for axis-singular parts, automatically chosen ?
    // Globatto<Real_wp> lobs_IE;
    // Matrix<Real_wp> Val_PhiR_IE, Grad_PhiR_IE, Val_PsiR_IE, Val_PhiL_IE,
    // Grad_PhiL_IE, Val_PsiL_IE;
    // arrays to interpolate basis functions to quadrature points
    Matrix<Real_wp> Val_PhiR_IE; 
    //!< interpolation of Jt basis functions on regular quadrature points
    Matrix<Real_wp> Grad_PhiR_IE; 
    //!< derivative of Jt basis functions on regular quadrature points
    Matrix<Real_wp> Val_PsiR_IE; 
    //!< interpolation of Jb basis functions on regular quadrature points
    Matrix<Real_wp> Val_PhiSx_IE; 
    //!< interpolation of Jt basis functions on singular quadrature points
    Matrix<Real_wp> Grad_PhiSx_IE; 
    //!< derivative of Jt basis functions on singular quadrature points
    Matrix<Real_wp> Val_PsiSx_IE; 
    //!< interpolation of Jb basis functions on singular quadrature points
    Matrix<Real_wp> Val_PhiSy_IE; 
    //!< interpolation of Jt basis functions on singular quadrature points
    Matrix<Real_wp> Grad_PhiSy_IE; 
    //!< derivative of Jt basis functions on singular quadrature points 
    Matrix<Real_wp> Val_PsiSy_IE; 
    //!< interpolation of Jb basis functions on singular quadrature points
    // arrays to interpolate basis functions to quadrature points
    Matrix<Real_wp> Val_PhiRx_IE;
    //!< interpolation of Jt basis functions on axis-singular quadrature points
    Matrix<Real_wp> Grad_PhiRx_IE; 
    //!< derivative of Jt basis functions on axis-singular quadrature points
    Matrix<Real_wp> Val_PsiRx_IE; 
    //!< interpolation of Jb basis functions on axis-singular quadrature points
    Matrix<Real_wp> Val_PhiRy_IE; 
    //!< interpolation of Jt basis functions on axis-singular quadrature points
    Matrix<Real_wp> Grad_PhiRy_IE; 
    //!< derivative of Jt basis functions on axis-singular quadrature points
    Matrix<Real_wp> Val_PsiRy_IE; 
    //!< interpolation of Jb basis functions on axis-singular quadrature points
    //! internal mesh
    Mesh<Dimension2> mesh_sing;
    //! 3-D integration points
    VectR3 Points3D_sing; VectReal_wp Weights3D_sing; //!< 3-D integration weights
    
    //! stiffness matrix for integral equation
    Matrix<Real_wp> matstiff_tt_IE;
    Matrix<Real_wp> matstiff_bb_IE; //!< stiffness matrix for integral equation
    Matrix<Real_wp> matstiff_tb_IE; //!< stiffness matrix for integral equation
    Matrix<Real_wp> matstiff_bt_IE; //!< stiffness matrix for integral equation
    
    // reference to EllipticProblem
    EllipticProblem<TypeElement, TypeEquation>& var_harmonic;
    
    int nb_edges_integral; //!< = nb_boundaries_ref
    int nodl_integral; //!< number of dofs for the integral equation
    int boundary_condition; //!< boundary condition related to the surface
    int nb_dof_Jt; //!< number of dofs on an edge for unknown Jt
    int nb_dof_Jb; //!< number of dofs on an edge for unknown Jb
    int nodl_Jt; //!< number of dofs for unknown Jt (J . t)
    int nodl_Jb; //!< number of dofs for unknown Jb (J . b)
    // matching array edge of integral equation <-> referenced edge on the 2-D mesh
    //                dof  of integral equation <-> dof of finite element method
    // and inverse of these two arrays
    IVect IndexEdgeInt_to_EdgeRef; 
    //!< matching array between edges of surfacic mesh and edges_ref of volumic mesh
    IVect Inverse_EI_ER; //!< reciprocal array of IndexEdgeInt_to_EdgeRef
    
    // matching array edge of integral equation <-> referenced edge on the 2-D mesh
    //                dof  of integral equation <-> dof of finite element method
    // and inverse of these two arrays
    IVect IndexDofInt_to_DofVol; //!< matching array between integral dofs and volumic dofs
    IVect Inverse_DI_DV; //!< reciprocal array of IndexDofInt_to_DofVol
    
    // variables used to display currents on the object
    IVect EdgeInterp_IE; //!< for each point of display grid, edge number
    VectReal_wp CoorInterp_IE; //!<  for each point of display grid, local coordinate in the edge
    int nbPoints_grid_boundary_ie; //!< number of points used to display current on the surface
    VectReal_wp RadiusInterp_IE; //!< radius for each point of display grid
    VectReal_wp ZiInterp_IE; //!< z-coordinate for each point of display grid
    VectBool SignVolDof_IE; //!< dof signs
    VectBool EdgeIE_On_Axe;  //!< edges on axe
    
    // Local to Global numbering dofs for integral equation
    // NodleInt(i,j) returns the global dof number of the local dof i on the edge j
    Matrix<int> NodleInt; //!< numbering array for dofs on the surface
    VectReal_wp Ri_lob_IE; //!< radius on Gauss-Lobatto points
    VectReal_wp Zi_lob_IE; //!< z-coordinate on Gauss-Lobatto points
    VectReal_wp Ds_lob_IE; //!< lineic element on Gauss-Lobatto points
    VectReal_wp Tx_lob_IE; //!< x-component of tangential vector on Gauss-Lobatto points
    VectReal_wp Tz_lob_IE; //!< z-component of tangential vector on Gauss-Lobatto points
#ifdef MONTJOIE_WITH_HIGH_CONDUCTIVITY_MODEL
    VectReal_wp K1_lob_IE; //!< interpolation of curvatures K2
    VectReal_wp K2_lob_IE; //!< interpolation of curvatures K1 
#endif
    // VectReal_wp RadiusDof_IE, ZiDof_IE, DsDof_IE, TxDof_IE, TzDof_IE;
    Real_wp space_step_IE; //!< space step

    VarIntegralEquationAxi(EllipticProblem<TypeElement, TypeEquation>& var);
    
    void SetInputData(const string& description_field, const VectString& parameters);
    
    void GetDofsBoundary(int ref, SurfacicMesh<Dimension2>& mesh_surf);
    void ComputeSectionGrid_IE(SurfacicMesh<Dimension2>& mesh_surf );
    void ConstructSurfacicMesh_IE(int ref, SurfacicMesh<Dimension2>& mesh_surf,
                                  Mesh<Dimension2>& mesh_refined);
    void WriteOutputFile_CurrentsJ(const VectComplex_wp& U0, const string& output_Jt,
				   const string& output_Jb, const string& output_Jnorme,
                                   bool IsUnknownJ);
    void AddInterpolated_Currents2D(const VectComplex_wp& U0, VectComplex_wp& output_Jt,
                                    VectComplex_wp& output_Jb,
				    const SurfacicMesh<Dimension2>& mesh_surf,
                                    bool IsUnknown_J) const;
    
    void WriteOutputFile_Currents2D(const VectComplex_wp& U0,
                                    const SurfacicMesh<Dimension2>& mesh_surf,
				    const string& name_Jt, const string& name_Jb, bool UnknownJ);
    
    void WriteOutputFile_Currents3D(const VectComplex_wp& U0, const Mesh<Dimension2>& mesh2d,
                                    const Mesh<Dimension3>& mesh3d, 
				    const string& name_Jt, const string& name_Jb,
                                    const string& name_Jnorme, bool IsUnknown_J);
    
    //! returns global dof number for Jt
    int GetNumberDofT(int num_edge, int node) const { return NodleInt(num_edge, node); }
    //! returns global dof number for Jb
    int GetNumberDofB(int num_edge, int node) const { return NodleInt(num_edge, nb_dof_Jt+node); }
    void ComputeInterpolatedVal(Real_wp& ri, Real_wp& zi, Real_wp& txi, Real_wp& tzi,
				Real_wp& dsi, const Globatto<Real_wp>& lob,
				const VectReal_wp& Ri_lob, const VectReal_wp& Zi_lob,
                                const VectReal_wp& Tx_lob, 
				const VectReal_wp& Tz_lob, const VectReal_wp& Ds_lob,
                                int offset1, const Real_wp& x) const;
    
    void ComputeG1_Gcos_Gsin_EFIE(Globatto<Real_wp>& gauss,
				  const Real_wp& phi0, const Real_wp& phi1,
				  Real_wp& dist, Real_wp& prodr, Complex_wp& G1,
				  Complex_wp& Gcos, Complex_wp& Gsin) const;
    
    void ComputeG1_Gcos_Gsin_MFIE(Globatto<Real_wp>& gauss,
				  const Real_wp& phi0, const Real_wp& phi1,
				  Real_wp& dist, Real_wp& prodr, Complex_wp& G1,
				  Complex_wp& Gcos, Complex_wp& Gsin) const;
    
    void ComputeG1_Gcos_Gsin_CFIE(const Globatto<Real_wp>& gauss,
				  const Real_wp& phi0, const Real_wp& phi1, 
                                  Real_wp& dist, Real_wp& prodr, 
				  Complex_wp& G1_efie,	Complex_wp& Gcos_efie,
                                  Complex_wp& Gsin_efie,
				  Complex_wp& G1_mfie,	Complex_wp& Gcos_mfie,
                                  Complex_wp& Gsin_mfie) const;
    
    void ComputeCoefficients_EFIE(Real_wp& ri, Real_wp& zi, Real_wp& txi,
                                  Real_wp& tzi, Real_wp& dsi, 
				  Real_wp& rj, Real_wp& zj, Real_wp& txj,
                                  Real_wp& tzj, Real_wp& dsj,
				  Real_wp& prodr, Real_wp& prod_ds, Complex_wp& im, Complex_wp& m2, 
				  Complex_wp& G1, Complex_wp& Gcos,
                                  Complex_wp& Gsin, Complex_wp& A1, Complex_wp& A2,
				  Complex_wp& A3, Complex_wp& A4,
                                  Complex_wp& A5, Complex_wp& A6,Complex_wp& A7, 
				  Complex_wp& A8, Complex_wp& A9) const;
    
    void ComputeCoefficients_MFIE(Real_wp& ri, Real_wp& zi, Real_wp& txi,
                                  Real_wp& tzi, Real_wp& dsi, 
				  Real_wp& rj, Real_wp& zj, Real_wp& txj,
                                  Real_wp& tzj, Real_wp& dsj,
				  Real_wp& diffz, Real_wp& prodr, Complex_wp& G1,
                                  Complex_wp& Gcos, Complex_wp& Gsin, 
				  Complex_wp& F11, Complex_wp& F12,
                                  Complex_wp& F21, Complex_wp& F22) const;
    
    void UpdateEFIEMatrices(int i, int j, const Complex_wp& A1, const Complex_wp& A2,
                            const Complex_wp& A3, const Complex_wp& A4,
			    const Complex_wp& A5, const Complex_wp& A6,
                            const Complex_wp& A7, const Complex_wp& A8, const Complex_wp& A9,
			    const Matrix<Real_wp>& Val_Phi1,
                            const Matrix<Real_wp>& Grad_Phi1, const Matrix<Real_wp>& Val_Psi1,
			    const Matrix<Real_wp>& Val_Phi2, const Matrix<Real_wp>& Grad_Phi2,
                            const Matrix<Real_wp>& Val_Psi2,
			    Matrix<Complex_wp>& mat_tt, Matrix<Complex_wp>& mat_bb,
			    Matrix<Complex_wp>& mat_tb, Matrix<Complex_wp>& mat_bt) const;
    
    void UpdateMFIEMatrices(int i, int j, const Complex_wp& F11, const Complex_wp& F12,
                            const Complex_wp& F21, const Complex_wp& F22,
			    const Matrix<Real_wp>& Val_Phi1, const Matrix<Real_wp>& Val_Psi1,
			    const Matrix<Real_wp>& Val_Phi2, const Matrix<Real_wp>& Val_Psi2,
			    Matrix<Complex_wp>& mat_tt, Matrix<Complex_wp>& mat_bb,
			    Matrix<Complex_wp>& mat_tb, Matrix<Complex_wp>& mat_bt) const;
    
    void UpdateEFIEMatrices(const Complex_wp& A1, const Complex_wp& A2,
                            const Complex_wp& A3, const Complex_wp& A4,
			    const Complex_wp& A5, const Complex_wp& A6,
                            const Complex_wp& A7, const Complex_wp& A8, const Complex_wp& A9,
			    const VectReal_wp& Val_Phi1, const VectReal_wp& Grad_Phi1,
                            const VectReal_wp& Val_Psi1,
			    const VectReal_wp& Val_Phi2, const VectReal_wp& Grad_Phi2,
                            const VectReal_wp& Val_Psi2,
			    Matrix<Complex_wp>& mat_tt, Matrix<Complex_wp>& mat_bb,
			    Matrix<Complex_wp>& mat_tb, Matrix<Complex_wp>& mat_bt) const;
    
    void UpdateMFIEMatrices(const Complex_wp& F11, const Complex_wp& F12,
                            const Complex_wp& F21, const Complex_wp& F22,
			    const VectReal_wp& Val_Phi1, const VectReal_wp& Val_Psi1,
			    const VectReal_wp& Val_Phi2, const VectReal_wp& Val_Psi2,
			    Matrix<Complex_wp>& mat_tt, Matrix<Complex_wp>& mat_bb,
			    Matrix<Complex_wp>& mat_tb, Matrix<Complex_wp>& mat_bt) const;

    Real_wp EvaluatePhi1(const Real_wp& r1, const Real_wp& z1, const Real_wp& r2,
                         const Real_wp& z2, const Real_wp& Rmax) const;
    
    void PerformRegularIntegration(const Real_wp& phi_init, Globatto<Real_wp>& gauss_X,
                                   Globatto<Real_wp>& gaussY,
				   Globatto<Real_wp>& gauss_phi, Matrix<Real_wp>& Val_PhiX,
                                   Matrix<Real_wp>& Val_PsiX, Matrix<Real_wp>& Grad_PhiX,
				   Matrix<Real_wp>& Val_PhiY, Matrix<Real_wp>& Val_PsiY,
                                   Matrix<Real_wp>& Grad_PhiY,
				   int offset1, int offset2, Matrix<Complex_wp>& mat_tt,
                                   Matrix<Complex_wp>& mat_bb,
				   Matrix<Complex_wp>& mat_tb, Matrix<Complex_wp>& mat_bt,
                                   Matrix<Complex_wp>& coup11, Matrix<Complex_wp>& coup22,
				   Matrix<Complex_wp>& coup12, Matrix<Complex_wp>& coup21);
    
    void PerformGaussSquaredIntegration(const Real_wp& phi_init, Globatto<Real_wp>& gauss_X,
                                        Globatto<Real_wp>& gaussY,
					Globatto<Real_wp>& gauss_phi, Matrix<Real_wp>& Val_PhiX,
                                        Matrix<Real_wp>& Val_PsiX, Matrix<Real_wp>& Grad_PhiX,
                                        int offset1, int offset2, Matrix<Complex_wp>& mat_tt,
                                        Matrix<Complex_wp>& mat_bb,
					Matrix<Complex_wp>& mat_tb, Matrix<Complex_wp>& mat_bt,
                                        Matrix<Complex_wp>& coup11, Matrix<Complex_wp>& coup22,
					Matrix<Complex_wp>& coup12, Matrix<Complex_wp>& coup21);
    
    void PerformSingularIntegration(const Real_wp& phi_init, Globatto<Real_wp>& gauss_X,
                                    Globatto<Real_wp>& gauss_Y,
				    VectR2& CoorQuadY, VectReal_wp& WeightsQuadY,
				    int offset1, int offset2, Matrix<Complex_wp>& mat_tt,
                                    Matrix<Complex_wp>& mat_bb,
				    Matrix<Complex_wp>& mat_tb, Matrix<Complex_wp>& mat_bt,
                                    Matrix<Complex_wp>& coup11, Matrix<Complex_wp>& coup22,
				    Matrix<Complex_wp>& coup12, Matrix<Complex_wp>& coup21);
    
    void InitComputation_IntegralEquation();
    
    template<class MatrixIntegral>
    void ComputeMatrixIntegralEquation(MatrixIntegral& mat_integral);
    
    template<class MatrixIntegral>
    void ComputeMatrixIntegralEquation(MatrixIntegral& mat_integral, Matrix<Complex_wp>& M11,
                                       Matrix<Complex_wp>& M12,
                                       Matrix<Complex_wp>& M21, Matrix<Complex_wp>& M22);

    void ComputeMatrixCoupling_FiniteElement(Matrix<Complex_wp>& mat_integral,
                                             const Matrix<Complex_wp, Symmetric, RowSymPacked>&
                                             schur_matrix);
    
    template<class MatrixIntegral>
    void ComputeMatrixCoupling_IntegralEquation(MatrixIntegral& mat_integral,
                                                const Matrix<Complex_wp>& M11,
                                                const Matrix<Complex_wp>& M12,
						const Matrix<Complex_wp>& M21,
                                                const Matrix<Complex_wp>& M22);
    
    void GetQuadratureFormulas_Singularity(const Real_wp& xs, const Real_wp& phi1,
                                           Mesh<Dimension2>& mesh_sing,
					   Globatto<Real_wp>& gauss_Y,
                                           VectR2& CoorQuadY, VectReal_wp& WeightsQuadY);
    
    void ComputeRightHandSideIntegralEquation(VectComplex_wp& rhs_integral, int num_angle);
    
    void ComputeQuadratureFormula_AxiSingularity(const Globatto<Real_wp>& gauss_reg,
                                                 VectR3& Points3D, VectReal_wp& Weights3D) const;
    
    void GetComponentsSolution_IntegralEquation(VectComplex_wp& full_sol, VectComplex_wp& J);
    void GetComponentsSolution_CouplingSystem(VectComplex_wp& full_sol, VectComplex_wp& J,
                                              VectComplex_wp& K);
    
    void ComputeSolution_MultiIncidence(Matrix<Complex_wp>& mat_integral,
                                        IVect& pivot_IE, VectComplex_wp& rhs_integral, 
					VectR3_Complex_wp& RCS_monostat_horizontal,
                                        VectR3_Complex_wp& RCS_monostat_vertical,
					VectReal_wp& time_process);
    
    void AddMassTerms(Globatto<Real_wp>& gauss_reg, Matrix<Real_wp>& Val_PhiR,
		      Matrix<Real_wp>& Val_PsiR, int offset1,
		      Matrix<Complex_wp>& mat_tt, Matrix<Complex_wp>& mat_bb,
		      Matrix<Complex_wp>& coup11, Matrix<Complex_wp>& coup22);
  
    void ComputeMassMatrix_IntegralEquation(int Nt, int Nb, const Globatto<Real_wp>& lob,
                                            const Globatto<Real_wp>& gauss_reg,
					    const VectReal_wp& Ri_lob, const VectReal_wp& Zi_lob,
                                            const VectReal_wp& Tx_lob, const VectReal_wp& Tz_lob,
					    const VectReal_wp& Ds_lob,
                                            const Matrix<Real_wp>& Val_PhiR,
                                            const Matrix<Real_wp>& Val_PsiR);
    
  };

}

#define MONTJOIE_FILE_INTEGRAL_EQUATION_MAXWELL_AXI_HXX
#endif

