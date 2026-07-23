#define MONTJOIE_WITH_ONE_DIM
#define MONTJOIE_WITH_TWO_DIM
// #define MONTJOIE_WITH_THREE_DIM

#define MONTJOIE_WITH_NODAL_H1
#define SELDON_WITH_PRECONDITIONING

#include "Elliptic/Helmholtz/MontjoieHelmholtz.hxx"
//#include "Elliptic/Helmholtz/Helmholtz1D.hxx"
//#include "Elliptic/Helmholtz/Helmholtz1D.cxx"

using namespace Montjoie;

// Variable globale pour la source
Real_wp cg=0.0;
Real_wp wg=0.0;
bool filtrage;
int condinterface=0; // 0 : dirichlet, 1 : impedance

template<class Dimension>
class SpaceFileSource
{};

//! source in space, whose values are read in a file
template<>
class SpaceFileSource<Dimension2>
{
public :
  VectReal_wp values;
  Globatto<Real_wp> lob; //!< interpolation data
  Real_wp xmin, xmax, dx;
    
  void Init(const string& input_file, const Real_wp& x0, const Real_wp& x1);
    
  Real_wp Evaluate(const Real_wp& x);
  
};


template<>
class SpaceFileSource<Dimension3>
{
public :
  
};

//! initialization of the source
void SpaceFileSource<Dimension2>::Init(const string& input_file,
				       const Real_wp& x0, const Real_wp& x1)
{
  values.ReadText(input_file);
  
  // cubic interpolation
  VectReal_wp Points(4); Points.Fill(); Mlt(1.0/3.0, Points);
  lob.AffectPoints(Points);
    
  // time steps
  xmin = x0;
  xmax = x1;
  dx = (xmax - xmin)/(values.GetM()-1);
}
  
//! evaluation of the source at position x
/*!
  \param[in] x fixed position where the source has to be computed
  \returns value of the source
*/
Real_wp SpaceFileSource<Dimension2>::Evaluate(const Real_wp& x)
{
  if ((x < xmin)||(x > xmax))
    return 0.0;
  
  // finding interval where x is
  int n0 = toInteger(floor((x-xmin)/(3.0*dx))); n0 = max(0, n0);
  int i = 3*n0, j = i + 3;
  if (j >= values.GetM())
    {
      j = values.GetM() - 1;
      i = j-3;
    }
    
  Real_wp lambda = (x - dx*i - xmin)/(3.0*dx);
  // cubic interpolation
  Real_wp f = 0;
  for (int k = 0; k < 4; k++)
    f += lob.EvaluatePhi(k, lambda)*values(i+k);
    
  return f;
}

//! classe pour definir la source
template<class Dimension>
class UserSourceDirichlet : public VirtualSourceFEM<Complex_wp, Dimension>
{
public :
  typedef typename Dimension::R_N R_N;
  
  SpaceFileSource<Dimension> src_file;
  bool source_on_file;
  
  //! reference de la frontiere ou on impose dirichlet
  int ref_dirichlet;
  
  //! constructor with given problem
  template<class TypeEquation>
  UserSourceDirichlet(EllipticProblem<TypeEquation>& var)
    : VirtualSourceFEM<Complex_wp, Dimension>(var)
  {
    source_on_file = false;
    if (Dimension::dim_N==2)
      ref_dirichlet = 3;
    else
      ref_dirichlet = 4;
  }
  
  bool IsNonNull_SourceDirichlet(int cond)
  {
    // return false;
    return (ref_dirichlet==this->ref_boundary_);
  }
  
  //! Dirichlet inhomogene
  void EvaluateFunction(int i, int j, const R_N& X, VectComplex_wp& f)
  {
    f.Fill(0);
    if (this->ref_boundary_ == ref_dirichlet)
      {
	if (source_on_file)
	  f(0) = src_file.Evaluate(X(0));
	else
          {
            f(0)= 0.0;
            if(X(1)==0)
              {
                if(X(0)!=0)
                  f(0)=sin(wg*X(0)/cg)/X(0);
                else
                  f(0)=wg/cg;
              }
            //f(0) = exp(-(2000/(2*12))*(2000/(2*12))*X(0)*X(0)); // on suppose xs=0          
          }
	 //f(0) = exp(-0.5*X(0)*X(0));
	// f(0) = exp(-0.5*(X(0)*X(0)+X(1)*X(1)));
      }
  }
  
  bool IsNonNull_SourceSurfacic(int cond_ref)
  {
    return false;
    // return (ref_dirichlet==this->ref_boundary_);
  }

  //! source de type du/dn + T u = f
  void EvaluateSourceSurfacic(int k, const SetPoints<Dimension>& PointsElem,
			      const SetMatrices<Dimension>& MatricesElem, VectComplex_wp& f)
  {
    typename Dimension::R_N point = PointsElem.GetPointQuadratureBoundary(k);
    f(0) = GetFunction(point(0));
  }
  
};


// on recupere les coordonnees et les numeros des ddls places a X(num_coor) = cte
template<class TypeEquation, class Vector1>
void GetDofFixedCoordinate(const EllipticProblem<TypeEquation>& var,
			   Vector1& PointsDof, IVect& NumDof, const Real_wp& cte, int num_coor)
{
  typedef typename TypeEquation::Dimension Dimension;
  typedef typename Dimension::R_N R_N;
  const Mesh<Dimension>& mesh = var.mesh;
  SetPoints<Dimension> PointsElem;
  Vector<R_N> s; R_N point;
  NumDof.Reallocate(var.GetNbDof());
  PointsDof.Reallocate(var.GetNbDof());
  int nb_dof_found = 0;
  // boucle sur les elements
  for (int i = 0; i < mesh.GetNbElt(); i++)
    {
      // sommets de l'element
      mesh.GetVerticesElement(i, s);
      
      // calcul des coordonnees des ddls
      var.GetReferenceElement(i).FjElemDof(s, PointsElem, mesh, i);
      
      // boucle sur les ddls
      int nb_dof_loc = var.mesh_num.GetNbLocalDof(i);
      for (int j = 0; j < nb_dof_loc; j++)
	{
	  point = PointsElem.GetPointDof(j);
	  if (abs(point(num_coor) - cte) < 1e-7)
	    {
	      // on a trouve un ddl sur l'hyperplan considere
	      NumDof(nb_dof_found) = var.mesh_num.Element(i).GetNumberDof(j);
	      PointsDof(nb_dof_found) = point;
	      nb_dof_found++;
	    }
	}
    }
  
  // on trie et enleve les doublons
  RemoveDuplicate(nb_dof_found, PointsDof, NumDof);
  
  // on ajuste la taille
  PointsDof.Resize(nb_dof_found); NumDof.Resize(nb_dof_found);
  
}


template<class TypeEquation, class Vector1>
void GetDuDy(EllipticProblem<TypeEquation>& var,
	     const Vector1& x_sol, const IVect& NumDof, Vector1& v, Vector1& dv_dy)
{
  typedef typename TypeEquation::Dimension Dimension;
  typedef typename TypeEquation::Complexe Complexe;
  
  Vector<Vector<TinyVector<Complexe, 1> > > Unodal;
  Vector<Vector<TinyVector<Complexe, Dimension::dim_N> > > GradNodal;
  
  var.ComputeNodalUgradU(x_sol, Unodal, GradNodal, false, true);
  
  IVect InverseDof(var.GetNbDof()); InverseDof.Fill(-1);
  for (int i = 0; i < NumDof.GetM(); i++)
    InverseDof(NumDof(i)) = i;
  
  dv_dy.Reallocate(NumDof.GetM());
  v.Reallocate(NumDof.GetM());
  // Correction -var.GetNbEltPML() pour calculer hors des pml
  for (int i = 0; i < var.mesh.GetNbElt()-var.GetNbEltPML(); i++)
    for (int j = 0; j < GradNodal(i).GetM(); j++)
      if (InverseDof(var.mesh_num.Element(i).GetNumberDof(j)) != -1)
	{
	  int num_dof = InverseDof(var.mesh_num.Element(i).GetNumberDof(j));
	  dv_dy(num_dof) = GradNodal(i)(j)(Dimension::dim_N-1);
	}
  
  for (int i = 0; i < NumDof.GetM(); i++)
    v(i) = x_sol(NumDof(i));
}


template<class TypeEquation>
void GetWeight(EllipticProblem<TypeEquation>& var, const IVect& NumDof,
	       const Real_wp& cte, int num_coor, VectReal_wp& weight)
{
  typedef typename TypeEquation::Dimension Dimension;
  IVect InverseDof(var.GetNbDof()); InverseDof.Fill(-1);
  for (int i = 0; i < NumDof.GetM(); i++)
    InverseDof(NumDof(i)) = i;
  
  weight.Reallocate(NumDof.GetM()); weight.Fill(0);
  for (int i = 0; i < var.mesh.GetNbBoundaryRef(); i++)
    {
      int nv1 = var.mesh.BoundaryRef(i).numVertex(0);
      int nv2 = var.mesh.BoundaryRef(i).numVertex(1);
      if ((abs(var.mesh.Vertex(nv1)(num_coor) - cte) < 1e-7)&&
	  (abs(var.mesh.Vertex(nv2)(num_coor) - cte) < 1e-7))
	{
	  int num_face = i;
	  int num_elem = var.mesh.Boundary(i).numElement(0);
	  int num_loc = var.mesh.Element(num_elem).GetPositionBoundary(num_face);
	  const ElementReference<Dimension, 1>& Fb = var.GetReferenceElement(num_elem);

	  for (int j = 0; j < Fb.GetNbDofBoundary(num_loc); j++)
	    {
	      int num_dof_loc = Fb.GetLocalNumber(num_loc, j);
	      int num_dof = var.mesh_num.Element(num_elem).GetNumberDof(num_dof_loc);
	      if (InverseDof(num_dof) != -1)
		{
		  int k = InverseDof(num_dof);
		  Real_wp length = var.mesh.Vertex(nv1).Distance(var.mesh.Vertex(nv2));
		  weight(k) += length*Fb.WeightsQuadratureBoundary(j, num_loc);
		}
	    }
	}


    }
}

template<class TypeEquation>
void GetInterpolate1D(EllipticProblem<TypeEquation>& var, const VectComplex_wp& V0,
                      VectComplex_wp& V0_interp, const IVect& ElementInterp, const VectReal_wp& CoorInterp)
{
  typedef typename TypeEquation::Dimension Dimension;
  int order = var.mesh.GetOrder();
  VectReal_wp phi;
  
  //  cout << "boucle jusqu'a : " << V0_interp.GetM() << endl;
  for (int ii = 0; ii < V0_interp.GetM(); ii++)
    {           
      int num_elem = ElementInterp(ii);
      const ElementReference<Dimension, 1>& Fb = var.GetReferenceElement(num_elem);
      //if(num_elem==65)
      //  cout << "   !!!!!!! ii = " << ii << endl;
      //cout << "ii = " << ii << "    et num_elem  == " << num_elem << endl;
      Real_wp xloc = CoorInterp(ii);
      
      Fb.ComputeValuesPhi(xloc, phi);
      Complex_wp valu(0.0);
      for (int k = 0; k < phi.GetM(); k++)
        {
          //cout << " k = " << k << endl;
          valu += phi(k)*V0(var.mesh.GetNumberDof(num_elem, k));
        }
      
      V0_interp(ii) = valu;
    }
}

void WriteGrid2D(const Real_wp& xmin, const Real_wp& xmax, const Real_wp& ymin,
                 const Real_wp& ymax, int nbPoints_x, int nbPoints_y,
                 const VectReal_wp& Voutput, const string& file_name)
{
  ofstream FileStream(file_name.data());
  int itmp = 1; FileStream.write(reinterpret_cast<char*>(&itmp),sizeof(int));
  typedef Real_wp real; real tmp;
  tmp = xmin; FileStream.write(reinterpret_cast<char*>(&tmp),sizeof(real));
  tmp = xmax; FileStream.write(reinterpret_cast<char*>(&tmp),sizeof(real));
  FileStream.write(reinterpret_cast<char*>(const_cast<int*>(&nbPoints_x)),sizeof(int));
  tmp = ymin; FileStream.write(reinterpret_cast<char*>(&tmp),sizeof(real));
  tmp = ymax; FileStream.write(reinterpret_cast<char*>(&tmp),sizeof(real));
  FileStream.write(reinterpret_cast<char*>(const_cast<int*>(&nbPoints_y)),sizeof(int));
  
  Voutput.Write(FileStream);
  
  FileStream.close();
}


class SourceDirichlet1D
{
public :
  void EvaluateFunction(const Real_wp& x, VectComplex_wp& f)
  {
    // f(0) = 100*exp(-(2000/(2*12))*(2000/(2*12))*x*x);             
    //    Real_wp radius = 4;
    if(filtrage)
      {
        if(x==0)
          f(0)=wg/cg;
        else
          f(0)=sin(wg*x/cg)/x;
      }
    else
      f(0) = 0.0; //exp(-7*x*x/16);
      //    f(0) = exp(-7.0*x*x/(radius*radius));
  }
};

class MyData : public InputDataProblem_Base
{
public :
  int nz;
  int nx;
  int nf;
  int ncouches;
  int nit; 
  int nmat;  
  Vector<bool> paraxial;
  bool filtrage; // filtrage de la source ?
  Real_wp vitesse;
  Real_wp rho;
  Real_wp mu;
  Real_wp zmax; // profondeur max
  Real_wp xmin, xmax; // largeur en x
  bool homogene;
  
  void SetInputData(const string& keyword, const Vector<string>& parameters)
  {
    if (!keyword.compare("NombrePointsZ"))
      {
        nz = to_num<int>(parameters(0));
      }
    
    if (!keyword.compare("NombrePointsX"))
      {
	nx = to_num<int>(parameters(0));
      }

    if(!keyword.compare("NbreFreq"))
      {
        nf = to_num<int>(parameters(0));
      }
    if(!keyword.compare("CondInterface"))
      {
        condinterface = to_num<int>(parameters(0));
      }
    if(!keyword.compare("Paraxial"))
      {
        int tmp;
        
        int np= to_num<int>(parameters(0));
        paraxial.Reallocate(np);
        for(int ii=0;ii<np;ii++)
          {
            tmp= to_num<int>(parameters(ii+1));
            paraxial(ii)=(tmp==1);
          }        
      }
    if(!keyword.compare("Xmin"))
      {
	xmin = to_num<Real_wp>(parameters(0));
      }
    if(!keyword.compare("Xmax"))
      {
	xmax = to_num<Real_wp>(parameters(0));
      }
    if(!keyword.compare("Zmax"))
      {
	zmax = to_num<Real_wp>(parameters(0));
      }
    if(!keyword.compare("Rho"))
      {
	rho = to_num<Real_wp>(parameters(0));
	if(mu!=0)
	  vitesse=sqrt(mu/rho);
      }
    if(!keyword.compare("Mu"))
      {
	mu = to_num<Real_wp>(parameters(0));
	if(rho!=0)
	  vitesse=sqrt(mu/rho);
      }
    if (!keyword.compare("NombreCouches"))
      {
	ncouches=to_num<int>(parameters(0));
      }
    if (!keyword.compare("NombreIterations"))
      {
	nit=to_num<int>(parameters(0));
      }
    if (!keyword.compare("Nmat"))
      {
	nmat=to_num<int>(parameters(0));
      }
    if(!keyword.compare("Filtrage"))
      {
        int tmp= to_num<int>(parameters(0));
        filtrage=(tmp==1);
      }
    if(!keyword.compare("Homogene"))
      {
        int tmp= to_num<int>(parameters(0));
        homogene=(tmp==1);
      }
  }
};


void SplineInterp(const VectReal_wp & Xe, const VectComplex_wp & Ue, const VectReal_wp & Xs, VectComplex_wp & Us) 
{
  int Ne=Xe.GetM();  // Nbre de points en entree
  int Ns=Xs.GetM();  // Nbre de points en sortie
  Matrix<Complex_wp, General, ArrayRowSparse> A(Ne-2,Ne-2);

  VectReal_wp hx(Ne-1);
  VectComplex_wp F(Ne-2);
  
  Us.Reallocate(Ns);  
  
  Vector<int> ElementInterp(Ns);

  bool dedans=false;   

  // Determination de la correspondant des elements
  for(int ii=0;ii<Ns-1;ii++)
    {
      dedans=false;
      int l=-1;
      while((dedans==false) && (l < Ne-1))
        {
          l++;
          if((Xs(ii) >= Xe(l)) && (Xs(ii) < Xe(l+1)))
            {
              dedans=true;
            }          
        }
        ElementInterp(ii)=l;
    }
  
  ElementInterp(Ns-1)=Ne-2;
    
  // Determination du pas du maillage
  for(int ii=0; ii<Ne-1;ii++)
    {
      hx(ii)=Xe(ii+1)-Xe(ii);
    }
  
  for(int ii=1; ii <Ne-2;ii++)
    {
      F(ii-1) = (6/hx(ii))*(Ue(ii+1)-Ue(ii)) - (6/hx(ii-1))*(Ue(ii)-Ue(ii-1));
      A.Get(ii-1,ii-1) = 2*(hx(ii)+hx(ii-1));
      if(ii-1 > 1)
        A.Get(ii-1,ii-2) = hx(ii-1);
      
      if(ii-1 < Ne-2)
        A.Get(ii-1,ii) = hx(ii);
    }
  
  VectComplex_wp M(Ne);

  MatrixMumps<Complex_wp> mat_lu;
  mat_lu.HideMessages();
  Seldon::GetLU(A, mat_lu);
  SolveLU(mat_lu, F); 

  M(0)=0.0;
  M(Ne-1)=0.0;
  for(int ii=0;ii<Ne-2;ii++)
    M(ii+1)=F(ii);
  
  // Us(Ns-1)=Ue(Ne-1);

  // Calcul de l'interpolation par spline cubique   
  for(int ii=0;ii<Ns;ii++)
    {
      int ind=ElementInterp(ii);
      Real_wp x=Xs(ii);
      Us(ii)=M(ind+1)*((x-Xe(ind))*(x-Xe(ind))*(x-Xe(ind)))/(6*hx(ind))+ (M(ind+1)*(Xe(ind+1)-x)*(Xe(ind+1)-x)*(Xe(ind+1)-x))/(6*hx(ind)) + Ue(ind) - (M(ind)*hx(ind)*hx(ind))/Real_wp(6.0)+(Ue(ind+1)-Ue(ind) - (hx(ind)*hx(ind)/6)*(M(ind+1)-M(ind)))*(x-Xe(ind))/hx(ind);
    }

}

void MatrixParaxial(Complex_wp Alpha, Complex_wp Beta, VectComplex_wp Damp, Real_wp hp, Real_wp wg, Real_wp dz, int Nx, VectReal_wp c, Matrix<Complex_wp, General, ArrayRowSparse> & Aplus,Matrix<Complex_wp, General, ArrayRowSparse> & Amoins)
{
  // Allocation memoire pour les 2 matrices utilisees
  Aplus.Reallocate(Nx+8,Nx+8);
  Aplus.Fill(0.0);
  Amoins.Reallocate(Nx+8,Nx+8);
  Amoins.Fill(0.0);
  
  for(int kk=0; kk<Nx+8; kk++)
    {
      Aplus.Get(kk,kk) = 1.0-(2.0/(hp*hp))*(c(kk)*c(kk)/(wg*wg)*Alpha - 0.5*c(kk)*Iwp*Beta*dz/wg); //1-Iwp*2*c*dz/(4*wg*hp*hp);
      Amoins.Get(kk,kk) = 1.0-(2.0/(hp*hp))*(c(kk)*c(kk)/(wg*wg)*Alpha + 0.5*c(kk)*Iwp*Beta*dz/wg); //1-Iwp*2*c*dz/(4*wg*hp*hp);
      
      if(kk>0)
        {
          Aplus.Get(kk,kk-1) = 1.0/(hp*hp)*(c(kk)*c(kk)/(wg*wg)*Alpha - 0.5*c(kk)*Iwp*Beta*dz/wg);  //Iwp*c*dz/(4*wg*hp*hp);
          Amoins.Get(kk,kk-1) = 1.0/(hp*hp)*(c(kk)*c(kk)/(wg*wg)*Alpha + 0.5*c(kk)*Iwp*Beta*dz/wg);  //Iwp*c*dz/(4*wg*hp*hp);                      
        }
      
      if(kk<Nx+7)
        {
          Aplus.Get(kk,kk+1) = 1.0/(hp*hp)*(c(kk)*c(kk)*Alpha/(wg*wg) - 0.5*c(kk)*Iwp*Beta*dz/wg); //Iwp*c*dz/(4*wg*hp*hp);
          Amoins.Get(kk,kk+1) = 1.0/(hp*hp)*(c(kk)*c(kk)*Alpha/(wg*wg) + 0.5*c(kk)*Iwp*Beta*dz/wg); //Iwp*c*dz/(4*wg*hp*hp);                      
        }      
    }

  // Prise en compte des PML
  for(int kk=0;kk<4;kk++)
    {
      Aplus.Get(kk,kk) = 0.5*((1.0/Damp(kk))+(1.0/Damp(kk+1))) - (Damp(kk)+Damp(kk+1))/(hp*hp)*(c(kk)*c(kk)*Alpha/(wg*wg) - 0.5*c(kk)*Iwp*Beta*dz/wg);                  
      Aplus.Get(kk,kk+1) = (Damp(kk+1)/(hp*hp)*(c(kk)*c(kk)*Alpha/(wg*wg) - 0.5*c(kk)*Iwp*Beta*dz/wg));
      Aplus.Get(kk+1,kk) = (Damp(kk+1)/(hp*hp)*(c(kk)*c(kk)*Alpha/(wg*wg) - 0.5*c(kk)*Iwp*Beta*dz/wg));
      
      Amoins.Get(kk,kk) = 0.5*((1.0/Damp(kk))+(1.0/Damp(kk+1))) - (Damp(kk)+Damp(kk+1))/(hp*hp)*(c(kk)*c(kk)*Alpha/(wg*wg) + 0.5*c(kk)*Iwp*Beta*dz/wg);
      Amoins.Get(kk,kk+1) = (Damp(kk+1)/(hp*hp)*(c(kk)*c(kk)*Alpha/(wg*wg) + 0.5*c(kk)*Iwp*Beta*dz/wg));                      
      Amoins.Get(kk+1,kk) = (Damp(kk+1)/(hp*hp)*(c(kk)*c(kk)*Alpha/(wg*wg) + 0.5*c(kk)*Iwp*Beta*dz/wg));
      
      Aplus.Get(Nx+7-kk,Nx+7-kk) = 0.5*((1.0/Damp(kk))+(1.0/Damp(kk+1))) - (Damp(kk)+Damp(kk+1))/(hp*hp)*(c(kk)*c(kk)*Alpha/(wg*wg) - 0.5*c(kk)*Iwp*Beta*dz/wg);                  
      Aplus.Get(Nx+7-kk,Nx+7-kk-1) = (Damp(kk+1)/(hp*hp)*(c(kk)*c(kk)*Alpha/(wg*wg) - 0.5*c(kk)*Iwp*Beta*dz/wg));
      Aplus.Get(Nx+7-kk-1,Nx+7-kk) = (Damp(kk+1)/(hp*hp)*(c(kk)*c(kk)*Alpha/(wg*wg) - 0.5*c(kk)*Iwp*Beta*dz/wg));
      
      Amoins.Get(Nx+7-kk,Nx+7-kk) = 0.5*((1.0/Damp(kk))+(1.0/Damp(kk+1))) - (Damp(kk)+Damp(kk+1))/(hp*hp)*(c(kk)*c(kk)*Alpha/(wg*wg) + 0.5*c(kk)*Iwp*Beta*dz/wg);
      Amoins.Get(Nx+7-kk,Nx+7-kk-1) = (Damp(kk+1)/(hp*hp)*(c(kk)*c(kk)*Alpha/(wg*wg) + 0.5*c(kk)*Iwp*Beta*dz/wg));                      
      Amoins.Get(Nx+7-kk-1,Nx+7-kk) = (Damp(kk+1)/(hp*hp)*(c(kk)*c(kk)*Alpha/(wg*wg) + 0.5*c(kk)*Iwp*Beta*dz/wg));
      
    }
  
  int kk=4;
  Aplus.Get(kk,kk) = 0.5*((1.0/Damp(kk))+(1.0/Damp(kk+1))) - (Damp(kk)+Damp(kk+1))/(hp*hp)*(c(kk)*c(kk)*Alpha/(wg*wg) - 0.5*c(kk)*Iwp*Beta*dz/wg);
  Aplus.Get(Nx+7-kk,Nx+7-kk) = 0.5*((1.0/Damp(kk))+(1.0/Damp(kk+1))) - (Damp(kk)+Damp(kk+1))/(hp*hp)*(c(kk)*c(kk)*Alpha/(wg*wg) - 0.5*c(kk)*Iwp*Beta*dz/wg);
  Amoins.Get(kk,kk) = 0.5*((1.0/Damp(kk))+(1.0/Damp(kk+1))) - (Damp(kk)+Damp(kk+1))/(hp*hp)*(c(kk)*c(kk)*Alpha/(wg*wg) + 0.5*c(kk)*Iwp*Beta*dz/wg);                  
  Amoins.Get(Nx+7-kk,Nx+7-kk) = 0.5*((1.0/Damp(kk))+(1.0/Damp(kk+1))) - (Damp(kk)+Damp(kk+1))/(hp*hp)*(c(kk)*c(kk)*Alpha/(wg*wg) + 0.5*c(kk)*Iwp*Beta*dz/wg);
              
} 

template<class TypeEquation>
void RunAll(EllipticProblem<TypeEquation>& var, const string& input_file, const string& name_element)
{
  // quelques typedef
  typedef typename TypeEquation::Dimension Dimension;
  typedef typename Dimension::R_N R_N;
  typedef typename TypeEquation::Complexe Complexe;
  GlobalGenericMatrix<Complexe> nat_mat;

  // initialisation proprietes physiques (rho, lambda, mu, etc)
  var.InitIndices(20);
  
  // on lit le fichier de donnees
  ReadInputFile(input_file, var);

  All_LinearSolver glob_solver(var);
  ReadInputFile(input_file, glob_solver); 

  // lecture des donnees pour le paraxial
  MyData my_data; 
  ReadInputFile(input_file, my_data);


    // liste des noms de fichiers de maillage dans mesh_data
  Vector<Vector<string> > mesh_data = var.mesh_data;
  var.mesh_data.Clear();

  Vector<Vector<R_N> > PointsDofBefore(my_data.ncouches), PointsDofAfter(my_data.ncouches);
  Vector<IVect > NumDofBefore(my_data.ncouches), NumDofAfter(my_data.ncouches);

  Vector<EllipticProblem<TypeEquation> > varT(my_data.ncouches);
  Vector<EllipticProblem<HelmholtzEquation1D> > varT1D(my_data.ncouches);
  
  int Ncouches=my_data.ncouches;
  BasicTimer cpt_tps;
  int ref_cpt = cpt_tps.GetNumber();
  Real_wp cumultps=0.0;
  Vector<Vector<Real_wp> > tps_resol(Ncouches);
  
  for(int ii=0; ii<Ncouches;ii++)
    tps_resol(ii).Reallocate(my_data.nit);

  // Donnees pour les DF paraxial
  VectReal_wp Xp(my_data.nx+8);
  VectReal_wp Xp2(my_data.nx);
  int Nxs=200;
  VectReal_wp Xs(Nxs);
  VectComplex_wp Us(Nxs);
  Real_wp hp;
  //int Np=8; //Nbre de points ajoutes en tout dans les pml (4 de chaque cote)
  VectReal_wp Sigma;
  VectComplex_wp Damp;
  Matrix<Complex_wp, General, ArrayRowSparse> Aplus;
  Matrix<Complex_wp, General, ArrayRowSparse> Amoins;
  
  
  // Definitions des problemes
  cpt_tps.Start(ref_cpt);
  for(int ii=0;ii<my_data.ncouches;ii++)
    {     
      varT(ii) = var;
      varT(ii).mesh_data.Clear();
      varT(ii).mesh_data.PushBack(mesh_data(ii));

      //if(!my_data.paraxial(ii)) // Si domaine(ii) Helmholtz
          //{
      varT(ii).mesh.SetAdditionPML(var.mesh.PML_BOTH_SIDES, var.mesh.PML_NEGATIVE_SIDE, -1, 1);
      varT(ii).mesh.SetThicknessPML(1.0);
      //(3.0*0.82*my_data.zmax)/80.0; //1.0;
      
      varT(ii).grid_to_be_computed = true;
      varT(ii).ComputeMeshAndFiniteElement(name_element);
      //varT(ii).mesh.Write("maillage.mesh");

      varT(ii).PerformOtherInitializations();
          
      GetDofFixedCoordinate(varT(ii), PointsDofAfter(ii), NumDofAfter(ii),varT(ii).mesh.GetYmax(), Dimension::dim_N-1);
      
      GetDofFixedCoordinate(varT(ii), PointsDofBefore(ii), NumDofBefore(ii),varT(ii).mesh.GetYmin()+varT(ii).mesh.GetThicknessPML(), Dimension::dim_N-1);
           
      //      DISP(varT(ii).GetNbEltPML());
      
          //  }
      if(my_data.paraxial(ii))           
        { // Si domaine(ii) Paraxial
          // indices physiques
          Xp.Reallocate(my_data.nx+8);
          hp=(my_data.xmax - my_data.xmin)/(my_data.nx-1);
          cout << "hp = " << hp << endl;
          
          for(int kk=0;kk<my_data.nx;kk++)
            {
              int pp=4+kk;
              Xp(pp)=my_data.xmin+kk*hp;
              Xp2(kk)=my_data.xmin+kk*hp;
            }
          
          for(int kk=0;kk<4;kk++)
            {
              int pp=3-kk;
              Xp(pp)=Xp(4)-(kk+1)*hp;
              Xp(my_data.nx+4+kk)=Xp(my_data.nx+4-1)+(kk+1)*hp;
            }
          
          Real_wp hs=hp=(my_data.xmax-my_data.xmin)/(Nxs-1);
          
          for(int kk=0; kk<Nxs;kk++)
            Xs(kk)=my_data.xmin+kk*hs;
          

          cout << " Xs = " << Xs(0) << "     " << Xs(Nxs-1) << endl;
          cout << "xmin = " << my_data.xmin << endl;
          cout << " Xp = " << Xp(0) << "     " << Xp(my_data.nx+7) << endl;

          
          Sigma.Reallocate(5);
          Sigma(4)=0.167;
          Sigma(3)=0.520;
          Sigma(2)=1.273;
          Sigma(1)=2.952;
          Sigma(0)=9.250;
          
          Damp.Reallocate(6);
          Damp.Fill(1.0);
          
          varT1D(ii).InitIndices(10);          
          varT1D(ii).ref_rho.Fill(1.0); varT1D(ii).ref_mu.Fill(1.0);         // !!!!!!!!!!!!!!!!!   A MODIFIER  !!!!!!!!!!!!!!!!!! 
          varT1D(ii).omega2 = -1.0;
          varT1D(ii).mesh.SetOrder(varT(ii).mesh_num.GetOrder());
          
          // maillage 1-D regulier
          varT1D(ii).mesh_data.Reallocate(1);
          varT1D(ii).mesh_data(0).Reallocate(4);
          varT1D(ii).mesh_data(0)(0) = string("REGULAR");
          varT1D(ii).mesh_data(0)(1) = varT(ii).mesh_data(0)(1); //to_str(my_data.nx); //varT(ii).mesh_data(0)(1); //to_str(my_data.nx);
          varT1D(ii).mesh_data(0)(2) = varT(ii).mesh_data(0)(2); //to_str(my_data.xmin);          
          varT1D(ii).mesh_data(0)(3) = varT(ii).mesh_data(0)(3); //to_str(my_data.xmax);
          

          varT1D(ii).ComputeMeshAndFiniteElement(string("EDGE_LOBATTO"));   
          
          // calcul de coefficients geometriques 
          varT1D(ii).ComputeMassMatrix();

        }
    }
  cpt_tps.Stop(ref_cpt);
  cout << "Temps boucle Initialisation = " << cpt_tps.GetSeconds(ref_cpt) << endl;
  cpt_tps.Stop(ref_cpt);
  
  
  // classe pour la source
  UserSourceDirichlet<Dimension> fct_source(varT(0));
  // fct_source.src_file.Init("toto.dat", -6, 6);
  fct_source.source_on_file = false;
  Vector<Complexe> source_rhs, x_sol, v, dv_dy; VectReal_wp weight;

  filtrage=my_data.filtrage; // sert pour la source
  
  int Nf=my_data.nf;

  int Nmat=my_data.nmat;   
  

  //  Vector<int> TypeMat(Ncouches);
  //for(int ii=0;ii<Ncouches;ii++)
  //  TypeMat(ii)=to_num<int>(varT(ii).mesh_data(0)(6));


  //for(int ii=0;ii<Ncouches;ii++)
  //  {
  //    cout << "Couche " << ii << endl;
  //    cout << " ref_rho = " << endl;
  //    varT(ii).ref_rho(0).Print();
  //    cout << endl << endl << " ref_mu = " << endl;
  //    varT(ii).ref_mu.Print();
  //    getchar();
  //  }

  // Lecture des parametres physiques
  Vector<Real_wp> rho(Ncouches);
  Vector<Real_wp> mu(Ncouches);
  Vector<Real_wp> c(Ncouches);
  
  for(int ii=0;ii<Nmat;ii++)
    {
      rho(ii)=abs(Complex_wp(var.ref_rho(ii+1)));
      mu(ii)=abs(Complex_wp(var.ref_mu(ii+1)(0,0)));
      c(ii)=sqrt(mu(ii)/rho(ii));        
    }
  c(0)=1.0;
  cg=c(0);
  if(my_data.ncouches > 1)
    {
      for(int ll=1; ll<my_data.ncouches; ll++)
        c(ll)=cg;
    }
  
  cout << "Domaines decompose en " << Ncouches << " couches." << endl;
  //  cout << "rho = " << rho(0) << " mu = " << mu(0) << " c = " << cg << endl;

  Vector<bool> paraxial(Ncouches);
  for(int ii=0;ii<Ncouches;ii++)
    paraxial(ii)=my_data.paraxial(ii);
  

  // Definition de la source
  Real_wp ts=0.82*my_data.zmax/cg;
  // !!!!!!!!!!!!!! RETOUCHER TS EN MULTICOUCHES !!!!!!!!!!!!!!!!
  // for(int ii=0;ii<Nmat;ii++)
  //  ts+=(Lz/Real_wp(Nmat))/c(ii);
  //  ts*=0.82;
  Real_wp fs=10/ts;
  Real_wp df=(8.0/3.0)*fs/((Real_wp)Nf);
  Real_wp ws=2*pi_wp*fs;
  Real_wp dw=2*pi_wp*df;

  
  Vector<Vector<Real_wp> > Image(Ncouches);
  VectComplex_wp V0;
  VectComplex_wp V0t, DuDz;

  cout << "Frequence : " << endl << Nf << " Frequences de " << dw << " a " << Nf*dw << endl << endl;
  cout << "ts = " << ts << " et fs = " << fs  << endl << "df = " << df << "et  ws = " << ws << endl << "dw = " << dw << endl;
  cout << "Nbre de domaines = " << mesh_data.GetM() << endl;
  cout << "           Appuyer sur une touche. " << endl;
  //getchar();

  // Pour le paraxial
  Real_wp xmin=my_data.xmin;
  Real_wp xmax=my_data.xmax;
  int Nz=my_data.nz;

  Real_wp dz;
  Real_wp hcouches=my_data.zmax/Real_wp(my_data.ncouches);
  
  // Boucle sur le frequence   
  VectComplex_wp gradv(my_data.nx);
  VectComplex_wp tracev(my_data.nx);
  
  for(int ww=1;ww<my_data.nit+1;ww++)  
    {      
      wg=2*pi_wp*df*ww;
      if(my_data.nit==1)
        wg=2*pi_wp*df*my_data.nf;
      var.SetOmega(wg);
      //wg=var.GetOmega();

      // Affichage d'informations
      cout << endl << endl << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl << "               iteration = " << ww << endl;
      cout << "               omega = " <<  var.GetOmega() << endl << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"  << endl << endl;

      // boucle sur tous les domaines a traiter
      for(int i = 0; i < mesh_data.GetM(); i++)
        {          

          varT(i).SetOmega(wg);
          Real_wp Lz=i*hcouches;


          cout << "Nbre de Elt = " << varT(i).mesh.GetNbElt() << endl;
          cout << "Nbre de Dof = " << varT(i).GetNbDof() << endl;
          cout << "Nbre Elts dans PML = " << varT(i).GetNbEltPML() << endl;
          cout << "NumDofBefore(ii).Size = " <<  NumDofBefore(i).GetM() << endl;
          cout << "MeshStep = " << varT(i).mesh.GetMeshSize() << endl;

          if(paraxial(i)) 
            {
              cout << "Domaine " << i << " : paraxial." << endl;              
                          

              dz=fabs(hcouches/Real_wp(Nz-1)); // pas en z             
              
              // Coefficient d'amortissement
              for(int kk=0;kk<5;kk++)
                Damp(kk)=1.0/(1.0 + Iwp*Sigma(kk)*(c(i)/(wg*hp)));
                        
              Damp(5)=1.0;
             
              // Les 2 matrices utilisees
              Aplus.Reallocate(my_data.nx+8,my_data.nx+8);
              Aplus.Fill(0.0);
              Amoins.Reallocate(my_data.nx+8,my_data.nx+8);
              Amoins.Fill(0.0);

              //  Pade 45 complexe
              Complex_wp Alpha=0.2-Iwp*0.1;
              Complex_wp Beta=0.5;

              // Parax 15
              //Complex_wp Alpha=0.0;
              //Complex_wp Beta=0.5;

              // Pade 45 normal
              //Complex_wp Alpha=0.25;
              //Complex_wp Beta=0.5;
                                          
              for(int kk=0; kk<my_data.nx+8; kk++)
                {
                  Aplus.Get(kk,kk) = 1.0-(2.0/(hp*hp))*(c(i)*c(i)/(wg*wg)*Alpha - 0.5*c(i)*Iwp*Beta*dz/wg); //1-Iwp*2*c*dz/(4*wg*hp*hp);
                  Amoins.Get(kk,kk) = 1.0-(2.0/(hp*hp))*(c(i)*c(i)/(wg*wg)*Alpha + 0.5*c(i)*Iwp*Beta*dz/wg); //1-Iwp*2*c*dz/(4*wg*hp*hp);
                  
                  if(kk>0)
                    {
                      Aplus.Get(kk,kk-1) = 1.0/(hp*hp)*(c(i)*c(i)/(wg*wg)*Alpha - 0.5*c(i)*Iwp*Beta*dz/wg);  //Iwp*c*dz/(4*wg*hp*hp);
                      Amoins.Get(kk,kk-1) = 1.0/(hp*hp)*(c(i)*c(i)/(wg*wg)*Alpha + 0.5*c(i)*Iwp*Beta*dz/wg);  //Iwp*c*dz/(4*wg*hp*hp);                      
                    }

                  if(kk<my_data.nx+7)
                    {
                      Aplus.Get(kk,kk+1) = 1.0/(hp*hp)*(c(i)*c(i)*Alpha/(wg*wg) - 0.5*c(i)*Iwp*Beta*dz/wg); //Iwp*c*dz/(4*wg*hp*hp);
                      Amoins.Get(kk,kk+1) = 1.0/(hp*hp)*(c(i)*c(i)*Alpha/(wg*wg) + 0.5*c(i)*Iwp*Beta*dz/wg); //Iwp*c*dz/(4*wg*hp*hp);                      
                    }
                  
                }
 
              // Prise en compte des PML
              for(int kk=0;kk<4;kk++)
                {
                  Aplus.Get(kk,kk) = 0.5*((1.0/Damp(kk))+(1.0/Damp(kk+1))) - (Damp(kk)+Damp(kk+1))/(hp*hp)*(c(i)*c(i)*Alpha/(wg*wg) - 0.5*c(i)*Iwp*Beta*dz/wg);                  
                  Aplus.Get(kk,kk+1) = (Damp(kk+1)/(hp*hp)*(c(i)*c(i)*Alpha/(wg*wg) - 0.5*c(i)*Iwp*Beta*dz/wg));
                  Aplus.Get(kk+1,kk) = (Damp(kk+1)/(hp*hp)*(c(i)*c(i)*Alpha/(wg*wg) - 0.5*c(i)*Iwp*Beta*dz/wg));
                  
                  Amoins.Get(kk,kk) = 0.5*((1.0/Damp(kk))+(1.0/Damp(kk+1))) - (Damp(kk)+Damp(kk+1))/(hp*hp)*(c(i)*c(i)*Alpha/(wg*wg) + 0.5*c(i)*Iwp*Beta*dz/wg);
                  Amoins.Get(kk,kk+1) = (Damp(kk+1)/(hp*hp)*(c(i)*c(i)*Alpha/(wg*wg) + 0.5*c(i)*Iwp*Beta*dz/wg));                      
                  Amoins.Get(kk+1,kk) = (Damp(kk+1)/(hp*hp)*(c(i)*c(i)*Alpha/(wg*wg) + 0.5*c(i)*Iwp*Beta*dz/wg));
                  
                  Aplus.Get(my_data.nx+7-kk,my_data.nx+7-kk) = 0.5*((1.0/Damp(kk))+(1.0/Damp(kk+1))) - (Damp(kk)+Damp(kk+1))/(hp*hp)*(c(i)*c(i)*Alpha/(wg*wg) - 0.5*c(i)*Iwp*Beta*dz/wg);                  
                  Aplus.Get(my_data.nx+7-kk,my_data.nx+7-kk-1) = (Damp(kk+1)/(hp*hp)*(c(i)*c(i)*Alpha/(wg*wg) - 0.5*c(i)*Iwp*Beta*dz/wg));
                  Aplus.Get(my_data.nx+7-kk-1,my_data.nx+7-kk) = (Damp(kk+1)/(hp*hp)*(c(i)*c(i)*Alpha/(wg*wg) - 0.5*c(i)*Iwp*Beta*dz/wg));
                  
                  Amoins.Get(my_data.nx+7-kk,my_data.nx+7-kk) = 0.5*((1.0/Damp(kk))+(1.0/Damp(kk+1))) - (Damp(kk)+Damp(kk+1))/(hp*hp)*(c(i)*c(i)*Alpha/(wg*wg) + 0.5*c(i)*Iwp*Beta*dz/wg);
                  Amoins.Get(my_data.nx+7-kk,my_data.nx+7-kk-1) = (Damp(kk+1)/(hp*hp)*(c(i)*c(i)*Alpha/(wg*wg) + 0.5*c(i)*Iwp*Beta*dz/wg));                      
                  Amoins.Get(my_data.nx+7-kk-1,my_data.nx+7-kk) = (Damp(kk+1)/(hp*hp)*(c(i)*c(i)*Alpha/(wg*wg) + 0.5*c(i)*Iwp*Beta*dz/wg));
                                                         
                }
                  
              int kk=4;
              Aplus.Get(kk,kk) = 0.5*((1.0/Damp(kk))+(1.0/Damp(kk+1))) - (Damp(kk)+Damp(kk+1))/(hp*hp)*(c(i)*c(i)*Alpha/(wg*wg) - 0.5*c(i)*Iwp*Beta*dz/wg);
              Aplus.Get(my_data.nx+7-kk,my_data.nx+7-kk) = 0.5*((1.0/Damp(kk))+(1.0/Damp(kk+1))) - (Damp(kk)+Damp(kk+1))/(hp*hp)*(c(i)*c(i)*Alpha/(wg*wg) - 0.5*c(i)*Iwp*Beta*dz/wg);
              Amoins.Get(kk,kk) = 0.5*((1.0/Damp(kk))+(1.0/Damp(kk+1))) - (Damp(kk)+Damp(kk+1))/(hp*hp)*(c(i)*c(i)*Alpha/(wg*wg) + 0.5*c(i)*Iwp*Beta*dz/wg);                  
              Amoins.Get(my_data.nx+7-kk,my_data.nx+7-kk) = 0.5*((1.0/Damp(kk))+(1.0/Damp(kk+1))) - (Damp(kk)+Damp(kk+1))/(hp*hp)*(c(i)*c(i)*Alpha/(wg*wg) + 0.5*c(i)*Iwp*Beta*dz/wg);
              
              //int Nx = varT1D(i).GetNbDof();
              // objet contenant alpha et beta
              //GlobalH1_Matrix<Complex_wp> nat_mat1D;

              //nat_mat1D.SetCoefMass(1.0);

              //  Pade complexe
              //Complex_wp Alpha=0.2-Iwp*0.1;
              //Complex_wp Beta=0.5;

              // Pade normal
              //Complex_wp Alpha=0.25;
              //Complex_wp Beta=0.5;

              //if(Parax45==true)
              //  nat_mat1D.SetCoefStiffness(-Alpha*c(i)*c(i)/(wg*wg) + 0.5*Beta*Iwp*c(i)*dz/wg);
              //else
              //  nat_mat1D.SetCoefStiffness(0.25*Iwp*c(i)*dz/wg);
                              
              //Matrix<Complex_wp, Symmetric, ArrayRowSymSparse> Dh_plus_betaK, Dh_minus_betaK;                          

              // varT1D(i).AddMatrixFEM(Dh_plus_betaK, nat_mat1D);

              
              //if(Parax45==true)
              //  nat_mat1D.SetCoefStiffness(-Alpha*c(i)*c(i)/(wg*wg) - 0.5*Beta*Iwp*c(i)*dz/wg);
              //else
              //  nat_mat1D.SetCoefStiffness(-0.25*Iwp*c(i)*dz/wg);
              
              //varT1D(i).AddMatrixFEM(Dh_minus_betaK, nat_mat1D);
                            

  
              MatrixMumps<Complex_wp> mat_lu;
              mat_lu.HideMessages();
              Seldon::GetLU(Aplus, mat_lu);

              
              // calcul de la grille d'interpolation 1-D -> image
              // IVect ElementInterp(my_data.nx);
              //VectReal_wp CoorInterp(my_data.nx);
              
              // Donnees servant a recuperer la source
              // VectReal_wp PosInterp(my_data.nx);
              // VectReal_wp s(2);
              //Real_wp step_x = (xmax - xmin)/(my_data.nx-1);
              //for (int ii = 0; ii < varT1D(i).mesh.GetNbElt(); ii++)
              // {
              //   varT1D(i).mesh.GetVerticesElement(ii, s);
              //   int ixmin = std::max(int(ceil((s(0)-xmin)/step_x)), 0);
              //    int ixmax = std::min(int(floor((s(1)-xmin)/step_x)), my_data.nx-1);
                  
              //    for (int ix = ixmin; ix <= ixmax; ix++)
              //      {
                      //PosInterp(ix)= xmin + ix*step_x;
              //        Real_wp pos_x = xmin + ix*step_x;
              //        ElementInterp(ix) = ii;
              //        CoorInterp(ix) = (pos_x - s(0)) / ( s(1) - s(0) );
              //      }
              //  }
              
              // Source
              VectComplex_wp V1(my_data.nx+8), Voutput(Nxs*my_data.nz);
              VectComplex_wp VoutReduit(my_data.nx);
              VectComplex_wp Vout(my_data.nx+8);
              if(i==0) //source
                {
                  
                  //SourceDirichlet1D fct_source1D;
                  V0.Reallocate(my_data.nx+8);
                  V0.Fill(0); V1.Fill(0);
                  //varT1D(i).AddVolumicProjection(Complex_wp(1), V0, fct_source1D);
                  for(int ll=4;ll<my_data.nx+4;ll++)
                    {                                            
                      if(Xp(ll)!=0)
                        V0(ll)=sin(wg*Xp(ll)/c(i))/Xp(ll);
                      else
                        V0(ll)=wg/c(i);
                    }
                  Complex_wp coef=Iwp*(wg/ws)*exp(-wg*wg/(ws*ws))*exp(-Iwp*wg*ts); // partie frequentielle
                  for(int ii=0;ii<V0.GetSize();ii++)
                    V0(ii)*=coef;
                }
              else
                {                                                                                                                                 
                  if(!my_data.paraxial(i-1)) // domain i-1 : Helmholtz                                                                                 
                    { 
                      cout << "HELMHOLTZ AVANT !!!!!!!!!!!!!" << endl;
                      V0.Reallocate(my_data.nx+8);
                      V0.Fill(0);
                      
                      VectReal_wp Xh(NumDofBefore(i).GetM()-2*varT(i).mesh_num.GetOrder());
                      VectComplex_wp Uh(NumDofBefore(i).GetM()-2*varT(i).mesh_num.GetOrder());

                      for(int kk=0;kk<NumDofBefore(i).GetM()-2*varT(i).mesh_num.GetOrder();kk++)
                        {
                          Xh(kk)=PointsDofBefore(i-1)(kk+varT(i).mesh_num.GetOrder())(0);
                          Uh(kk)=v(kk+varT(i).mesh_num.GetOrder());
                        }

                      
                      SplineInterp(Xh,Uh,Xp2,VoutReduit);


                      for(int kk=0; kk<my_data.nx; kk++)
                        V0(kk+4)=VoutReduit(kk)*exp(-Iwp*wg*Lz/c(i));
                    }
                  else // domaine i-1 : Parax
                    {
                      // V0=V0t;
                      // rien a faire V0 = ok
                    }
                }
            
              
              // solution a la surface du domaine
              // GetInterpolate1D(varT1D(i), V0, V0_interp, ElementInterp, CoorInterp); 
              
              Vout=V0;
              if(i!=0)
                Mlt(exp(Iwp*wg*Lz/c(i)), Vout);

              VoutReduit.Fill(0);
              for(int kk=0; kk<my_data.nx;kk++)
                VoutReduit(kk)=Vout(kk+4);

              SplineInterp(Xp2,VoutReduit,Xs,Us);
              
              for (int ii = 0; ii < Nxs; ii++)
                Voutput(ii) = Us(ii);                                
                //Voutput(ii) = V0_interp(ii);                                
                
              for (int k = 1; k < Nz; k++)
                {
                  
                  V0t=V0;
                  //Mlt(Dh_minus_betaK, V0, V1);    // V1 =  Dh_minus_betaK * V0             
                  Mlt(Amoins, V0, V1);    // V1 =  Amoins * V0                               
                  SolveLU(mat_lu, V1); // V1 = inv(Aplus) * V1
                  
                  //GetInterpolate1D(varT1D(i), V1, V0_interp, ElementInterp, CoorInterp); 
                  
                  
                  Vout=V1;
                  Complexe coefp = exp(Iwp*wg*(Real_wp(k)*dz + Lz)/c(i));
                  Mlt(coefp, Vout);

                  for(int kk=0; kk<my_data.nx;kk++)
                    VoutReduit(kk)=Vout(kk+4);

                  SplineInterp(Xp2,VoutReduit,Xs,Us);
                  
                  
                  for (int ii = 0; ii < Nxs; ii++)
                    Voutput(Nxs*k + ii) = Us(ii); //*coef;
                  
                  V0 = V1;
                  //v.Reallocate(V1.GetM());
                  //v = V1;
                }  
              
              if(i+1<my_data.ncouches)
                {
                  if(!paraxial(i+1))
                    {
                      Complexe coefp = exp(Iwp*wg*Real_wp((Nz-1)*dz + Lz)/c(i));
                      
                      // V0t = trace de u avant le bord
                      // V1 = trace de u sur le bord
                      
                      DuDz=V1;
                      Add(Complex_wp(-1, 0), V0t, DuDz);
                      Mlt(1/dz,DuDz); // DvDz
                      Mlt(coefp,DuDz);


                      V0t=V1;
                      Mlt(coefp,V0t); // V0t = u

                      Add(-Iwp*wg/cg,V0t,DuDz);  // DuDz 
                      
                      for(int kk=0;kk<my_data.nx; kk++)
                        {
                          tracev(kk)=V0t(kk+4);
                          gradv(kk)=DuDz(kk+4);
                        }
                    }
                }
              
              if(ww==1)
                {        
                  Image(i).Reallocate(Nxs*my_data.nz);
                  for(int kk=0; kk<Nxs*my_data.nz; kk++)
                    Image(i)(kk)=Voutput(kk).real();
                }
              else
                {                  
                  for(int ii=0;ii<Voutput.GetSize();ii++)
                    Image(i)(ii)+=Voutput(ii).real();
                }
              
              if(ww==10)
                {
                  std::string fichiersortie="Migration10freqP_";
                  std::string final="dom";
                  std::ostringstream oss1;
                  std::ostringstream oss2;
                  oss1 << i;
                  oss2 << Ncouches;
                  fichiersortie+=oss2.str();
                  fichiersortie+=final;
                  fichiersortie+=oss1.str();
                  fichiersortie+=".dat";
                  //WriteGrid2D(xmin, xmax, 0.0, my_data.zmax,my_data.nx,my_data.nz, Image(i), fichiersortie);
                  WriteGrid2D(xmin, xmax,to_num<Real_wp>(varT(i).mesh_data(0)(5)), to_num<Real_wp>(varT(i).mesh_data(0)(4)),Nxs,my_data.nz, Image(i), fichiersortie);
                }

              if(ww==40)
                {
                  std::string fichiersortie="Migration40freqP_";
                  std::string final="dom";
                  std::ostringstream oss1;
                  std::ostringstream oss2;
                  oss1 << i;
                  oss2 << Ncouches;
                  fichiersortie+=oss2.str();
                  fichiersortie+=final;
                  fichiersortie+=oss1.str();
                  fichiersortie+=".dat";
                  //varT(i).WriteOutputFile(Image(i),fichiersortie);
                  //WriteGrid2D(xmin, xmax, 0.0, my_data.zmax,my_data.nx,my_data.nz, Image(i), fichiersortie);
                  WriteGrid2D(xmin, xmax, to_num<Real_wp>(varT(i).mesh_data(0)(5)), to_num<Real_wp>(varT(i).mesh_data(0)(4)),Nxs,my_data.nz, Image(i), fichiersortie);
                  //WriteGrid2D(xmin, xmax, -10.0, 0.0, my_data.nx, my_data.nz, Image(i), fichiersortie);
                  //WriteGrid2D(xmin, xmax, 0, my_data.zmax, my_data.nx, Nz, Image(i), "MigrationParax40.dat");
                }

              if(ww==my_data.nit)
                {
                  //Image(i).WriteText("toto.dat");
                  std::string fichiersortie="MigrationP_";
                  std::string final="dom";
                  std::ostringstream oss1;
                  std::ostringstream oss2;
                  oss1 << i;
                  oss2 << Ncouches;
                  fichiersortie+=oss2.str();
                  fichiersortie+=final;
                  fichiersortie+=oss1.str();
                  fichiersortie+=".dat";
                  //varT(i).WriteOutputFile(Image(i),fichiersortie);
                  //WriteGrid2D(xmin, xmax, 0, my_data.zmax, my_data.nx, Nz, Image(i), "MigrationParax.dat");
                  //WriteGrid2D(xmin, xmax, -10.0, 0.0, my_data.nx, my_data.nz, Image(i), fichiersortie);
                  //WriteGrid2D(xmin, xmax, 0.0, my_data.zmax,my_data.nx,my_data.nz, Image(i), fichiersortie);
                  //WriteGrid2D(xmin, xmax, to_num<Real_wp>(varT(i).mesh_data(0)(4)),to_num<Real_wp>(varT(i).mesh_data(0)(5)), Nxs,my_data.nz, Image(i), fichiersortie);
                  WriteGrid2D(xmin, xmax,to_num<Real_wp>(varT(i).mesh_data(0)(5)),to_num<Real_wp>(varT(i).mesh_data(0)(4)),Nxs,my_data.nz, Image(i), fichiersortie);
                }

        }
      else  // ######################  HELMHOLTZ #############################
            { 
              cout << "Domaine " << i << " : helmholtz." << endl;
              cpt_tps.Start(ref_cpt);
              // calcul des donnees geometriques (matrices jacobiennes DF_i)          
              varT(i).ComputeMassMatrix();
              All_LinearSolver solver_tmp(varT(i));
                       
              cpt_tps.Stop(ref_cpt);
              cout << "Domaine " << i << " ww = " << ww << " Temps ComputeMass = " << cpt_tps.GetSeconds(ref_cpt) << endl;
              cpt_tps.Reset(ref_cpt);
              // factorisation du systeme lineaire a resoudre (si resolution directe)
              // calcul de la matrice et du preconditionneur (si resolution iterative)
              
              //glob_solverT(i).PerformFactorizationStep(nat_mat);

              // on calcule la source
              int nb_ddl = varT(i).GetNbDof();
              source_rhs.Reallocate(nb_ddl);
              source_rhs.Zero();
              if (i == 0)
                {
		  Vector<Vector<Complexe> > rhs_vec;
		  rhs_vec.SetData(1, &source_rhs);
		  
                  // premier domaine, on calcule la source
                  varT(i).ComputeGenericSource(rhs_vec, fct_source);

                  Complexe coef=Iwp*(wg/ws)*exp(-wg*wg/(ws*ws))*exp(-Iwp*wg*ts);
                  Mlt(coef, source_rhs);

                  // Multiplication par le terme freqentiel
                  //                  for(int ii=0;ii<source_rhs.GetSize();ii++)
                  //  source_rhs(ii)*=Iwp*(wg/ws)*exp(-wg*wg/(ws*ws))*exp(-Iwp*wg*ts);
		  rhs_vec.Nullify();
                }
              else
                {
                  if(my_data.paraxial(i-1))
                    {
                      VectReal_wp weight2;
                      cout << "                 INTERFACE P - H " << endl;
                      VectComplex_wp tracevh;
                      VectComplex_wp gradvh;
                                            
                      //   VectReal_wp Xs(PointsDofAfter(i).GetM()-2*int(varT(i).mesh.thickness_PML)*varT(i).mesh.GetOrder());
                      //for(int kk=0; kk< PointsDofAfter(i).GetM()-2*varT(i).mesh.thickness_PML*varT(i).mesh.GetOrder(); kk++)
                      //  Xs(kk) = PointsDofAfter(i)(kk+varT(i).mesh.GetOrder())(0);
                      
                      //source_rhs.Zero();
                      GetWeight(varT(i), NumDofAfter(i), varT(i).mesh.GetYmax(), Dimension::dim_N-1, weight2);

                      VectReal_wp Xh(NumDofAfter(i).GetM()-2*varT(i).mesh_num.GetOrder()); // Maillage helmholtz

                      for(int kk=0;kk<NumDofAfter(i).GetM()-2*varT(i).mesh_num.GetOrder();kk++)
                        Xh(kk)=PointsDofAfter(i)(kk+varT(i).mesh_num.GetOrder())(0);
                      
                      // interpolation de la trace et du gradient Paraxial vers Helmholtz
                      SplineInterp(Xp2,tracev,Xh,tracevh);
                      SplineInterp(Xp2,gradv,Xh,gradvh);
                      
                      if(condinterface==0) // Transmission par dirichlet
                        { 
                          for(int kk=0; kk<NumDofAfter(i).GetM()-2*varT(i).mesh_num.GetOrder();kk++)
                            source_rhs(NumDofAfter(i)(kk+varT(i).mesh_num.GetOrder()))=tracevh(kk)*weight2(kk+varT(i).mesh_num.GetOrder());
  
                        }
                      else 
                        {
                          for(int kk=0; kk<NumDofAfter(i).GetM()-2*varT(i).mesh_num.GetOrder();kk++)
                            source_rhs(NumDofAfter(i)(kk+varT(i).mesh_num.GetOrder()))=(gradvh(kk) - Iwp*varT(i).GetOmega()*tracevh(kk))*weight2(kk+varT(i).mesh_num.GetOrder());
                        }
                      
                      
                      //for(int kk=0; kk<V0.GetM();kk++)
                      //  source_rhs(NumDofAfter(i)(kk+varT(i).mesh_num.GetOrder()))=V0(kk); //(DuDz(kk) - Iwp*varT(i).GetOmega()*V0(kk))*weight(kk+varT(i).mesh_num.GetOrder());
                    }
                  else
                    {
                      
                      GetWeight(varT(i), NumDofAfter(i), varT(i).mesh.GetYmax(), Dimension::dim_N-1, weight);

                      // on injecte la solution du domaine d'avant
                      if(condinterface==0)
                        {         
                          cout << "Dirichlet a l'interface" << endl;
                          for (int j = 0; j < NumDofAfter(i).GetM(); j++)
                            source_rhs(NumDofAfter(i)(j)) = v(j)*weight(j);
                        }                                            
                      else
                        {
                          cout << "Impedance a l'interface" << endl;                          
                          //GetWeight(varT(i), NumDofAfter(i), varT(i).mesh.GetYmin(), Dimension::dim_N-1, weight);
                          for (int j = 0; j < PointsDofBefore(i).GetM(); j++)
                            source_rhs(NumDofAfter(i)(j)) = (dv_dy(j) - Iwp*varT(i).GetOmega()*v(j))*weight(j);
                        }
                    }
                }
              
              // Multiplication par le terme freqentiel
              //Real_wp tomega=var.GetOmega();
              //for(int ii=0;ii<source_rhs.GetSize();ii++)
              // source_rhs(ii)*=Iwp*tomega/ws*exp(-tomega*tomega/(ws*ws))*exp(Iwp*tomega*ts);

              
              // on calcule la solution
              //glob_solverT(i).ComputeSolution(source_rhs, x_sol, nat_mat);
              //time_t tstart, tend;
              //time(&tstart);
              cpt_tps.Start(ref_cpt);
              solver_tmp.PerformFactorizationStep(nat_mat);
              solver_tmp.ComputeSolution(source_rhs, x_sol, nat_mat);
              cpt_tps.Stop(ref_cpt);
              Real_wp tps_solve=cpt_tps.GetSeconds(ref_cpt);
              tps_resol(i)(ww-1)=tps_solve;
              cumultps+=tps_solve;
              //cout << "Domaine " << i << " : Temps factorisation+resolution = " << tps_solve << " et cumul = " << cumultps << endl;
              cpt_tps.Reset(ref_cpt);
              if(ww==1) 
                {
                  Image(i).Reallocate(x_sol.GetM());
                  for(int kk=0; kk<x_sol.GetM(); kk++)
                    Image(i)(kk)=x_sol(kk).real();          
                }
              else
                {                            
                  for(int ii=0;ii<x_sol.GetM();ii++)
                    Image(i)(ii)+=x_sol(ii).real();
                }
              
              
              //      int test_input; cout<<"Appuyez sur une touche pour passer au domaine suivant "<<endl; cin>>test_input;
              //GetDofFixedCoordinate(var, PointsDofBefore, NumDofBefore,
              //                     var.mesh.GetYmin()+var.mesh.thickness_PML, Dimension::dim_N-1);
              
              // on calcule le gradient u sur les points nodaux
              //GetDuDy(varT(i), x_sol, NumDofBefore(i), v, dv_dy);
              
              GetDuDy(varT(i), x_sol, NumDofBefore(i), v, dv_dy);
              
              if(ww==10)
                {
                  std::string fichiersortie="Migration10freq_";
                  std::string final="dom";
                  std::ostringstream oss1;
                  std::ostringstream oss2;
                  oss1 << i;
                  oss2 << Ncouches;
                  fichiersortie+=oss2.str();
                  fichiersortie+=final;
                  fichiersortie+=oss1.str();
                  fichiersortie+=".dat";
                  VectComplex_wp Imagetmp(Image(i).GetM());                  
                  for(int kk=0; kk<Image(i).GetM(); kk++)
                    Imagetmp(kk)=Complexe(Image(i)(kk));

                  //varT(i).WriteOutputFile(VectComplex_wp(Image(i)),fichiersortie);  
                  varT(i).WriteOutputFile(Imagetmp,fichiersortie);  
                } 

              if(ww==40)
                {
                  std::string fichiersortie="Migration40freq_";
                  std::string final="dom";
                  std::ostringstream oss1;
                  std::ostringstream oss2;
                  oss1 << i;
                  oss2 << Ncouches;
                  fichiersortie+=oss2.str();
                  fichiersortie+=final;
                  fichiersortie+=oss1.str();
                  fichiersortie+=".dat";

                  VectComplex_wp Imagetmp(Image(i).GetM());                  
                  for(int kk=0; kk<Image(i).GetM(); kk++)
                    Imagetmp(kk)=Complexe(Image(i)(kk));
                  varT(i).WriteOutputFile(Imagetmp,fichiersortie);  
                } 

              if(ww==my_data.nit)
                {
                  std::string fichiersortie="Migration_";
                  std::string final="dom";
                  std::ostringstream oss1;
                  std::ostringstream oss2;
                  oss1 << i;
                  oss2 << Ncouches;
                  fichiersortie+=oss2.str();
                  fichiersortie+=final;
                  fichiersortie+=oss1.str();
                  fichiersortie+=".dat";

                  VectComplex_wp Imagetmp(Image(i).GetM());                  
                  for(int kk=0; kk<Image(i).GetM(); kk++)
                    Imagetmp(kk)=Complexe(Image(i)(kk));
                  varT(i).WriteOutputFile(Imagetmp,fichiersortie);                  
                  }
            } // fin else parax
        } // fin boucle i
    } // fin boucle ww
  
  VectReal_wp tps_moy(Ncouches);
  Real_wp cumul_moy=0.0;
  for(int ii=0; ii<Ncouches;ii++)
    {
      tps_moy(ii)=0.0;
      for(int kk=0;kk<my_data.nit;kk++)
        tps_moy(ii)+=tps_resol(ii)(kk);

      tps_moy(ii)*=1/Real_wp(my_data.nit);
      cumul_moy+=tps_moy(ii);
      cout << "Domaine " << ii << " temps moyen = " << tps_moy(ii) << " et cumul = " << cumul_moy << endl;
    }
}



int main(int argc, char **argv) 
{

  time_t tstart,tend; // mesure du temps de calcul
  time (&tstart);
  InitMontjoie(argc, argv);
    
  int Ne=200;
  int Ns=100;

  Real_wp he=1.0/Real_wp(Ne-1);
  Real_wp hs=1.0/Real_wp(Ns-1);
  VectReal_wp Xe(Ne);
  VectComplex_wp Ue(Ne);
  VectReal_wp Xs(Ns);
  VectComplex_wp Us(Ns);

  for(int ii=0;ii<Ne;ii++)
    {
      Xe(ii)=ii*he;
      Ue(ii)=Xe(ii)*Xe(ii)*sqrt(2*(1+tanh(Xe(ii))))-5*Xe(ii);
    }
  
  for(int ii=0; ii<Ns;ii++)
    Xs(ii)=ii*hs;

  cout << "Avant Interp" << endl;
  SplineInterp(Xe, Ue, Xs,Us);

  cout << "Apres Interp" << endl;

  VectReal_wp Uer(Ne);
  VectReal_wp Usr(Ns);

  for(int ii=0;ii<Ne;ii++)
    Uer(ii)=Ue(ii).real();

  for(int ii=0;ii<Ns;ii++)
    Usr(ii)=Us(ii).real();
  
  for(int ii=0;ii<Ne;ii++)
    Uer(ii)=Ue(ii).real();

  Usr.WriteText("Usr.dat");
  Uer.WriteText("Uer.dat");
  

  
  
  if (argc>1)
    {
      EllipticProblem<HelmholtzEquation<Dimension2> > var;
      RunAll(var, string(argv[1]), string("TRIANGLE_LOBATTO"));
    }
  else
    {
      cout<<"This code needs a data file in argument"<<endl;
      cout<<"helmholtz2D.x nom_fichier"<<endl;
      cout<<"is a good syntax"<<endl;
    }
  time(&tend);
  
  cout << "Temps d'execution = " << difftime(tend,tstart) << endl;
  cout<<"End of the program"<<endl; 

#ifdef SELDON_WITH_MPI  
  MPI_Finalize();
#endif
    
  return 0;
}

