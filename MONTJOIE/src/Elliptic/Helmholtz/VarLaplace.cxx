#ifndef MONTJOIE_FILE_VAR_LAPLACE_CXX

#include "VarLaplaceInline.cxx"

namespace Montjoie
{
  
  /*****************************
   * VarComputationRCS_Laplace *
   *****************************/
  
  
  //! Sets values of attributes to default values
  template<class Dimension>
  void VarComputationRCS_Laplace<Dimension>::InitDefaultValues()
  {
    store_values_on_surface = false;
    //store_values_on_surface = true;
    n0_current = 0;
    t0 = 0.0;
    num_t0 = 0;
    t0_start = 0.0;
    coef_u = 1.0;
    
    compute_L2_norm = false;
    file_name_normL2 = string("norm.dat");
    t_begin_calculL2 = 0.0;
    
    size_buffer = 100;
    first_output_ext = true;
  }

  
  //! initialization before computation of the far field
  template<class Dimension>
  void VarComputationRCS_Laplace<Dimension>::InitComputationRCS(bool assemble)
  {
    VarComputationRCS_Base<Dimension>::InitComputationRCS(assemble);
    
    if (this->PointsOutside.GetM() <= 0)
      return;
    
    // computation of tmin, tmax
    Real_wp rmin(1e300), rmax(0);
    for (int i = 0; i < this->PointsOutside.GetM(); i++)
      for (int j = 0; j < this->var_mesh.GetNbAllQuadraturePoints(); j++)
        {
          Real_wp r = this->PointsOutside(i).Distance(this->var_mesh.GetQuadraturePoint(j));
          rmin = min(r, rmin);
          rmax = max(r, rmax);
        }
    
    int ref = -1;
    if (this->var_mesh.GetNbBoundary() > 0)
      {
        int num_elem = this->var_mesh.GetElementNumberOfSurface(0);
        ref = var_problem.mesh.Element(num_elem).GetReference();
      }
    
    Real_wp c0 = 0;
    if (ref >= 0)
      c0 = var_problem.GetVelocityOfMedia(ref);
    
#ifdef SELDON_WITH_MPI
    Real_wp rmin2 = rmin; Vector<int64_t> xtmp;
    MpiAllreduce(var_problem.comm_group_mode, &rmin2, xtmp,
                 &rmin, 1, MPI_MIN);
    
    Real_wp rmax2 = rmax;
    MpiAllreduce(var_problem.comm_group_mode, &rmax2, xtmp,
                 &rmax, 1, MPI_MAX);
    
    Real_wp c0_tmp = c0;
    MpiAllreduce(var_problem.comm_group_mode, &c0_tmp, xtmp,
                 &c0, 1, MPI_MAX);        
#endif
    
    invCinfinity = 1.0/c0;
    
    tmin = rmin * invCinfinity;
    int nmin = toInteger(floor(tmin/this->deltat));
    tmin = nmin*this->deltat;
    tmax = rmax * invCinfinity;
    int nmax = toInteger(ceil(tmax/this->deltat));
    tmax = nmax*this->deltat;
  }
  
  
  //! datas are read from disk to resume the simulation and the computation of far field
  template<class Dimension>
  void VarComputationRCS_Laplace<Dimension>
  ::LoadDatas(ParameterOutputReprise& output_reprise)
  {
    //typedef typename TypeElement::Dimension Dimension;
    if (this->PointsOutside.GetM() <= 0)
      return;
    
    coef_u = 1.0/(4.0*pi_wp);
    output_reprise.Read(this->PointsOutside);
    // this->PointsOutside.Read(prefix + "PotentielRetardePoints" + suffix);
    
    if (store_values_on_surface)
      {
        // reprise non implementee dans ce cas
        cout << " Not implemented" << endl;
        abort();
      }
    else
      {
        //string file_name = prefix + "PotentielRetardeDatas" + suffix;        
        //ifstream file_in(file_name.data());
        
        // on lit tmin, tmax, invCinfinity, t0
        VectReal_wp param(4);
        output_reprise.Read(param);
        tmin = param(0); tmax = param(1); invCinfinity = param(2); t0 = param(3);
        num_t0 = 0; t0_start = t0;
        
        //file_in.read(reinterpret_cast<char*>(&tmin), sizeof(Real_wp));
        //file_in.read(reinterpret_cast<char*>(&tmax), sizeof(Real_wp));
        //file_in.read(reinterpret_cast<char*>(&invCinfinity), sizeof(Real_wp));
        //file_in.read(reinterpret_cast<char*>(&t0), sizeof(Real_wp));
        
        // on lit n0_current et la taille de valU_outside
        int nb_iterate;
        output_reprise.Read(nb_iterate);
        output_reprise.Read(n0_current);
        
        //file_in.read(reinterpret_cast<char*>(&n0_current), sizeof(int));
        //file_in.read(reinterpret_cast<char*>(&nb_iterate), sizeof(int));
        this->valU_outside.Reallocate(nb_iterate);
        
        // puis les donnees
        for (int i = 0; i < this->valU_outside.GetM(); i++)
          output_reprise.Read(this->valU_outside(i));
        
        output_reprise.Read(normL2);
        //this->valU_outside(i).Read(file_in);        
        //normL2.Read(file_in);
        
        //file_in.close();
      }
  }
  
  
  //! variables needed to resume the computation of far field are stored in the disk
  template<class Dimension>
  void VarComputationRCS_Laplace<Dimension>
  ::SaveDatas(ParameterOutputReprise& output_reprise)
  {
    //typedef typename TypeElement::Dimension Dimension;
    if (this->PointsOutside.GetM() <= 0)
      return;
    
    // on ecrit les points
    output_reprise.Write(this->PointsOutside);
    
    //this->PointsOutside.Write(prefix + "PotentielRetardePoints" + suffix);
    
    if (store_values_on_surface)
      {
        // reprise non implementee dans ce cas
        cout << " Not implemented" << endl;
        abort();
      }
    else
      {
        // on ecrit tmin, tmax, invCinfinity, t0
        VectReal_wp param(4);
        param(0) = tmin; param(1) = tmax; param(2) = invCinfinity; param(3) = t0;
        output_reprise.Write(param);
        
        //file_out.write(reinterpret_cast<char*>(&tmin), sizeof(Real_wp));
        //file_out.write(reinterpret_cast<char*>(&tmax), sizeof(Real_wp));
        //file_out.write(reinterpret_cast<char*>(&invCinfinity), sizeof(Real_wp));
        //file_out.write(reinterpret_cast<char*>(&t0), sizeof(Real_wp));
        
        // on ecrit n0_current et la taille de valU_outside
        int nb_iterate = this->valU_outside.GetM();
        output_reprise.Write(nb_iterate);
        output_reprise.Write(n0_current);
        
        // file_out.write(reinterpret_cast<char*>(&n0_current), sizeof(int));
        // file_out.write(reinterpret_cast<char*>(&nb_iterate), sizeof(int));
        
        // puis les donnees
        for (int i = 0; i < this->valU_outside.GetM(); i++)
          output_reprise.Write(this->valU_outside(i));
        
        output_reprise.Write(normL2);
        //this->valU_outside(i).Write(file_out);        
        //normL2.Write(file_out);
        
        //file_out.close();
      }
  }
  
  
  //! computation of the far field is advanced, output is written if needed
  template<class Dimension>
  void VarComputationRCS_Laplace<Dimension>
  ::WriteOutput(int nt, const Real_wp& t, const Real_wp& dt,
		const VectReal_wp& Pn, const VectReal_wp& PnPoint)
  {
    if (this->PointsOutside.GetM() <= 0)
      return;
    
    bool double_prec = var_output.OutputWrittenInDoublePrecision();
    
    if (nt == 0)
      {
        int nb_points = this->PointsOutside.GetM();
        if (store_values_on_surface)
          {
            int N = toInteger(ceil((tmax - tmin)/dt)) + 1;
            delayed_trace_Pn.Reallocate(N);
            delayed_trace_PnPoint.Reallocate(N);
            delayed_trace_dPnDn.Reallocate(N);
            delayed_trace_time.Reallocate(N);
            delayed_trace_time.Fill(0);
            
            int Nall = this->var_mesh.GetNbAllQuadraturePoints();
            for (int i = 0; i < N; i++)
              {
                delayed_trace_Pn(i).Reallocate(Nall);
                delayed_trace_PnPoint(i).Reallocate(Nall);
                delayed_trace_dPnDn(i).Reallocate(Nall);
              }
            
            FillZero(delayed_trace_Pn);
            FillZero(delayed_trace_PnPoint);
            FillZero(delayed_trace_dPnDn);
          }
        else
          {
            int N = toInteger(ceil((tmax - tmin)/this->deltat)) + 1;
            valU_outside.Reallocate(nb_points);
            for (int i = 0; i < nb_points; i++)
              valU_outside(i).Reallocate(N);
            
            FillZero(valU_outside);
          }
        
        t0 = t + tmin; n0_current = 0; coef_u = 1.0/(4.0*pi_wp);
        num_t0 = 0; t0_start = t0;
        int Np = toInteger(floor((tmin-t) / this->deltat));
#ifdef SELDON_WITH_MPI
        int rank_proc; MPI_Comm_rank(var_problem.comm_group_mode, &rank_proc);
        if (rank_proc == 0)
          {
#endif
            string file_name = var_output.output_points_outside_param.GetTotalFieldFile();
            ofstream file_out(file_name.data());            
            VectReal_wp zeros(nb_points+1); zeros.Fill(0);
            // on ecrit la precision
            WriteBinaryTypeData_DoubleOrFloat(zeros, file_out, double_prec);
            
            VectReal_wp points(Dimension::dim_N*nb_points);
            points.Fill(0);
            for (int i = 0; i < nb_points; i++)
              CopyVector(this->PointsOutside(i), i, points);
            
            // on ecrit tous les points
            WriteBinaryDoubleOrFloat(points, file_out, double_prec);
            
            // puis on ecrit les 0 pour t=0 a tmin
            DISP(nb_points);
            for (int i = 0; i < Np; i++)
              {
                zeros(0) = t + i*this->deltat;
                WriteBinaryDoubleOrFloat(zeros, file_out, double_prec);
              }
            
            file_out.close();
            
#ifdef SELDON_WITH_MPI
          }
#endif
        
        if (compute_L2_norm)
          {
            normL2.Reallocate(nb_points);
            normL2.Fill(0);
          }
        /*Matrix<Real_wp> sis;
          sis.ReadText("SismoCubeRefExact.dat");
          
          tmax_ref = sis(sis.GetM()-1, 0);
          deltat_ref = sis(1, 0) - sis(0, 0);
        
          int Ns = sis.GetM();
          u_time_ref.Reallocate(Ns);
          du_time_ref.Reallocate(Ns);
          uP_time_ref.Reallocate(Ns);
          uP_time_ref.Fill(0); u_time_ref.Fill(0);
          du_time_ref.Fill(0);
          for (int i = 0; i < Ns; i++)
          {
          if ((i > 0) && (i < Ns-1))
          uP_time_ref(i) = (sis(i+1, 1) - sis(i-1, 1))/(2.0*deltat_ref);
        
          u_time_ref(i) = sis(i, 1);
          du_time_ref(i) = -sis(i, 4);
          }
        */
      }
    
    if (store_values_on_surface)
      {
        int n = nt%delayed_trace_Pn.GetM();
        delayed_trace_time(n) = t;
        
        this->var_mesh.ComputeEnHnOnBoundary(var_problem, Pn,
                                             delayed_trace_Pn(n),
                                             delayed_trace_dPnDn(n), false, true);
        
        this->var_mesh.ComputeEnHnOnBoundary(var_problem, PnPoint,
                                             delayed_trace_PnPoint(n),
                                             delayed_trace_dPnDn(n), false, false);
      }
    else
      {
        this->var_mesh.ComputeEnHnOnBoundary(var_problem, Pn,
                                             Pn_current,
                                             dPnDn_current, false, true);
        
        this->var_mesh.ComputeEnHnOnBoundary(var_problem, PnPoint,
                                             PnPoint_current,
                                             dPnDn_current, false, false);
        
        int N = this->PointsOutside.GetM();
        for (int i = 0; i < N; i++)
          {
            this->UpdateTimeIntegralRepresentation(nt, t, dt, this->var_mesh,
                                                   this->PointsOutside(i),
                                                   valU_outside(i));
          }
      }
    
    if ((t+tmin) >= (t0 - epsilon_machine))
      {
        int N = this->PointsOutside.GetM();
        VectReal_wp sol(N);
        Real_wp scal_u;
        if (store_values_on_surface)
          for (int i = 0; i < N; i++)
            {
              this->ComputeTimeIntegralRepresentation(nt, t, dt, this->var_mesh,
                                                      this->PointsOutside(i), scal_u);
              
              sol(i) = scal_u;
            }
        else
          {
            for (int i = 0; i < N; i++)
              {
                sol(i) = coef_u*valU_outside(i)(n0_current);
                valU_outside(i)(n0_current) = 0.0;
              }
            
            n0_current = (n0_current+1)%valU_outside(0).GetM();
          }
        
#ifdef SELDON_WITH_MPI
        int rank_proc; MPI_Comm_rank(var_problem.comm_group_mode, &rank_proc);
        int nb_proc; MPI_Comm_size(var_problem.comm_group_mode, &nb_proc);
        if (nb_proc > 1)
          {
            VectReal_wp tmp(sol); Vector<int64_t> xtmp;
            MpiReduce(var_problem.comm_group_mode, tmp, xtmp,
                      sol, N, MPI_SUM, 0);
          }
#else
        int rank_proc(0);
#endif
        
        if (rank_proc == 0)
          {
            VectReal_wp vec(N+1);
            vec(0) = t + tmin;
            for (int i = 0; i < N; i++)
              vec(i+1) = sol(i);
            
            //            DISP(first_output_ext);
            if (first_output_ext)
              {
                output_buffer.SetBinary();
                output_buffer.Init(var_output.output_points_outside_param.GetTotalFieldFile(),
                                   size_buffer, false, false);
                
                first_output_ext = false;
              }
            
            //            DISP(vec);
            output_buffer.AddVect(vec);
            
            if (compute_L2_norm && (t > t_begin_calculL2) )
              {
                for (int i = 0; i < N; i++)
                  normL2(i) += this->deltat*sol(i)*sol(i);
                
                ofstream file_norm(file_name_normL2.data());
                for (int i = 0; i < N; i++)
                  vec(i+1) = sqrt(normL2(i) / (t+tmin - t_begin_calculL2));
                
                WriteBinaryDoubleOrFloat(vec, file_norm, double_prec);
                file_norm.close();
              }
          }
        
        num_t0++;
        t0 = t0_start + num_t0*this->deltat;
      }
  }
  
  
  //! not implemented in 2-D
  template<>
  void VarComputationRCS_Laplace<Dimension2>
  ::ComputeTimeIntegralRepresentation(int nt, const Real_wp& t, const Real_wp& dt,
				      const MeshInterpolationFEM<Dimension2>& mesh_, const R2& pointX,
				      Real_wp& scal_u) const
  {
    abort();
  }
  
  
  //! not implemented in 2-D
  template<>
  void VarComputationRCS_Laplace<Dimension2>
  ::UpdateTimeIntegralRepresentation(int nt, const Real_wp& t, const Real_wp& dt,
				     const MeshInterpolationFEM<Dimension2>& mesh_, const R2& pointX,
				     VectReal_wp& scal_u) const
  {
    abort();
  }
  
  
#ifdef MONTJOIE_WITH_THREE_DIM
  //! computation of u^{far} from u, du/dt, du/dn on a surface of integration
  template<>
  void VarComputationRCS_Laplace<Dimension3>
  ::ComputeTimeIntegralRepresentation(int nt, const Real_wp& t, const Real_wp& dt,
				      const MeshInterpolationFEM<Dimension3>& mesh_, const R3& pointX,
				      Real_wp& scal_u) const
  {
    scal_u = 0.0;
    R3 pointY, normaleY, XmY;
    Real_wp R, invR, invR2, ny_dot_XmY, tau, tn, invDt = 1.0/dt, L;
    int nb_iter = delayed_trace_Pn.GetM();
    Real_wp Pn, Pn_point, dPn_dn(0), poids;
    //Real_wp PnRef(0), PnPointRef(0), dPn_dn_ref(0);
    for (int k = 0; k < mesh_.GetNbAllQuadraturePoints(); k++)
      {
        // on recupere l'oppose de la normale pour qu'elle soit sortante
        normaleY(0) = -mesh_.GetQuadratureNormale(k)(0);
        normaleY(1) = -mesh_.GetQuadratureNormale(k)(1);
        normaleY(2) = -mesh_.GetQuadratureNormale(k)(2);
        
        pointY =  mesh_.GetQuadraturePoint(k);
        poids = mesh_.GetQuadratureWeight(k);
        XmY = pointX - pointY;
        ny_dot_XmY = DotProd(normaleY, XmY);
        
        R = Norm2(XmY);
        invR = 1.0/R; invR2 = invR*invR;
        
        // temps retarde
        tau = t+tmin - R*invCinfinity;
        int num_iter  = toInteger(ceil((t-tau)/dt));
        if ((num_iter < 0) || (num_iter > nb_iter))
          {
            cout << "Cas impossible" << endl;
            abort();
          }
        
        int n = (nt - num_iter)%nb_iter;
        if (n < 0 )
          n += nb_iter;
        
        int np1 = n + 1;
        if (np1 == nb_iter)
          np1 = 0;
        
        // on recupere P, dP/dt et dP/dn a ce temps la
        tn = t - num_iter*dt;
        L = (tau-tn)*invDt;
        Pn = delayed_trace_Pn(n)(k)*(1.0-L) + delayed_trace_Pn(np1)(k)*L;
        Pn_point = delayed_trace_PnPoint(n)(k)*(1.0-L) + delayed_trace_PnPoint(np1)(k)*L;
        dPn_dn = delayed_trace_dPnDn(n)(k)*(1.0-L) + delayed_trace_dPnDn(np1)(k)*L;
        
        /*if (tau > 0)
          if (abs(tn - delayed_trace_time(n)) > 1e-12)
          {
          DISP(tn); DISP(delayed_trace_time(n));
          abort();
          }
        */  
        /*
          if ((tau > 0) && (tau < tmax_ref))
          {
          int nref = toInteger(floor(tau/deltat_ref));
          Real_wp t0_ref = nref*deltat_ref;
          Real_wp lambda = (tau - t0_ref) / deltat_ref;
          PnRef = (1.0-lambda)*u_time_ref(nref) + lambda*u_time_ref(nref+1);
          PnPointRef = (1.0-lambda)*uP_time_ref(nref) + lambda*uP_time_ref(nref+1);
          dPn_dn_ref = (1.0-lambda)*du_time_ref(nref) + lambda*du_time_ref(nref+1);
            
          if (abs(tn - delayed_trace_time(n)) > 1e-12)
          {
          DISP(tn); DISP(delayed_trace_time(n));
          abort();
          }
            
          if ((abs(Pn-PnRef) > 1e-8) || (abs(PnPointRef - Pn_point) > 1e-8)
          || (abs(dPn_dn_ref - dPn_dn) > 1e-8))
          {
          DISP(Pn); DISP(PnRef);
          DISP(Pn_point); DISP(PnPointRef);
          DISP(dPn_dn); DISP(dPn_dn_ref);       
          DISP(tau); DISP(tn); DISP(L);
          DISP(delayed_trace_Pn(n)(k)(0));
          DISP(delayed_trace_Pn(np1)(k)(0));
          DISP(delayed_trace_PnPoint(n)(k)(0));
          DISP(delayed_trace_PnPoint(np1)(k)(0));
          DISP(delayed_trace_dPnDn(n)(k)(0));
          DISP(delayed_trace_dPnDn(np1)(k)(0));
          abort();
          }            
          }
          else
          {
          PnRef = 0; PnPointRef = 0; dPn_dn_ref = 0;
          }
        
          Pn = PnRef; 
          Pn_point = PnPointRef;
          dPn_dn = dPn_dn_ref;
        */
        scal_u += poids*(ny_dot_XmY*invR2*(invR*Pn + invCinfinity*Pn_point)
                         + invR*dPn_dn);
      }
    
    scal_u *= 1.0/(4.0*pi_wp);
  }
  
  
  //! updating values of u^{far} from u, du/dt, du/dn on a surface of integration
  template<>
  void VarComputationRCS_Laplace<Dimension3>
  ::UpdateTimeIntegralRepresentation(int nt, const Real_wp& t, const Real_wp& dt,
				     const MeshInterpolationFEM<Dimension3>& mesh_, const R3& pointX,
				     VectReal_wp& scal_u) const
  {
    R3 pointY, normaleY, XmY;
    Real_wp R, invR(0), invR2(0), ny_dot_XmY(0), tau(0), invDt = 1.0/dt, L;
    Real_wp Pn(0), Pn_point(0), dPn_dn(0), poids(0); bool precompute;
    Real_wp t_advance, tm;
    for (int k = 0; k < mesh_.GetNbAllQuadraturePoints(); k++)
      {
        // calcul de R = |X - Y|
        pointY =  mesh_.GetQuadraturePoint(k);
        XmY = pointX - pointY;
        
        R = Norm2(XmY);
        
        // temps 'avance' : t + |x-y|/c
        t_advance = t + R*invCinfinity;
        
        // on cherche tm tel que tm >= t + |x-y|/c
        int n = toInteger(ceil((t_advance - t0) / this->deltat));
        tm = t0 + n*this->deltat;
        
        // si tm < t + |x-y|/c + dt, on met a jour u(tm)
        precompute = false;
        if (tm < t_advance + dt)
          {
            precompute = true;
            // on recupere l'oppose de la normale pour qu'elle soit sortante
            normaleY(0) = -mesh_.GetQuadratureNormale(k)(0);
            normaleY(1) = -mesh_.GetQuadratureNormale(k)(1);
            normaleY(2) = -mesh_.GetQuadratureNormale(k)(2);
            
            ny_dot_XmY = DotProd(normaleY, XmY);
            invR = 1.0/R; invR2 = invR*invR;
            poids = mesh_.GetQuadratureWeight(k);
            
            // calcul de tau = tm - |x-y|/c
            tau = tm - R*invCinfinity;
            L = (tau - t)*invDt;
            
            // on ajoute la contribution
            Pn = Pn_current(k)*(1.0-L);
            Pn_point = PnPoint_current(k)*(1.0-L);
            dPn_dn = dPnDn_current(k)*(1.0-L);
            int np = (n + n0_current)%scal_u.GetM();
            
            scal_u(np) += poids*(ny_dot_XmY*invR2*(invR*Pn + invCinfinity*Pn_point)
                                 + invR*dPn_dn);
          }
        
        // on cherche tm tel que tm <= t + |x-y|/c
        if (tm != t_advance)
          {
            n = toInteger(floor((t_advance-t0) / this->deltat));
            tm = t0 + n*this->deltat;
            
            // si tm > t + |x-y|/c - dt, on met a jour u(tm)
            if (tm > t_advance - dt)
              {
                if (!precompute)
                  {
                    normaleY(0) = -mesh_.GetQuadratureNormale(k)(0);
                    normaleY(1) = -mesh_.GetQuadratureNormale(k)(1);
                    normaleY(2) = -mesh_.GetQuadratureNormale(k)(2);
                    
                    ny_dot_XmY = DotProd(normaleY, XmY);
                    invR = 1.0/R; invR2 = invR*invR;
                    poids = mesh_.GetQuadratureWeight(k);                    
                  }
                
                // calcul de tau = tm - |x-y|/c
                tau = tm - R*invCinfinity;
                L = (tau - t + dt)*invDt;
                // on ajoute la contribution
                Pn = Pn_current(k)*L;
                Pn_point = PnPoint_current(k)*L;
                dPn_dn = dPnDn_current(k)*L;
                int np = (n + n0_current)%scal_u.GetM();
                
                scal_u(np) += poids*(ny_dot_XmY*invR2*(invR*Pn + invCinfinity*Pn_point)
                                     + invR*dPn_dn);                
              }
          }
      }
  }
#endif
  
  
  /*********************************
   * IncidentWaveProjector_Laplace *
   *********************************/
  
  
  //! initialisation of source with time t0
  template<class Dimension>
  void IncidentWaveProjector_Laplace<Dimension>::Init(const Real_wp& t0)
  {
    wave_pulse = this->incident_wave.GetTimeSource();
    
    if (var_source.GetNbParameterSource(0) <= 0)
      {
        cout << "Parameters needed for source" << endl;
        abort();
      }
    
    this->t = t0;
    const Vector<string>& param = var_source.GetParameterSource(0, 0);
    c0 = var_boundary.GetMaximumVelocityPML();    
    k_wave = var_problem.GetWaveVector();
    Mlt(1.0/Norm2(k_wave), k_wave);
    
    if (param(1) == "AUTO")
      {    
        offset_arg = GetMinimumDotProd(k_wave, var_problem.GetXmin(), var_problem.GetXmax(),
                                       var_problem.GetYmin(), var_problem.GetYmax(),
                                       var_problem.GetZmin(), var_problem.GetZmax());        
      }
    else
      offset_arg = to_num<Real_wp>(param(1));
  }
  
  
  //! Evaluates incident field
  template<class Dimension>
  void IncidentWaveProjector_Laplace<Dimension>
  ::EvaluateFunction(int i, int j, const R_N& x, VectReal_wp& f)
  {
    if (wave_pulse == NULL)
      return;

    Real_wp arg = this->t + (offset_arg - DotProd(k_wave, x))/c0;
    Real_wp fs = wave_pulse->Evaluate(arg);
    f(0) = fs;    
  }
  
  
  //! gradient of incident field
  template<class Dimension>
  void IncidentWaveProjector_Laplace<Dimension>
  ::EvaluateFunctionGradient(int i, int j, const R_N& x, VectReal_wp& f, VectReal_wp& df)
  {
    if (wave_pulse == NULL)
      return;
        
    Real_wp arg = this->t + (offset_arg - DotProd(k_wave, x))/c0;
    Real_wp fs = wave_pulse->Evaluate(arg);
    Real_wp dfs = -wave_pulse->EvaluateDerivative(arg)/c0;
    
    f(0) = fs;
    for (int k = 0; k < Dimension::dim_N; k++)
      df(k) = dfs*k_wave(k);
  }


  /********************************
   * DiffractedWaveSource_Laplace *
   ********************************/
  

  //! initialisation of source with time t0
  template<class Dimension>
  Real_wp DiffractedWaveSource_Laplace<Dimension>
  ::Init(const Real_wp& t0, const Real_wp& dt, int print_level, int n, bool scal)
  {
#ifdef SELDON_WITH_MPI
    int rank_proc; MPI_Comm_rank(var_problem.comm_group_mode, &rank_proc);
#else
    int rank_proc(0);
#endif

    Real_wp tlimit(0);
    wave_pulse = incident_wave.GetTimeSource();
    
    if (var_source.GetNbParameterSource(0) <= 0)
      {
        abort();
      }
    
    t = t0;
    scalar_eq = scal;
    
    const Vector<string>& param = var_source.GetParameterSource(0, 0);
    c0 = var_boundary.GetMaximumVelocityPML();    
    k_wave = var_problem.GetWaveVector();
    Mlt(1.0/Norm2(k_wave), k_wave);
    
    if (param(1) == "AUTO")
      {    
        offset_arg = GetMinimumDotProd(k_wave, var_problem.GetXmin(), var_problem.GetXmax(),
                                       var_problem.GetYmin(), var_problem.GetYmax(),
                                       var_problem.GetZmin(), var_problem.GetZmax());        

        Real_wp off_max = GetMaximumDotProd(k_wave, var_problem.GetXmin(), var_problem.GetXmax(),
					    var_problem.GetYmin(), var_problem.GetYmax(),
					    var_problem.GetZmin(), var_problem.GetZmax());        
	
	tlimit = off_max - offset_arg;
      }
    else
      {
	offset_arg = to_num<Real_wp>(param(1));
	tlimit = to_num<Real_wp>(param(2));
      }
    
    if ((t0 <= dt+20*epsilon_machine) && (print_level >= 2) && (rank_proc == 0))
      cout << "Offset in expression of plane wave " << offset_arg << endl;
    
    num_deriv = n;
    
    if (n > 1)
      {
        cout << "Higher-order derivatives not implemented" << endl;
        abort();
      }
    
    return tlimit;
  }
  

  //! returns true for Dirichlet condition
  template<class Dimension>
  bool DiffractedWaveSource_Laplace<Dimension>::PresenceDirichlet() const
  {
    return dirichlet_cond;
  }

  
  //! For Dirichlet condition
  template<class Dimension>
  void DiffractedWaveSource_Laplace<Dimension>
  ::EvaluateFunction(int i, int j, const R_N& x, VectReal_wp& f)
  {
    if (wave_pulse == NULL)
      return;

    Real_wp arg = t + (offset_arg - DotProd(k_wave, x))/c0;
    Real_wp fs(0);
    
    if (num_deriv == 0)
      fs = wave_pulse->Evaluate(arg);
    else
      fs = wave_pulse->EvaluateDerivative(arg);
    
    f(0) = -fs;    
  }
  
  
  //! initialisation of element i for volume integrals
  template<class Dimension>
  void DiffractedWaveSource_Laplace<Dimension>
  ::InitElement(int num_elem, const VectR_N& s)
  {
    // not implemented
  }
  
  
  //! returns true if there is a volume integral
  template<class Dimension>
  bool DiffractedWaveSource_Laplace<Dimension>
  ::IsNonNullVolumetricSource(const VectR_N& s)
  {
    return false;
  }
    
  
  //! fills f for volume integral \int f \varphi
  template<class Dimension>
  void DiffractedWaveSource_Laplace<Dimension>
  ::EvaluateVolumetricSource(int i, int j, const R_N& x, VectReal_wp& f)
  {
  }
  
  
  //! returns true if there a volume integral \int f \nabla \varphi
  template<class Dimension>
  bool DiffractedWaveSource_Laplace<Dimension>
  ::IsNonNullGradientSource(const VectR_N& s)
  {
    return false;
  }
    
  
  //! fills f for volume integral \int f \nabla \varphi
  template<class Dimension>
  void DiffractedWaveSource_Laplace<Dimension>
  ::EvaluateGradientSource(int i, int j, const R_N& x, VectReal_wp& f)
  {
  }
    
  
  //! initialisation of surface i
  template<class Dimension>
  void DiffractedWaveSource_Laplace<Dimension>
  ::InitSurface(int i, int num_face, int num_elem, int num_loc)
  {
    VirtualSourceFEM<Real_wp, Dimension>::InitSurface(i, num_face, num_elem, num_loc);
  }
  
  
  //! fills f for surface integral \int f \varphi
  template<class Dimension>
  void DiffractedWaveSource_Laplace<Dimension>
  ::EvaluateSurfacicSource(int k, const SetPoints<Dimension>& PointsElem,
                           const SetMatrices<Dimension>& MatricesElem, VectReal_wp& f)
  {
    if (wave_pulse == NULL)
      return;

    if ((!scalar_eq) || var_problem.InsidePML(this->num_elem_))
      {
        f.Fill(0);
        return;
      }

    // coef = k . n
    Real_wp coef = DotProd(k_wave, MatricesElem.GetNormaleQuadratureBoundary(k));
    int ref = var_problem.mesh.Element(this->num_elem_).GetReference();
    coef *= var_laplace.ref_mu(ref).GetConstant()(0, 0);
    
    // on evalue ct - k.x + offset_arg
    Real_wp arg = t + (offset_arg - DotProd(k_wave, PointsElem.GetPointQuadratureBoundary(k)))/c0;
    
    if (num_deriv == 0)
      {
        // primitive en temps
        Real_wp fs = wave_pulse->Evaluate(arg)/c0;
        f(0) = coef*fs;
      }
    else
      {
        // on evalue f' en ce point
	Real_wp fprime = wave_pulse->EvaluateDerivative(arg)/c0; 
	
        // la source vaut -du/dn = k.n f'(ct - k.x + offset_arg)/c0
        f(0) = coef*fprime;
      }
  }
  
  
  //! returns true if there is a surface integral \int f \varphi
  template<class Dimension>
  bool DiffractedWaveSource_Laplace<Dimension>
  ::IsNonNullSurfacicSource(int ref)
  {
    if (var_problem.mesh.GetBoundaryCondition(ref) == BoundaryConditionEnum::LINE_NEUMANN)
      return true;
    
    return false;
  }

  
  /***************************
   * TotalWaveSource_Laplace *
   ***************************/
  
  
  //! initialisation of source with time t0  
  template<class Dimension>
  Real_wp TotalWaveSource_Laplace<Dimension>
  ::Init(const Real_wp& t0, const Real_wp& dt, int print_level, int n, bool scal)
  {
#ifdef SELDON_WITH_MPI
    int rank_proc; MPI_Comm_rank(var_problem.comm_group_mode, &rank_proc);
#else
    int rank_proc(0);
#endif

    Real_wp tlimit(0);
    wave_pulse = incident_wave.GetTimeSource();
    
    if (var_source.GetNbParameterSource(0) <= 0)
      {
        abort();
      }
    
    t = t0;
    scalar_eq = scal;
    
    const Vector<string>& param = var_source.GetParameterSource(0, 0);
    c0 = var_boundary.GetMaximumVelocityPML();
    k_wave = var_problem.GetWaveVector();
    Mlt(1.0/Norm2(k_wave), k_wave);
    
    if (param(1) == "AUTO")
      {
        offset_arg = GetMinimumDotProd(k_wave, var_problem.GetXmin(), var_problem.GetXmax(),
                                       var_problem.GetYmin(), var_problem.GetYmax(),
                                       var_problem.GetZmin(), var_problem.GetZmax());        

        Real_wp off_max = GetMaximumDotProd(k_wave, var_problem.GetXmin(), var_problem.GetXmax(),
					    var_problem.GetYmin(), var_problem.GetYmax(),
					    var_problem.GetZmin(), var_problem.GetZmax());        
	
	tlimit = off_max - offset_arg;
      }
    else
      {
	offset_arg = to_num<Real_wp>(param(1));
	tlimit = to_num<Real_wp>(param(2));
      }
    
    if ((t0 <= dt+20*epsilon_machine) && (print_level >= 2) && (rank_proc == 0))
      cout << "Offset in expression of plane wave " << offset_arg << endl;
    
    num_deriv = n;
    if (n > 1)
      {
        cout << "Higher-order derivatives not implemented" << endl;
        abort();
      }

    return tlimit;
  }
  

  //! For Dirichlet condition
  template<class Dimension>
  void TotalWaveSource_Laplace<Dimension>
  ::EvaluateFunction(int i, int j, const R_N& x, VectReal_wp& f)
  {
  }
  
  
  //! initialisation of element i for volume integrals
  template<class Dimension>
  void TotalWaveSource_Laplace<Dimension>::InitElement(int num_elem, const VectR_N& s)
  {
  }
  
  
  //! returns true if there is a volume integral
  template<class Dimension>
  bool TotalWaveSource_Laplace<Dimension>
  ::IsNonNullVolumetricSource(const VectR_N& s)
  {
    return false;
  }
    

  //! fills f for volume integral \int f \varphi  
  template<class Dimension>
  void TotalWaveSource_Laplace<Dimension>
  ::EvaluateVolumetricSource(int i, int j, const R_N& x, VectReal_wp& f)
  {
  }
  
  
  //! returns true if there a volume integral \int f \nabla \varphi
  template<class Dimension>
  bool TotalWaveSource_Laplace<Dimension>
  ::IsNonNullGradientSource(const VectR_N& s)
  {
    return false;
  }
    

  //! fills f for volume integral \int f \nabla \varphi  
  template<class Dimension>
  void TotalWaveSource_Laplace<Dimension>
  ::EvaluateGradientSource(int i, int j, const R_N& x, VectReal_wp& f)
  {
  }
    
  
  //! initialisation of surface i
  template<class Dimension>
  void TotalWaveSource_Laplace<Dimension>
  ::InitSurface(int i, int num_face, int num_elem, int num_loc)
  {
    VirtualSourceFEM<Real_wp, Dimension>::
      InitSurface(i, num_face, num_elem, num_loc);
  }
  
  
  //! fills f for surface integral \int f \varphi  
  template<class Dimension>
  void TotalWaveSource_Laplace<Dimension>
  ::EvaluateSurfacicSource(int k, const SetPoints<Dimension>& PointsElem,
                           const SetMatrices<Dimension>& MatricesElem, VectReal_wp& f)
  {
    if (wave_pulse == NULL)
      return;
    
    if (!scalar_eq)
      {
        f.Fill(0);
        return;
      }

    // coef = 1 - k . n
    Real_wp coef = 1.0 - DotProd(k_wave, MatricesElem.GetNormaleQuadratureBoundary(k));
    int ref = var_problem.mesh.Element(this->num_elem_).GetReference();
    coef *= var_laplace.ref_mu(ref).GetConstant()(0, 0);
    
    // on evalue ct - k.x + offset_arg
    Real_wp arg = t + (offset_arg - DotProd(k_wave, PointsElem.GetPointQuadratureBoundary(k)))/c0;
    
    //DISP(num_deriv);
    if (num_deriv == 0)
      {
        Real_wp fs = wave_pulse->Evaluate(arg)/c0;
	f(0) = coef*fs;
      }
    else
      {
        // on evalue f' en ce point
        Real_wp fprime = wave_pulse->EvaluateDerivative(arg)/c0; 
	        
        // la source vaut du/dn + 1/c du/dt = (1 - k.n) f'(ct - k.x + offset_arg)
        f(0) = coef*fprime;
      }
  }
  
  
  //! returns true if there is a surface integral \int f \varphi
  template<class Dimension>
  bool TotalWaveSource_Laplace<Dimension>::IsNonNullSurfacicSource(int ref)
  {
    if (var_problem.mesh.GetBoundaryCondition(ref) == BoundaryConditionEnum::LINE_ABSORBING)
      return true;
    
    return false;
  }
  

  /*********************
   * LaplaceEquationDG *
   *********************/  
  

  //! returns the matrix M, in the integral \f$ \int_K M u v \f$ 
  /*!
    \param[in] vars given problem
    \param[in] i number of the element where M needs to be evaluated
    \param[in] j number of the local quadrature point in the element
    \param[in] nat_mat mass and stiffness coefficients
    \param[in] ref reference of the physical domain
    \param[out] mass matrix M
  */  
  template<class Dimension>
  template<class TypeEquation, class T, class MatMass>
  void LaplaceEquationDG<Dimension>::
  GetTensorMass(const EllipticProblem<TypeEquation>& var,
                int i, int j, const GlobalGenericMatrix<T>& nat_mat, int ref, MatMass& mass)
  {
    FillZero(mass);
    mass(0, 0) = var.ref_rho(ref).GetCoefficient(var, i, j)*nat_mat.GetCoefMass();
    typename Dimension::MatrixN_Nsym mu = var.ref_invMu(ref).GetCoefficient(var, i, j);
    Mlt(-nat_mat.GetCoefMass(), mu);
    for (int p = 0; p < Dimension::dim_N; p++)
      for (int q = 0; q < Dimension::dim_N; q++)
        mass(p+1, q+1) = mu(p, q);    
  }
  

  //! Applies matrix M to a vector
  /*!
    \param[in] var given problem
    \param[in] i number of the element where M needs to be evaluated
    \param[in] j number of the local quadrature point in the element
    \param[in] nat_mat mass and stiffness coefficients
    \param[in] ref reference of the physical domain
    \param[in] Un vector to be multiplied by M
    \param[out] Vn result vector Vn = M Un
  */    
  template<class Dimension> template<class TypeEquation, class T, class Vector1>
  void LaplaceEquationDG<Dimension>::
  ApplyTensorMass(const EllipticProblem<TypeEquation>& var, int i, int k,
                  const GlobalGenericMatrix<T>& nat_mat, int ref, Vector1& Un, Vector1& Vn)
  {
    Vn.Fill(0);
    Vn(0) = Un(0)*var.ref_rho(ref).GetCoefficient(var, i, k);

    typename Dimension::MatrixN_Nsym mu = var.ref_invMu(ref).GetCoefficient(var, i, k);
    for (int n = 0; n < Dimension::dim_N; n++)
      for (int k = 0; k < Dimension:: dim_N; k++)
        Vn(n+1) -= Un(k+1)*mu(n, k);
    
    Vn *= nat_mat.GetCoefMass();
  }
  
  
  //! computation of matrix N associated to the boundary condition
  /*!
    \param[out] Nabc matrix N
    \param[in] normale outward normale
    \param[in] ref reference of the boundary
    \param[in] iquad element number
    \param[in] k local quadrature point number
    \param[in] nat_mat object containing mass and stiffness coefficients
    \param[in] ref_d reference of the element
    \param[in] vars given problem
    \param[in] Fb finite element associated with the element
  */
  template<class Dimension> template<class T, class Matrix1, class TypeEquation>
  void LaplaceEquationDG<Dimension>::
  GetNabc(Matrix1& Nabc, typename Dimension::R_N& normale,
	  int ref, int iquad, int k, const GlobalGenericMatrix<T>& nat_mat, int ref_d,
	  const EllipticProblem<TypeEquation>& vars, const ElementReference<Dimension, 1>& Fb)
  {
    int cond = vars.mesh.GetBoundaryCondition(ref);
    Nabc.Fill(0);
    Real_wp coef = 1.0;
    if (vars.FirstOrderFormulation() && (!vars.GetSymmetrizationUse()))
      coef = -1.0;
    
    if (cond == BoundaryConditionEnum::LINE_DIRICHLET)
      {
	// condition Dirichlet
	// N = (0, n; -n^t, 0)
	for (int i = 0; i < Dimension::dim_N; i++)
	  {
	    Nabc(0, i+1) = -normale(i);
	    Nabc(i+1, 0) = -coef*normale(i);
	  }
	
        Nabc *= nat_mat.GetCoefStiffness();
      }
    else if (cond == BoundaryConditionEnum::LINE_NEUMANN)
      {
	for (int i = 0; i < Dimension::dim_N; i++)
	  {
	    Nabc(0, i+1) = normale(i);
	    Nabc(i+1, 0) = coef*normale(i);
	  }
	
        Nabc *= nat_mat.GetCoefStiffness();
      }
    else if (cond == BoundaryConditionEnum::LINE_IMPEDANCE)
      {
	Real_wp mu = vars.ref_mu(ref_d).GetCoefficient(vars, iquad, k)(0, 0);
	// we multiply by the impedance present in the data file
	vars.MltParamCondition(ref, 0, mu);
	
	// then we form the matrix
	Nabc(0, 0) = mu*nat_mat.GetCoefStiffness();	
	T mu0 = nat_mat.GetCoefStiffness()/mu;
	if (vars.FirstOrderFormulation() && (!vars.GetSymmetrizationUse()))
	  mu0 = -mu0;
	
	for (int i = 0; i < Dimension::dim_N; i++)
	  for (int j = 0; j < Dimension::dim_N; j++)
	    Nabc(i+1, j+1) = -mu0*normale(i)*normale(j);
      }
    else if (cond == BoundaryConditionEnum::LINE_ABSORBING)
      {
        coef = vars.coefficient_impedance_absorbing(ref_d);
	Nabc(0,0) = coef; coef = 1.0/coef;
	if ((!vars.FirstOrderFormulation()) || vars.GetSymmetrizationUse())
	  coef = -coef;
	
	for (int i = 0; i < Dimension::dim_N; i++)
	  for (int j = 0; j < Dimension::dim_N; j++)
	    Nabc(i+1, j+1) = coef*normale(i)*normale(j);
	
        Nabc *= nat_mat.GetCoefStiffness();
      }
  }
  
  
  //! Computation of "penalization" matrices C
  /*!
    \param[out] Nabc penalization matrix C
    \param[in] normale outward normale
    \param[in] iquad element number
    \param[in] k local quadrature point number
    \param[in] nat_mat object containing mass and stiffness coefficients
    \param[in] ref reference of the boundary
    \param[in] vars given problem
    \param[in] Fb finite element associated with the element
   */
  template<class Dimension> template<class Matrix1, class TypeEquation, class T>
  void LaplaceEquationDG<Dimension>::
  GetPenalDG(Matrix1& Nabc, typename Dimension::R_N& normale, int iquad,
             int k, int nf, const GlobalGenericMatrix<T>& nat_mat, int ref, int ref2,
	     const EllipticProblem<TypeEquation>& vars,
	     const ElementReference<Dimension, 1>& Fb)
  {
    Nabc.Zero();
    T coef = nat_mat.GetCoefStiffness();
    
    if (vars.upwind_fluxes)
      {        
        // upwind fluxes : see Hesthaven paper for Maxwell's equations
        Real_wp Ym = vars.coefficient_impedance_absorbing(ref);
        Real_wp Yp = vars.coefficient_impedance_absorbing(ref2);
        Real_wp Zm = 1.0/Ym;
        Real_wp Zp = 1.0/Yp;
        Real_wp inv_Zbar = 1.0/(Zm + Zp); 
        Real_wp inv_Ybar = 1.0/(Ym + Yp);
        T coef_u = (Zp - Zm)*inv_Zbar*nat_mat.GetCoefStiffness();
        T coef_v = (Yp - Ym)*inv_Ybar*nat_mat.GetCoefStiffness();

	// multiplication by two because of PoidsFlux which is divided by two
	coef *= Real_wp(2);

	Nabc(0, 0) = -inv_Zbar*coef;
	if ((!vars.FirstOrderFormulation()) || vars.GetSymmetrizationUse())
	  {
	    coef = -coef;
	    coef_v = -coef_v;
	  }

	for (int i = 0; i < Dimension::dim_N; i++)
	  {
            Nabc(0, i+1) = -coef_u*normale(i);
            Nabc(i+1, 0) = -coef_v*normale(i);
            for (int j = 0; j < Dimension::dim_N; j++)
              Nabc(i+1, j+1) = -inv_Ybar*coef*normale(i)*normale(j);
          }

        return;
      }
    
    Nabc(0,0) = vars.alpha_penalization;
    coef = vars.delta_penalization;
    if ((!vars.FirstOrderFormulation()) || vars.GetSymmetrizationUse())
      coef = -coef;

    for (int i = 0; i < Dimension::dim_N; i++)
      for (int j = 0; j < Dimension::dim_N; j++)
	Nabc(i+1,j+1) = coef*normale(i)*normale(j);
    
    Nabc *= nat_mat.GetCoefStiffness();
  }
  
  
  //! Multiplication by penalization matrices
  /*!
    \param[in] normale outward normale
    \param[in] Vn vector to multiply
    \param[out] Un result vector Un = C*Vn
    \param[in] iquad element number
    \param[in] k local quadrature point number    
    \param[in] nat_mat object containing mass and stiffness coefficients
    \param[in] ref reference of the boundary
    \param[in] vars given problem
    \param[in] Fb finite element associated with the element
   */
  template<class Dimension> template<class Vector1, class T, class TypeEquation>
  void LaplaceEquationDG<Dimension>::
  MltPenalDG(const typename Dimension::R_N& normale, const Vector1& Vn, Vector1& Un,
	     int i, int n, int nf, const GlobalGenericMatrix<T>& nat_mat,
             int ref, int ref2, const EllipticProblem<TypeEquation>& vars,
	     const ElementReference<Dimension, 1>& Fb)
  {
    // int ref_elem = vars.mesh.elements(i).GetReference();
    // condition absorbante sur les sauts
    Un(0) = Vn(0)*vars.alpha_penalization;
    T vloc(0);
    for (int i = 0; i < Dimension::dim_N; i++)
      vloc += Vn(i+1)*normale(i);
    
    vloc *= vars.delta_penalization;
    for (int i = 0; i < Dimension::dim_N; i++)
      Un(i+1) = vloc*normale(i);
    
    Un *= nat_mat.GetCoefStiffness();
  }

  
  //! mutliplication by matrix N associated to the boundary condition
  /*!
    \param[in] normale outward normale
    \param[in] ref reference of the boundary
    \param[in] Vn vector to multiply
    \param[out] Un result Un = N Vn
    \param[in] num_elem1 element number
    \param[in] num_point local quadrature point number
    \param[in] nat_mat object containing mass and stiffness coefficients
    \param[in] ref2 reference of the element
    \param[in] vars given problem
    \param[in] Fb finite element associated with the element
   */  
  template<class Dimension> template<class Vector1, class T, class TypeEquation>
  void LaplaceEquationDG<Dimension>::
  MltNabc(typename Dimension::R_N& normale, int ref, const Vector1& Vn, Vector1& Un,
          int num_elem1, int num_point, const GlobalGenericMatrix<T>& nat_mat, int ref_d,
          const EllipticProblem<TypeEquation>& vars,
	  const ElementReference<Dimension, 1>& Fb)
  {
    // typedef typename TypeElement::Dimension Dimension;
    // typedef typename Dimension::R_N R_N;
    int cond = vars.mesh.GetBoundaryCondition(ref); 
    Un.Fill(0);
    Real_wp coef = 1.0;
    if (vars.FirstOrderFormulation())
      coef = -1.0;

    // DISP(cond); DISP(GlobalVariables::LINE_NEUMANN); exit(0);
    if (cond == BoundaryConditionEnum::LINE_DIRICHLET)
      {
	// condition Dirichlet
	// N = (0, n; -n^t, 0)
	for (int i = 0; i < Dimension::dim_N; i++)
	  {
	    Un(0) += coef*Vn(i+1)*normale(i);
	    Un(i+1) = Vn(0)*normale(i);
	  }

        Un *= nat_mat.GetCoefStiffness();
      }
    else if (cond == BoundaryConditionEnum::LINE_NEUMANN)
      {
	for (int i = 0; i < Dimension::dim_N; i++)
	  {
	    Un(0) -= coef*Vn(i+1)*normale(i);
	    Un(i+1) = -Vn(0)*normale(i);
	  }
        Un *= nat_mat.GetCoefStiffness();
      }
    else if (cond == BoundaryConditionEnum::LINE_ABSORBING)
      {
	// condition absorbante
	Un(0) = -Vn(0)*vars.coefficient_impedance_absorbing(ref_d);
	T vloc(0);
	for (int i = 0; i < Dimension::dim_N; i++)
	  vloc += Vn(i+1)*normale(i);
	
	vloc /= vars.coefficient_impedance_absorbing(ref_d);
	for (int i = 0; i < Dimension::dim_N; i++)
	  Un(i+1) = vloc*normale(i);
	
        Un *= nat_mat.GetCoefStiffness();
      }
    
    // DISP(Un); DISP(Vn);
  }    
  

  //! Evaluates the incident field
  template<class Dimension>
  void IncidentWaveProjector_LaplaceDG<Dimension>
  ::EvaluateFunction(int i, int j, const R_N& x, VectReal_wp& f)
  {
    if (this->wave_pulse == NULL)
      return;

    Real_wp arg = this->t + (this->offset_arg - DotProd(this->k_wave, x))/this->c0;
    Real_wp fs = this->wave_pulse->Evaluate(arg);
    f(0) = fs;    
  }


  //! gradient of incident field
  template<class Dimension>
  void IncidentWaveProjector_LaplaceDG<Dimension>
  ::EvaluateFunctionGradient(int i, int j, const R_N& x, VectReal_wp& f, VectReal_wp& df)
  {
    if (this->wave_pulse == NULL)
      return;

    Real_wp arg = this->t + (this->offset_arg - DotProd(this->k_wave, x))/this->c0;
    Real_wp fs = this->wave_pulse->Evaluate(arg);
    Real_wp dfs = -this->wave_pulse->EvaluateDerivative(arg)/this->c0;
    
    f(0) = fs;
    for (int k = 0; k < Dimension::dim_N; k++)
      df(k) = dfs*this->k_wave(k);
  }

  
  //! initialisation of source with time t0
  template<class Dimension>
  Real_wp DiffractedWaveSource_LaplaceDG<Dimension>
  ::Init(const Real_wp& t0, const Real_wp& dt, int print_level, int n, bool scal)
  {
#ifdef SELDON_WITH_MPI
    int rank_proc; MPI_Comm_rank(var_problem.comm_group_mode, &rank_proc);
#else
    int rank_proc(0);
#endif

    Real_wp tlimit(0);
    wave_pulse = incident_wave.GetTimeSource();
    
    if (var_source.GetNbParameterSource(0) <= 0)
      {
        abort();
      }
    
    t = t0;
    scalar_eq = scal;
        
    const Vector<string>& param = var_source.GetParameterSource(0, 0);
    c0 = var_boundary.GetMaximumVelocityPML();
    k_wave = var_problem.GetWaveVector();
    Mlt(1.0/Norm2(k_wave), k_wave);
    
    if (param(1) == "AUTO")
      {
        offset_arg = GetMinimumDotProd(k_wave, var_problem.GetXmin(), var_problem.GetXmax(),
                                       var_problem.GetYmin(), var_problem.GetYmax(),
                                       var_problem.GetZmin(), var_problem.GetZmax());

        Real_wp off_max = GetMaximumDotProd(k_wave, var_problem.GetXmin(), var_problem.GetXmax(),
					    var_problem.GetYmin(), var_problem.GetYmax(),
					    var_problem.GetZmin(), var_problem.GetZmax());        
	
	tlimit = off_max - offset_arg;        
      }
    else
      {
	offset_arg = to_num<Real_wp>(param(1));
	tlimit = to_num<Real_wp>(param(2));
      }
    
    if ((t0 <= dt+20*epsilon_machine) && (print_level >= 2) && (rank_proc == 0))
          cout << "Offset in expression of plane wave " << offset_arg << endl;    
    
    num_deriv = n;
    if (n > 1)
      {
        cout << "Higher-order derivatives not implemented" << endl;
        abort();
      }

    return tlimit;
  }
    
  
  //! initialisation of element i for volume integrals
  template<class Dimension>
  void DiffractedWaveSource_LaplaceDG<Dimension>
  ::InitElement(int num_elem, const VectR_N& s)
  {
  }
  
  
  //! returns true if there is a volume integral
  template<class Dimension>
  bool DiffractedWaveSource_LaplaceDG<Dimension>
  ::IsNonNullVolumetricSource(const VectR_N& s)
  {
    return false;
  }
    
  
  //! fills f for volume integral \int f \varphi
  template<class Dimension>
  void DiffractedWaveSource_LaplaceDG<Dimension>
  ::EvaluateVolumetricSource(int i, int j, const R_N& x, VectReal_wp& f)
  {
  }
  
  
  //! returns true if there a volume integral \int f \nabla \varphi
  template<class Dimension>
  bool DiffractedWaveSource_LaplaceDG<Dimension>
  ::IsNonNullGradientSource(const VectR_N& s)
  {
    return false;
  }
    
  
  //! fills f for volume integral \int f \nabla \varphi
  template<class Dimension>
  void DiffractedWaveSource_LaplaceDG<Dimension>
  ::EvaluateGradientSource(int i, int j, const R_N& x, VectReal_wp& f)
  {
  }
    
  
  //! initialisation of surface i
  template<class Dimension>
  void DiffractedWaveSource_LaplaceDG<Dimension>
  ::InitSurface(int i, int num_face, int num_elem, int num_loc)
  {
    VirtualSourceFEM<Real_wp, Dimension>::InitSurface(i, num_face, num_elem, num_loc);
  }
  
  
  //! fills f for surface integral \int f \varphi
  template<class Dimension>
  void DiffractedWaveSource_LaplaceDG<Dimension>
  ::EvaluateSurfacicSource(int k, const SetPoints<Dimension>& PointsElem,
                           const SetMatrices<Dimension>& MatricesElem, VectReal_wp& f)
  {
    if (wave_pulse == NULL)
      return;
    
    if (var_problem.InsidePML(this->num_elem_))
      {
        f.Fill(0);
        return;
      }

    // coef = k . n
    typename Dimension::R_N normale = MatricesElem.GetNormaleQuadratureBoundary(k);
    Real_wp coef = DotProd(k_wave, normale);
    
    // on evalue t - (k.x + offset_arg)/c0
    Real_wp arg = t + (offset_arg - DotProd(k_wave, PointsElem.GetPointQuadratureBoundary(k))) / this->c0;
    
    if (num_deriv != 0)
      {
	cout << "Not implemented" << endl;
        abort();
      }

    if (var_problem.FormulationDG() == ElementReference_Base::HDG)
      {
	Real_wp fs = wave_pulse->Evaluate(arg);
	if (reference_condition == BoundaryConditionEnum::LINE_DIRICHLET)
	  {
	    for (int i = 1; i <= Dimension::dim_N; i++)
	      f(i) = -normale(i-1)*fs;
	  }
	else
	  {
	    int ref = var_problem.mesh.Element(this->num_elem_).GetReference();
            coef *= var_laplace.ref_mu(ref).GetConstant()(0, 0)/this->c0;
	    f(0) = -coef*fs;
	  }

	return;
      }
    
    Real_wp fs = wave_pulse->Evaluate(arg);
    if (reference_condition == BoundaryConditionEnum::LINE_DIRICHLET)
      {
	if (!scalar_eq)
	  {
	    for (int i = 0; i < Dimension::dim_N; i++)
	      f(i) = -normale(i)*fs;
	  }
      }
    else
      {
	int ref = var_problem.mesh.Element(this->num_elem_).GetReference();
	coef *= var_laplace.ref_mu(ref).GetConstant()(0, 0)/c0;
	if (scalar_eq)
	  f(0) = coef*fs;
      }
  }
  
  
  //! returns true if there is a surface integral \int f \varphi
  template<class Dimension>
  bool DiffractedWaveSource_LaplaceDG<Dimension>
  ::IsNonNullSurfacicSource(int ref)
  {
    reference_condition = var_problem.mesh.GetBoundaryCondition(ref);
    if ((reference_condition == BoundaryConditionEnum::LINE_DIRICHLET)
        || (reference_condition == BoundaryConditionEnum::LINE_NEUMANN))
      return true;
    
    return false;
  }


  //! initialisation of source with time t0  
  template<class Dimension>
  Real_wp TotalWaveSource_LaplaceDG<Dimension>
  ::Init(const Real_wp& t0, const Real_wp& dt, int print_level, int n, bool scal)
  {
#ifdef SELDON_WITH_MPI
    int rank_proc; MPI_Comm_rank(var_problem.comm_group_mode, &rank_proc);
#else
    int rank_proc(0);
#endif

    Real_wp tlimit(0);
    wave_pulse = incident_wave.GetTimeSource();
    
    if (var_source.GetNbParameterSource(0) <= 0)
      {
        abort();
      }
    
    t = t0;
    scalar_eq = scal;
        
    const Vector<string>& param = var_source.GetParameterSource(0, 0);
    c0 = var_boundary.GetMaximumVelocityPML();
    k_wave = var_problem.GetWaveVector();
    Mlt(1.0/Norm2(k_wave), k_wave);
    
    if (param(1) == "AUTO")
      {
        offset_arg = GetMinimumDotProd(k_wave, var_problem.GetXmin(), var_problem.GetXmax(),
                                       var_problem.GetYmin(), var_problem.GetYmax(),
                                       var_problem.GetZmin(), var_problem.GetZmax());

        Real_wp off_max = GetMaximumDotProd(k_wave, var_problem.GetXmin(), var_problem.GetXmax(),
					    var_problem.GetYmin(), var_problem.GetYmax(),
					    var_problem.GetZmin(), var_problem.GetZmax());        
	
	tlimit = off_max - offset_arg;
      }
    else
      {
	offset_arg = to_num<Real_wp>(param(1));
	tlimit = to_num<Real_wp>(param(2));
      }
    
    if ((t0 <= dt+20*epsilon_machine) && (print_level >= 2) && (rank_proc == 0))
      cout << "Offset in expression of plane wave " << offset_arg << endl;
    
    num_deriv = n;
    if (n > 1)
      {
        cout << "Higher-order derivatives not implemented" << endl;
        abort();
      }

    return tlimit;
  }
  

  //! For Dirichlet condition
  template<class Dimension>
  void TotalWaveSource_LaplaceDG<Dimension>
  ::EvaluateFunction(int i, int j, const R_N& x, VectReal_wp& f)
  {
  }
  
  
  //! initialisation of element i for volume integrals
  template<class Dimension>
  void TotalWaveSource_LaplaceDG<Dimension>::InitElement(int num_elem, const VectR_N& s)
  {
  }
  
  
  //! returns true if there is a volume integral
  template<class Dimension>
  bool TotalWaveSource_LaplaceDG<Dimension>
  ::IsNonNullVolumetricSource(const VectR_N& s)
  {
    return false;
  }
    

  //! fills f for volume integral \int f \varphi  
  template<class Dimension>
  void TotalWaveSource_LaplaceDG<Dimension>
  ::EvaluateVolumetricSource(int i, int j, const R_N& x, VectReal_wp& f)
  {
  }
  
  
  //! returns true if there a volume integral \int f \nabla \varphi
  template<class Dimension>
  bool TotalWaveSource_LaplaceDG<Dimension>
  ::IsNonNullGradientSource(const VectR_N& s)
  {
    return false;
  }
    

  //! fills f for volume integral \int f \nabla \varphi  
  template<class Dimension>
  void TotalWaveSource_LaplaceDG<Dimension>
  ::EvaluateGradientSource(int i, int j, const R_N& x, VectReal_wp& f)
  {
  }
    
  
  //! initialisation of surface i
  template<class Dimension>
  void TotalWaveSource_LaplaceDG<Dimension>
  ::InitSurface(int i, int num_face, int num_elem, int num_loc)
  {
    VirtualSourceFEM<Real_wp, Dimension>::InitSurface(i, num_face, num_elem, num_loc);
  }
  
  
  //! fills f for surface integral \int f \varphi  
  template<class Dimension>
  void TotalWaveSource_LaplaceDG<Dimension>
  ::EvaluateSurfacicSource(int k, const SetPoints<Dimension>& PointsElem,
                           const SetMatrices<Dimension>& MatricesElem, VectReal_wp& f)
  {
    if (wave_pulse == NULL)
      return;
    
    // coef = 1 - k . n
    typename Dimension::R_N normale = MatricesElem.GetNormaleQuadratureBoundary(k);
    Real_wp coef = 1.0 - DotProd(k_wave, normale);
    
    // on evalue ct - k.x + offset_arg
    Real_wp arg = t + (offset_arg - DotProd(k_wave, PointsElem.GetPointQuadratureBoundary(k)))/c0;

    if (num_deriv != 0)
      {
	cout << "Not implemented" << endl;
        abort();
      }

    Real_wp fs = wave_pulse->Evaluate(arg);
    
    if (var_problem.FormulationDG() == ElementReference_Base::HDG)
      {
	int ref = var_problem.mesh.Element(this->num_elem_).GetReference();
	coef *= var_laplace.coefficient_impedance_absorbing(ref);
	f(0) = -coef*fs;
	
	return;
      }

    coef *= 0.5;
    if (scalar_eq)
      {
	int ref = var_problem.mesh.Element(this->num_elem_).GetReference();
	coef *= var_laplace.ref_mu(ref).GetConstant()(0, 0);
	f(0) = coef*fs;
      }
    else
      {            
	fs /= c0;
	for (int i = 0; i < Dimension::dim_N; i++)
	  f(i) = coef*fs*normale(i);
      }
  }
  
  
  //! returns true if there is a surface integral \int f \varphi
  template<class Dimension>
  bool TotalWaveSource_LaplaceDG<Dimension>::IsNonNullSurfacicSource(int ref)
  {
    if (var_problem.mesh.GetBoundaryCondition(ref) == BoundaryConditionEnum::LINE_ABSORBING)
      return true;
    
    return false;
  }

  
  //! impedance
  template<class Complexe, class Dimension>
  void ImpedanceABC_Laplace<Complexe, Dimension>
  ::EvaluateImpedancePhi(int i, int num_elem, int num_edge, int num_loc, int k,
                         const GlobalGenericMatrix<Complexe>& nat_mat, int ref,
                         const SetPoints<Dimension>& Pts, const SetMatrices<Dimension>& Mat)
  {
    this->read_param_condition = false;
    this->coef_grad = 0;
    if (var_helm.coefficient_impedance_absorbing.GetM() != 0)
      {
        // for acoustic equation
	this->coef_phi = var_helm.coefficient_impedance_absorbing(ref)*nat_mat.GetCoefDamping();
        if (var_problem.FormulationDG() == ElementReference_Base::HDG)
	  this->coef_phi = -var_helm.coefficient_impedance_absorbing(ref)*nat_mat.GetCoefStiffness();
      }
    else
      {
        abort();
      }
  }
  
}

#define MONTJOIE_FILE_VAR_LAPLACE_CXX
#endif
