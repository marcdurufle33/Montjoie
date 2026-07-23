subroutine mlt_add_test(alpha, level, X, beta, B, nodl, nb_elt, nb_dof_elt, &
     Nodle, PtrRh, IndRh, ValRh, PtrRht, IndRht, ValRht, dfjm1, invJacobian)
  
  implicit none
  
  real(kind=8) :: alpha, beta
  integer :: level
  real(kind=8) :: X(*), B(*), ValRh(*), ValRht(*)
  integer :: nodl, nb_elt, nb_dof_elt
  real(kind = 8) :: dfjm1(*), invJacobian(nb_elt,*)
  integer :: Nodle(nb_elt, *)
  integer :: PtrRh(*), IndRh(*), PtrRht(*), IndRht(*)
  
  integer :: i, j, nb_coef, k, num_dof, nnz
  real(kind=8) :: Uh(nb_dof_elt), Vh(2*nb_dof_elt), Prod_Uh(nb_dof_elt)
  real(kind=8) :: tmp(2), vh_loc(2), val, mass
  
  do i = 1, nodl
     B(i) = B(i)*beta
  end do
  
  nb_coef = 4*nb_dof_elt
  !print*, "X", X(1:nodl)
  !print*, "Ptr", PtrRh(1:2*nb_dof_elt+1)
  nnz = PtrRh(2*nb_dof_elt+1)
  !print*, "Ind", IndRh(1:nnz)
  !print*, "Val", ValRh(1:nnz)
  do i = 1, nb_elt
     
     do j = 1, nb_dof_elt
        num_dof = Nodle(i, j)+1
        Uh(j) = X(num_dof)
     end do
     !print*, "i", i
     !print*, "Uh", Uh
     
     do j = 1, 2*nb_dof_elt
        val = 0.0d0
        do k = PtrRh(j)+1, PtrRh(j+1)
           val = val + ValRh(k)*Uh(IndRh(k)+1)
        end do
        Vh(j) = val
     end do
     !print*, "Vh", Vh
     
     ! multiplication by J_i / omega_k DF_i^{-1} \mu DF_i^{*-1}
     do j = 1, nb_dof_elt
        
        tmp(1) = Vh(2*j-1)
        tmp(2) = Vh(2*j)
        
        vh_loc(1) = dfjm1((i-1)*nb_coef+4*j-3)*tmp(1) + dfjm1((i-1)*nb_coef+4*j-1)*tmp(2)
        vh_loc(2) = dfjm1((i-1)*nb_coef+4*j-2)*tmp(1) + dfjm1((i-1)*nb_coef+4*j)*tmp(2)
        
        mass = invJacobian(i, j)
        vh_loc(1) = vh_loc(1)*mass
        vh_loc(2) = vh_loc(2)*mass
        
        tmp(1) = dfjm1((i-1)*nb_coef+4*j-3)*vh_loc(1) + dfjm1((i-1)*nb_coef+4*j-2)*vh_loc(2)
        tmp(2) = dfjm1((i-1)*nb_coef+4*j-1)*vh_loc(1) + dfjm1((i-1)*nb_coef+4*j)*vh_loc(2)
	  
        Vh(2*j-1) = tmp(1)
        Vh(2*j) = tmp(2)
        
     end do
     !print*, "Vh", Vh
     
     do j = 1, nb_dof_elt
        val = 0.0d0
        do k = PtrRht(j)+1, PtrRht(j+1)
           val = val + ValRht(k)*Vh(IndRht(k)+1)
        end do
        Prod_Uh(j) = val
     end do
     !print*, "Prod_Uh", Prod_Uh
     
     do j = 1, nb_dof_elt
        num_dof = Nodle(i, j)+1
        B(num_dof) = B(num_dof) - alpha*Prod_Uh(j);
     end do
  end do
  !print*, "B", B(1:nodl)
  
end subroutine mlt_add_test

subroutine mlt_add_sparse(m, n, Ptr, Ind, Val, X, Y)
  
  integer :: m, n, Ptr(*), Ind(*)
  real(8) :: Val(*), X(*), Y(*), eval
  
  do j = 1, m
     eval = 0.0d0
     do k = Ptr(j)+1, Ptr(j+1)
        eval = eval + Val(k)*X(Ind(k)+1)
     end do
     Y(j) = eval
  end do
  
end subroutine mlt_add_sparse
