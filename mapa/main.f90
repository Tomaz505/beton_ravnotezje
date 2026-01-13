program main
    use rutine


    integer, parameter :: dp = selected_real_kind(15,307)
    real(dp), parameter :: pi = 3.141592654
    integer :: n_c, n_a, n_p, u, analiza  !st vozlisc, stevilo armaturnih palic, stevilo kablov

    real(dp) :: f_c, e_c, phi_cr,x_area
    real(dp) :: f_s, e_s
    real(dp) :: f_p, e_p
    real(dp) ::  a=0.0_dp, ix=0.0_dp
    real(dp) :: med, ned

    real(dp), allocatable :: xy_c(:,:), xy_s(:,:), xy_p(:,:), r_s(:), r_p(:), p_p(:)
    real(dp) :: def_pl(3), eps_sh(3)

    open(newunit = u, file = "in.txt",status = "old")
    read(u,*) analiza, med,ned,n_c,n_s,n_p
    close(u)

    if (np == 0) then
        allocate(xy_c(2,n_c), xy_s(2,n_s), r_s(n_s))
        open(newunit = u, file = "in.txt",status = "old")
            read(u,*) analiza,med,ned,n_c,n_s,n_p,f_c,e_c,phi_cr,f_s,e_s,xy_c,eps_sh,xy_s,r_s
        close(u)
    else
        allocate(xy_c(2,n_c), xy_s(2,n_s),xy_p(2,n_p), r_s(n_s), r_p(n_p),p_p(n_p))
        open(newunit = u, file = "in.txt",status = "old")
            read(u,*) analiza, med,ned,n_c,n_s,n_p,f_c,e_c,phi_cr,f_s,e_s,xy_c,eps_sh,xy_s,r_s,f_p,e_p,xy_p,r_p,p_p
        close(u)
    end if

    !Skladnost koordinatnega sistema
    eps_sh(2) = -eps_sh(2)




    !RAČUN KARAKTERISTIK IN PREMIK TEZISCA V IZODISCE
    block
        real(dp) :: sx=0.0_dp

        call area_n(xy_c,n_c,a,0)
        call area_n(xy_c,n_c,sx,1)

        xy_c(2,:) = xy_c(2,:)-sx/a
        xy_s(2,:) = xy_s(2,:)-sx/a

        if (n_p /= 0) then
            xy_p(2,:) = xy_p(2,:) - sx/a
        end if

        call area_n(xy_c,n_c,ix,2)

    end block



    goto 10



    !PRINTANJE REZULTATOV
100 block
        real(dp) ::  def_plc(3), z_extr(2)
        def_plc(:) = 0
        z_extr(1) = minval(xy_c(2,:))
        z_extr(2) = maxval(xy_c(2,:))

        def_plc = (def_pl-eps_sh)/(1+phi_cr)

        call write_js(xy_c,n_c,xy_s,n_s,r_s,def_plc, def_pl, z_extr,analiza,e_c,f_c,e_s,f_s)

        goto 110
    end block !-> 110




    !RAČUN NERAZPOKANEGA PREREZA
 10 block
        real(dp) :: c_mat(2,2), f_vec(2)
        f_vec = (/ned, med/)
        c_mat(:,:) = 0

        do i =1,n_s
            c_mat(1,2) = c_mat(1,2) - 2*e_s*xy_s(2,i)*pi/4 *r_s(i)**2
            c_mat(1,1) = c_mat(1,1) + 2*e_s*pi/4 *r_s(i)**2
            c_mat(2,2) = c_mat(2,2) + 2*e_s* xy_s(2,i)**2 *pi/4 *r_s(i)**2
        end do

        c_mat(1,1) = c_mat(1,1) +  a*e_c/(1+phi_cr)
        c_mat(2,2) = c_mat(2,2) + ix*e_c/(1+phi_cr)
        c_mat(2,1) = c_mat(1,2)

        f_vec(1) = f_vec(1) +  (eps_sh(1)*a+eps_sh(3)*ix)*e_c/(1+phi_cr)
        f_vec(2) = f_vec(2) + eps_sh(3)*ix*e_c/(1+phi_cr)


        def_pl(1:2) = (/c_mat(2,2)*f_vec(1)-c_mat(1,2)*f_vec(2) , -c_mat(1,2)*f_vec(1)+c_mat(1,1)*f_vec(2)/) /(c_mat(1,1)*c_mat(2,2)-c_mat(1,2)**2)
        def_pl(3) = 0

        goto (100,20,30), analiza !write_js(xy_c,n_c,xy_s,n_s,r_s,def_plc, def_pl, z_extr,analiza)
    end block !-> 100




    !RACUN RAZPOKANEGA LINEARNEGA MATERIALA
20  block
        real(dp) :: xy_eff(2,2*n_c), z_crack(2), def_pl_crac(3), z_extr(2), jac_eq(2,2), f_eq(2),z_crac(2),dlta_eps(2)
        real(dp) :: i0,i1,i2,i3,eps_h,kapa_h

        !korak diference
        eps_h  = 2.0_dp**(-25.0_dp)
        kapa_h = 2.0_dp**(-25.0_dp)

        z_extr(1) = minval(xy_c(2,:))
        z_extr(2) = maxval(xy_c(2,:))

        def_pl_crac = def_pl
        if (eps_sh(3) == 0) then

            do k1 = 1,10

            f_eq(:) = 0.0_dp
            jac_eq(:,:) = 0.0_dp


            call crac_linmat_linsh(xy_c,n_c,xy_s,n_s,def_pl_crac,z_extr,eps_sh,e_s,e_c,f_eq,r_s,phi_cr)
            call crac_linmat_linsh(xy_c,n_c,xy_s,n_s,def_pl_crac+(/eps_h,0.0_dp,0.0_dp/),z_extr,eps_sh,e_s,e_c,jac_eq(:,1),r_s,phi_cr)
            call crac_linmat_linsh(xy_c,n_c,xy_s,n_s,def_pl_crac+(/0.0_dp,kapa_h,0.0_dp/),z_extr,eps_sh,e_s,e_c,jac_eq(:,2),r_s,phi_cr)


            jac_eq(:,1) = (jac_eq(:,1)-f_eq)/(eps_h)
            jac_eq(:,2) = (jac_eq(:,2)-f_eq)/(kapa_h)

            f_eq = f_eq-(/ned,med/)

            dlta_eps(1) = (f_eq(1)*jac_eq(2,2) - f_eq(2)*jac_eq(1,2))
            dlta_eps(2) = (-f_eq(1)*jac_eq(2,1) + f_eq(2)*jac_eq(1,1))
            dlta_eps = dlta_eps/(jac_eq(1,1)*jac_eq(2,2) - jac_eq(1,2)*jac_eq(2,1))


            def_pl_crac(1:2) = def_pl_crac(1:2) - dlta_eps


            end do
            def_pl = def_pl_crac

        else


!             do k1 = 1,10
!
!             f_eq(:) = 0.0_dp
!             jac_eq(:,:) = 0.0_dp
!
!
!             call crac_linmat_nlsh(xy_c,n_c,xy_s,n_s,def_pl_crac,z_extr,eps_sh,e_s,e_c,f_eq,r_s,phi_cr)
!             call crac_linmat_nlsh(xy_c,n_c,xy_s,n_s,def_pl_crac+(/eps_h,0.0_dp,0.0_dp/),z_extr,eps_sh,e_s,e_c,jac_eq(:,1),r_s,phi_cr)
!             call crac_linmat_nlsh(xy_c,n_c,xy_s,n_s,def_pl_crac+(/0.0_dp,kapa_h,0.0_dp/),z_extr,eps_sh,e_s,e_c,jac_eq(:,2),r_s,phi_cr)
!
!
!             jac_eq(:,1) = (jac_eq(:,1)-f_eq)/(eps_h)
!             jac_eq(:,2) = (jac_eq(:,2)-f_eq)/(kapa_h)
!
!             f_eq = f_eq-(/ned,med/)
!
!             dlta_eps(1) = (f_eq(1)*jac_eq(2,2) - f_eq(2)*jac_eq(1,2))
!             dlta_eps(2) = (-f_eq(1)*jac_eq(2,1) + f_eq(2)*jac_eq(1,1))
!             dlta_eps = dlta_eps/(jac_eq(1,1)*jac_eq(2,2) - jac_eq(1,2)*jac_eq(2,1))
!
!
!             def_pl_crac(1:2) = def_pl_crac(1:2) - dlta_eps
!
!
!             end do
!             def_pl = def_pl_crac
!
!
         end if
!                                                            |
        goto (100,100,30),analiza  !write_js(xy_c,n_c,xy_s,n_s,r_s,def_plc, def_pl, z_extr,analiza)
    end block !-> 100



    !RACUN RAZPOKANEGA NELINEARNEGA MATERIALA
30  block
        real(dp) :: xy_eff(2,2*n_c), z_crack(2), def_pl_crac(3), z_extr(2), jac_eq(2,2), f_eq(2),z_crac(2),dlta_eps(2)
        real(dp) :: i0,i1,i2,i3,eps_h,kapa_h

        !korak diference
        eps_h  = 2.0_dp**(-25.0_dp)
        kapa_h = 2.0_dp**(-25.0_dp)

        z_extr(1) = minval(xy_c(2,:))
        z_extr(2) = maxval(xy_c(2,:))
        def_pl_crac = def_pl
        !if (eps_sh(3) == 0) then

            do k1 = 1,1

            f_eq(:) = 0.0_dp
            jac_eq(:,:) = 0.0_dp

            call crac_nlmat_linsh(xy_c,n_c,xy_s,n_s,def_pl_crac,z_extr,eps_sh,f_s,e_s,f_c,f_eq,r_s)
            call crac_nlmat_linsh(xy_c,n_c,xy_s,n_s,def_pl_crac+(/eps_h,0.0_dp,0.0_dp/),z_extr,eps_sh,f_s,e_s,f_c,jac_eq(:,1),r_s)
            call crac_nlmat_linsh(xy_c,n_c,xy_s,n_s,def_pl_crac+(/0.0_dp,kapa_h,0.0_dp/),z_extr,eps_sh,f_s,e_s,f_c,jac_eq(:,2),r_s)


            jac_eq(:,1) = (jac_eq(:,1)-f_eq)/(eps_h)
            jac_eq(:,2) = (jac_eq(:,2)-f_eq)/(kapa_h)

            f_eq = f_eq-(/ned,med/)

            dlta_eps(1) = (f_eq(1)*jac_eq(2,2) - f_eq(2)*jac_eq(1,2))
            dlta_eps(2) = (-f_eq(1)*jac_eq(2,1) + f_eq(2)*jac_eq(1,1))
            dlta_eps = dlta_eps/(jac_eq(1,1)*jac_eq(2,2) - jac_eq(1,2)*jac_eq(2,1))


            def_pl_crac(1:2) = def_pl_crac(1:2) - dlta_eps


            end do
            def_pl = def_pl_crac

            print *, def_pl
        !end if
        goto 100
    end block


    !INTERAKCIJSKI DIAGRAM
40  block
    end block


110     print*," "
        print*,"    Račun končan."
        print*," "








end program
