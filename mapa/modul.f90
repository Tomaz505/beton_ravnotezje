module rutine
    private dp
    public area_n
    integer, parameter :: dp = selected_real_kind(15,307)

    contains

    subroutine area_n(xy,n,a,deg)
        !n-stevilo vozlisc
        !deg-stopnja potence y^deg v integralu
        !a - vhodna količina
        !xy - koordinate vozlisc
        integer :: n,deg
        real(dp) :: a, xy(:,:)

        do i=1,n-1
            x = 0.0
            do j= 0,deg
                x = x + xy(2,i)**(deg-j)*xy(2,i+1)**(j)
            end do
            a = a + x*(xy(1,i)*xy(2,i+1)-xy(1,i+1)*xy(2,i))*2/(1+deg)/(2+deg)
        end do

        return
    end subroutine


    subroutine write_js(xy,n,xys,ns,ds,def_c,def_s,h)
        !n-stevilo vozlisc
        !deg-stopnja potence y^deg v integralu
        !a - vhodna količina
        !xy - koordinate vozlisc
        integer :: n,ns
        real(dp) :: xy(:,:),def_c(:),def_s(:),h(:),xys(:,:),ds(:)

        open (unit = 1,file = "mapa/out.js",status = "old")
        write(1,*) "const sec_coor = `"
        do i=1,n
            write(1,*) xy(:,i)
        end do
        do i=n-1,1,-1
            write(1,*) (/ -xy(1,i),xy(2,i)/)
        end do
        write(1,*) "`;"

        write(1,*) "const a_coor = `"
        do i=1,ns
            write(1,*) xys(:,i) , ds(i)
        end do
        do i=ns,1,-1
            write(1,*) (/ -xys(1,i),xys(2,i),  ds(i)/)
        end do
        write(1,*) "`;"

        write(1,*) "const eps_c = `"
        do i=1,2
            write(1,*) def_c(i),h(i)
        end do
        write(1,*) "`;"

        write(1,*) "const eps_s = `"
        do i=1,2
            write(1,*) def_s(i),h(i)
        end do
        write(1,*) "`;"


        close(1)
    end subroutine



end module
