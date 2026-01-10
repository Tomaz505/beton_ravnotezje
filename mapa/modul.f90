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



end module
