module m_dmpcld

contains

subroutine s_dmpcld( ni, nj, nk, val, cf, sf, dt )

  use m_comblk

  implicit none

  integer, intent(in) :: ni
  integer, intent(in) :: nj
  integer, intent(in) :: nk
  real, intent(inout) :: val(0:ni+1,0:nj+1,1:nk)
  character(5), intent(in) :: cf
  character(1), intent(in) :: sf
  real, intent(in) :: dt

  integer :: i, j, k
  real :: dti

  dti=1.0e0/dt

!$omp parallel default(shared)

  if(cf(1:5)=='dqadt')then

    if(sf(1:1)=='o')then

!$omp do schedule(runtime) private(i,j,k)

       do k=1,nk-1
       do j=1,nj-1
       do i=1,ni-1
         dqadt(i,j,k)=val(i,j,k)*dti
       end do
       end do
       end do

!$omp end do

     else if(sf(1:1)=='m')then

!$omp do schedule(runtime) private(i,j,k)

       do k=1,nk-1
       do j=1,nj-1
       do i=1,ni-1
         dqadt(i,j,k)=val(i,j,k)-dqadt(i,j,k)
         dqadt(i,j,k)=dqadt(i,j,k)*dti
       end do
       end do
       end do

!$omp end do

     else if(sf(1:1)=='p')then

!$omp do schedule(runtime) private(i,j,k)

       do k=1,nk-1
       do j=1,nj-1
       do i=1,ni-1
         dqadt(i,j,k)=val(i,j,k)+dqadt(i,j,k)
         dqadt(i,j,k)=dqadt(i,j,k)*dti
       end do
       end do
       end do

!$omp end do

     else if(sf(1:1)=='r')then

!$omp do schedule(runtime) private(i,j,k)

       do k=1,nk-1
       do j=1,nj-1
       do i=1,ni-1
         val(i,j,k)=dqadt(i,j,k)
       end do
       end do
       end do

!$omp end do

     end if

  else if(cf(1:4)=='dqdt') then

    if(sf(1:1)=='o')then

!$omp do schedule(runtime) private(i,j,k)

      do k=1,nk-1
      do j=1,nj-1
      do i=1,ni-1
        dqdt(i,j,k)=val(i,j,k)*dti
      end do
      end do
      end do

!$omp end do

    else if(sf(1:1)=='m')then

!$omp do schedule(runtime) private(i,j,k)

      do k=1,nk-1
      do j=1,nj-1
      do i=1,ni-1
        dqdt(i,j,k)=val(i,j,k)-dqdt(i,j,k)
        dqdt(i,j,k)=dqdt(i,j,k)*dti
      end do
      end do
      end do

!$omp end do

     else if(sf(1:1)=='p')then

!$omp do schedule(runtime) private(i,j,k)

       do k=1,nk-1
       do j=1,nj-1
       do i=1,ni-1
         dqdt(i,j,k)=val(i,j,k)+dqdt(i,j,k)
         dqdt(i,j,k)=dqdt(i,j,k)*dti
       end do
       end do
       end do

!$omp end do

     else if(sf(1:1)=='r')then

!$omp do schedule(runtime) private(i,j,k)

       do k=1,nk-1
       do j=1,nj-1
       do i=1,ni-1
         val(i,j,k)=dqdt(i,j,k)
       end do
       end do
       end do

!$omp end do

    end if

  else if(cf(1:5)=='vdvc ') then

    if(sf(1:1)=='o')then

!$omp do schedule(runtime) private(i,j,k)

      do k=1,nk-1
      do j=1,nj-1
      do i=1,ni-1
        vdvc(i,j,k)=val(i,j,k)*dti
      end do
      end do
      end do

!$omp end do

    else if(sf(1:1)=='m')then

!$omp do schedule(runtime) private(i,j,k)

      do k=1,nk-1
      do j=1,nj-1
      do i=1,ni-1
        vdvc(i,j,k)=val(i,j,k)-vdvc(i,j,k)
        vdvc(i,j,k)=vdvc(i,j,k)*dti
      end do
      end do
      end do

!$omp end do

    else if(sf(1:1)=='r')then

!$omp do schedule(runtime) private(i,j,k)

       do k=1,nk-1
       do j=1,nj-1
       do i=1,ni-1
         val(i,j,k)=vdvc(i,j,k)
       end do
       end do
       end do

!$omp end do

    end if

  else if(cf(1:5)=='vdvc2') then

    if(sf(1:1)=='o')then

!$omp do schedule(runtime) private(i,j,k)

      do k=1,nk-1
      do j=1,nj-1
      do i=1,ni-1
        vdvc2(i,j,k)=val(i,j,k)*dti
      end do
      end do
      end do

!$omp end do

    else if(sf(1:1)=='m')then

!$omp do schedule(runtime) private(i,j,k)

      do k=1,nk-1
      do j=1,nj-1
      do i=1,ni-1
        vdvc2(i,j,k)=val(i,j,k)-vdvc2(i,j,k)
        vdvc2(i,j,k)=vdvc2(i,j,k)*dti
      end do
      end do
      end do

!$omp end do

    else if(sf(1:1)=='r')then

!$omp do schedule(runtime) private(i,j,k)

       do k=1,nk-1
       do j=1,nj-1
       do i=1,ni-1
         val(i,j,k)=vdvc2(i,j,k)
       end do
       end do
       end do

!$omp end do

    end if

  end if

!$omp end parallel

end subroutine s_dmpcld

end module m_dmpcld
