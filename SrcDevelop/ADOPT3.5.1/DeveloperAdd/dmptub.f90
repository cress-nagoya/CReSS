module m_dmptub

contains

subroutine s_dmptub( ni, nj, nk, val, cf, sf, dt )

  use m_comtub

  implicit none

  integer, intent(in) :: ni
  integer, intent(in) :: nj
  integer, intent(in) :: nk
  real, intent(inout) :: val(0:ni+1,0:nj+1,1:nk)
  character(6), intent(in) :: cf
  character(1), intent(in) :: sf
  real, intent(in) :: dt

  integer :: i, j, k
  real :: dti

  dti=1.0e0/dt

!$omp parallel default(shared)

  if(cf(1:6)=='numdu ')then

    if(sf(1:1)=='o')then

!$omp do schedule(runtime) private(i,j,k)

       do k=1,nk-1
       do j=1,nj-1
       do i=1,ni-1
         numdu(i,j,k)=val(i,j,k)*dti
       end do
       end do
       end do

!$omp end do

     else if(sf(1:1)=='m')then

!$omp do schedule(runtime) private(i,j,k)

       do k=1,nk-1
       do j=1,nj-1
       do i=1,ni-1
         numdu(i,j,k)=val(i,j,k)-numdu(i,j,k)
         numdu(i,j,k)=numdu(i,j,k)*dti
       end do
       end do
       end do

!$omp end do

     else if(sf(1:1)=='p')then

!$omp do schedule(runtime) private(i,j,k)

       do k=1,nk-1
       do j=1,nj-1
       do i=1,ni-1
         numdu(i,j,k)=val(i,j,k)+numdu(i,j,k)
         numdu(i,j,k)=numdu(i,j,k)*dti
       end do
       end do
       end do

!$omp end do

     else if(sf(1:1)=='r')then

!$omp do schedule(runtime) private(i,j,k)

       do k=1,nk-1
       do j=1,nj-1
       do i=1,ni-1
         val(i,j,k)=numdu(i,j,k)
       end do
       end do
       end do

!$omp end do

     end if

  else if(cf(1:6)=='numdv ')then

    if(sf(1:1)=='o')then

!$omp do schedule(runtime) private(i,j,k)

       do k=1,nk-1
       do j=1,nj-1
       do i=1,ni-1
         numdv(i,j,k)=val(i,j,k)*dti
       end do
       end do
       end do

!$omp end do

     else if(sf(1:1)=='m')then

!$omp do schedule(runtime) private(i,j,k)

       do k=1,nk-1
       do j=1,nj-1
       do i=1,ni-1
         numdv(i,j,k)=val(i,j,k)-numdv(i,j,k)
         numdv(i,j,k)=numdv(i,j,k)*dti
       end do
       end do
       end do

!$omp end do

     else if(sf(1:1)=='p')then

!$omp do schedule(runtime) private(i,j,k)

       do k=1,nk-1
       do j=1,nj-1
       do i=1,ni-1
         numdv(i,j,k)=val(i,j,k)+numdv(i,j,k)
         numdv(i,j,k)=numdv(i,j,k)*dti
       end do
       end do
       end do

!$omp end do

     else if(sf(1:1)=='r')then

!$omp do schedule(runtime) private(i,j,k)

       do k=1,nk-1
       do j=1,nj-1
       do i=1,ni-1
         val(i,j,k)=numdv(i,j,k)
       end do
       end do
       end do

!$omp end do

     end if

  else if(cf(1:6)=='numdw ')then

    if(sf(1:1)=='o')then

!$omp do schedule(runtime) private(i,j,k)

       do k=1,nk-1
       do j=1,nj-1
       do i=1,ni-1
         numdw(i,j,k)=val(i,j,k)*dti
       end do
       end do
       end do

!$omp end do

     else if(sf(1:1)=='m')then

!$omp do schedule(runtime) private(i,j,k)

       do k=1,nk-1
       do j=1,nj-1
       do i=1,ni-1
         numdw(i,j,k)=val(i,j,k)-numdw(i,j,k)
         numdw(i,j,k)=numdw(i,j,k)*dti
       end do
       end do
       end do

!$omp end do

     else if(sf(1:1)=='p')then

!$omp do schedule(runtime) private(i,j,k)

       do k=1,nk-1
       do j=1,nj-1
       do i=1,ni-1
         numdw(i,j,k)=val(i,j,k)+numdw(i,j,k)
         numdw(i,j,k)=numdw(i,j,k)*dti
       end do
       end do
       end do

!$omp end do

     else if(sf(1:1)=='r')then

!$omp do schedule(runtime) private(i,j,k)

       do k=1,nk-1
       do j=1,nj-1
       do i=1,ni-1
         val(i,j,k)=numdw(i,j,k)
       end do
       end do
       end do

!$omp end do

     end if

  else if(cf(1:6)=='numdpt')then

    if(sf(1:1)=='o')then

!$omp do schedule(runtime) private(i,j,k)

       do k=1,nk-1
       do j=1,nj-1
       do i=1,ni-1
         numdpt(i,j,k)=val(i,j,k)*dti
       end do
       end do
       end do

!$omp end do

     else if(sf(1:1)=='m')then

!$omp do schedule(runtime) private(i,j,k)

       do k=1,nk-1
       do j=1,nj-1
       do i=1,ni-1
         numdpt(i,j,k)=val(i,j,k)-numdpt(i,j,k)
         numdpt(i,j,k)=numdpt(i,j,k)*dti
       end do
       end do
       end do

!$omp end do

     else if(sf(1:1)=='p')then

!$omp do schedule(runtime) private(i,j,k)

       do k=1,nk-1
       do j=1,nj-1
       do i=1,ni-1
         numdpt(i,j,k)=val(i,j,k)+numdpt(i,j,k)
         numdpt(i,j,k)=numdpt(i,j,k)*dti
       end do
       end do
       end do

!$omp end do

     else if(sf(1:1)=='r')then

!$omp do schedule(runtime) private(i,j,k)

       do k=1,nk-1
       do j=1,nj-1
       do i=1,ni-1
         val(i,j,k)=numdpt(i,j,k)
       end do
       end do
       end do

!$omp end do

     end if

  else if(cf(1:6)=='numdqv')then

    if(sf(1:1)=='o')then

!$omp do schedule(runtime) private(i,j,k)

       do k=1,nk-1
       do j=1,nj-1
       do i=1,ni-1
         numdqv(i,j,k)=val(i,j,k)*dti
       end do
       end do
       end do

!$omp end do

     else if(sf(1:1)=='m')then

!$omp do schedule(runtime) private(i,j,k)

       do k=1,nk-1
       do j=1,nj-1
       do i=1,ni-1
         numdqv(i,j,k)=val(i,j,k)-numdqv(i,j,k)
         numdqv(i,j,k)=numdqv(i,j,k)*dti
       end do
       end do
       end do

!$omp end do

     else if(sf(1:1)=='p')then

!$omp do schedule(runtime) private(i,j,k)

       do k=1,nk-1
       do j=1,nj-1
       do i=1,ni-1
         numdqv(i,j,k)=val(i,j,k)+numdqv(i,j,k)
         numdqv(i,j,k)=numdqv(i,j,k)*dti
       end do
       end do
       end do

!$omp end do

     else if(sf(1:1)=='r')then

!$omp do schedule(runtime) private(i,j,k)

       do k=1,nk-1
       do j=1,nj-1
       do i=1,ni-1
         val(i,j,k)=numdqv(i,j,k)
       end do
       end do
       end do

!$omp end do

     end if

  else if(cf(1:6)=='turbu ')then

    if(sf(1:1)=='o')then

!$omp do schedule(runtime) private(i,j,k)

       do k=1,nk-1
       do j=1,nj-1
       do i=1,ni-1
         turbu(i,j,k)=val(i,j,k)*dti
       end do
       end do
       end do

!$omp end do

     else if(sf(1:1)=='m')then

!$omp do schedule(runtime) private(i,j,k)

       do k=1,nk-1
       do j=1,nj-1
       do i=1,ni-1
         turbu(i,j,k)=val(i,j,k)-turbu(i,j,k)
         turbu(i,j,k)=turbu(i,j,k)*dti
       end do
       end do
       end do

!$omp end do

     else if(sf(1:1)=='p')then

!$omp do schedule(runtime) private(i,j,k)

       do k=1,nk-1
       do j=1,nj-1
       do i=1,ni-1
         turbu(i,j,k)=val(i,j,k)+turbu(i,j,k)
         turbu(i,j,k)=turbu(i,j,k)*dti
       end do
       end do
       end do

!$omp end do

     else if(sf(1:1)=='r')then

!$omp do schedule(runtime) private(i,j,k)

       do k=1,nk-1
       do j=1,nj-1
       do i=1,ni-1
         val(i,j,k)=turbu(i,j,k)
       end do
       end do
       end do

!$omp end do

     end if

  else if(cf(1:6)=='turbv ')then

    if(sf(1:1)=='o')then

!$omp do schedule(runtime) private(i,j,k)

       do k=1,nk-1
       do j=1,nj-1
       do i=1,ni-1
         turbv(i,j,k)=val(i,j,k)*dti
       end do
       end do
       end do

!$omp end do

     else if(sf(1:1)=='m')then

!$omp do schedule(runtime) private(i,j,k)

       do k=1,nk-1
       do j=1,nj-1
       do i=1,ni-1
         turbv(i,j,k)=val(i,j,k)-turbv(i,j,k)
         turbv(i,j,k)=turbv(i,j,k)*dti
       end do
       end do
       end do

!$omp end do

     else if(sf(1:1)=='p')then

!$omp do schedule(runtime) private(i,j,k)

       do k=1,nk-1
       do j=1,nj-1
       do i=1,ni-1
         turbv(i,j,k)=val(i,j,k)+turbv(i,j,k)
         turbv(i,j,k)=turbv(i,j,k)*dti
       end do
       end do
       end do

!$omp end do

     else if(sf(1:1)=='r')then

!$omp do schedule(runtime) private(i,j,k)

       do k=1,nk-1
       do j=1,nj-1
       do i=1,ni-1
         val(i,j,k)=turbv(i,j,k)
       end do
       end do
       end do

!$omp end do

     end if

  else if(cf(1:6)=='turbpt')then

    if(sf(1:1)=='o')then

!$omp do schedule(runtime) private(i,j,k)

       do k=1,nk-1
       do j=1,nj-1
       do i=1,ni-1
         turbpt(i,j,k)=val(i,j,k)*dti
       end do
       end do
       end do

!$omp end do

     else if(sf(1:1)=='m')then

!$omp do schedule(runtime) private(i,j,k)

       do k=1,nk-1
       do j=1,nj-1
       do i=1,ni-1
         turbpt(i,j,k)=val(i,j,k)-turbpt(i,j,k)
         turbpt(i,j,k)=turbpt(i,j,k)*dti
       end do
       end do
       end do

!$omp end do

     else if(sf(1:1)=='p')then

!$omp do schedule(runtime) private(i,j,k)

       do k=1,nk-1
       do j=1,nj-1
       do i=1,ni-1
         turbpt(i,j,k)=val(i,j,k)+turbpt(i,j,k)
         turbpt(i,j,k)=turbpt(i,j,k)*dti
       end do
       end do
       end do

!$omp end do

     else if(sf(1:1)=='r')then

!$omp do schedule(runtime) private(i,j,k)

       do k=1,nk-1
       do j=1,nj-1
       do i=1,ni-1
         val(i,j,k)=turbpt(i,j,k)
       end do
       end do
       end do

!$omp end do

     end if

  else if(cf(1:6)=='turbqv')then

    if(sf(1:1)=='o')then

!$omp do schedule(runtime) private(i,j,k)

       do k=1,nk-1
       do j=1,nj-1
       do i=1,ni-1
         turbqv(i,j,k)=val(i,j,k)*dti
       end do
       end do
       end do

!$omp end do

     else if(sf(1:1)=='m')then

!$omp do schedule(runtime) private(i,j,k)

       do k=1,nk-1
       do j=1,nj-1
       do i=1,ni-1
         turbqv(i,j,k)=val(i,j,k)-turbqv(i,j,k)
         turbqv(i,j,k)=turbqv(i,j,k)*dti
       end do
       end do
       end do

!$omp end do

     else if(sf(1:1)=='p')then

!$omp do schedule(runtime) private(i,j,k)

       do k=1,nk-1
       do j=1,nj-1
       do i=1,ni-1
         turbqv(i,j,k)=val(i,j,k)+turbqv(i,j,k)
         turbqv(i,j,k)=turbqv(i,j,k)*dti
       end do
       end do
       end do

!$omp end do

     else if(sf(1:1)=='r')then

!$omp do schedule(runtime) private(i,j,k)

       do k=1,nk-1
       do j=1,nj-1
       do i=1,ni-1
         val(i,j,k)=turbqv(i,j,k)
       end do
       end do
       end do

!$omp end do

     end if

  else if(cf(1:6)=='vish  ')then

!$omp do schedule(runtime) private(i,j,k)

     do k=1,nk-1
     do j=1,nj-1
     do i=1,ni-1
       vish(i,j,k)=val(i,j,k)*dti
     end do
     end do
     end do

!$omp end do

  else if(cf(1:6)=='visv  ')then

!$omp do schedule(runtime) private(i,j,k)

     do k=1,nk-1
     do j=1,nj-1
     do i=1,ni-1
       visv(i,j,k)=val(i,j,k)*dti
     end do
     end do
     end do

!$omp end do

  else if(cf(1:6)=='difh  ')then

!$omp do schedule(runtime) private(i,j,k)

     do k=1,nk-1
     do j=1,nj-1
     do i=1,ni-1
       difh(i,j,k)=val(i,j,k)*dti
     end do
     end do
     end do

!$omp end do

  else if(cf(1:6)=='difv  ')then

!$omp do schedule(runtime) private(i,j,k)

     do k=1,nk-1
     do j=1,nj-1
     do i=1,ni-1
       difv(i,j,k)=val(i,j,k)*dti
     end do
     end do
     end do

!$omp end do

  end if

!$omp end parallel

end subroutine s_dmptub

end module m_dmptub
