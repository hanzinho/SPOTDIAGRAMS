program spot
  use spotmod
  implicit none

  integer :: nb=50000,ib,nc,i
  logical :: lok

  integer, parameter  :: nm=100 ! allow up to 100 components
  type (comp)         :: c(nm)
  type (lray),pointer :: r,rr(:) => NULL()
  type (bundle)       :: b
  character (len=500) :: fout,cd
  character           :: colour='g'

  call read_design(c,nc,nm)

!-Open output file 
  call get_option(fout,"-o")
  write(*,'(2a)')"Open output file: ",trim(fout)
  open(8,file=trim(fout),form='formatted')

  call get_option(cd,"-c") ; if(cd/="")read(cd,'(a)')colour
  call get_option(cd,"-n") ; if(cd/="")read(cd,*)nb
  allocate(rr(nb))

  call set_color(c,nc,colour)
  call ini_bundle(b,c(1)%a,nb,0.0078125)
! call ini_bundle(b,c(1)%a,nb,0.01)
! call ini_bundle(b,c(1)%a,nb,0.0)
  ib=1
  do 
    r => rr(ib)
    call shoot_ray(r,b,lok)  ; if (.not.lok) exit
    do i=1,nc-1
       call cut(r,c(i),lok)  ; if (.not.lok) exit
       call ref(r,c(i),lok)  ; if (.not.lok) exit
    enddo
    if (lok) ib=ib+1
  enddo
  nb=ib-1
  write(*,*)'nb: ',nb

! Find location in focal plane
  call find_focus(rr,nb,c(nc))
  do ib=1,nb
     r => rr(ib)
     call cut(r,c(nc),lok)
     if (lok) write(8,'(2f12.6)')r%b(:2)
  enddo
     

!-Clean up
  deallocate (rr)

end program spot
