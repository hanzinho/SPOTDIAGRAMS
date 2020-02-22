subroutine read_design(c,ic,nc)
  use spotmod
  implicit none
  integer           :: nc
  type (comp)       :: c(nc)

  character (len=500) :: dfile,cline
  integer             :: iu=12,i,j
  integer             :: n,ic

!-Substrings, allow upto kn of them
  integer,parameter   :: kn=100
  integer             :: ks(kn),is
  character (len=100) :: ch(kn),cc

!-Find design file from program arguments
  call get_option(dfile,"-i")
  write(*,'(2a)')"Open design file: ",trim(dfile)

!-Scan the design file
  ic=-1
  open(iu,file=trim(dfile),form='formatted')
  do
    read(iu,'(A200)',end=9)cline

    i=scan(cline,'#')  ! strip comment and empty lines
      if (i/=0) cline=trim(cline(:i-1))
      if (len_trim(cline)==0) cycle
    write(*,'(a)')trim(cline)
    ic=ic+1

    call split_string(cline,ks,kn,is)
    if (ic==0) then
      do i=1,is
         ch(i)=trim(cline(ks(i):ks(i+1)-1))
         call to_lower(ch(i))
      enddo
    else
      do i=1,is
         cc=trim(cline(ks(i):ks(i+1)-1))
         call to_lower(cc)
         if (index(cc,'-',.true.)==len_trim(cc)) cycle
         if (ch(i)=="type"    ) c(ic)%type=cc(:1)
         if (ch(i)=="surface" ) c(ic)%surf=cc(:1)
         if (ch(i)=="medium"  ) c(ic)%medium=cc
         if (ch(i)=="aperture") read(cc,*)c(ic)%a
         if (ch(i)=="radius"  ) read(cc,*)c(ic)%r
         if (ch(i)=="scp"     ) read(cc,*)c(ic)%scp
         if (ch(i)=="e"       ) read(cc,*)c(ic)%scp
         if (ch(i)=="e"       ) c(ic)%scp=1.-c(ic)%scp**2
         if (ch(i)=="distance") read(cc,*)c(ic)%d
      enddo
    endif
  enddo
9 close(iu)
  write(*,*)'Number of components: ',ic

!-Add focal plane and post-process components
  ic=ic+1
  c(ic)%type='v'
  c(ic)%surf='f'
  
  do i=ic,1,-1
     if(i>1 )c(i)%d=c(i-1)%d
     if(i==1)c(i)%d=0.
     call find_glass_code(c(i))
  enddo

end subroutine read_design

!------------------------------------------------------------------

subroutine find_glass_code(c)
  use spotmod
  implicit none
  type (comp) :: c

  character*20 :: m
  real         :: xld=587.56        ! wave length d
  real         :: fmc=486.13-656.27 ! wave length m - c
  real         :: xnd,dn
  logical      :: lcode

  m=trim(c%medium) ; call to_lower(m)
  if (lcode(m)) then
     read(m(1:3),*)xnd ; xnd=0.001*xnd ! (n - 1)
     read(m(4:6),*)dn  ; dn =0.100*dn  ! Abbe number (n-1)/(nf-nc)
     dn=xnd/dn
     xnd=1.0+xnd
     c%ri%xn(:)=xnd + dn*(c%ri%xl(:)-xld)/fmc
  else
     if (m=='air') c%ri%xn(:)=1.0
  endif
end subroutine find_glass_code

!------------------------------------------------------------------
subroutine get_option(cstr,o)
  implicit none
  character(len=*) :: cstr,o 

  integer :: i,j
  j=-1
  do i=1, iargc()
     call getarg(i,cstr)
     if(trim(cstr)==trim(o)) j=i
     if (i==j+1) exit
  enddo
  if (j==-1) cstr=""

end subroutine get_option
!------------------------------------------------------------------

subroutine to_lower(cstring)
  implicit none
  character(len=*) :: cstring
  integer :: i,j
  integer :: iu=iachar('A'), il=iachar('a')

  do i=1,len(cstring)
     j=iachar(cstring(i:i))-iu
     if (j>=0 .and. j<=26) cstring(i:i)=achar(il+j)
  enddo

end subroutine to_lower
!------------------------------------------------------------------

function lcode(cstring)
  implicit none
  logical lcode
  character(len=*) :: cstring

  integer :: i,j,i0=iachar('0')

  lcode=(len_trim(cstring)==6)
  do i=1,6
     j=iachar(cstring(i:i))-i0
     if (j<0 .or. j>9)exit 
  enddo
  if (i<=6) lcode=.false.

end function lcode
!------------------------------------------------------------------

subroutine split_string(cline,kstart,n,k)
  implicit none
  character (len=*) :: cline
  integer :: n, kstart(n),k

  integer :: i,j

  j=0; k=0
  do i=1,len_trim(cline)
     if (cline(i:i)==' ') then
        j=i
     else
        if (i==j+1) then
           k=k+1
           kstart(k)=i
        endif
     endif
  enddo
  kstart(k+1)=len_trim(cline)+1

end subroutine split_string
