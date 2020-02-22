module spotmod
  implicit none

  type lray
    real :: a(3), b(3)
  end type lray

  type rindex  ! indices of refraction
    character(len=5) :: cl="rdgbv"
    real :: xn(5)
    real :: xl(5) = (/656.27,&  ! r (red   ); C
                    & 587.56,&  ! d (yellow): d
                    & 546.07,&  ! g (green ): e
                    & 486.13,&  ! b (blue  ); F
                    & 404.66/)  ! v (violet); h
  end type rindex

  type comp
    real      :: n(3), s(3)
    character :: type     ! 'm' (mirror) 
                          ! 'l' (lens) 
                          ! 'v' (vignet) 
    character :: surf     ! 'c' (cone) 
                          ! 'f' (flat) 

    real      :: np2m     ! (n2/n1)^2 - 1, ratio breaking indices
    real      :: a = 1.E20! aperture
    real      :: r = 1.E20! radius
    real      :: d = 0.   ! position on optical axis
    real      :: scp=1.   ! Swarzschild Constant + 1 = 1-e^2

    character(len=20) :: medium = "air"
    character(len= 6) :: glcode
    type (rindex)     :: ri
  end type comp

  type bundle
    character :: type ! 'r' (random)
                      ! 'c' (circular)
    integer   :: n,i  ! number of rays
    real      :: a    ! aperture of bundle
    real      :: t=0. ! tangent of incidence angle
  end type bundle

contains
!------------------------------------------------------------------

function in(a,b)            ! Inner product
  real :: in, a(3), b(3)
  in=sum(a*b)
end function in
!------------------------------------------------------------------

subroutine renorm_ray(r,rn)
  type (lray)   :: r,rn

  rn%a=r%a/r%a(3)
  rn%b=r%b-r%b(3)*rn%a
end subroutine renorm_ray
!------------------------------------------------------------------

subroutine find_focus(rr,nr,c)
  integer     :: nr
  type (lray) :: rr(nr)
  type (comp) :: c

  real        :: as,dm
  type (lray) :: r
  integer     :: i
  real        :: a(3)=0.,b(3)=0.,ab=0.,aa=0.,bb=0.

  do i=1,nr
     call renorm_ray(rr(i),r) ; r%a(3)=0.
     a=a+r%a/nr
     b=b+r%b/nr
     aa=aa+in(r%a,r%a)/nr
     ab=ab+in(r%a,r%b)/nr
     bb=bb+in(r%b,r%b)/nr
  enddo
  as=aa-in(a,a) 
  c%d=(in(a,b)-ab)/as
  dm=bb-in(b,b)-as*c%d*c%d
  dm=sqrt(max(0.,dm))
  write(*,'(a,f10.4)')"Focus distance from last component (mm): ",c%d
  write(*,'(a,f10.4)')"Size of bundle      in focal plane (mm): ",dm
end subroutine find_focus
!------------------------------------------------------------------

subroutine set_color(c,nc,cl)
  integer    :: nc 
  type (comp):: c(nc)
  character  :: cl

  integer :: i,j
  real :: x1,x2=1.0

  do i=1,nc-1
     j=index(c(i)%ri%cl,cl) ! Find wave length
     if(j==0) call abort()
     x1=c(i)%ri%xn(j)
     c(i)%np2m=(x1/x2)**2. - 1.0
     x2=x1
  enddo

end subroutine set_color
!------------------------------------------------------------------
subroutine ini_bundle(b,a,n,t)
  type (bundle) :: b
  real    :: a,t
  integer :: n

  integer*4,parameter :: seed = 86456
  call srand(seed)

  b%a=a
  b%t=t
  b%i=0
  b%n=n
end subroutine ini_bundle
  
!------------------------------------------------------------------

subroutine shoot_ray(r,b,lok)
  type (bundle) :: b
  type (lray)   :: r
  logical       :: lok

  real :: x,y

  lok=(b%i<b%n) ; if (.not.lok) return

  r%a=(/0.,0.,1./) ; r%a(1)=b%t
  do
    x=2.*rand()-1.
    y=2.*rand()-1.
    if (x*x+y*y <= 1.0) then
       r%b(1)=x*0.5*b%a
       r%b(2)=y*0.5*b%a
       r%b(3)=0
       exit
    endif
  enddo
  b%i=b%i+1

end subroutine shoot_ray
!------------------------------------------------------------------
subroutine shape(z,dz,x,c)  ! Component Shape
! Cone surface
  real        :: z,dz,x
  type (comp) :: c

  real :: r,d

  select case (c%surf)
    case ('c') ! Cone surface
      r =abs(c%r)
      d =sqrt(r*r-c%scp*x*x)
      z =x*x/(r+d)
      dz=x/d
      if (c%r<0.) then
         z =- z
         dz=-dz
      endif
    case default ! flat surface
      z=0; dz=0
  end select
end subroutine shape
!------------------------------------------------------------------

subroutine cut(r,c,lok)
! Find intersection between light ray and component
  type (lray) :: r
  type (comp) :: c
  logical     :: lok

  real :: lam, flam, h,z,dz
  real, parameter :: eps=1.E-12 ! thats 10^-5 times the size of an atom
  integer :: iter,niter=10

  lok=.true.

!-First shift coordinate system towards component
  lam=(c%d-r%b(3))/r%a(3)
  r%b(:)=r%b(:)+lam*r%a(:)
  r%b(3)=0.
  
  if (c%type == 'm') r%a(3)=-r%a(3) ! Mirror z-axis

!-First guess
  lam=0.  ! Flat surface

!-Newton-Raphson loop
  do iter=1,niter
  !-Surface location
    c%s=lam*r%a+r%b 

  !-Distance from axis
    h=sqrt(sum(c%s(:2)*c%s(:2)))
    call shape(z,dz,h,c)

  !-Surface gradient
    if (h==0.) c%n(:2) = 0.
    if (h/=0.) c%n(:2) = dz*c%s(:2)/h 
               c%n( 3) =-1.

  !-Update estimate 
    flam=z-c%s(3) 
    lam=lam-flam/in(r%a,c%n) 
    if (abs(flam) < eps) exit  ! Done
  enddo
  if (iter>niter)then
     write(*,*)'Warning: too many iterations'
     write(*,*)"c%s: ",c%s
     write(*,*)"r%a: ",r%a
     lok=.false.
  endif

!-Update light ray
  r%b=c%s

!- Check aperture
   if(sum(c%s(:2)*c%s(:2))>0.25*c%a*c%a)lok=.false.

end subroutine cut
!------------------------------------------------------------------

subroutine ref(r,c,lok)
! Ref-lect or ref-ract light ray
  type (lray) :: r
  type (comp) :: c
  logical     :: lok

  integer :: i
  real    :: an(3),p,a

  lok=.true.
  p=in(r%a,c%n)/in(c%n,c%n)
  an=p*c%n   

  select case (c%type)
    case ('m')    
      a=-2
    case ('l')    
      a= in(r%a,r%a)/in(an,an)
      if (a<-c%np2m) then
         write(*,*)"Light ray cannot exit medium"
         lok=.false. 
         return
      endif
      a= sqrt(1.+c%np2m*a)-1.
    case default ; a= 0
  end select

! Update light ray direction
  r%a = r%a + a*an
end subroutine ref

end module spotmod
