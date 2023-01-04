program div_vor
implicit none
real,dimension(20,131,105)::u,v,du_dx,dv_dx,du_dy,dv_dy,div,rel_vor,abs_vor
real::t,f(131),dx(131)
integer::i,j,k,p(20),dy

open(300,file='pres_levels.dat')
open(100,file='u_20_levels.dat')
open(200,file='v_20_levels.dat')
open(400,file='read_data.dat')
open(500,file='abs_vor.csv')
open(600,file='divergence.csv')

do k=1,20
read(300,*)p(k)
end do

do k=1,20
do i=1,131
read(100,*)(u(k,i,j),j=1,105)
read(200,*)(v(k,i,j),j=1,105)
end do
end do

dy=0.5*110000.0
i=1
do t=-20.0,45.0,0.5
dx(i)=0.625*cos(t*(3.14/180.0))*110000.0
if(i.le.131)i=i+1
end do
i=1
do t=-20.0,45.0,0.5
f(i)=2.0*7.292*(10.0**(-5.0))*sin(t*(3.14/180.0))
if(i.le.131)i=i+1
end do
write(*,*)f

call derivative(u,du_dx,du_dy)
call derivative(v,dv_dx,dv_dy)

do k=1,20
do i=1,131
do j=1,105
div(k,i,j)=du_dx(k,i,j)+dv_dy(k,i,j)
rel_vor(k,i,j)=dv_dx(k,i,j)-du_dy(k,i,j)
abs_vor(k,i,j)=rel_vor(k,i,j)+f(i)
end do 
end do
end do

do k=1,20
do i=1,131
write(400,*)f(i)
write(500,10)(abs_vor(k,i,j),j=1,105)
write(600,10)(div(k,i,j),j=1,105)
end do
end do

10 format(105(f12.7,1x))

contains
subroutine derivative(p,dp_dx,dp_dy)
implicit none
integer::k,i,j
real,dimension(20,131,105)::p,dp_dx,dp_dy

do k=1,20
do i=1,131
do j=1,105

if(j.gt.1.and.j.lt.105)dp_dx(k,i,j)=(p(k,i,j+1)-p(k,i,j-1))/(2*dx(i))
if(i.gt.1.and.i.lt.131)dp_dy(k,i,j)=(p(k,i+1,j)-p(k,i-1,j))/(2*dy)

if(i.eq.1)dp_dy(k,i,j)=(p(k,i+1,j)-p(k,i,j))/(dy)
if(j.eq.1)dp_dx(k,i,j)=(p(k,i,j+1)-p(k,i,j))/(dx(i))

if(i.eq.131)dp_dy(k,i,j)=(p(k,i,j)-p(k,i-1,j))/(dy)
if(j.eq.105)dp_dx(k,i,j)=(p(k,i,j)-p(k,i,j-1))/(dx(i))
 
end do 
end do
end do

end subroutine derivative
end program div_vor




