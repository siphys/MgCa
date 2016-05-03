module constants
implicit none
integer,parameter:: dp=selected_real_kind(15,300)
integer,parameter:: l=20 !size of plane lxl
integer,parameter:: iter=1 !number of iterations
end module




module montecarlo
use constants
implicit none
contains

subroutine MC(lattice)
integer,parameter:: dp=selected_real_kind(15,300)
!integer,intent(in)::
integer,dimension(:,:),allocatable,intent(inout)::lattice
real(kind=dp)::r,energy,energy_trial,diff
integer::i,j,k,x_1,y_1,x_0,y_0,nn_old0,nn_new0,nn_old1,nn_new1,m
i=1

do while(i<=iter) !do for a number of iterations 
    do while(m<=l**2)
        !!!!!!!switch both a Mg and Ca at two random points on the lattice!!!!
        call RANDOM_NUMBER(r)
        x_0=(l)*r+1 !random x pos on lattice
        call RANDOM_NUMBER(r)
        y_0=l*r+1!random y
        if( lattice(x_0,y_0) .eq. 0)then
            do while(lattice(x_0,y_0) .eq. 0)!do until atom is not oxygen
                call RANDOM_NUMBER(r)
                x_0=(l)*r+1 
                call RANDOM_NUMBER(r)
                y_0=l*r+1
            end do 
        end if
        nn_old0=nn(x_0,y_0)!number of atoms of same type neighbouring before switch
     
        call RANDOM_NUMBER(r)
        x_1=(l)*r+1 !second atom to be switched
        call RANDOM_NUMBER(r)
        y_1=l*r+1
        if( (lattice(x_1,y_1) .eq. 0 ).or. (lattice(x_1,y_1) .eq. lattice(x_0,y_0)) )then!if atom is not opposite the first chosen
            do while((lattice(x_1,y_1) .eq. 0 ).or. (lattice(x_1,y_1) .eq. lattice(x_0,y_0)) )!do until atom opposite the first switch is found
                call RANDOM_NUMBER(r)
                x_1=(l)*r+1 !random x pos on lattice
                call RANDOM_NUMBER(r)
                y_1=l*r+1
            end do 
        end if
        nn_old1=nn(x_1,y_1)!number of neighbours of the same type as atom 2
        lattice(x_0,y_0)=-lattice(x_0,y_0)!switch position note oxygen remains fixed
        nn_new0=nn(x_0,y_0)!nn after switch
        lattice(x_1,y_1)=-lattice(x_1,y_1)!switch second
        nn_new1=nn(x_1,y_1)
       
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111
        diff=nn_new1-nn_old1 + nn_new0-nn_old0 !diff in nn, more nn of same type aprox more energy
        if(diff > 0)then !energy of new state higher energy
             call RANDOM_NUMBER(r)
             if( r > exp(diff) )then
                  lattice(x_0,y_0)=-lattice(x_0,y_0)!revert atoms back
                  lattice(x_1,y_1)=-lattice(x_1,y_1)
             end if
        end if

        m=m+1
   end do
   
   j=1
   k=1
   do while(j<=l)
      do while(k<=l)
         if(lattice(k,j) == 1)then
              write(10,*)'Mg'
          else if(lattice(k,j) == 0)then
              write(10,*)'O'
          else
              write(10,*)'Ca'
          end if
          k=k+1
       end do
       k=1
       j=j+1
   end do
   write(10,*)
   i=i+1
end do



contains
function nn(a,b)
integer::a,b,k,nn
nn=0

k=a+1
if (k > l)then
  k=1
end if
do while(lattice(k,b) == 0)!make sure to skip oxygen
   k=k+1
   if (k > l)then
      k=1
   end if
end do
if( lattice(a,b) == lattice(k,b))then
  nn=nn+1
end if
k=a-1
if (k < 0)then
  k=l
end if
do while(lattice(k,b) == 0)
   k=k-1
   if (k < 0)then
      k=l
   end if
end do
if( lattice(a,b) == lattice(k,b))then
  nn=nn+1
end if

k=b+1
if (k > l)then
  k=1
end if
do while(lattice(a,k) == 0)
   k=k+1
   if (k > l)then
      k=1
   end if
end do
if( lattice(a,b) == lattice(a,k))then
  nn=nn+1
end if
k=b-1
if (k < 0)then
  k=l
end if
do while(lattice(a,k) == 0)
   k=k-1
   if (k < 0)then
      k=l
   end if
end do
if( lattice(a,b) == lattice(a,k))then
  nn=nn+1
end if

end function

end subroutine




end module

program monte_carlo
use montecarlo
use constants
implicit none
integer :: i,j,k!xyz
character(2)::atom
real(kind=dp)::r
integer,dimension(:,:),allocatable :: lattice

allocate (lattice(l,l))

open (unit=15,file='atoms.dat')
OPEN (UNIT=12,FILE='pos.txt')
k=0
i=1
j=1!!Set up array for the given surface 
do while(j<=l)
   do while(i<=l)   
      read(15,*)atom
      if (atom == 'Mg')then
         lattice(i,j)=1
      else if (atom == 'O')then
         lattice(i,j)=0
      else
         lattice(i,j)=-1
      end if
      write(12,*)i,j,k
      i=i+1
   end do
   i=1
   j=j+1
end do



OPEN (UNIT=10,FILE='permutations.txt')
call MC(lattice)





deallocate (lattice)
close(10)
close(12)
END PROGRAM
