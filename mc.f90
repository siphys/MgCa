module constants
implicit none
integer,parameter:: dp=selected_real_kind(15,300)
integer,parameter:: l=3 !size of plane lxl
integer,parameter:: iter=1e3 !number of iterations
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
integer::i,j,x_1,y_1,x_0,y_0,nn_old0,nn_new0,nn_old1,nn_new1
i=1

do while(i<=iter) !do for a number of iterations 
    !!!!!!!switch both a Mg and Ca at two random points on the lattice!!!!
    call RANDOM_NUMBER(r)
    x_0=(2*l)*r+1 !random x pos on lattice
    call RANDOM_NUMBER(r)
    y_0=l*r+1!random y
    if( lattice(x_0,y_0) .eq. 0)then
        do while(lattice(x_0,y_0) .eq. 0)!do until atom is not oxygen
            call RANDOM_NUMBER(r)
            x_0=(2*l)*r+1 
            call RANDOM_NUMBER(r)
            y_0=l*r+1
        end do 
    end if
    nn_old0=nn(x_0,y_0)!number of atoms of same type neighbouring before switch
 
    call RANDOM_NUMBER(r)
    x_1=(2*l)*r+1 !second atom to be switched
    call RANDOM_NUMBER(r)
    y_1=l*r+1
    if( (lattice(x_1,y_1) .eq. 0 ).or. (lattice(x_1,y_1) .eq. lattice(x_0,y_0)) )then!if atom is not opposite the first chosen
        do while((lattice(x_1,y_1) .eq. 0 ).or. (lattice(x_1,y_1) .eq. lattice(x_0,y_0)) )!do until atom opposite the first switch is found
            call RANDOM_NUMBER(r)
            x_1=(2*l)*r+1 !random x pos on lattice
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

    j=1
    do while(j<=l)
          write(10,*)lattice(1,j),lattice(2,j),lattice(3,j),lattice(4,j),lattice(5,j),lattice(6,j)!write out permutations
          j=j+1
    end do
    write(10,*)"#########"
    i=i+1
end do


contains
function nn(a,b)
integer::a,b,nn
nn=0
if( (a+1) > 2*l)then
   if ( lattice(a,b) == lattice(1,b))then
        nn=nn+1
   end if
else
   if ( lattice(a,b) == lattice(a+1,b))then
      nn=nn+1
   end if
end if
if( (a-1) < 0)then
   if ( lattice(a,b) == lattice(2*l,b))then
        nn=nn+1
   end if
else
   if ( lattice(a,b) == lattice(a-1,b))then
      nn=nn+1
   end if
end if
if( (b+1) > l)then
   if ( lattice(a,b) == lattice(a,1))then
        nn=nn+1
   end if
else
   if ( lattice(a,b) == lattice(a,b+1))then
      nn=nn+1
   end if
end if
if( (b-1) < 0)then
   if ( lattice(a,b) == lattice(a,l))then
        nn=nn+1
   end if
else
   if ( lattice(a,b) == lattice(a,b-1))then
      nn=nn+1
   end if
end if
end function

end subroutine




end module

program monte_carlo
use montecarlo
use constants
implicit none
integer :: i,j
character(2)::atom
real(kind=dp)::r
integer,dimension(:,:),allocatable :: lattice

allocate (lattice(2*l,l))

open (unit=15,file='atoms.dat')
i=1
j=1!!Set up array for the given surface 
do while(i<=l)
   do while(j<=2*l)   
      read(15,*)atom
      if (atom == 'Mg')then
         lattice(j,i)=1
      else if (atom == 'O')then
         lattice(j,i)=0
      else
         lattice(j,i)=-1
      end if
      j=j+1
   end do
   j=1
   i=i+1
end do



OPEN (UNIT=10,FILE='interface.txt')
call MC(lattice)




deallocate (lattice)
close(10)
END PROGRAM
