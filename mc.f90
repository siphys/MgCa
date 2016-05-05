module constants
implicit none
integer,parameter:: dp=selected_real_kind(15,300)
integer,parameter:: l=10 !size of plane lxl
integer,parameter:: iter=10000 !number of iterations
integer :: count = 0
integer :: perms_to_print = 1
end module




module montecarlo
use constants
implicit none
contains

subroutine MC(lattice)
integer,parameter:: dp=selected_real_kind(15,300)
!integer,intent(in)::
integer,dimension(:,:),allocatable,intent(inout)::lattice
real(kind=dp)::r,diff
integer::i,j,k,x_1,y_1,x_0,y_0,E_old,E_new,m
i=1

do i=1,iter !do for a number of iterations 
    m=0
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
        E_old=E()!number of neighbours of the same type as atom 2
        lattice(x_0,y_0)=-lattice(x_0,y_0)!switch position note oxygen remains fixed
        lattice(x_1,y_1)=-lattice(x_1,y_1)!switch second
        E_new=E()
       
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111
        diff=E_new-E_old !diff in nn, more nn of same type aprox more energy
        if(diff > 0)then !energy of new state higher energy
             call RANDOM_NUMBER(r)
             if( r*l > exp(diff/30) )then
  
                  lattice(x_0,y_0)=-lattice(x_0,y_0)!revert atoms back
                  lattice(x_1,y_1)=-lattice(x_1,y_1)
             end if
        end if
        m=m+1
   end do


   if(i>iter-perms_to_print)then
     j=1
     k=1

     do while(j<=l)
        do while(k<=l)
           if(lattice(k,j) == 1)then
                write(10,*)'Mg'
            else if(lattice(k,j) == 0)then
                write(10,*)'O'
            else
            end if
            k=k+1
         end do
         k=1
         j=j+1
     end do
     write(10,*)
   end if
   
end do




contains
function E()
integer::x,y
real(kind=dp)::E
E=0
x=1
y=1

do while(y<=l)
   x=1
   do while(x<=l)
      j=1
      k=1
      do while(k<=l)
         j=1
         do while(j<=l)
            if( (j .ne. x) .and. (k .ne. y))then
                E=E+2*lattice(x,y)*2*lattice(j,k)!energy calc using coulomb
            end if
            j=j+1
         end do
         k=k+1
      end do
      x=x+1
   end do
   y=y+1
end do


end function

end subroutine




end module

program monte_carlo
use montecarlo
use constants
implicit none
integer :: i,j,k!xyz
real(kind=dp)::r
integer,dimension(:,:),allocatable :: lattice

allocate (lattice(l,l))
lattice = 0

open (unit=15,file='atoms.dat')
OPEN (UNIT=12,FILE='pos.txt')
k=0
i=1
j=1!!Set up array for the given surface 



count = 0
do j=1, l 
    do i=1, l !CREATES 50% RATIO OF O to 25% Ca 25% Mg 
        if (mod(i+j, 2)==0) then    
            lattice(i, j) = 0 !oxygen
        else 
            if (mod(count, 2)==0) then
                lattice(i, j)=1 !Mg
            else
                lattice(i, j) = -1!Ca
            end if
            count = count + 1
        end if
        write(12,*)i,j,0
    end do
end do


OPEN (UNIT=10,FILE='permutations.txt')
call MC(lattice)





deallocate (lattice)
close(10)
close(12)
END PROGRAM
