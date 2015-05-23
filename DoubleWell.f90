!-------------------------------------------------------------------------
!> @brief  double well benchmark for monte carlo simulation
!>          
!> @date  22-May-2015
!> @author HY
!-------------------------------------------------------------------------
Program DoubleWell
implicit none
integer 	::Ninit                          !<steps to  equilibrate initially
integer 	::Nmeasure                       !<steps for measurement
real*8		::RefT				 !<reference   temperature
real*8		::MeasureT			 !<measurement temperature
real*8  	::PotentialB			 !<potential barrier in -1.5<= x <= 1.5
real*8          ::Coordinate			 !<position of the particle 
real*8          ::Potential			 !<potential of the particle
real*8          ::MeanCoord			 !<mean coordinate
real*8		::MeanPotential			 !<mean potential
!real*8		::statCoord			 !<expectation of coord of multipule MC simulation
!real*8		::statPotential			 !<expectation of potential of multipule MC simulation
!integer	::Nloop				 !<cycle of MC simulation
!integer	::i				 !<cycle counter
real*8	        ::errorCoord,errorPotential
!-----default----------------------------------------------
Ninit=1E3
Nmeasure=1E5
MeasureT=0.1
!Nloop=1E3
write(6,'(A30,X,I6)') 'steps for inital equilibrium:',Ninit
write(6,'(A30,X,I6)') 'steps for measurement:',       Nmeasure
write(6,'(A30,X,E20.12)') 'measurement temperature',  MeasureT
!write(6,'(A30,X,E20.12)') 'num of cycle of MC simulation',  Nloop

!-----read initial condition--------------------------------
write(6,'(A)')  'please input the potential B between(-1.5) and (-0.5):'
read(5,*) PotentialB
write(6,'(A)')      'please input the initial coordinate:'
read(5,*) Coordinate
write(6,'(A)')      'please input the reference temperature:'
read(5,*) RefT

write(6,'(A30)') 	  'finish reading input.'
write(6,'(A30,X,E20.12)') 'potential B :',	    PotentialB
write(6,'(A30,X,E20.12)') 'inital coordinate :',    Coordinate
write(6,'(A30,X,E20.12)') 'reference temperature :',RefT

!---starting MC simulation--------------------------------

!----inital steps to equilibrate---------------------------
call ETMC(Ninit,Coordinate,PotentialB,MeasureT,RefT,MeanCoord,MeanPotential,errorCoord,errorPotential)

!---measurement--------------------------------------------
Coordinate=MeanCoord
call ETMC(Nmeasure,Coordinate,PotentialB,MeasureT,RefT,MeanCoord,MeanPotential,errorCoord,errorPotential)
write(6,'(2(A,X,E20.12,X))') 'MeanCoord=',MeanCoord,'errorCoord',errorCoord
write(6,'(2(A,X,E20.12,X))') 'MeanPotential=',MeanPotential,'errorPotential',errorPotential

end Program DoubleWell

subroutine getPotential(Coordinate,PotentialB,Potential)
!---------------------------------------------------
!get inital Potential energy according to coordinate
!---------------------------------------------------
implicit none
real*8,intent(inout) 	::Coordinate
real*8,intent(in) 	::PotentialB
real*8,intent(out)	::Potential

if    (Coordinate >=  1.5d0 ) then	 
      	 	 Potential=1.d9		 !use this to simulate infinity
elseif(Coordinate >  0.5d0 .and. Coordinate <  1.5d0) then 
      	 	Potential= 0.d0
elseif(Coordinate >= -0.5d0 .and. Coordinate <=  0.5d0) then
       	 	Potential= 1.d0
elseif(Coordinate >  -1.5d0 .and. Coordinate < -0.5d0) then
         	Potential=PotentialB
elseif(Coordinate <= -1.5d0) then
         	Potential=1.d9		!use this to simulate infinity
end if

return

end subroutine getPotential


subroutine ETMC(Nstep,Coordinate,PotentialB,MeasureT,RefT,MeanCoord,MeanPotential,errorCoord,errorPotential)
!---------------------------------------------------                              
!Effective Temperature Monte Carlo
!
!	 <A e^-(beta-beta*)H>|RefT
!<A>|T = --------------------------
!	 <  e^-(beta-beta*)H>|RefT
!---------------------------------------------------
implicit none
integer,intent(in)	 ::Nstep
real*8, intent(in)	 ::Coordinate,PotentialB,MeasureT,RefT
real*8, intent(out)	 ::MeanPotential,MeanCoord,errorCoord,errorPotential

integer			 ::i
real*8		         ::Coord,CoordOld,Pot,PotOld,dX,rand,criteria
real*8			 ::const,coef,sumcoef
real*8			 ::Potential,Mean2Coord,Mean2Potential

const=1.0/MeasureT-1.0/RefT
Coord=Coordinate
!-----get inital potential energy--------------------------
call getPotential(Coord,PotentialB,Potential)
!write(6,'(A,X,E20.12)') 'initial Potential at this position is:',Potential
Pot  =Potential

sumcoef	      =0.d0
MeanPotential =0.d0
MeanCoord     =0.d0
Mean2Potential=0.d0
Mean2Coord    =0.d0
errorPotential=0.d0
errorCoord    =0.d0

coef=dexp(-const*Pot)
sumcoef=sumcoef+coef
MeanCoord=MeanCoord+Coord*coef
MeanPotential=MeanPotential+Pot*coef
Mean2Coord=Mean2Coord+Coord*Coord*coef
Mean2Potential=Mean2Potential+Pot*Pot*coef

dX=0.5d0

call random_seed()

do i=1,Nstep   
    call random_number(rand) !0<=rand<1
    !update coord
    CoordOld=Coord
    Coord=Coord + (2.d0*rand-1.d0)*dX
  
    !update potential
    PotOld=Pot
    call getPotential(Coord,PotentialB,Pot)
   

    !decide whether to step forward
    criteria=dexp(-(Pot-PotOld)/RefT)
    call random_number(rand)

    if ( Pot <= PotOld ) then
        Pot =Pot
    elseif(Pot>PotOld) then
     if( criteria >  rand ) then
        Pot =Pot
     else
      Pot=PotOld
      Coord=CoordOld
     endif
    end if
    
    !info sampling
    coef=dexp(-const*Pot)
    sumcoef=sumcoef+coef
    MeanCoord=MeanCoord+Coord*coef
    MeanPotential=MeanPotential+Pot*coef
    Mean2Coord=Mean2Coord+Coord*Coord*coef
    Mean2Potential=Mean2Potential+Pot*Pot*coef
end do
  
MeanCoord=MeanCoord/sumcoef
MeanPotential=MeanPotential/sumcoef
Mean2Coord=Mean2Coord/sumcoef
Mean2Potential=Mean2Potential/sumcoef

errorCoord    =dsqrt((Mean2Coord-MeanCoord*MeanCoord)/Nstep)
errorPotential=dsqrt((Mean2Potential-MeanPotential*MeanPotential)/Nstep)
return
end subroutine ETMC

        
