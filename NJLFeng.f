       subroutine getnjlFeng(ixj)
      PARAMETER (MAXPTN=40001,MAXR=1000)
      PARAMETER (Muln=3000,Nevn=1000,ncon=101,nxg=101,nyg=101,nzg=101)	  
      implicit real*8 (a-h,o-z)
      common/particle/ par(9,Muln,Nevn) !, ntotal(Nevn) 

      COMMON /prec1/GX0(MAXPTN),GY0(MAXPTN),GZ0(MAXPTN),FT0(MAXPTN),
     &       PX0(MAXPTN), PY0(MAXPTN), PZ0(MAXPTN), E0(MAXPTN),
     &       XMASS0(MAXPTN), ITYP0(MAXPTN)
      COMMON /ilist7/ LSTRG0(MAXPTN), LPART0(MAXPTN)
      COMMON /ilistnjl/ LSTRGnjl(MAXPTN,MAXR), LPARTnjl(MAXPTN,MAXR)
c      COMMON /PARA1/ MUL
      COMMON /NJLMUL/ NJLMUL(MAXR)
      COMMON /PARA1/ MUL
      common /blastwave/parBL(9,MAXPTN)

            call system('rm NJLMUL.dat')

            call system('touch NJLMUL.dat')
      open(unit=32, file='NJLMUL.dat', status='unknown')
! number of partitions (max=10)
         write(*,*)"getnjl, events=",ixj
        npart = 1 
!      do npara=1, nevent !<----------------------------------------------------- new event
        npara = ixj
!       iparthd(npara)=0 !<======this is the array recording the # of hadronized partons in each event.
       !<=========Initially, all the partons aren't hadronized, so set all the components to zero.
       k=0
       do npartition=1, npart
!          read(17,*) iev, i2, ipart(npara), bim, i5, i6, i7, i8 
!		  read(17,*) nsg(npara) 

!          do i=1, ipart(npara)
            write(32,*)ixj,MUL
            close(32)
            call system('./NJLballini.sh')
            !call NJLball_ini(mul)

             do i=1,mul
!             read(17,*) id,lstring,lpart,xx,yy,zz,qx,qy,qz,xmass,time 
             k=k+1 !k=i
                !par(1,k,npara)=parBL(1,I) !GX0(I)
                !  par(2,k,npara)=parBL(2,I) !GY0(I)
               !par(3,k,npara)=parBL(3,I) !GZ0(I)
               !par(4,k,npara)=FT0(I)
               ! par(5,k,npara)=parBL(5,I) !PX0(I)
               !par(6,k,npara)=parBL(6,I) !PY0(I)
               !par(7,k,npara)=parBL(7,I) !PZ0(I)
               !par(8,k,npara)=parBL(8,I) !XMASS0(I)
               !par(9,k,npara)=parBL(9,I) !float(ITYP0(I))
!c			 ls(k,npara) = LSTRG0(I)
!c			 lp(k,npara) = LPART0(I) 
!propagate to formation time
       !E = SQRT(par(5,k,npara)**2+par(6,k,npara)**2+par(7,k,npara)**2
       !&+par(8,k,npara)**2)
        !par(1,k,npara) = par(1,k,npara) + par(5,k,npara)/E*par(4,k,npara)
       !par(2,k,npara) = par(2,k,npara) + par(6,k,npara)/E*par(4,k,npara)
        !par(3,k,npara) = par(3,k,npara) + par(7,k,npara)/E*par(4,k,npara)

        LSTRGnjl(i,ixj) = LSTRG0(i)
        LPARTnjl(i,ixj) = LPART0(i)

        enddo ! for particles in one partition
       enddo ! for particles in one event
!    itotal(npara)=ntotal
!       ntotal(npara)=k
!	enddo ! for events
         NJLMUL(npara) = mul
         write(*,*)"getnjl, quark number=",mul
        RETURN
        END

       subroutine getnjlFengback(nevent)
      PARAMETER (MAXPTN=40001,MAXR=1000)
      implicit real*8 (a-h,o-z)
      COMMON /precnjl/GXnjl(MAXPTN,MAXR),GYnjl(MAXPTN,MAXR)
     &,GZnjl(MAXPTN,MAXR),FTnjl(MAXPTN,MAXR)
     &,PXnjl(MAXPTN,MAXR), PYnjl(MAXPTN,MAXR), PZnjl(MAXPTN,MAXR)
     &,Enjl(MAXPTN,MAXR),XMASSnjl(MAXPTN,MAXR), ITYPnjl(MAXPTN,MAXR)
      COMMON /NJLMUL/ NJLMUL(MAXR)	 
        open (unit=99, file='zpcBW.dat', status='unknown')	  
        !open(unit=106, file='NJLEVE.dat', status='unknown')	
        !read(106,*)nevent
		
	  do npara=1, nevent

          read(99,*) iev, i2, ipart, r1, i4, i5, i6, i7,i8  	 
       !write(*,*)"data back:",nevent,npara,ipart		  
	  do i=1,ipart !NJLMUL(npara)
          read(99,*) qx, qy, qz, ispc, xmass, xx, yy, zz, time	
	      Eqs=sqrt(qx**2.+qy**2.+qz**2.)	
       !write(*,*)"data back:",i,qx, qy, qz, ispc	  
	    GXnjl(i,npara) = xx
	   GYnjl(i,npara) = yy
	   GZnjl(i,npara) = zz
           FTnjl(i,npara) = time
	  PXnjl(i,npara) = qx
	   PYnjl(i,npara) = qy
	  PZnjl(i,npara) = qz
	   XMASSnjl(i,npara) = xmass
	   Enjl(i,npara) = Eqs
	
       !write(*,*)"data back:",xx,yy,zz,time,xmass,Eqs	   
	   ITYPnjl(i,npara) = int(ispc)
	   enddo
       enddo
        RETURN
        END	 
		



       subroutine HadReSamp(npara0,ity)
        PARAMETER (MAXSTR=150001,MAXSTR1=10000)
      PARAMETER (MAXPTN=40001,MAXR=1000)
      implicit real*8 (a-h,o-z)		  
      COMMON /precnjl/GXnjl(MAXPTN,MAXR),GYnjl(MAXPTN,MAXR)
     &,GZnjl(MAXPTN,MAXR),FTnjl(MAXPTN,MAXR)
     &,PXnjl(MAXPTN,MAXR), PYnjl(MAXPTN,MAXR), PZnjl(MAXPTN,MAXR)
     &,Enjl(MAXPTN,MAXR),XMASSnjl(MAXPTN,MAXR), ITYPnjl(MAXPTN,MAXR) 
        common /tmpquark/ ITYu(MAXSTR1),
     &       GXu(MAXSTR1), GYu(MAXSTR1), GZu(MAXSTR1), FTu(MAXSTR1),
     &       PXu(MAXSTR1), PYu(MAXSTR1), PZu(MAXSTR1), PEu(MAXSTR1),
     &       XMu(MAXSTR1),ITYd(MAXSTR1),
     &       GXd(MAXSTR1), GYd(MAXSTR1), GZd(MAXSTR1), FTd(MAXSTR1),
     &       PXd(MAXSTR1), PYd(MAXSTR1), PZd(MAXSTR1), PEd(MAXSTR1),
     &       XMd(MAXSTR1)
      COMMON /NJLMUL/ NJLMUL(MAXR)	 
      common /kjquark/kjgx,kjgy,kjgz
        !open (unit=99, file='zpcBW.dat', status='unknown')	  
        !open(unit=106, file='NJLEVE.dat', status='unknown')	
        !read(106,*)nevent
		
              print *,"HadReSamp input:",npara0,ity,gx,gy,gz  		
	
       numqu = 0
       numqd = 0
	  do npara=npara0, npara0	  
	  do i=1,NJLMUL(npara)	 
         ityq = ITYPnjl(i,npara) 
          if (ityq.eq.1)then
              numqd=numqd+1
			  
              Gxd(numqd) = GXnjl(i,npara)
              Gyd(numqd) = Gynjl(i,npara)
              Gzd(numqd) = Gznjl(i,npara)
              ftd(numqd) = ftnjl(i,npara)			  
           endif

          if (ityq.eq.2)then
              numqu=numqu+1
			  
              Gxu(numqu) = GXnjl(i,npara)
              Gyu(numqu) = Gynjl(i,npara)
              Gzu(numqu) = Gznjl(i,npara)
              ftu(numqu) = ftnjl(i,npara)			  
           endif		   
		   
	   enddo
       enddo
	   
        !print *,"hadron resample begin, events=",npara0  

              Is = 0
              if (abs(ity).eq.2112)then
               Is=1
                 r = rand()
				 k = int( r*float(numqd) )+1
                  if (k.gt.numqd)k=numqd
                 print *,r,k,numqd,gxd(k),gyd(k),gzd(k)
                 kjgx = gxd(k)
                 kjgy = gyd(k)
                 kjgz = gzd(k)
                 !ft = ftd(k)
               elseif  (abs(ity).eq.2212)then
               Is=1
			   
                 r = rand()
				 k = int( r*float(numqu) )+1
                  if (k.gt.numqu)k=numqu
                  print *,r,k,numqu,gxu(k),gyu(k),gzu(k)
                 kjgx = gxu(k)
                 kjgy = gyu(k)
                 kjgz = gzu(k)
                 !ft = ftu(k)
				 
              endif
			  
              !if (Is.eq.1)then

              !print *,"before:",gx,gy,gz               
              !gx=gx1
              !gy=gy1
              !gz=gz1
              !ftar(I)=ft
              !print *, "after:",gx,gy,gz
              !endif
			  
         !print *,"hadron resample done!, events=",npara0		   
        RETURN
        END			
