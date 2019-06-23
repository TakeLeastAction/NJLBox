	program NJLBoxGv

	implicit real*8 (a-h,o-z)
    integer,parameter :: NumPar= 8000, Ntest=1000, NCenX=0,NCenY=0,NCenZ=0                         !,NCenZ=3  NumPar: maximum number of particles, Ntest: maximum number of test particle
                                                                                          ! 7.7GeV NumPar=2000,Ntest=1000
    integer,parameter :: Ngx = 31, Ngy = 31, Ngz = 31, NgN = 150                                              ! must be larger than nx,ny,nz
!********************************* particle ******************************
    dimension par(9,NumPar,Ntest), q(9,NumPar,Ntest), q0(3,NumPar,Ntest), ihd(NumPar,Ntest), qq(9,NumPar,Ntest)
    dimension icoll(NumPar,Ntest)                                                                                       ! ihd: <======It's an array denoting whether a parton is hadronized or not. ihd=1 means hadronized.
                                                                                           ! par: x,y,z,t,px,py,pz,m,id
    dimension spin(3,NumPar,Ntest)                                                           ! spin                                                                                    
!**************************** spectators **************************************                                                                                           
    dimension spec(11,NumPar,Ntest),nospec(Ntest)                                                                                       
 
!********************************* cells  ********************************     
    common/grid1/nx, ny, nz                                                                ! grid number
    common/grid2/dx, dy, dz                                                                ! size of grids
    common/scalar/ Qmu(Ngx,Ngy,Ngz), Qmd(Ngx,Ngy,Ngz), Qms(Ngx,Ngy,Ngz)                    ! Quark mass in each cell; 
    common/BoundMass/ BoundQmu(Ngx+2,Ngy+2,Ngz+2), BoundQmd(Ngx+2,Ngy+2,Ngz+2), BoundQms(Ngx+2,Ngy+2,Ngz+2)  ! Quark mass in each cell and one beyond boundary; 
	common/numberb/ numbd(Ngx,Ngy,Ngz),numbu(Ngx,Ngy,Ngz),numbs(Ngx,Ngy,Ngz),numb(Ngx,Ngy,Ngz)               ! numbd: d and anti-d number in each cell
    common/chargeb/ rhobd(Ngx,Ngy,Ngz),rhobu(Ngx,Ngy,Ngz),rhobs(Ngx,Ngy,Ngz)               ! charged baryon density of u quark, d quark ,s quark
	common/force/ forceux(Ngx,Ngy,Ngz),forceuy(Ngx,Ngy,Ngz),forceuz(Ngx,Ngy,Ngz)  
	
    common/vector/ At(Ngx,Ngy,Ngz), Ax(Ngx,Ngy,Ngz), Ay(Ngx,Ngy,Ngz), Az(Ngx,Ngy,Ngz)      ! total vector poential in each cell
    common/Uvector/ Atu(Ngx,Ngy,Ngz), Axu(Ngx,Ngy,Ngz), Ayu(Ngx,Ngy,Ngz), Azu(Ngx,Ngy,Ngz) ! vector potential of (anti-)u quark in each cell
    common/Dvector/ Atd(Ngx,Ngy,Ngz), Axd(Ngx,Ngy,Ngz), Ayd(Ngx,Ngy,Ngz), Azd(Ngx,Ngy,Ngz) ! 
    common/Svector/ Ats(Ngx,Ngy,Ngz), Axs(Ngx,Ngy,Ngz), Ays(Ngx,Ngy,Ngz), Azs(Ngx,Ngy,Ngz) !  
    common/pUvector/ pAxu(Ngx,Ngy,Ngz), pAyu(Ngx,Ngy,Ngz), pAzu(Ngx,Ngy,Ngz)                  ! previous vector potential of u(anti-u) quark
    common/pDvector/ pAxd(Ngx,Ngy,Ngz), pAyd(Ngx,Ngy,Ngz), pAzd(Ngx,Ngy,Ngz)                  ! 
    common/pSvector/ pAxs(Ngx,Ngy,Ngz), pAys(Ngx,Ngy,Ngz), pAzs(Ngx,Ngy,Ngz)                  !    

    common/BoundUvector/ BoundAtu(Ngx+2,Ngy+2,Ngz+2), BoundAxu(Ngx+2,Ngy+2,Ngz+2), BoundAyu(Ngx+2,Ngy+2,Ngz+2), BoundAzu(Ngx+2,Ngy+2,Ngz+2) ! vector potential of (anti-)u quark in each cell
    common/BoundDvector/ BoundAtd(Ngx+2,Ngy+2,Ngz+2), BoundAxd(Ngx+2,Ngy+2,Ngz+2), BoundAyd(Ngx+2,Ngy+2,Ngz+2), BoundAzd(Ngx+2,Ngy+2,Ngz+2) ! 
    common/BoundSvector/ BoundAts(Ngx+2,Ngy+2,Ngz+2), BoundAxs(Ngx+2,Ngy+2,Ngz+2), BoundAys(Ngx+2,Ngy+2,Ngz+2), BoundAzs(Ngx+2,Ngy+2,Ngz+2) !  
    common/BoundpUvector/ BoundpAtu(Ngx+2,Ngy+2,Ngz+2), BoundpAxu(Ngx+2,Ngy+2,Ngz+2), BoundpAyu(Ngx+2,Ngy+2,Ngz+2), BoundpAzu(Ngx+2,Ngy+2,Ngz+2) ! vector potential of (anti-)u quark in each cell
    common/BoundpDvector/ BoundpAtd(Ngx+2,Ngy+2,Ngz+2), BoundpAxd(Ngx+2,Ngy+2,Ngz+2), BoundpAyd(Ngx+2,Ngy+2,Ngz+2), BoundpAzd(Ngx+2,Ngy+2,Ngz+2) ! 
    common/BoundpSvector/ BoundpAts(Ngx+2,Ngy+2,Ngz+2), BoundpAxs(Ngx+2,Ngy+2,Ngz+2), BoundpAys(Ngx+2,Ngy+2,Ngz+2), BoundpAzs(Ngx+2,Ngy+2,Ngz+2) !  
 
    common/vectorgviso/Atgv(Ngx,Ngy,Ngz),Axgv(Ngx,Ngy,Ngz),Aygv(Ngx,Ngy,Ngz),Azgv(Ngx,Ngy,Ngz) 
    dimension pxs(Ngx,Ngy,Ngz), pys(Ngx,Ngy,Ngz), pzs(Ngx,Ngy,Ngz)                            ! totoal momentum in each cell
    dimension vxs(Ngx,Ngy,Ngz), vys(Ngx,Ngy,Ngz), vzs(Ngx,Ngy,Ngz)                            ! velocity of each cell

    dimension rhobxu(Ngx,Ngy,Ngz), rhobyu(Ngx,Ngy,Ngz), rhobzu(Ngx,Ngy,Ngz)                   ! charged baryon vector density of u quark
    dimension rhobxd(Ngx,Ngy,Ngz), rhobyd(Ngx,Ngy,Ngz), rhobzd(Ngx,Ngy,Ngz)                   ! charged baryon vector density of d quark
    dimension rhobxs(Ngx,Ngy,Ngz), rhobys(Ngx,Ngy,Ngz), rhobzs(Ngx,Ngy,Ngz)                   ! charged baryon vector density of s quark
    dimension gamLV(Ngx,Ngy,Ngz)                                                              ! lorenz factor for each cell
!-------------------scalar part ------------------------------
    dimension dd0(Ngx,Ngy,Ngz), uu0(Ngx,Ngy,Ngz), ss0(Ngx,Ngy,Ngz)                            ! sum 1/E/dV in each cell of u,d,s quark
    dimension scalarU(Ngx,Ngy,Ngz), scalarD(Ngx,Ngy,Ngz), scalarS(Ngx,Ngy,Ngz)                ! condensation in each cell    
!-------------------vector part ------------------------------	    
    dimension EdnS(Ngx,Ngy,Ngz),EdnV(Ngx,Ngy,Ngz),EdnT(Ngx,Ngy,Ngz),EdnK(Ngx,Ngy,Ngz)         ! energy density of scalar, vector field.and total energy density in each cell
!-------------------total part  ------------------------------ 
    dimension Edd(Ngx,Ngy,Ngz), Euu(Ngx,Ngy,Ngz), Ess(Ngx,Ngy,Ngz)                            ! total energy of u,d,s quark in each cell;
    dimension EddL(Ngx,Ngy,Ngz), EuuL(Ngx,Ngy,Ngz), EssL(Ngx,Ngy,Ngz),EdnTL(Ngx,Ngy,Ngz),EdnKL(Ngx,Ngy,Ngz)! total energy o u, d, s quark in each cell boost to Local center of mass frame
    dimension EdnEM(Ngx,Ngy,Ngz),EdnEMP(Ngx,Ngy,Ngz)
    dimension EMcrt(NgN,Ngx,Ngy,Ngz),EMcrx(NgN,Ngx,Ngy,Ngz),EMcry(NgN,Ngx,Ngy,Ngz),EMcrz(NgN,Ngx,Ngy,Ngz) ! record charge density and  Electron current density 
    dimension NEMnu(NgN,Ngx,Ngy,Ngz),NEMnd(NgN,Ngx,Ngy,Ngz),NEMns(NgN,Ngx,Ngy,Ngz)
    dimension NEMnub(NgN,Ngx,Ngy,Ngz),NEMndb(NgN,Ngx,Ngy,Ngz),NEMnsb(NgN,Ngx,Ngy,Ngz)
    
    dimension Exf(Ngx,Ngy,Ngz),Eyf(Ngx,Ngy,Ngz),Ezf(Ngx,Ngy,Ngz)
    dimension Bxf(Ngx,Ngy,Ngz),Byf(Ngx,Ngy,Ngz),Bzf(Ngx,Ngy,Ngz)
    
    dimension specEx(Ngx,Ngy,Ngz),specEy(Ngx,Ngy,Ngz),specEz(Ngx,Ngy,Ngz)
    dimension specBx(Ngx,Ngy,Ngz),specBy(Ngx,Ngy,Ngz),specBz(Ngx,Ngy,Ngz)    
    !-----------------Energy-momentum tensor ---------------------
    dimension T00(Ngx,Ngy,Ngz),T11(Ngx,Ngy,Ngz),T22(Ngx,Ngy,Ngz),T33(Ngx,Ngy,Ngz)             ! Energy-momentum tensor in each cell
	dimension T01(Ngx,Ngy,Ngz),T02(Ngx,Ngy,Ngz),T03(Ngx,Ngy,Ngz)
	dimension T12(Ngx,Ngy,Ngz),T13(Ngx,Ngy,Ngz),T23(Ngx,Ngy,Ngz)
	dimension ETens(4,4),ETens1(4,4)                                                         ! Energy-momentum tensor
    dimension xLorMat(4,4),xLorMatinv(4,4)                                                   ! Lorentz transformation Matrix
       
    
	dimension ipart(Ntest),nsg(Ntest),ls(NumPar,Ntest),lp(NumPar,Ntest),iparthd(Ntest)        ! ipart: number of quarks in each event
    dimension itotal(10000), ntotal(10000), TimeN(NgN),dzN(NgN)
!************************************ evolution *************************************************** 
    character(4) :: datatype1
    character(6) :: datatype2
	character(9) :: datatype3
	character(9) :: datatype4
    character(11) :: datatype5
    character(4) :: tit          ! title
    integer :: Nstep        

    
    
!******************************* saving ***************************************  
    open (unit=17, file='box.dat', status='unknown')                                         !changed by sunkj  
    open (unit=100, file='zpcBW.dat', status='unknown')    !------------------------------- Final state quark information   
    open (unit=108, file='zpc_Freeze.dat', status='unknown')    !------------------------------- Freeze-out data    
    open (unit=109, file='zpc_Freeze_Ed.dat', status='unknown')    !------------------------------- Freeze-out data       
    open (unit=101, file='ana/central_NJL.dat',status='unknown') !------------------------------- Central cell evolution 
    open (unit=111, file='ana/central_mass.dat',status='unknown')
    open (unit=102, file='total_NJL.dat',status='unknown')  !------------------------------- Total energy Evolution  
    open (unit=103, file='cen_sden_NJL.dat', status='unknown')
    open (unit=104, file='cen_rhob_NJL.dat', status='unknown')
    open (unit=105, file='cen_rhobx_NJL.dat', status='unknown')
    open (unit=106, file='cen_rhoby_NJL.dat', status='unknown')
    open (unit=107, file='cen_rhobz_NJL.dat', status='unknown')
    open (unit=110, file='ana/collision.dat', status='unknown')
	open (unit=115, file='ana/energykinetic.dat', status='unknown')
    open (unit=1001, file='NumFt.dat',status='unknown')
    open (unit=1002, file='cen_EM.dat',status='unknown')
    !open(5,file='NJLixj.txt')             !number of runs
    !read(5,*)Nixj
    !close(5)
    open(6,file='NJLEVE.dat')          !number of test particle
    read(6,*)nevent
    close(6)
    Nixj=1

!***************************** Control Parameter ********************************
	qLen = 10.0  ! box length  [-qLen, qLen]  
	IperiodicBound = 1  !1 periodic boundary
	iglobal = 4                        !-1:global (central energy density)hadronization 0: local energydensity hadronization
                                        !1: global mass hadronization,2: local quark mass hadronization 
                                        !4: until tend.										
    ievolution = 0                                                                          !1: record evolution information, 0: do not record evolution information
    icollision = 1  ! 0 no collisions
    inewton = 1  ! 0 no hamiltonian equation 
    inewtonM = 1 ! 0 no mass difference contribution
    inewtonPA = 0 ! 0 no contribution from \partial A/\partical t
    iEM = 0                                                                                  !1: including electro-magnetic force	
    ispec = 0                                                                               ! spectators
    ispin = 0                                                                               ! spin 
    ichiral = 1   ! 0 vacuum  1 chiral   2 current
    !print *, "IEM",iem
    thd = 0.3                                                                               !time after which hadronization is allowed
	tend = 30.0                                                                               !end time      
    if(ispec.eq.1)then
        open (unit=18, file='spec.dat', status='unknown')   
    endif
!*************************** Parameters  ****************************************
	! vector coupling strength (gv/Gs)
    ! 1./3. for Fierz transformation
    ! 2.2/3. for Weise's paper
	ratio = 1./3. 

    ! number of events(number of test particle) (max=1000)
    !	nevent = 1000

    ! number of partitions (max=10)
	npart = 1 

    ! parton cross section (mb)
	cs = 2. 
    sigma=cs/10./npart !<------------------------------------------------------------------ cross section is 10 mb = 1 fm^2

    ! seed # for random numbers
	nseed = irun*10+1
    

	do ixj = 1,Nixj !<-------------------------------------------------------------------- runs begin
	!write(101,*) ixj
	write(*,*) ixj 

	
	
	
	pi=3.14159265
	nc=3

	nt=40                                                                                   ! iteration number for quark mass

    tau=0. !0.05                                                                                ! initial time

    nx=21                                                                                   !changed by sunkj
	ny=21                                                                                   !changed by sunkj 
	nz=21                                                                                   !changed by sunkj 
    dx=qLen*2.0/nx !1.0!5.d-1
	dy=qLen*2.0/ny !1.0!5.d-1
	dz=qLen*2.0/nz !1.0 !5.d-1 


    nxh=int(nx/2.+1.d-8)+1
	nyh=int(ny/2.+1.d-8)+1
	nzh=int(nz/2.+1.d-8)+1

!*********************************************** Lagrangian for NJL Model *************************************************************
! Lagrangian = L_0 + L_sym + L_det + L_Gv + L_GA + L_Gis + L_Giv
! \psi = (u; d; s )
! L_0 = \psi_bar (i\gamma^\mu\partial_\mu - M) \psi
! L_sym = Gs\sum_a=0^8 [(\psi_bar \lambda^a \psi)^2 +(\psi_bar i \gamma_5 \lambda^a \psi)^2]
! L_det = - K [det(\psi_bar(1+\gamma_5)\psi) + det(\psi_bar(1-\gamma_5)\psi)]
! L_Gv = Gv\sum_a=0^8 (\psi_bar \gamma_\mu \lambda^a \psi)^2
! L_GA = GA \sum_a=0^8 (\psi_bar \gamma_\mu \gamma_5 \lambda^a \psi)^2
! L_gviso = gviso(\psi_bar \gamma^\mu \psi)^2
! L_Gis = Gis [(\psi_bar \tau_vec \psi)^2+(\psi \tau_vec \gamma_5 \psi)^2]
! L_Giv = Giv [(\psi_bar \gamma_\mu \tau_vec \psi)^2 +(\psi_bar \gamma_\mu \gamma_5\tau_vec\psi)^2]
! Mean Field Approximation
! Mu* = Mu - 4Gs*uu + 2K*dd*ss -2Gis*(uu-dd)
! Md* = Md - 4Gs*dd + 2K*uu*ss +2Gis*(uu-dd)
! Ms* = Ms - 4Gs*ss + 2K*uu*dd
! i\partial_u\mu* = i\partial_u\mu + 4Gv*\rho_u\mu + 2Giv(\rho_u\mu - rho_d\mu)
! i\partial_d\mu* = i\partial_d\mu + 4Gv*\rho_d\mu - 2Giv(\rho_u\mu - rho_d\mu)
! i\partial_s\mu* = i\partial_s\mu + 4Gv*\rho_s\mu
! L_den = -2Gs(uu^2 + dd^2+ss^2)+4*K*uu*dd*ss -2Gv*(\rho_u\mu^2 +\rho_d\mu^2+\rho_s\mu^2) - Gis*(uu-dd)^2 -Giv*(\rho_u\mu-\rho_d\mu)^2
! uu,dd,ss is the local condensation which are calculated in comoving reference! 
       xm0=5.5d-3 !3.6d-3        ! current mass for u(anti-u),d(anti-d) quarks; GeV
	   xms0=140.7d-3 !87.d-3       ! current mass for s,anti-s quarks ; GeV
       xlam=0.6023 !0.75         ! momentum cut;unit is GeV
       gg= 1.835/xlam**2.  !1.82/xlam**2.   ! gg=Gs     unit = GeV^-2
	   gk=12.36/xlam**5 !8.9/xlam**5.   ! gk = K
       gis=gg/2.0           ! gis = Gis
	   gv= 0.0 !gg !gg/4.  !ratio*(2.*gg)  ! gv = 2/3 Gs = - 2 gv   ;  be careful!!!!
       giv =-5.0*gg !gg/2.          ! giv = -Giv             ;  be careful!!!!    

       gviso= 0.0 !0.5*gg !gg/2.         ! gviso: 

! case Weise
!       xm0=3.6d-3
!	   xms0=87.d-3
!       xlam=0.75 ! unit is GeV
!       gg=1.82/xlam**2. ! in Weise's original paper, gg=3.6
!	   gk=8.9/xlam**5.
!	   gv=ratio*(2.*gg)
!       gis=0.0 !gg !gv !3.*gg !5.*gg        ! by Sunkj
!       giv =0.0 !gg ! gv !2.*gv !gv !gis         ! by sunkj
!       gv = 0.000001

!***************************** this is for vacuum *********************************    
	xmu=xm0
    xmd=xm0
	xms=xms0
    xlam2=xlam**2.
    do i=1, 100                                                         ! 100 iterations to obtain effective quark mass in vacuum 
       alpha=log(xlam/xmu +sqrt((xlam/xmu)**2.+1.))
       xmu2=xmu**2.
	   uu=-xmu*((2.*nc)/(2.*pi)**3.)*pi*xmu2*(sinh(2.*alpha)-2.*alpha)
       
       alpha=log(xlam/xmd +sqrt((xlam/xmd)**2.+1.))
       xmd2=xmd**2.
	   dd=-xmd*((2.*nc)/(2.*pi)**3.)*pi*xmd2*(sinh(2.*alpha)-2.*alpha)       
	 

       alpha=log(xlam/xms +sqrt((xlam/xms)**2.+1.))
       xms2=xms**2.
	   ss=-xms*((2.*nc)/(2.*pi)**3.)*pi*xms2*(sinh(2.*alpha)-2.*alpha)

	   xmu=xm0 -4.*gg*uu +2.*gk*dd*ss-2.*gis*(uu-dd)  ! changed by sunkj
       xmd=xm0 -4.*gg*dd +2.*gk*uu*ss+2.*gis*(uu-dd)  ! changed by sunkj
	   xms=xms0-4.*gg*ss +2.*gk*uu*dd   
    enddo
	uuv=uu ! uu condensate in vacuum
    ddv=dd
    ssv=ss ! ss condensate in vacuum
    
    xmuv=xmu!<========== record the constitute quark mass in vacuum
    xmdv=xmd
    xmsv=xms
    !<==========
    alpha=log(xlam/xm0 +sqrt((xlam/xm0)**2.+1.))
    ubottom=((2.*nc)/(2.*pi)**3.)*pi*xm0*xm0*(sinh(2.*alpha)-2.*alpha) ! for  d

    alpha=log(xlam/xm0 +sqrt((xlam/xm0)**2.+1.))
    dbottom=((2.*nc)/(2.*pi)**3.)*pi*xm0*xm0*(sinh(2.*alpha)-2.*alpha) ! for  d    
    
    alpha=log(xlam/xms0 +sqrt((xlam/xms0)**2.+1.))
    sbottom=((2.*nc)/(2.*pi)**3.)*pi*xms0*xms0*(sinh(2.*alpha)-2.*alpha) ! for s


!************************************ this is for vacuum energy density  ***************************************
    deg=2.*nc
    xmu2=xmu**2.
    xmd2=xmd**2.                                        ! changed by sunkj
    xms2=xms**2.
	Eu0=sqrt(xmu2 +xlam2)
    Ed0=sqrt(xmd2 +xlam2)                                ! changed by sunkj
    Es0=sqrt(xms2 +xlam2)
	veu1=-deg/(16.*pi*pi)*(xlam*(2.*xlam2+xmu2)*Eu0 -xmu2*xmu2*log((Eu0+xlam)/xmu))      ! unit = GeV^4
	ved1=-deg/(16.*pi*pi)*(xlam*(2.*xlam2+xmd2)*Ed0 -xmd2*xmd2*log((Ed0+xlam)/xmd))      !  
	ves1=-deg/(16.*pi*pi)*(xlam*(2.*xlam2+xms2)*Es0 -xms2*xms2*log((Es0+xlam)/xms))

    ve2=2.*gg*(uuv*uuv +ddv*ddv +ssv*ssv) -4.*gk*uuv*ddv*ssv + gis*(uu-dd)**2.            ! unit = GeV^4

    vedn=(veu1 + ved1 +ves1) +ve2                                                         ! total vacuum energy density unit = GeV^4

    call SRAND(nseed) 
	kk=0
    do npara=1, nevent !<----------------------------------------------------- new event
       print *, "event = ", npara     

       iparthd(npara)=0 !<======this is the array recording the # of hadronized partons in each event.
       !<=========Initially, all the partons aren't hadronized, so set all the components to zero.
       k=0
       do npartition=1, npart

          read(17,*) iev, i2,ipart(npara), r1, i4, i5, i6, i7

          do i=1, ipart(npara)

             read(17,*) qx,qy,qz,id,xmass,xx,yy,zz,time
             !if (npara.eq.100)then
             ! print *,id,xx,yy,zz,qx,qy,qz,xmass,time 
             !endif
                k=k+1
	         par(1,k,npara)=xx
	         par(2,k,npara)=yy
	         par(3,k,npara)=zz
             par(4,k,npara)=time
	         par(5,k,npara)=qx
	         par(6,k,npara)=qy
	         par(7,k,npara)=qz
             par(8,k,npara)=xmass
             par(9,k,npara)=float(id)
            !if (npara.eq.100)then
            !  print *,k,npara,par(1,k,npara),par(9,k,npara)
            ! endif

	      enddo ! for particles in one partition
       enddo ! for particles in one event
!    itotal(npara)=ntotal
       ntotal(npara)=k
    enddo ! for events 

        !print *,"last particle:", par(1,4642,100),par(9,4642,100)



 
! *******************************************  Spectators  ************************************************

    if(ispec.eq.1)then
    do nspec =1,nevent
       read(18,*)NNOZPC,ITYPN, GXN, GYN, GZN, FTN,PXN, PYN, PZN, EEN, XMN  ! X Y Z T PX PY PZ E M
       nospec(nspec) = NNOZPC
       spec(1,1,nspec)=ITYPN
       spec(2,1,nspec)=GXN
       spec(3,1,nspec)=GYN
       spec(4,1,nspec)=GZN
       spec(5,1,nspec)=FTN
       spec(6,1,nspec)=PXN
       spec(7,1,nspec)=PYN
       spec(8,1,nspec)=PZN
       spec(9,1,nspec)=EEN
       spec(10,1,nspec)=XMN
       do k=1,NNOZPC-1
           read(18,*)NNOZPC,ITYPN, GXN, GYN, GZN, FTN,PXN, PYN, PZN, EEN, XMN  ! X Y Z T PX PY PZ E M
           spec(1,k+1,nspec)=ITYPN
           spec(2,k+1,nspec)=GXN
           spec(3,k+1,nspec)=GYN
           spec(4,k+1,nspec)=GZN
           spec(5,k+1,nspec)=FTN
           spec(6,k+1,nspec)=PXN
           spec(7,k+1,nspec)=PYN
           spec(8,k+1,nspec)=PZN
           spec(9,k+1,nspec)=EEN
           spec(10,k+1,nspec)=XMN           
       enddo
               
    enddo
    endif

        !print *,"last particle:", par(1,4642,100),par(9,4642,100)
!****** changed, important********************
        itotal2 = 0
        do i=1, nevent  ! changed, important
           itotal(i)=0
        enddo

        do npara=1, nevent
       !mp=itotal(npara)
       np=0
           do i=1, ntotal(npara)
           icoll(i,npara)=0               ! record collision number
           
              mp = i
                         np=np+1
                     q(1,mp,npara)=par(1,i,npara)
                     q(2,mp,npara)=par(2,i,npara)
                     q(3,mp,npara)=par(3,i,npara)
             q(4,mp,npara)=par(4,i,npara)
                     q(5,mp,npara)=par(5,i,npara)
                     q(6,mp,npara)=par(6,i,npara)
                     q(7,mp,npara)=par(7,i,npara)
             q(8,mp,npara)=par(8,i,npara)
                     q(9,mp,npara)=par(9,i,npara)
                     ihd(mp,npara)=0  !<======It's an array denoting whether a parton is hadronized or not. ihd=1 means hadronized.
              ! endif ! if Iglobal = 4
           enddo
        !print *, npara,q(1,mp,npara),par(1,mp,npara),ihd(mp,npara)
       itotal(npara)=itotal(npara)+np
       itotal2=itotal2+np
    enddo

        !print *,"last particle:", par(1,4642,100),par(9,4642,100)
! *************************************** Begin Time Evolution *******************************************
    t=tau                                       ! initial time 
    NumF = 0
    NumFt = 0
    NumOut = 0
    do i=1, nx
	do j=1, ny
	do k=1, nz
       Qmu(i,j,k)=xmu                           ! u quark mass in each cell(in vacuum)
       Qmd(i,j,k)=xmd                           ! d quark mass in each cell(in cacuum)
       Qms(i,j,k)=xms                           ! s quark mass in each cell(in vacuum)
	   At(i,j,k)=0.d0 !<--- phi
	   Ax(i,j,k)=0.d0 !<--- Ax
	   Ay(i,j,k)=0.d0 !<--- Ay
	   Az(i,j,k)=0.d0 !<--- Az
	   Atgv(i,j,k)=0.d0 !<--- phi
	   Axgv(i,j,k)=0.d0 !<--- Ax
	   Aygv(i,j,k)=0.d0 !<--- Ay
	   Azgv(i,j,k)=0.d0 !<--- Az
       
	   Atu(i,j,k)=0.d0 !<--- phi of u              
	   Axu(i,j,k)=0.d0 !<--- Ax  of u             
	   Ayu(i,j,k)=0.d0 !<--- Ay  of u            
	   Azu(i,j,k)=0.d0 !<--- Az  of u           
       Atd(i,j,k)=0.d0 !<--- phi of d
	   Axd(i,j,k)=0.d0 !<--- Ax  of d
	   Ayd(i,j,k)=0.d0 !<--- Ay  of d
	   Azd(i,j,k)=0.d0 !<--- Az  of d
       Ats(i,j,k)=0.d0 !<--- phi of s
	   Axs(i,j,k)=0.d0 !<--- Ax  of s
	   Ays(i,j,k)=0.d0 !<--- Ay  of s
	   Azs(i,j,k)=0.d0 !<--- Az  of s                   


	   pAxu(i,j,k)=0.d0 !<--- previous Ax       ! changed by sunkj
	   pAyu(i,j,k)=0.d0 !<--- Ay                ! changed by sunkj
	   pAzu(i,j,k)=0.d0 !<--- Az                ! changed by sunkj
	   pAxd(i,j,k)=0.d0 !<--- Ax
	   pAyd(i,j,k)=0.d0 !<--- Ay
	   pAzd(i,j,k)=0.d0 !<--- Az
	   pAxs(i,j,k)=0.d0 !<--- Ax
	   pAys(i,j,k)=0.d0 !<--- Ay
	   pAzs(i,j,k)=0.d0 !<--- Az                       
       
	   T00(i,j,k)=0.d0 
	   T01(i,j,k)=0.d0 
	   T02(i,j,k)=0.d0 
	   T03(i,j,k)=0.d0 
	   T11(i,j,k)=0.d0 
	   T12(i,j,k)=0.d0 
	   T13(i,j,k)=0.d0 
	   T22(i,j,k)=0.d0 
	   T23(i,j,k)=0.d0 
	   T33(i,j,k)=0.d0 
       
       scalarU(i,j,k)=0.d0
       scalarD(i,j,k)=0.d0 
       scalarS(i,j,k)=0.d0              
	enddo
	enddo
	enddo    
    !do i=1, nevent
	!   itotal(i)=0
	!enddo

    !itotal2=0
!   dt=0.02*exp(t/1.9)
    !dt=2.d-2                                                        !initial time step
    Nstep = 1                                                       ! record time step
    datatype1='.dat'
    datatype2='EM.dat'
    datatype3='force.dat'
	datatype4='force.dat'
    datatype5='density.dat'
 30 continue !<----------------------------------------------------- new time step
!************************************** New time step *********************************************************

    !call CPU_TIME(t0)
    print *, "t = ", t
    print *, "Nstep = ", Nstep

	!if (t.ge.25)then
	!    icollision=0
	!endif

   !dz = 2*t/(nz-1)
   dzN(Nstep) = dz
   !NcenZ =3 ! floor(t/dz*0.5)                                          ! central cell in z direction
   ncoll=0                                                             ! number of collisions
   dV=dx*dy*dz*5.07**3.                                                                    ! 5.07 = 1/0.19733; GeV^-3   
   TimeN(Nstep)=t     ! record time
!   if (NcenZ<1)then
!        NcenZ = 1
!   endif
!   NcenZ=3                                                          ! need to be delete
     write(tit,"(I4)")Nstep   
     if (ievolution .eq. 1)then 
        open (unit=31, file=tit//datatype1, status='unknown')
        if(iEM.eq.1)then
        open (unit=32, file=tit//datatype2, status='unknown')
        endif
     endif

     !open (unit=41, file='ana/'//tit//datatype3, status='unknown')
	 open (unit=42, file=tit//datatype4, status='unknown')
     open (unit=43, file=tit//datatype5, status='unknown') 
!    do i=1, nevent  ! changed, important
!           itotal(i)=0
!        enddo
!    do npara=1, nevent
!       !mp=itotal(npara)
!       np=0
!	   do i=1, ntotal(npara)
!           icoll(i,npara)=0               ! record collision number
!          ! if ((par(4,i,npara) .gt. t-dt) .and. (par(4,i,npara) .le. t)) then !
!          ! if Iglobal = 4
!              !print *, "np = ", np
!             !mp=mp+1
!              mp = i
!			 np=np+1
!		     q(1,mp,npara)=par(1,i,npara)
!		     q(2,mp,npara)=par(2,i,npara)
!		     q(3,mp,npara)=par(3,i,npara)
!             q(4,mp,npara)=par(4,i,npara)
!		     q(5,mp,npara)=par(5,i,npara)
!		     q(6,mp,npara)=par(6,i,npara)
!		     q(7,mp,npara)=par(7,i,npara)
!            q(8,mp,npara)=par(8,i,npara)
!		     q(9,mp,npara)=par(9,i,npara)
!		     ihd(mp,npara)=0  !<======It's an array denoting whether a parton is hadronized or not. ihd=1 means hadronized.
!	      ! endif ! if Iglobal = 4
!	   enddo
!       itotal(npara)=itotal(npara)+np
!       itotal2=itotal2+np
!    enddo
! ***********************************spin**************************************
     if (ispin.eq.1)then
         do npara=1, nevent
             do i=1, ntotal(npara)
55 continue                 
                 u=(rand()-0.5)*2.
                 v=(rand()-0.5)*2.
                 r2=u**2.+v**2.
                 if (r2.gt.1.)goto 55
                 spin(1,i,npara)=2*u*sqrt(1-r2)                    ! uniform distribution in 3d direction
                 spin(2,i,npara)=2*v*sqrt(1-r2)
                 spin(3,i,npara)=1-2*r2
                 if(abs(par(9,i,npara)).le.1.5)then
                     spin(1,i,npara)=spin(1,i,npara)*(-0.1)    ! unit = fm magnetic moment for d quark
                     spin(2,i,npara)=spin(2,i,npara)*(-0.1)
                     spin(3,i,npara)=spin(3,i,npara)*(-0.1)
                 else if (abs(par(9,i,npara)).le.2.5)then
                     spin(1,i,npara)=spin(1,i,npara)*(0.2)    ! unit = fm  magnetic moment for u quark
                     spin(2,i,npara)=spin(2,i,npara)*(0.2)
                     spin(3,i,npara)=spin(3,i,npara)*(0.2)
                 else
                     spin(1,i,npara)=spin(1,i,npara)*(-0.066)    ! unit = fm  magnetic moment for s quark
                     spin(2,i,npara)=spin(2,i,npara)*(-0.066)
                     spin(3,i,npara)=spin(3,i,npara)*(-0.066)
                 endif
                 
                     
             enddo
         enddo
         
         endif

	dt=0.2 !*exp(t/1.9)!--------------------------------------------- time step---------------                  
    if (iglobal.lt.4) then
    if (itotal2 .lt. 1) then    ! if there is no particle freezed out 
       t=t+dt
    	   goto 30
    	endif
    endif


! goto 77
!********************************************************************** potential

	do i=1, nx
	   do j=1, ny
	      do k=1, nz
             pxs(i,j,k)=0.d0
             pys(i,j,k)=0.d0
             pzs(i,j,k)=0.d0
          enddo
       enddo
	enddo

!    print *,"testing:",ihd(1,1)
! step 0 : sum momentum in each fluid cell
    do npara=1, nevent
	   do i=1, itotal(npara) 
          if(ihd(i,npara).eq.1) goto 32 !xjfeng
	      mx=int((q(1,i,npara)+nx*dx/2.)/dx) +1      ! x[-0.5 0.5), mx =22
	      my=int((q(2,i,npara)+ny*dy/2.)/dy) +1      ! int(1.1)=1,int(1.5)=1,int(-1.1)=-1,int(-1.5)=-1
	      mz=int((q(3,i,npara)+nz*dz/2.)/dz) +1
          if ((mx .lt. 1) .or. (mx .gt. nx)) goto 32
          if ((my .lt. 1) .or. (my .gt. ny)) goto 32
          if ((mz .lt. 1) .or. (mz .gt. nz)) goto 32
          pxs(mx,my,mz)= pxs(mx,my,mz) +q(5,i,npara)    ! totoal momentum in each cell : unit = GeV
          pys(mx,my,mz)= pys(mx,my,mz) +q(6,i,npara)
          pzs(mx,my,mz)= pzs(mx,my,mz) +q(7,i,npara)
 32 continue
       enddo
    enddo    
  

! this is for scalar mean-fields

    it=0
 31 continue

    !print *, "time = ", t, "  Nstep = ", Nstep
    !print *, "interation: ", " mass in central cell = ", xmudcen
! step 1 : initializing
	do i=1, nx
	   do j=1, ny
	      do k=1, nz

			 Edd(i,j,k)=0.d0
			 Euu(i,j,k)=0.d0
			 Ess(i,j,k)=0.d0
             EdnS(i,j,k)=0.d0
             EdnV(i,j,k)=0.d0
             EdnK(i,j,k)=0.d0
             EdnT(i,j,k)=0.d0
			 EddL(i,j,k)=0.d0
			 EuuL(i,j,k)=0.d0
			 EssL(i,j,k)=0.d0
             EdnKL(i,j,k)=0.d0
             EdnTL(i,j,k)=0.d0             
             EdnEM(i,j,k)=0.d0  ! Electro-magnetic energy density
             EdnEMP(i,j,k)=0.d0 ! interaction energy density
             Exf(i,j,k)=0.d0
             Eyf(i,j,k)=0.d0
             Ezf(i,j,k)=0.d0
             Bxf(i,j,k)=0.d0
             Byf(i,j,k)=0.d0
             Bzf(i,j,k)=0.d0
             
             specEx(i,j,k)=0.d0
             specEy(i,j,k)=0.d0
             specEz(i,j,k)=0.d0
             specBx(i,j,k)=0.d0
             specBy(i,j,k)=0.d0
             specBz(i,j,k)=0.d0
             
             dd0(i,j,k)=0.d0
             uu0(i,j,k)=0.d0
             ss0(i,j,k)=0.d0

          enddo
       enddo
	enddo


    !print *, "time = ", t, "  Nstep = ", Nstep
    !print *, "step 1: ", " mass in central cell = ", xmudcen
! step 2 : sum energy in each fluid cell
    do npara=1, nevent
	   do i=1, itotal(npara)
             ! print *, "Ihadronization = ", ihd(i,npara)
          if(ihd(i,npara).eq.1) goto 33 !<===== if a parton is hadronized, it won't contribute to the mean field.
	      !print *, "mx = ", mx, " my = ", my, " mz = ", mz
              mx=int((q(1,i,npara)+nx*dx/2.)/dx) +1
	      my=int((q(2,i,npara)+ny*dy/2.)/dy) +1
	      mz=int((q(3,i,npara)+nz*dz/2.)/dz) +1
          if ((mx .lt. 1) .or. (mx .gt. nx)) goto 33
          if ((my .lt. 1) .or. (my .gt. ny)) goto 33
          if ((mz .lt. 1) .or. (mz .gt. nz)) goto 33
          q2=q(5,i,npara)**2. +q(6,i,npara)**2. +q(7,i,npara)**2.
          if (abs(q(9,i,npara)) .lt. 1.5) then		       ! for u, d quarks
             Edd(mx,my,mz)= Edd(mx,my,mz) +sqrt(q(8,i,npara)**2.+q2)            ! total energy : unit = GeV
          else if (abs(q(9,i,npara)) .lt. 2.5) then		       ! for u, d quarks
             Euu(mx,my,mz)= Euu(mx,my,mz) +sqrt(q(8,i,npara)**2.+q2)
		  else                                             ! for s quarks
             Ess(mx,my,mz)= Ess(mx,my,mz) +sqrt(q(8,i,npara)**2.+q2)
          endif
 33 continue
       enddo
	enddo

! step 3 : find fluid velocity
	do i=1, nx
	   do j=1, ny
	      do k=1, nz
             if (Edd(i,j,k)+Euu(i,j,k)+Ess(i,j,k) .gt. 1.d-9) then
                vxs(i,j,k)= pxs(i,j,k)/(Edd(i,j,k)+Euu(i,j,k)+Ess(i,j,k))      ! cell velocity : unit = 1
                vys(i,j,k)= pys(i,j,k)/(Edd(i,j,k)+Euu(i,j,k)+Ess(i,j,k))
                vzs(i,j,k)= pzs(i,j,k)/(Edd(i,j,k)+Euu(i,j,k)+Ess(i,j,k))
             endif
          enddo
       enddo
	enddo    
    
    
    
! step 4 : move to CM frame
    do npara=1, nevent
	   do i=1, itotal(npara)
         if(ihd(i,npara).eq.1) goto 34 !<=====if a parton is hadronized, it won't contribute to the mean field.
	      mx=int((q(1,i,npara)+nx*dx/2.)/dx) +1
	      my=int((q(2,i,npara)+ny*dy/2.)/dy) +1
	      mz=int((q(3,i,npara)+nz*dz/2.)/dz) +1
          if ((mx .lt. 1) .or. (mx .gt. nx)) goto 34
          if ((my .lt. 1) .or. (my .gt. ny)) goto 34
          if ((mz .lt. 1) .or. (mz .gt. nz)) goto 34
          q2=q(5,i,npara)**2. +q(6,i,npara)**2. +q(7,i,npara)**2.
          E1=sqrt(q(8,i,npara)**2.+q2)
          bx=vxs(mx,my,mz)
		  by=vys(mx,my,mz)
		  bz=vzs(mx,my,mz)
          gam = 1./sqrt(1.-bx**2.-by**2.-bz**2.)
          call xLorentz(E1,q(5,i,npara),q(6,i,npara),q(7,i,npara), bx,by,bz, E2,px,py,pz)

          q2= px*px +py*py +pz*pz
		  if (q2 .gt. xlam2) goto 34

          if (abs(q(9,i,npara)) .lt. 1.5) then
             dd0(mx,my,mz)= dd0(mx,my,mz) +1./E1/dV                 ! unit = GeV ^2
          else if (abs(q(9,i,npara)) .lt. 2.5) then
             uu0(mx,my,mz)= uu0(mx,my,mz) +1./E1/dV
          else
             ss0(mx,my,mz)= ss0(mx,my,mz) +1./E1/dV
		  endif
34        continue
       enddo
	enddo    

	xmuucen = 0.d0                   ! average quark mass in central cell
	xmudcen = 0.d0 
	xmuscen = 0.d0 
    totalS=0.
	centerS=0.
! step 5 : determine condensate and scalar mean-fields
	do i=1, nx
	   do j=1, ny
	      do k=1, nz

             if (Edd(i,j,k) .gt. 1.d-9) then
			    alpha=log(xlam/Qmd(i,j,k) +sqrt((xlam/Qmd(i,j,k))**2.+1.))
                dmass2=Qmd(i,j,k)**2.
			    constI=((2.*nc)/(2.*pi)**3.)*pi*dmass2*(sinh(2.*alpha)-2.*alpha)
		        dd= -Qmd(i,j,k) *(constI -dd0(i,j,k)/npart/nevent)
                if (dd .gt. 0.d0) dd=0.d0
			 else
			    dd=ddv                              !  quark mass in vacuum
			 endif

             if (Euu(i,j,k) .gt. 1.d-9) then
			    alpha=log(xlam/Qmu(i,j,k) +sqrt((xlam/Qmu(i,j,k))**2.+1.))
                umass2=Qmu(i,j,k)**2.
			    constI=((2.*nc)/(2.*pi)**3.)*pi*umass2*(sinh(2.*alpha)-2.*alpha)
		        uu= -Qmu(i,j,k) *(constI -uu0(i,j,k)/npart/nevent)                 ! unit = GeV^3
                if (uu .gt. 0.d0) uu=0.d0
            else
			    uu=uuv                            !  quark mass in vacuum
            endif

             if (Ess(i,j,k) .gt. 1.d-9) then
			    alpha=log(xlam/Qms(i,j,k) +sqrt((xlam/Qms(i,j,k))**2.+1.))
                smass2=Qms(i,j,k)**2.
			    constI=((2.*nc)/(2.*pi)**3.)*pi*smass2*(sinh(2.*alpha)-2.*alpha)
		        ss= -Qms(i,j,k) *(constI -ss0(i,j,k)/npart/nevent)
                if (ss .gt. 0.d0) ss=0.d0
			 else
			    ss=ssv                    !  quark mass in vacuum
			 endif

	         Qmd(i,j,k)= xm0 -4.*gg*dd +2.*gk*uu*ss + 2.*gis*(uu-dd)   ! effective mass of d quark: unit =  GeV
	         Qmu(i,j,k)= xm0 -4.*gg*uu +2.*gk*dd*ss - 2.*gis*(uu-dd)   ! 
	         Qms(i,j,k)= xms0-4.*gg*ss +2.*gk*uu*dd                    !
                 !Qmu(i,j,k) = Qms(i,j,k)  !test
                 !Qmd(i,j,k) = Qms(i,j,k)  ! test
                 !Qms(i,j,k) = Qmu(i,j,k)    !test
				 
             if (ichiral.eq.0)then
                 Qmu(i,j,k) = xmuv-0.1  !test
                 Qmd(i,j,k) = xmdv-0.1  ! test
                 Qms(i,j,k) = xmsv-0.1    !test
                 elseif (ichiral.eq.2)then
                 Qmu(i,j,k) = xm0
                 Qmd(i,j,k) = xm0
                 Qms(i,j,k) = xms0
              endif
			  
             if (Qmd(i,j,k)<xm0) then                                  ! Added by sunkj, suggested by Liewen Chen
                 Qmd(i,j,k) = xm0                                      ! which is not consistent~~
             endif
             
             if (Qmu(i,j,k)<xm0) then                                ! Added by sunkj
                 Qmu(i,j,k) = xm0
             endif
             
             if (Qms(i,j,k)<xms0) then                                ! Added by sunkj
                 Qms(i,j,k) = xms0
             endif
             
             
             scalarU(i,j,k) = uu                                         ! condensation of u and anti-u: unit=GeV^3 
             scalarD(i,j,k) = dd                                         ! 
             scalarS(i,j,k) = ss
         if ((i .ge. nxh-NCenX) .and. (i .le. nxh+NCenX)) then
	      if ((j .ge. nyh-NCenY) .and. (j .le. nyh+NCenY)) then
	      if ((k .ge. nzh-NCenZ) .and. (k .le. nzh+NCenZ)) then
			 xmudcen = xmudcen + Qmd(i,j,k) !d quark mass in central cells : unit = GeV
			 xmuucen = xmuucen + Qmu(i,j,k) !u quark mass in central cells
			 xmuscen = xmuscen + Qms(i,j,k) !s quark mass in central cells
	      endif
	      endif
	      endif

! this is to calculate energy density from scalar mean fields

             deg=2.*nc
			 xmd2=Qmd(i,j,k)**2.
			 xmu2=Qmu(i,j,k)**2.
			 xms2=Qms(i,j,k)**2.
	         Ed0=sqrt(xmd2 +xlam2)
	         Eu0=sqrt(xmu2 +xlam2)
             Es0=sqrt(xms2 +xlam2)
	         ved1=-deg/(16.*pi*pi)*(xlam*(2.*xlam2+xmd2)*Ed0 -xmd2*xmd2*log((Ed0+xlam)/sqrt(xmd2)))  !GeV^4
	         veu1=-deg/(16.*pi*pi)*(xlam*(2.*xlam2+xmu2)*Eu0 -xmu2*xmu2*log((Eu0+xlam)/sqrt(xmu2)))
	         ves1=-deg/(16.*pi*pi)*(xlam*(2.*xlam2+xms2)*Es0 -xms2*xms2*log((Es0+xlam)/sqrt(xms2)))

             ve2=2.*gg*(uu*uu +dd*dd +ss*ss) -4.*gk*uu*dd*ss + gis*(uu-dd)**2.     ! GeV^4

             pot=(ved1 +veu1 +ves1)  -vedn  +ve2
             totalS=totalS +pot                                                    ! unit = GeV^4
             EdnS(i,j,k) = pot*5.07**3.                                            ! unit = GeV/fm^3     
             if ((i .ge. nxh-NCenX) .and. (i .le. nxh+NCenX)) then        ! changed by sunkj
	         if ((j .ge. nyh-NCenY) .and. (j .le. nyh+NCenY)) then       ! changed by sunkj
	         if ((k .ge. nzh-NCenZ) .and. (k .le. nzh+NCenZ)) then       ! changed by sunkj
                centerS=centerS +pot
	         endif
	         endif
	         endif

          enddo
       enddo
	enddo

	xmudcen = xmudcen/((2*NCenX+1.)*(2*NCenY+1.0)*(2*NCenZ+1.)) !d quark mass in central cells! changed by sunkj  :unit = GeV
	xmuucen = xmuucen/((2*NCenX+1.)*(2*NCenY+1.0)*(2*NCenZ+1.)) !u quark mass in central cells! changed by sunkj
	xmuscen = xmuscen/((2*NCenX+1.)*(2*NCenY+1.0)*(2*NCenZ+1.)) !s quark mass in central cells! changed by sunkj

! step 6 : update particle mass
    summas=0.d0
	sumnum=0.d0
    do npara=1, nevent
	   do i=1, itotal(npara)
          
          if(ihd(i,npara).eq.1) goto 35 !<=====if a parton is hadronized, its mass won't be affected by mean field.
	      mx=int((q(1,i,npara)+nx*dx/2.)/dx) +1
	      my=int((q(2,i,npara)+ny*dy/2.)/dy) +1
	      mz=int((q(3,i,npara)+nz*dz/2.)/dz) +1
          if ((mx .lt. 1) .or. (mx .gt. nx)) goto 35
          if ((my .lt. 1) .or. (my .gt. ny)) goto 35 
          if ((mz .lt. 1) .or. (mz .gt. nz)) goto 35            
          if (abs(q(9,i,npara)) .lt. 1.5) then
		     q(8,i,npara)= Qmd(mx,my,mz)
			 summas=summas+Qmd(mx,my,mz)
			 sumnum=sumnum+1.d0
          else if (abs(q(9,i,npara)) .lt. 2.5) then
		     q(8,i,npara)= Qmu(mx,my,mz)
			 summas=summas+Qmu(mx,my,mz)
			 sumnum=sumnum+1.d0
		  else
		     q(8,i,npara)= Qms(mx,my,mz)
          endif       
 35 continue
       enddo
	enddo

    it=it+1
	if (it .le. nt) goto 31
!-------------------------------------******************iteration ends**********************--------------------------------!    
! step 7 : update total momentum & energy in each cell AFTER interation

	do i=1, nx
	   do j=1, ny
	      do k=1, nz
			 Edd(i,j,k)=0.d0
			 Euu(i,j,k)=0.d0
			 Ess(i,j,k)=0.d0
          enddo
       enddo
    enddo
    
    do npara=1, nevent
	   do i=1, itotal(npara)

          if(ihd(i,npara).eq.1) goto 39 !<=====if a parton is hadronized, it won't contribute to the mean field.
	      mx=int((q(1,i,npara)+nx*dx/2.)/dx) +1
	      my=int((q(2,i,npara)+ny*dy/2.)/dy) +1
	      mz=int((q(3,i,npara)+nz*dz/2.)/dz) +1

          if ((mx .lt. 1) .or. (mx .gt. nx)) goto 39
          if ((my .lt. 1) .or. (my .gt. ny)) goto 39
          if ((mz .lt. 1) .or. (mz .gt. nz)) goto 39

          q2=q(5,i,npara)**2. +q(6,i,npara)**2. +q(7,i,npara)**2.

          if (abs(q(9,i,npara)) .lt. 1.5) then		       ! for u, d quarks
             Edd(mx,my,mz)= Edd(mx,my,mz) +sqrt(q(8,i,npara)**2.+q2)
          else if (abs(q(9,i,npara)) .lt. 2.5) then		   ! for u, d quarks
             Euu(mx,my,mz)= Euu(mx,my,mz) +sqrt(q(8,i,npara)**2.+q2)
		  else                                             ! for s quarks
             Ess(mx,my,mz)= Ess(mx,my,mz) +sqrt(q(8,i,npara)**2.+q2)
          endif

 39 continue

       enddo
    enddo    
    

    
! this is for vector mean-fields
!------------------------------ central cells---------------------------------------
	dentotal = 0.d0                           !total number density in central cells
	dennet = 0.d0                             !net number density in central cells
    scalD = 0.d0
    scalU = 0.d0
    scalS = 0.d0
    vecDt = 0.d0
    vecDx = 0.d0
    vecDy = 0.d0
    vecDz = 0.d0
    vecUt = 0.d0
    vecUx = 0.d0
    vecUy = 0.d0
    vecUz = 0.d0
    vecSt = 0.d0
    vecSx = 0.d0
    vecSy = 0.d0
    vecSz = 0.d0    
    
    totalE=0.d0
    energyden=0.d0

    
	do i=1, nx
	   do j=1, ny
	      do k=1, nz
			 rhobd(i,j,k)=0.d0           ! baryon density of d quark
			 rhobu(i,j,k)=0.d0 
             rhobs(i,j,k)=0.d0
			 rhobxu(i,j,k)=0.d0          ! current baryon density of u quark
			 rhobyu(i,j,k)=0.d0
             rhobzu(i,j,k)=0.d0
			 rhobxd(i,j,k)=0.d0
			 rhobyd(i,j,k)=0.d0
             rhobzd(i,j,k)=0.d0 
			 rhobxs(i,j,k)=0.d0
			 rhobys(i,j,k)=0.d0
             rhobzs(i,j,k)=0.d0  
             numbd(i,j,k)=0
             numbu(i,j,k)=0
             numbs(i,j,k)=0
             numb(i,j,k)=0
             gamLV(i,j,k)=1.d0
			 forceux(i,j,k) = 0.d0
			 forceuy(i,j,k) = 0.d0
			 forceuz(i,j,k) = 0.d0
          enddo
       enddo
    enddo 
    if (iEM.eq.1)then
	    do i=1, nx
	       do j=1, ny
	          do k=1, nz
                 EMcrt(Nstep,i,j,k)=0.0
                 EMcrx(Nstep,i,j,k)=0.0
                 EMcry(Nstep,i,j,k)=0.0
                 EMcrz(Nstep,i,j,k)=0.0
                 NEMnu(Nstep,i,j,k)=0.0
                 NEMnd(Nstep,i,j,k)=0.0
                 NEMns(Nstep,i,j,k)=0.0
                 NEMnub(Nstep,i,j,k)=0.0
                 NEMndb(Nstep,i,j,k)=0.0
                 NEMnsb(Nstep,i,j,k)=0.0
              enddo
           enddo
        enddo   
     endif
! step 2 : move to no current frame
    do npara=1, nevent
	   do i=1, itotal(npara)
          if(ihd(i,npara).eq.1) goto 36 !<=====if a parton is hadronized, it won't contribute to the mean field.

	      mx=int((q(1,i,npara)+nx*dx/2.)/dx) +1
	      my=int((q(2,i,npara)+ny*dy/2.)/dy) +1
	      mz=int((q(3,i,npara)+nz*dz/2.)/dz) +1
          if ((mx .lt. 1) .or. (mx .gt. nx)) then
              print *,"out of range mx"
             goto 36
             endif
          if ((my .lt. 1) .or. (my .gt. ny)) goto 36 
          if ((mz .lt. 1) .or. (mz .gt. nz)) goto 36  

            if (iEM.eq.1)then       
		    id1=q(9,i,npara)
		    x1=q(1,i,npara)
		    y1=q(2,i,npara)
		    z1=q(3,i,npara)
		    px1=q(5,i,npara)
		    py1=q(6,i,npara)
		    pz1=q(7,i,npara)
		    pm1=q(8,i,npara)
            Em1=pm1**2.+px1**2.+py1**2.+pz1**2.
            Em1=sqrt(Em1)
		    vx1=px1/Em1
		    vy1=py1/Em1
		    vz1=pz1/Em1
            
            EMcrt(Nstep,mx,my,mz)=EMcrt(Nstep,mx,my,mz) + id1/abs(id1)*1/3.
            EMcrx(Nstep,mx,my,mz)=EMcrx(Nstep,mx,my,mz) + id1/abs(id1)*1/3.*vx1
            EMcry(Nstep,mx,my,mz)=EMcry(Nstep,mx,my,mz) + id1/abs(id1)*1/3.*vy1
            EMcrz(Nstep,mx,my,mz)=EMcrz(Nstep,mx,my,mz) + id1/abs(id1)*1/3.*vz1

		    endif             
          
          
          
          q2=q(5,i,npara)**2. +q(6,i,npara)**2. +q(7,i,npara)**2.
          E1=sqrt(q(8,i,npara)**2.+q2)

          bx=vxs(mx,my,mz) 
		  by=vys(mx,my,mz) 
		  bz=vzs(mx,my,mz) 
          gamLV(mx,my,mz) = 1./sqrt(1.-bx**2.-by**2.-bz**2.)
          call xLorentz(E1,q(5,i,npara),q(6,i,npara),q(7,i,npara), bx,by,bz, E2,px,py,pz)

          q2= px*px +py*py +pz*pz
		  if (q2 .gt. xlam2) then
                print *, "momentum out of range:"  
                goto 36

             endif

!<======count the charge current density in the fireball frame.
!          chab(mx,my,mz)= chab(mx,my,mz) +q(9,i,npara)/abs(q(9,i,npara))
!          pxv(mx,my,mz)= pxv(mx,my,mz) +q(5,i,npara)/E1*q(9,i,npara)/abs(q(9,i,npara))
!          pyv(mx,my,mz)= pyv(mx,my,mz) +q(6,i,npara)/E1*q(9,i,npara)/abs(q(9,i,npara))
!          pzv(mx,my,mz)= pzv(mx,my,mz) +q(7,i,npara)/E1*q(9,i,npara)/abs(q(9,i,npara))








          if (abs(q(9,i,npara)) .lt. 1.5) then
  
              rhobd(mx,my,mz)= rhobd(mx,my,mz) +q(9,i,npara)/abs(q(9,i,npara))
              rhobxd(mx,my,mz)= rhobxd(mx,my,mz) +q(5,i,npara)/E1*q(9,i,npara)/abs(q(9,i,npara))
              rhobyd(mx,my,mz)= rhobyd(mx,my,mz) +q(6,i,npara)/E1*q(9,i,npara)/abs(q(9,i,npara))
              rhobzd(mx,my,mz)= rhobzd(mx,my,mz) +q(7,i,npara)/E1*q(9,i,npara)/abs(q(9,i,npara))    
              numbd(mx,my,mz)= numbd(mx,my,mz) +1  
              if (iEM.eq.1)then
              if(q(9,i,npara).gt.0.)then
                  NEMnd(Nstep,mx,my,mz)=NEMnd(Nstep,mx,my,mz)+1
              else 
                  NEMndb(Nstep,mx,my,mz)=NEMndb(Nstep,mx,my,mz)+1
              endif
              endif
          else if (abs(q(9,i,npara)) .lt. 2.5) then
              rhobu(mx,my,mz)= rhobu(mx,my,mz) +q(9,i,npara)/abs(q(9,i,npara))    
              rhobxu(mx,my,mz)= rhobxu(mx,my,mz) +q(5,i,npara)/E1*q(9,i,npara)/abs(q(9,i,npara))
              rhobyu(mx,my,mz)= rhobyu(mx,my,mz) +q(6,i,npara)/E1*q(9,i,npara)/abs(q(9,i,npara))
              rhobzu(mx,my,mz)= rhobzu(mx,my,mz) +q(7,i,npara)/E1*q(9,i,npara)/abs(q(9,i,npara))
              numbu(mx,my,mz)= numbu(mx,my,mz) +1 
              if(iEM.eq.1)then
              if(q(9,i,npara).gt.0.)then
                  NEMnu(Nstep,mx,my,mz)=NEMnu(Nstep,mx,my,mz)+1
              else 
                  NEMnub(Nstep,mx,my,mz)=NEMnub(Nstep,mx,my,mz)+1
              endif          
              endif
          else
              rhobs(mx,my,mz)= rhobs(mx,my,mz) +q(9,i,npara)/abs(q(9,i,npara))
              rhobxs(mx,my,mz)= rhobxs(mx,my,mz) +q(5,i,npara)/E1*q(9,i,npara)/abs(q(9,i,npara))
              rhobys(mx,my,mz)= rhobys(mx,my,mz) +q(6,i,npara)/E1*q(9,i,npara)/abs(q(9,i,npara))
              rhobzs(mx,my,mz)= rhobzs(mx,my,mz) +q(7,i,npara)/E1*q(9,i,npara)/abs(q(9,i,npara))
              numbs(mx,my,mz)= numbs(mx,my,mz) +1  
              if(iEM.eq.1)then
              if(q(9,i,npara).gt.0.)then
                  NEMns(Nstep,mx,my,mz)=NEMns(Nstep,mx,my,mz)+1
              else 
                  NEMnsb(Nstep,mx,my,mz)=NEMnsb(Nstep,mx,my,mz)+1
              endif      
              endif
          endif
          numb(mx,my,mz)= numb(mx,my,mz) +1

          q2=q(5,i,npara)**2. +q(6,i,npara)**2. +q(7,i,npara)**2.
          Eqst=sqrt(q2 +q(8,i,npara)**2.)
          if (abs(q(9,i,npara)) .lt. 1.5) then		       ! for u, d quarks
             EddL(mx,my,mz)= EddL(mx,my,mz) + E2   !Eqst
          else if (abs(q(9,i,npara)) .lt. 2.5) then		   ! for u, d quarks
             EuuL(mx,my,mz)= EuuL(mx,my,mz) + E2   !Eqst
		  else                                             ! for s quarks
             EssL(mx,my,mz)= EssL(mx,my,mz) + E2   !Eqst
          endif          
          
	      
          totalE=totalE +Eqst

          if ((mx .ge. nxh-NCenX) .and. (mx .le. nxh+NCenX)) then
	      if ((my .ge. nyh-NCenY) .and. (my .le. nyh+NCenY)) then
	      if ((mz .ge. nzh-NCenZ) .and. (mz .le. nzh+NCenZ)) then
             ! print *,i,npara,q(1,i,npara),q(2,i,npara),q(3,i,npara),q(4,i,npara),q(5,i,npara)
               energyden=energyden+Eqst
			 dentotal = dentotal + 1. 
			 dennet = dennet + q(9,i,npara)/abs(q(9,i,npara))
              !print *,i,npara,q(1,i,npara),q(9,i,npara)
              !print *,i,npara,par(1,i,npara),par(9,i,npara)
	      endif
	      endif
	      endif          
          
          
 36 continue

       enddo
    enddo
    if (iEM.eq.1)then
        
	    do i=1, nx
	       do j=1, ny
	          do k=1, nz
                 EMcrt(Nstep,i,j,k)=EMcrt(Nstep,i,j,k)/(Nevent+0.)/(dx*dy*dz)
                 EMcrx(Nstep,i,j,k)=EMcrx(Nstep,i,j,k)/(Nevent+0.)/(dx*dy*dz)
                 EMcry(Nstep,i,j,k)=EMcry(Nstep,i,j,k)/(Nevent+0.)/(dx*dy*dz)
                 EMcrz(Nstep,i,j,k)=EMcrz(Nstep,i,j,k)/(Nevent+0.)/(dx*dy*dz)
              enddo
           enddo
        enddo   

    
    

    do i=1, nx
	   do j=1, ny
	      do k=1, nz
              do mx=1,nx
                  do my=1,ny
                      do mz=1,nz
                          r1x = (i-nxh)*dx
                          r1y = (j-nyh)*dy
                          r1z = (k-nzh)*dz
                          r2x = (mx-nxh)*dx
                          r2y = (my-nyh)*dy
                          r2z = (mz-nzh)*dz
                          
                          rx= r2x-r1x
                          ry= r2y-r1y
                          rz= r2z-r1z
                          r2=rx**2.+ry**2.+rz**2.
                          a=(2*mz/(nz-1.))**2.-1.
                          b=4*mz*rz/(nz-1.)
                          c=r2
                          detm = b**2.-4.*a*c
                          if(detm.ge.0.)then
                             if(a.eq.0.)then
                                 if(b.eq.0.)then
                                     NstepI=0
                                 else
                                     dtm = -c/b
                                 endif
                              else   
                             dtm1 = (-b+sqrt(b**2.-4*a*c))/(2*a)
                             dtm2 = (-b-sqrt(b**2.-4*a*c))/(2*a)
                             if((dtm1.ge.dtm2).and.(dtm2.ge.0.))dtm=dtm2
                             if((dtm2.ge.dtm1).and.(dtm1.ge.0.))dtm=dtm1
                             if((dtm1.ge.0.).and.(dtm2.le.0.))dtm=dtm1
                             if((dtm1.le.0.).and.(dtm2.ge.0.))dtm=dtm2
                              endif
                              
                             do nt=1,Nstep
                                if (t-dtm.lt.TimeN(nt))then
                                    NstepI=nt
                                    goto 37
                                endif
                            enddo                              
                              
                          else
                          NstepI=0
                          endif
                          

37                        continue           
                       if(NstepI.ge.1)then   
                          rx= r2x-r1x
                          ry= r2y-r1y
                          r2z = (mz-nzh)*dzN(NstepI)
                          rz=r2z-r1z                        
                          r2=rx**2.+ry**2.+rz**2.
                                         
          if((abs(rx).gt.dx*1.5).and.(abs(ry).gt.dy*1.5).and.(abs(rz).gt.dzN(NstepI)*1.5))then
                             
                          EMJt=EMcrt(NstepI,mx,my,mz)
                          EMJx=EMcrx(NstepI,mx,my,mz)
                          EMJy=EMcry(NstepI,mx,my,mz)
                          EMJz=EMcrz(NstepI,mx,my,mz)

           call EMFieldcell(rx,ry,rz,dx*dy*dzN(NstepI),EMJt,EMJx,EMJy,EMJz,Ex,Ey,Ez,Bx,By,Bz)                   
                            Exf(i,j,k)=Exf(i,j,k)+Ex
                            Eyf(i,j,k)=Eyf(i,j,k)+Ey
                            Ezf(i,j,k)=Ezf(i,j,k)+Ez
                            Bxf(i,j,k)=Bxf(i,j,k)+Bx
                            Byf(i,j,k)=Byf(i,j,k)+By
                            Bzf(i,j,k)=Bzf(i,j,k)+Bz
          endif
          endif
                      enddo
                  enddo
              enddo
          enddo
       enddo
    enddo 
    
    
    call CPU_TIME(t0)
    print *, "t0 = ", t0
    if(ispec.eq.1)then
        do i=1, nx
	       do j=1, ny
	          do k=1, nz     
                    ndx=i-nxh
                    ndy=j-nyh
                    ndz=k-nzh
                    r1x = ndx*dx
                    r1y = ndy*dy
                    r1z = ndz*dz
                     do nspec =1,1   !nevent   in order to speed up
                        do k1=1,nospec(nspec)
                 if(abs(spec(1,k1,nspec)).eq.2212)then                 ! proton without neutron
                     sign = spec(1,k1,nspec)/abs(spec(1,k1,nspec))
                     
                         
                     
                     r2x = spec(2,k1,nspec)
                     r2y = spec(3,k1,nspec)
                     
                     p2x = spec(6,k1,nspec)
                     p2y = spec(7,k1,nspec)
                     p2z = spec(8,k1,nspec)
                     E2  = spec(9,k1,nspec)
                     v2x = p2x/E2
                     v2y = p2y/E2
                     v2z = p2z/E2
                     gamma1=1/sqrt(1-v2z**2.)
                     gamma2=E2/spec(10,k1,nspec)
                     r2z = spec(4,k1,nspec) + v2z*t
                     Rx = r1x-r2x
                     Ry = r1y-r2y
                     Rz = r1z-r2z

                     call EMFieldspec(Rx,Ry,Rz,v2x,v2y,v2z,sign,Ex,Ey,Ez,Bx,By,Bz)
                           specEx(i,j,k)=specEx(i,j,k)+Ex
                           specEy(i,j,k)=specEy(i,j,k)+Ey
                           specEz(i,j,k)=specEz(i,j,k)+Ez
                           specBx(i,j,k)=specBx(i,j,k)+Bx
                           specBy(i,j,k)=specBy(i,j,k)+By
                           specBz(i,j,k)=specBz(i,j,k)+Bz
                 endif          
                        enddo
                     enddo
                     
               enddo
          enddo
       enddo         
        
    endif
    

    call CPU_TIME(t1)
    print *,"time elpsed:",t1-t0
    if (ispec.eq.1)then
        do i=1, nx
	       do j=1, ny
	          do k=1, nz 
                  Exf(i,j,k)=Exf(i,j,k)+specEx(i,j,k) !/nevent
                  Eyf(i,j,k)=Eyf(i,j,k)+specEy(i,j,k) !/nevent
                  Ezf(i,j,k)=Ezf(i,j,k)+specEz(i,j,k) !/nevent
                  Bxf(i,j,k)=Bxf(i,j,k)+specBx(i,j,k) !/nevent
                  Byf(i,j,k)=Byf(i,j,k)+specBy(i,j,k) !/nevent
                  Bzf(i,j,k)=Bzf(i,j,k)+specBz(i,j,k) !/nevent
               enddo
          enddo
       enddo         
    endif
    
	do i=1, nx
	   do j=1, ny
	      do k=1, nz
           !EdnEMP(i,j,k)=EdnEMP(i,j,k)/(dx*dy*dz)                                   ! GeV/fm^3
           EdnEM(i,j,k)=1/(8.*pi)*(1./137.)*(Exf(i,j,k)**2.+Eyf(i,j,k)**2. &
                +Ezf(i,j,k)**2.+Bxf(i,j,k)**2.+Byf(i,j,k)**2.+Bzf(i,j,k)**2.)*0.19733     ! GeV/fm^3

           ealpha=1/137.
           hbarc=0.19733
           Exf(i,j,k)=Exf(i,j,k)*ealpha*hbarc                           !GeV/fm  eE
           Eyf(i,j,k)=Eyf(i,j,k)*ealpha*hbarc                           !GeV/fm  eE
           Ezf(i,j,k)=Ezf(i,j,k)*ealpha*hbarc                            !GeV/fm  eE
           Bxf(i,j,k)=Bxf(i,j,k)*ealpha*hbarc                            !GeV/fm  eB
           Byf(i,j,k)=Byf(i,j,k)*ealpha*hbarc                            !GeV/fm  eB
           Bzf(i,j,k)=Bzf(i,j,k)*ealpha*hbarc                           !GeV/fm  eB
          enddo
       enddo
    enddo     

    endif
              
     
    totalV=0.
	centerV=0
    totalVu=0.            ! changed by sunkj
	centerVu=0.             ! changed by sunkj
    totalVd=0.              ! changed by sunkj
	centerVd=0.             ! changed by sunkj
    totalVs=0.              ! changed by sunkj
	centerVs=0.             ! changed by sunkj

! step 3 : go back to Lab. frame and determine vector-mean fields
	do i=1, nx
	   do j=1, ny
	      do k=1, nz
              mx=i
              my=j
              mz=k
              if ((mx .ge. nxh-NCenX) .and. (mx .le. nxh+NCenX)) then
	          if ((my .ge. nyh-NCenY) .and. (my .le. nyh+NCenY)) then
	          if ((mz .ge. nzh-NCenZ) .and. (mz .le. nzh+NCenZ)) then
                 scalD = scalD + scalarD(mx,my,mz)
                 scalU = scalU + scalarU(mx,my,mz)
                 scalS = scalS + scalarS(mx,my,mz)                 
                 vecDt = vecDt + rhobd(mx,my,mz)
                 vecDx = vecDx + rhobxd(mx,my,mz)
                 vecDy = vecDy + rhobyd(mx,my,mz)
                 vecDz = vecDz + rhobzd(mx,my,mz)
                 vecUt = vecUt + rhobu(mx,my,mz)
                 vecUx = vecUx + rhobxu(mx,my,mz)
                 vecUy = vecUy + rhobyu(mx,my,mz)
                 vecUz = vecUz + rhobzu(mx,my,mz) 
                 vecSt = vecSt + rhobs(mx,my,mz)
                 vecSx = vecSx + rhobxs(mx,my,mz)
                 vecSy = vecSy + rhobys(mx,my,mz)
                 vecSz = vecSz + rhobzs(mx,my,mz)                  
	          endif
	          endif
              endif              
 
              
           
              

              
              
              
              
              
              
              
              
             rhobu(i,j,k) = rhobu(i,j,k)/dV/npart/nevent                                         ! current baryon density
             rhobxu(i,j,k) = rhobxu(i,j,k)/dV/npart/nevent                                       ! unit = GeV^3
             rhobyu(i,j,k) = rhobyu(i,j,k)/dV/npart/nevent
             rhobzu(i,j,k) = rhobzu(i,j,k)/dV/npart/nevent  
             rhobd(i,j,k) = rhobd(i,j,k)/dV/npart/nevent
             rhobxd(i,j,k) = rhobxd(i,j,k)/dV/npart/nevent
             rhobyd(i,j,k) = rhobyd(i,j,k)/dV/npart/nevent
             rhobzd(i,j,k) = rhobzd(i,j,k)/dV/npart/nevent                   
             rhobs(i,j,k) = rhobs(i,j,k)/dV/npart/nevent
             rhobxs(i,j,k) = rhobxs(i,j,k)/dV/npart/nevent
             rhobys(i,j,k) = rhobys(i,j,k)/dV/npart/nevent
             rhobzs(i,j,k) = rhobzs(i,j,k)/dV/npart/nevent 
             
             
             
			 comp0=(rhobd(i,j,k)+rhobu(i,j,k)+rhobs(i,j,k))                         ! GeV^3
			 comp1= (rhobxd(i,j,k)+rhobxu(i,j,k)+rhobxs(i,j,k))
			 comp2= (rhobyd(i,j,k)+rhobyu(i,j,k)+rhobys(i,j,k))
			 comp3= (rhobzd(i,j,k)+rhobzu(i,j,k)+rhobzs(i,j,k))

             At(i,j,k)= 4.*gv*comp0                                                     ! GeV
             Ax(i,j,k)= 4.*gv*comp1
             Ay(i,j,k)= 4.*gv*comp2
             Az(i,j,k)= 4.*gv*comp3
             Atgv(i,j,k)=2.*gviso*comp0
             Axgv(i,j,k)=2.*gviso*comp1
             Aygv(i,j,k)=2.*gviso*comp2
             Azgv(i,j,k)=2.*gviso*comp3
             Atu(i,j,k)= (4.*gv*rhobu(i,j,k) + 2.*giv*(rhobu(i,j,k)-rhobd(i,j,k)))  +Atgv(i,j,k)      ! be careful. the sign , and coefficient
             Axu(i,j,k)= (4.*gv*rhobxu(i,j,k) +2.*giv*(rhobxu(i,j,k)-rhobxd(i,j,k)))+Axgv(i,j,k)       
             Ayu(i,j,k)= (4.*gv*rhobyu(i,j,k)+2.*giv*(rhobyu(i,j,k)-rhobyd(i,j,k))) +Aygv(i,j,k)     
             Azu(i,j,k)= (4.*gv*rhobzu(i,j,k)+2.*giv*(rhobzu(i,j,k)-rhobzd(i,j,k))) +Azgv(i,j,k)      
             Atd(i,j,k)= (4.*gv*rhobd(i,j,k) - 2.*giv*(rhobu(i,j,k)-rhobd(i,j,k)))  +Atgv(i,j,k)    
             Axd(i,j,k)= (4.*gv*rhobxd(i,j,k) -2.*giv*(rhobxu(i,j,k)-rhobxd(i,j,k)))+Axgv(i,j,k)       
             Ayd(i,j,k)= (4.*gv*rhobyd(i,j,k)-2.*giv*(rhobyu(i,j,k)-rhobyd(i,j,k))) +Aygv(i,j,k)       
             Azd(i,j,k)= (4.*gv*rhobzd(i,j,k)-2.*giv*(rhobzu(i,j,k)-rhobzd(i,j,k))) +Azgv(i,j,k)          
             Ats(i,j,k)= 4.*gv*rhobs(i,j,k) +Atgv(i,j,k)                                          
             Axs(i,j,k)= 4.*gv*rhobxs(i,j,k)+Axgv(i,j,k)                                 
             Ays(i,j,k)= 4.*gv*rhobys(i,j,k)+Aygv(i,j,k)                                           
             Azs(i,j,k)= 4.*gv*rhobzs(i,j,k)+Azgv(i,j,k)                                          


 
! this is to calculate energy density from vector mean fields

             !totalV=totalV +Ut(i,j,k)**2./(2.*gv) 
             vednU = -2.*gv*(rhobu(i,j,k)**2.-rhobxu(i,j,k)**2.-rhobyu(i,j,k)**2.-rhobzu(i,j,k)**2.)       ! unit GeV^4/N^2              
             vednD = -2.*gv*(rhobd(i,j,k)**2.-rhobxd(i,j,k)**2.-rhobyd(i,j,k)**2.-rhobzd(i,j,k)**2.)       ! unit GeV^4/N^2
             vednS = -2.*gv*(rhobs(i,j,k)**2.-rhobxs(i,j,k)**2.-rhobys(i,j,k)**2.-rhobzs(i,j,k)**2.)       ! unit GeV^4/N^2
             
             cellednVu =  Atu(i,j,k)*rhobu(i,j,k)                                                        ! unit: GeV^4/N^2 A0 
             cellednVd =  Atd(i,j,k)*rhobd(i,j,k)                                                        ! unit: GeV^4/N^2
             cellednVs =  Ats(i,j,k)*rhobs(i,j,k)                                                        ! unit: GeV^4/N^2
             cellednVu = cellednVu  -(Axu(i,j,k)*rhobxu(i,j,k)+Ayu(i,j,k)*rhobyu(i,j,k)+Azu(i,j,k)*rhobzu(i,j,k))
             cellednVd = cellednVd  -(Axd(i,j,k)*rhobxd(i,j,k)+Ayd(i,j,k)*rhobyd(i,j,k)+Azd(i,j,k)*rhobzd(i,j,k))
             cellednVs = cellednVs  -(Axs(i,j,k)*rhobxs(i,j,k)+Ays(i,j,k)*rhobys(i,j,k)+Azs(i,j,k)*rhobzs(i,j,k))
             
             cellednVgviso=Atgv(i,j,k)*comp0
             cellednVgviso=cellednVgviso-(Axgv(i,j,k)*comp1+Aygv(i,j,k)*comp2+Azgv(i,j,k)*comp3)
             totalV = totalV+vednU + vednD+vednS+cellednVu+cellednVd+cellednVs+ &                        ! total energy density of vector part
             -giv*((rhobu(i,j,k)-rhobd(i,j,k))**2. -(rhobxu(i,j,k)-rhobxd(i,j,k))**2.- &
                (rhobyu(i,j,k)-rhobyd(i,j,k))**2.-(rhobzu(i,j,k)-rhobzd(i,j,k))**2. )-&
             gviso*(comp0**2.-comp1**2.-comp2**2.-comp3**2.)
             
             EdnV(i,j,k) =  vednU + vednD+vednS+cellednVu+cellednVd+cellednVs+ &                         ! energy density of vector part in each cell
             -giv*((rhobu(i,j,k)-rhobd(i,j,k))**2. -(rhobxu(i,j,k)-rhobxd(i,j,k))**2.- &                  ! unit = GeV^4
                (rhobyu(i,j,k)-rhobyd(i,j,k))**2.-(rhobzu(i,j,k)-rhobzd(i,j,k))**2. )-&
             gviso*(comp0**2.-comp1**2.-comp2**2.-comp3**2.)
             EdnV(i,j,k) = EdnV(i,j,k)*5.07**3.0                                                         ! unit = GeV/fm^3
             EdnK(i,j,k) = (Edd(i,j,k)+ Euu(i,j,k) +Ess(i,j,k))/(dx*dy*dz)/npart/nevent 
             
             EdnKL(i,j,k) = (EddL(i,j,k)+ EuuL(i,j,k) +EssL(i,j,k))/(dx*dy*dz*gamLV(i,j,k))/npart/nevent  ! Lorentz factor gamLV: unit = GeV/fm^3
             EdnT(i,j,k) = EdnS(i,j,k) + EdnV(i,j,k) + EdnK(i,j,k)
             EdnTL(i,j,k) = EdnS(i,j,k) + EdnV(i,j,k) + EdnKL(i,j,k)
             if(iEM.eq.1)then
                 EdnTL(i,j,k)=EdnTL(i,j,k)+EdnEM(i,j,k)+EdnEMP(i,j,k)
             endif
             
             if ((i .ge. nxh-NCenX) .and. (i .le. nxh+NCenX)) then
	         if ((j .ge. nyh-NCenY) .and. (j .le. nyh+NCenY)) then
	         if ((k .ge. nzh-NCenZ) .and. (k .le. nzh+NCenZ)) then
                centerV=centerV + vednU + vednD+vednS+cellednVu+cellednVd+cellednVs+ &
             -giv*((rhobu(i,j,k)-rhobd(i,j,k))**2. -(rhobxu(i,j,k)-rhobxd(i,j,k))**2.- &
                (rhobyu(i,j,k)-rhobyd(i,j,k))**2.-(rhobzu(i,j,k)-rhobzd(i,j,k))**2. )-&
             gviso*(comp0**2.-comp1**2.-comp2**2.-comp3**2.)
              if(iEM.eq.1)then
                 centerV=centerV+EdnEMP(i,j,k)+EdnEM(i,j,k)
             endif               
                
	         endif
	         endif
             endif



          enddo
       enddo
    enddo

    
    if((ievolution.eq.1).and.(iEM.eq.1))then
	    do i=1, nx
	       do j=1, ny
	          do k=1, nz   
                  write(32,2010)t, Exf(i,j,k),Eyf(i,j,k),Ezf(i,j,k),Bxf(i,j,k),Byf(i,j,k),Bzf(i,j,k), &
                                 EdnTL(i,j,k)
              enddo
           enddo
        enddo       
        
    endif
    close(32)
2010 format(1x,8(1x,F14.6))    
    
	   ! do k=1, ny   
           !write(41,*)numbu(nxh,k,:)/float(nevent)
        !enddo 
	!close(41)
	 !   do k=1, nx  
           !write(42,2010)t,dx,numbu(k,nyh,nzh)/float(nevent),numbd(k,nyh,nzh)/float(nevent), &
	!   numbs(k,nyh,nzh)/float(nevent)
        !enddo
        !close (42)
        !print *, "recording: ",numbu(11,11,:)
        !print *, "recording: ",rhobu(11,11,:)
        !print *, "recording: ",numbd(11,11,:)
        !print *, "recording:",rhobd(11,11,:)
       ! do k=1, nx 
          !print *, numbu(k,:,nzh) 
                   !do ky = 1,ny
        !   write(43,*)numbu(k,:,nzh)/float(nevent)
                   !enddo
        !enddo
        do i=1, nx
           do j=1, ny
              do k=1, nz

                write(43,17) rhobd(i,j,k),rhobu(i,j,k)
				
				!write(42,19) forceux(i,j,k),forceuy,forceuz(i,j,k)
                 enddo
       enddo
       enddo         

        close(43)

        print *, "boundary check",Nstep,rhobu(1,1,1)*dV,rhobu(1,2,1)*dV,rhobu(1,2,3)*dV
        print *, "boundary check",Nstep,rhobd(1,1,1)*dV,rhobd(1,2,1)*dV,rhobd(1,2,3)*dV
        print *, "boundary check",Nstep,Atu(1,1,1),Atu(1,2,1),Atu(1,2,3)
        print *, "boundary check",Nstep,Atd(1,1,1),Atd(1,2,1),Atd(1,2,3)


        ntotald = 0
        rhototald = 0.
        do k=1,nx
           do k1=1,ny
              do k2=1,nz

                ntotald = ntotald + numbd(k,k1,k2)
                rhototald = rhototald + rhobd(k,k1,k2)
                enddo
            enddo
        enddo
        print *,"total number of d: ",ntotald
        print *, "total density of d: ",rhototald*5.07**3.
 
	
!********************************************************************** cascade
77 continue    
    
!------------------------------------------Freeze Out---------------------------------    
    if(iglobal.eq.0)then
    Nfreeze = 0
	do npara=1, nevent
       do i=1, itotal(npara)-1
          if(ihd(i,npara).eq.1) goto 81  
	      mx=int((q(1,i,npara)+nx*dx/2.)/dx) +1
	      my=int((q(2,i,npara)+ny*dy/2.)/dy) +1
	      mz=int((q(3,i,npara)+nz*dz/2.)/dz) +1
          if ((mx .lt. 1) .or. (my .lt. 1) .or. (mz .lt. 1)) then
             
             print*,"out of range: ",q(1,i,npara)
             Nfreeze =1
	         goto 80
		  endif

	      if ((mx .gt. nx) .or. (my .gt. ny) .or. (mz .gt. nz)) then
             Nfreeze =1
			 goto 80
          endif
 
          
 !<============Determine the hadronization mass by flavor============
          if (abs(q(9,i,npara)).lt.2.5) then
            xmhd= 0.3   !0.5*xmuv
          else
            xmhd= 0.45  !0.5*xmsv
          end if
!<==================================================================

!<=============Once a parton is hadronized, record its information in array qq=====
		
        if(iglobal.eq.0) then 
          !if((q(8,i,npara).gt.xmhd).and.(ihd(i,npara).eq.0).and.t.gt.thd) then
          xNtot = (numbd(mx,my,mz)+numbu(mx,my,mz)+numbs(mx,my,mz))/(dx*dy*dz)/nevent
          if (((EdnTL(mx,my,mz).lt.0.6)).and.(ihd(i,npara).eq.0).and.t.gt.thd) then   !.and.(xNtot<2.)
               NumFt=NumFt+1
              write(109,*)q(9,i,npara),q(8,i,npara),EdnTL(mx,my,mz),xNtot,EdnTL(mx,my,mz)/xNtot
            if (Qmd(mx,my,mz)<0.05)then
                !if (Nfreeze.eq.0)then
                    NumF = NumF +1
!                    print *,"mass of d quark= ",EdnTL(mx,my,mz),Vd(mx,my,mz),numbd(mx,my,mz),numbu(mx,my,mz)
!                    print *,"condesation is : ",scalarD(mx,my,mz),scalarU(mx,my,mz),xm0 -4.*gg*scalarD(mx,my,mz)
!                    print *,"number density:" , numb(mx,my,mz)/(dx*dy*dz)/nevent
                !endif
            endif
80          continue	 
           
            iparthd(npara)=iparthd(npara)+1 !<====# of hadronized partons +1            
            do idx=1, 9
                qq(idx,i,npara)=q(idx,i,npara)
            end do
            qq(4,i,npara)=t-dt
!            if (qq(8,i,npara)<0.05)then
!            print *,"id ",qq(9,i,npara),"mass: ",qq(8,i,npara),par(8,i,npara)
!            print *,Vd(mx,my,mz),Vu(mx,my,mz),Vs(mx,my,mz)
!            print *, mx,my,mz,i,npara
!            print *, q(1,i,npara),q(2,i,npara),q(3,i,npara)
!            print *,int((q(1,i,npara)+nx*dx/2.)/dx) +1,int((q(2,i,npara)+ny*dy/2.)/dy) +1, &
!                    int((q(3,i,npara)+nz*dz/2.)/dz) +1
!            endif

            ihd(i,npara)=1
            
          end if
        endif
        Nfreeze = 0         
          
81          continue	          
 !<==============================================================    
       enddo
    enddo    
    endif
!-------------------------------------------------------Freeze Out end----------------------------------------------------    
   

    
	do i=1, nx
	   do j=1, ny
	      do k=1, nz

!			 if (numb(i,j,k) .lt. 1) numb(i,j,k)=1
             if (numbd(i,j,k) .lt. 1) numbd(i,j,k)=1  ! added by sunkj  Position is crucial important
             if (numbu(i,j,k) .lt. 1) numbu(i,j,k)=1  ! added by sunkj
             if (numbs(i,j,k) .lt. 1) numbs(i,j,k)=1  ! added by sunkj
          enddo
       enddo
    enddo
        if (icollision.eq.0)goto 71

	do npara=1, nevent
       do i=1, itotal(npara)-1
	           if(ihd(i,npara).eq.1) goto 70       
	      E1=sqrt(q(8,i,npara)**2. +q(5,i,npara)**2. +q(6,i,npara)**2. +q(7,i,npara)**2.)

	      do iq=i+1, itotal(npara)
		          if(ihd(iq,npara).eq.1) goto 60 
             dis=sqrt((q(1,i,npara)-q(1,iq,npara))**2.+(q(2,i,npara)-q(2,iq,npara))**2.+(q(3,i,npara)-q(3,iq,npara))**2.)
             if (dis .ge. 3.*sqrt(sigma/pi)) goto 60
             if (dis .lt. 1.d-9) goto 60 !<--- to avoid the successive collision of two particles
	         E2=sqrt(q(8,iq,npara)**2. +q(5,iq,npara)**2. +q(6,iq,npara)**2. +q(7,iq,npara)**2.)


             E=E1+E2
	         qx=q(5,i,npara)+q(5,iq,npara)
	         qy=q(6,i,npara)+q(6,iq,npara)
	         qz=q(7,i,npara)+q(7,iq,npara)

             bx=qx/E
	         by=qy/E
	         bz=qz/E

             beta2=bx*bx +by*by +bz*bz
	         gamma=1./sqrt(1.-beta2)

             call xLorentz(E1,q(5,i,npara), q(6,i,npara), q(7,i,npara),  bx,by,bz, Ei, px, py, pz)
	         call xLorentz(t, q(1,i,npara), q(2,i,npara), q(3,i,npara),  bx,by,bz, ti, xi, yi, zi )
             call xLorentz(t, q(1,iq,npara),q(2,iq,npara),q(3,iq,npara), bx,by,bz, tiq,xiq,yiq,ziq)


             s=E*E -qx*qx -qy*qy -qz*qz
	         pr2=(s-(q(8,i,npara)+q(8,iq,npara))**2.)*(s-(q(8,i,npara)-q(8,iq,npara))**2.)/(4.*s)
             pr2=max(pr2,1.d-9)
	         pr=sqrt(pr2)

             vxi=px/Ei
	         vyi=py/Ei
	         vzi=pz/Ei

	         Eiq=sqrt(q(8,iq,npara)**2.+px*px +py*py +pz*pz)
	         vxiq=-px/Eiq
	         vyiq=-py/Eiq
	         vziq=-pz/Eiq

	         xi0=xi-vxi*ti
	         yi0=yi-vyi*ti
	         zi0=zi-vzi*ti

	         xiq0=xiq-vxiq*tiq
	         yiq0=yiq-vyiq*tiq
	         ziq0=ziq-vziq*tiq

	         xnum=(vxi-vxiq)*(xi0-xiq0) +(vyi-vyiq)*(yi0-yiq0) +(vzi-vziq)*(zi0-ziq0)
	         deno=(vxi-vxiq)**2. +(vyi-vyiq)**2. +(vzi-vziq)**2.
             tc=-xnum/deno

             d2=(xi0-xiq0)**2.+(yi0-yiq0)**2.+(zi0-ziq0)**2. -xnum**2./deno

	         xic=xi0+vxi*tc
	         yic=yi0+vyi*tc
	         zic=zi0+vzi*tc

	         xiqc=xiq0+vxiq*tc
	         yiqc=yiq0+vyiq*tc
	         ziqc=ziq0+vziq*tc

             xc=(xic+xiqc)/2.
             yc=(yic+yiqc)/2.
             zc=(zic+ziqc)/2.

	         call xLorentz(tc,xc,yc,zc, -bx,-by,-bz, tc0,xc0,yc0,zc0)

             if ((d2 .lt. sigma/pi) .and. (tc0 .lt. t+dt) .and. (tc0 .gt. t)) then

	            t1=rand()
		        t2=rand()

		        cs=(t1-0.5)*2.
		        sn=sqrt(1.-cs*cs)
		        ph=t2*(2.*pi)

                Ecm1=sqrt(q(8,i,npara)**2. +pr2)
                Ecm2=sqrt(q(8,iq,npara)**2. +pr2)
                px=pr*sn*cos(ph)
		        py=pr*sn*sin(ph)
		        pz=pr*cs

		        call xLorentz(Ecm1,px,py,pz, -bx,-by,-bz, E1,qx1,qy1,qz1)

                

		        call xLorentz(Ecm2,-px,-py,-pz, -bx,-by,-bz, E2,qx2,qy2,qz2)

				 mx1=int((xc0+nx*dx/2.)/dx) +1
				 my1=int((yc0+ny*dy/2.)/dy) +1
				 mz1=int((zc0+nz*dz/2.)/dz) +1

				if ((mx1 .ge. 1) .and. (mx1 .le. nx).and.(my1 .ge. 1) .and. (my1 .le. ny).and.(mz1.ge. 1) .and. (mz1 .le. nz))then
				
				rho = (rhobd(mx1,my1,mz1)+rhobu(mx1,my1,mz1))/2.
				p1 = (qx1**2.+qy1**2.+qz1**2.)**(1./2.)
				pfermi = (pi**2.*rho)**(1./3.)
			        !print *, p1,pfermi	
				! if (p1.le.pfermi*3./4.)goto 60
				endif
				
				if ((mx2 .ge. 1) .and. (mx2 .le. nx).and.(my2 .ge. 1) .and. (my2 .le. ny).and.(mz2.ge. 1) .and. (mz2 .le. nz))then
				
				rho = (rhobd(mx2,my2,mz2)+rhobu(mx2,my2,mz2))/2.
				p2 = (qx2**2.+qy2**2.+qz2**2.)**(1./2.)
				pfermi = (pi**2.*rho)**(1./3.)
				
				!if (p2.le.pfermi*3./4.)goto 60
				endif
				
				
				!write(*,*)"p1  before cll:",q(1,i,npara),q(2,i,npara),q(3,i,npara)
				!write(*,*)"p1  after cll:",xc0,yc0,zc0
				!write(*,*)"p2  before cll:",q(1,iq,npara),q(2,iq,npara),q(3,iq,npara)
				!write(*,*)"p2  after cll:",xc0,yc0,zc0
				
                        !qx1 = q(5,iq,npara) ! by kaijia sun, just exchange p
                        !qy1 = q(6,iq,npara) ! by kaijia sun, just exchange p
                        !qz1 = q(7,iq,npara) ! by kaijia sun, just exchange p

                        !qx2 = q(5,i,npara) ! by kaijia sun, just exchange p
                        !qy2 = q(6,i,npara) ! by kaijia sun, just exchange p
                        !qz2 = q(7,i,npara) ! by kaijia sun, just exchange p

                        
                        q(5,i,npara)=qx1
 	                q(6,i,npara)=qy1
		        q(7,i,npara)=qz1
		        !q(1,i,npara)=xc0 
		        !q(2,i,npara)=yc0 
		        !q(3,i,npara)=zc0 
				
		        q(5,iq,npara)=qx2
		        q(6,iq,npara)=qy2
		        q(7,iq,npara)=qz2
		        !q(1,iq,npara)=xc0 
		        !q(2,iq,npara)=yc0 
		        !q(3,iq,npara)=zc0 

		        ncoll=ncoll+1
                icoll(i,npara)=icoll(i,npara)+1
                icoll(iq,npara)=icoll(iq,npara)+1
		        goto 70

             endif
 60 continue
	      enddo
70        continue           
       enddo
    enddo

71        continue
!********************************************************************** next time step

!*********************************************************************Boundary mass matrix***********************
	call xBound(Qmu,nx,ny,nz,BoundQmu)     
	call xBound(Qmd,nx,ny,nz,BoundQmd)
	call xBound(Qms,nx,ny,nz,BoundQms)
	
	call xBound(Atu,nx,ny,nz,BoundAtu)
	call xBound(Atd,nx,ny,nz,BoundAtd)
	call xBound(Ats,nx,ny,nz,BoundAts)
	
	call xBound(Axd,nx,ny,nz,BoundAxd)
	call xBound(Ayd,nx,ny,nz,BoundAyd)
	call xBound(Azd,nx,ny,nz,BoundAzd)
	
	call xBound(Axu,nx,ny,nz,BoundAxu)
	call xBound(Ayu,nx,ny,nz,BoundAyu)
	call xBound(Azu,nx,ny,nz,BoundAzu)
	call xBound(Axs,nx,ny,nz,BoundAxs)
	call xBound(Ays,nx,ny,nz,BoundAys)
	call xBound(Azs,nx,ny,nz,BoundAzs)
	
	call xBound(pAxd,nx,ny,nz,BoundpAxd)
	call xBound(pAyd,nx,ny,nz,BoundpAyd)
	call xBound(pAzd,nx,ny,nz,BoundpAzd)
	call xBound(pAxu,nx,ny,nz,BoundpAxu)
	call xBound(pAyu,nx,ny,nz,BoundpAyu)
	call xBound(pAzu,nx,ny,nz,BoundpAzu)
	call xBound(pAxs,nx,ny,nz,BoundpAxs)
	call xBound(pAys,nx,ny,nz,BoundpAys)
	call xBound(pAzs,nx,ny,nz,BoundpAzs)
        print *, "boundary matrix", Qmu(nxh,nyh,nzh),BoundQmu(nxh+1,nyh+1,nzh+1)
        print *, "boundary matrix", Qmd(nxh,nyh,nzh),BoundQmd(nxh+1,nyh+1,nzh+1)
        print *, "boundary matrix", Qms(nxh,nyh,nzh),BoundQms(nxh+1,nyh+1,nzh+1)
        print *, "boundary matrix", Axd(nxh,nyh,nzh),BoundAxd(nxh+1,nyh+1,nzh+1)
        print *, "boundary matrix", Ayd(nxh,nyh,nzh),BoundAyd(nxh+1,nyh+1,nzh+1)
        print *, "boundary matrix", Azd(nxh,nyh,nzh),BoundAzd(nxh+1,nyh+1,nzh+1)
        !print *, "boundary check",Nstep,rhobu(10,11,11)*dV,rhobu(12,11,11)*dV,Atu(10,11,11),Atu(12,11,11)
        !print *, "boundary check",Nstep,rhobd(10,11,11)*dV,rhobd(12,11,11)*dV,Atd(10,11,11),Atd(12,11,11)
        !print *, "boundary check",Nstep,BoundAtu(1+1,1+1,1+1),BoundAtu(1+1,2+1,1+1),BoundAtu(1+1,2+1,3+1)
        !print *, "boundary check",Nstep,BoundAtd(1+1,1+1,1+1),BoundAtd(1+1,2+1,1+1),BoundAtd(1+1,2+1,3+1)

	
	nxMaxBound = 0
	nyMaxBound = 0
	nxMinBound = 0
	nyMinBound = 0
	

        do mx=-1, nx+1
           do my=-1, ny+1
              do mz=-1, nz+1
                          
                          if (inewtonPA.eq.0) then
                          !BoundpAtd(mx+1,my+1,mz+1) = BoundAtd(mx+1,my+1,mz+1)
                          !BoundpAtu(mx+1,my+1,mz+1) = BoundAtu(mx+1,my+1,mz+1)
                          !BoundpAts(mx+1,my+1,mz+1) = BoundAts(mx+1,my+1,mz+1)                          

                          BoundpAxd(mx+1,my+1,mz+1) = BoundAxd(mx+1,my+1,mz+1)
                          BoundpAyd(mx+1,my+1,mz+1) = BoundAyd(mx+1,my+1,mz+1)
                          BoundpAzd(mx+1,my+1,mz+1) = BoundAzd(mx+1,my+1,mz+1)
                          BoundpAxu(mx+1,my+1,mz+1) = BoundAxu(mx+1,my+1,mz+1)
                          BoundpAyu(mx+1,my+1,mz+1) = BoundAyu(mx+1,my+1,mz+1)
                          BoundpAzu(mx+1,my+1,mz+1) = BoundAzu(mx+1,my+1,mz+1)
                          BoundpAxs(mx+1,my+1,mz+1) = BoundAxs(mx+1,my+1,mz+1)
                          BoundpAys(mx+1,my+1,mz+1) = BoundAys(mx+1,my+1,mz+1)
                          BoundpAzs(mx+1,my+1,mz+1) = BoundAzs(mx+1,my+1,mz+1)
               endif
             enddo
           enddo
        enddo


!print *, "boundary check:",Nstep,InewtonPA,(BoundAtd(1+1,2+1,3+1)-BoundAtd(1+1,2+1,1+1))/2./dz,forceuz(1,2,2),BoundAzu(1+1,2+1,2+1),BoundpAzu(2+1,2+1,2+1)

	
	do mx=1, nx
           do my=1, ny
              do mz=1, nz
		
						  
				Fux = -(BoundAxu(mx+1,my+1,mz+1)-BoundpAxu(mx+1,my+1,mz+1))/dt - (BoundAtu(mx+1+1,my+1,mz+1)-BoundAtu(mx-1+1,my+1,mz+1))/(2.*dx)
				Fuy = -(BoundAyu(mx+1,my+1,mz+1)-BoundpAyu(mx+1,my+1,mz+1))/dt - (BoundAtu(mx+1,my+1+1,mz+1)-BoundAtu(mx+1,my-1+1,mz+1))/(2.*dy)
				Fuz = -(BoundAzu(mx+1,my+1,mz+1)-BoundpAzu(mx+1,my+1,mz+1))/dt - (BoundAtu(mx+1,my+1,mz+1+1)-BoundAtu(mx+1,my+1,mz-1+1))/(2.*dz)
				Fuxs =-(BoundQmd(mx+1+1,my+1,mz+1)-BoundQmd(mx-1+1,my+1,mz+1))/(2.*dx) !scalar
				Fuys =-(BoundQmd(mx+1,my+1+1,mz+1)-BoundQmd(mx+1,my-1+1,mz+1))/(2.*dy)
				Fuzs =-(BoundQmd(mx+1,my+1,mz+1+1)-BoundQmd(mx+1,my+1,mz-1+1))/(2.*dz)
				
				forceux(mx,my,mz) = Fux
				forceuy(mx,my,mz) = Fuy
				forceuz(mx,my,mz) = Fuz
				write(42,19) Fux,Fuy,Fuz,Fuxs,Fuys,Fuzs
              enddo
       enddo
    enddo 
	!print *, "boundary check",Nstep,InewtonPA,(BoundAtd(1+1,2+1,3+1)-BoundAtd(1+1,2+1,1+1))/2./dz,forceuz(1,2,2),BoundAzu(1+1,2+1,2+1),BoundpAzu(1+1,2+1,2+1)
	!do i=1, nx
    !       do j=1, ny
    !          do k=1, nz
    !				write(42,19) forceux(i,j,k),forceuy(i,j,k),forceuz(i,j,k)
    !          enddo
    !   enddo
    !enddo         

        close(42)
!**********************************************end*************************
	
    do npara=1, nevent
          if(ievolution.eq.1)then
            write(31,15) npara,itotal(npara),bim,i2,i3,ls(i,npara),lp(i,npara)
          endif        
	   do i=1, itotal(npara)
          if(ihd(i,npara).eq.1) goto 78 !<=====if a parton is hadronized, it won't be affected by mean field.

          px2=q(5,i,npara)
	      py2=q(6,i,npara)
	      pz2=q(7,i,npara)

	      mx=int((q(1,i,npara)+nx*dx/2.)/dx) +1
	      my=int((q(2,i,npara)+ny*dy/2.)/dy) +1
	      mz=int((q(3,i,npara)+nz*dz/2.)/dz) +1

!goto 78

          !if ((mx .le. 1) .or. (my .le. 1) .or. (mz .le. 1)) then
!		     print *, "lower"
!print *, i, q(1,i,npara), q(2,i,npara), q(3,i,npara)
            ! EdnPar = 0
            ! Nfreeze =1
             if ((mx .lt. 1) .or. (my .lt. 1) .or. (mz .lt. 1)) then
             NumOut = NumOut +1
             endif
	         !goto 78
			 
		  !endif

	     ! if ((mx .ge. nx) .or. (my .ge. ny) .or. (mz .ge. nz)) then
!		     print *, "upper"
!print *, i, q(1,i,npara), q(2,i,npara), q(3,i,npara)
             !EdnPar = 0  ! total energy density of the cell.
             !Nfreeze =1
	      if ((mx .gt. nx) .or. (my .gt. ny) .or. (mz .gt. nz)) then             
             NumOut = NumOut +1
             endif
			 !goto 78
          !endif

  
          EdnPar=EdnTL(mx,my,mz)
          
	      Eqst=sqrt(q(5,i,npara)**2. +q(6,i,npara)**2. +q(7,i,npara)**2. +q(8,i,npara)**2.)

! momentum update

             bx=vxs(mx,my,mz)
             by=vys(mx,my,mz)
			 bz=vzs(mx,my,mz)


		  call xLorentz(Eqst,px2,py2,pz2, bx,by,bz, s0,s1,s2,s3)

		  q2=s1**2.+s2**2.+s3**2.
		  if (q2 .gt. xlam**2.) then
                     print *,'momentum out of range'  
                     goto 78
                  endif 

                  if (inewton.eq.0) goto 78
			  velx = px2/Eqst    ! quark velocity
			  vely = py2/Eqst 
			  velz = pz2/Eqst
                          
			  Eldx = -(BoundAxd(mx+1,my+1,mz+1)-BoundpAxd(mx+1,my+1,mz+1))/dt - (BoundAtd(mx+1+1,my+1,mz+1)-BoundAtd(mx-1+1,my+1,mz+1))/(2.*dx)
			  Eldy = -(BoundAyd(mx+1,my+1,mz+1)-BoundpAyd(mx+1,my+1,mz+1))/dt - (BoundAtd(mx+1,my+1+1,mz+1)-BoundAtd(mx+1,my-1+1,mz+1))/(2.*dy)
			  Eldz = -(BoundAzd(mx+1,my+1,mz+1)-BoundpAzd(mx+1,my+1,mz+1))/dt - (BoundAtd(mx+1,my+1,mz+1+1)-BoundAtd(mx+1,my+1,mz-1+1))/(2.*dz)
			  Elux = -(BoundAxu(mx+1,my+1,mz+1)-BoundpAxu(mx+1,my+1,mz+1))/dt - (BoundAtu(mx+1+1,my+1,mz+1)-BoundAtu(mx-1+1,my+1,mz+1))/(2.*dx)
			  Eluy = -(BoundAyu(mx+1,my+1,mz+1)-BoundpAyu(mx+1,my+1,mz+1))/dt - (BoundAtu(mx+1,my+1+1,mz+1)-BoundAtu(mx+1,my-1+1,mz+1))/(2.*dy)
			  Eluz = -(BoundAzu(mx+1,my+1,mz+1)-BoundpAzu(mx+1,my+1,mz+1))/dt - (BoundAtu(mx+1,my+1,mz+1+1)-BoundAtu(mx+1,my+1,mz-1+1))/(2.*dz)
			  Elsx = -(BoundAxs(mx+1,my+1,mz+1)-BoundpAxs(mx+1,my+1,mz+1))/dt - (BoundAts(mx+1+1,my+1,mz+1)-BoundAts(mx-1+1,my+1,mz+1))/(2.*dx)
			  Elsy = -(BoundAys(mx+1,my+1,mz+1)-BoundpAys(mx+1,my+1,mz+1))/dt - (BoundAts(mx+1,my+1+1,mz+1)-BoundAts(mx+1,my-1+1,mz+1))/(2.*dy)
			  Elsz = -(BoundAzs(mx+1,my+1,mz+1)-BoundpAzs(mx+1,my+1,mz+1))/dt - (BoundAts(mx+1,my+1,mz+1+1)-BoundAts(mx+1,my+1,mz-1+1))/(2.*dz)
			  Bmdx = 0. !(BoundAzd(mx+1,my+1+1,mz+1)-BoundAzd(mx+1,my-1+1,mz+1))/(2.*dy)-(BoundAyd(mx+1,my+1,mz+1+1)-BoundAyd(mx+1,my+1,mz-1+1))/(2.*dz)
			  Bmdy = 0. !(BoundAxd(mx+1,my+1,mz+1+1)-BoundAxd(mx+1,my+1,mz-1+1))/(2.*dz)-(BoundAzd(mx+1+1,my+1,mz+1)-BoundAzd(mx-1+1,my+1,mz+1))/(2.*dx)
			  Bmdz = 0. !(BoundAyd(mx+1+1,my+1,mz+1)-BoundAyd(mx-1+1,my+1,mz+1))/(2.*dx)-(BoundAxd(mx+1,my+1+1,mz+1)-BoundAxd(mx+1,my-1+1,mz+1))/(2.*dy)
			  Bmux = 0. !(BoundAzu(mx+1,my+1+1,mz+1)-BoundAzu(mx+1,my-1+1,mz+1))/(2.*dy)-(BoundAyu(mx+1,my+1,mz+1+1)-BoundAyu(mx+1,my+1,mz-1+1))/(2.*dz)
			  Bmuy = 0. !(BoundAxu(mx+1,my+1,mz+1+1)-BoundAxu(mx+1,my+1,mz-1+1))/(2.*dz)-(BoundAzu(mx+1+1,my+1,mz+1)-BoundAzu(mx-1+1,my+1,mz+1))/(2.*dx)
			  Bmuz = 0. !(BoundAyu(mx+1+1,my+1,mz+1)-BoundAyu(mx-1+1,my+1,mz+1))/(2.*dx)-(BoundAxu(mx+1,my+1+1,mz+1)-BoundAxu(mx+1,my-1+1,mz+1))/(2.*dy)
			  Bmsx = 0. !(BoundAzs(mx+1,my+1+1,mz+1)-BoundAzs(mx+1,my-1+1,mz+1))/(2.*dy)-(BoundAys(mx+1,my+1,mz+1+1)-BoundAys(mx+1,my+1,mz-1+1))/(2.*dz)
			  Bmsy = 0. !(BoundAxs(mx+1,my+1,mz+1+1)-BoundAxs(mx+1,my+1,mz-1+1))/(2.*dz)-(BoundAzs(mx+1+1,my+1,mz+1)-BoundAzs(mx-1+1,my+1,mz+1))/(2.*dx)
			  Bmsz = 0. !(BoundAys(mx+1+1,my+1,mz+1)-BoundAys(mx-1+1,my+1,mz+1))/(2.*dx)-(BoundAxs(mx+1,my+1+1,mz+1)-BoundAxs(mx+1,my-1+1,mz+1))/(2.*dy)
			  !Eldx=Eldx *(1.-1./numbd(mx,my,mz))
			  idq = q(9,i,npara)
          if (abs(q(9,i,npara)) .lt. 1.5) then
                 if (inewtonM.eq.1) then
	         q(5,i,npara)=q(5,i,npara)-dt*(BoundQmd(mx+1+1,my+1,mz+1)-BoundQmd(mx-1+1,my+1,mz+1))/(2.*dx)*q(8,i,npara)/Eqst
	         q(6,i,npara)=q(6,i,npara)-dt*(BoundQmd(mx+1,my+1+1,mz+1)-BoundQmd(mx+1,my-1+1,mz+1))/(2.*dy)*q(8,i,npara)/Eqst
	         q(7,i,npara)=q(7,i,npara)-dt*(BoundQmd(mx+1,my+1,mz+1+1)-BoundQmd(mx+1,my+1,mz-1+1))/(2.*dz)*q(8,i,npara)/Eqst
	         endif
			 q(5,i,npara)=q(5,i,npara)+dt*(Eldx + vely*Bmdz - velz*Bmdy)*idq/abs(idq)
	         q(6,i,npara)=q(6,i,npara)+dt*(Eldy + velz*Bmdx - velx*Bmdz)*idq/abs(idq)
	         q(7,i,npara)=q(7,i,npara)+dt*(Eldz + velx*Bmdy - vely*Bmdx)*idq/abs(idq)
                 xxx = dt*(Eldx + vely*Bmdz - velz*Bmdy)*idq/abs(idq)
                  !print *, "force: ",Bmdx,Bmdy,Bmdz
          else if (abs(q(9,i,npara)) .lt. 2.5) then
                 if (inewtonM.eq.1)then
	         q(5,i,npara)=q(5,i,npara)-dt*(BoundQmu(mx+1+1,my+1,mz+1)-BoundQmu(mx-1+1,my+1,mz+1))/(2.*dx)*q(8,i,npara)/Eqst
	         q(6,i,npara)=q(6,i,npara)-dt*(BoundQmu(mx+1,my+1+1,mz+1)-BoundQmu(mx+1,my-1+1,mz+1))/(2.*dy)*q(8,i,npara)/Eqst
	         q(7,i,npara)=q(7,i,npara)-dt*(BoundQmu(mx+1,my+1,mz+1+1)-BoundQmu(mx+1,my+1,mz-1+1))/(2.*dz)*q(8,i,npara)/Eqst
                 endif
			 q(5,i,npara)=q(5,i,npara)+dt*(Elux + vely*Bmuz - velz*Bmuy)*idq/abs(idq)
	         q(6,i,npara)=q(6,i,npara)+dt*(Eluy + velz*Bmux - velx*Bmuz)*idq/abs(idq)
	         q(7,i,npara)=q(7,i,npara)+dt*(Eluz + velx*Bmuy - vely*Bmux)*idq/abs(idq) 
          else
                 if (inewtonM.eq.1)then
	         q(5,i,npara)=q(5,i,npara)-dt*(BoundQms(mx+1+1,my+1,mz+1)-BoundQms(mx-1+1,my+1,mz+1))/(2.*dx)*q(8,i,npara)/Eqst
	         q(6,i,npara)=q(6,i,npara)-dt*(BoundQms(mx+1,my+1+1,mz+1)-BoundQms(mx+1,my-1+1,mz+1))/(2.*dy)*q(8,i,npara)/Eqst
	         q(7,i,npara)=q(7,i,npara)-dt*(BoundQms(mx+1,my+1,mz+1+1)-BoundQms(mx+1,my+1,mz-1+1))/(2.*dz)*q(8,i,npara)/Eqst
                 endif
			 q(5,i,npara)=q(5,i,npara)+dt*(Elsx + vely*Bmsz - velz*Bmsy)*idq/abs(idq)
	         q(6,i,npara)=q(6,i,npara)+dt*(Elsy + velz*Bmsx - velx*Bmsz)*idq/abs(idq)
	         q(7,i,npara)=q(7,i,npara)+dt*(Elsz + velx*Bmsy - vely*Bmsx)*idq/abs(idq)       
          endif

!*************************Electro-magnetic force*************************************
        if (iEM.eq.1)then
            
		id1=q(9,i,npara)
		x1=q(1,i,npara)
		y1=q(2,i,npara)
		z1=q(3,i,npara)
		px1=q(5,i,npara)
		py1=q(6,i,npara)
		pz1=q(7,i,npara)
		pm1=q(8,i,npara)
        Em1=pm1**2.+px1**2.+py1**2.+pz1**2.
        Em1=sqrt(Em1)       
		vx1=px1/Em1
		vy1=py1/Em1
		vz1=pz1/Em1
        if ((id1.gt.1.5).and.(id1.lt.2.5))then
            charge = 2./3.
        else
            charge = -1./3.
        endif
		!charge=1/3.
        Ftotx=id1/abs(id1)*charge*(Exf(mx,my,mz)+vy1*Bzf(mx,my,mz)-vz1*Byf(mx,my,mz))
		Ftoty=id1/abs(id1)*charge*(Eyf(mx,my,mz)+vz1*Bxf(mx,my,mz)-vx1*Bzf(mx,my,mz))
		Ftotz=id1/abs(id1)*charge*(Ezf(mx,my,mz)+vx1*Byf(mx,my,mz)-vy1*Bxf(mx,my,mz))
        
        q(5,i,npara)=q(5,i,npara)+dt*Ftotx
		q(6,i,npara)=q(6,i,npara)+dt*Ftoty
		q(7,i,npara)=q(7,i,npara)+dt*Ftotz
        !if ((mx.eq.25).and.(my.eq.25).and.(mz.eq.20))then
        !print *, Ftotx,Ftoty,Ftotz
        !endif

        if(ispin.eq.1)then
            Fxsp = spin(1,i,npara)*(Bxf(mx+1,my,mz)-Bxf(mx,my,mz))/dx + &
                   spin(2,i,npara)*(Byf(mx+1,my,mz)-Byf(mx,my,mz))/dx + &
                   spin(3,i,npara)*(Bzf(mx+1,my,mz)-Bzf(mx,my,mz))/dx
            Fysp = spin(1,i,npara)*(Bxf(mx,my+1,mz)-Bxf(mx,my,mz))/dy + &
                   spin(2,i,npara)*(Byf(mx,my+1,mz)-Byf(mx,my,mz))/dy + &
                   spin(3,i,npara)*(Bzf(mx,my+1,mz)-Bzf(mx,my,mz))/dy
            Fzsp = spin(1,i,npara)*(Bxf(mx,my,mz+1)-Bxf(mx,my,mz))/dz + &
                   spin(2,i,npara)*(Byf(mx,my,mz+1)-Byf(mx,my,mz))/dz + &
                   spin(3,i,npara)*(Bzf(mx,my,mz+1)-Bzf(mx,my,mz))/dz   
        q(5,i,npara)=q(5,i,npara)+dt*Fxsp
		q(6,i,npara)=q(6,i,npara)+dt*Fysp
		q(7,i,npara)=q(7,i,npara)+dt*Fzsp            
            
        endif
		endif



!***************************Electro-magnectic force end*************************************          
          
          
78        continue    
          
		  

    
		  
		  
		  
! position update
	      Eqst=sqrt(q(5,i,npara)**2. +q(6,i,npara)**2. +q(7,i,npara)**2. +q(8,i,npara)**2.)
	      vx=q(5,i,npara)/Eqst
	      vy=q(6,i,npara)/Eqst
	      vz=q(7,i,npara)/Eqst
	      q(1,i,npara)=q(1,i,npara) +dt*vx
          q(2,i,npara)=q(2,i,npara) +dt*vy
          q(3,i,npara)=q(3,i,npara) +dt*vz
          
		  qxx = q(1,i,npara)
		  qxy = q(2,i,npara)
		  qxz = q(3,i,npara)
		  if (IperiodicBound.eq.1) then
		     if(qxx.gt.qLen) then 
                        !print *,"Out of range:", i,npara, q(1,i,npara)
						nxMaxBound = nxMaxBound +1
			!print *, "Out of range x-max",nxMaxBound,rhobd(nx,nyh,nzh),Qmd(nx,nyh,nzh),BoundQmd(nx+2,nyh+1,nzh+1)
                        q(1,i,npara) = mod(qxx,qLen)-qLen
                        !print *, "In range:",i,npara,q(1,i,npara)
                        endif

			 if(qxx.lt.-qLen) then 
			                nxMinBound = nxMinBound +1
			!print *, "Out of range x-min",nxMinBound,rhobd(1,nyh,nzh),Qmd(1,nyh,nzh),BoundQmd(1,nyh+1,nzh+1)
                            q(1,i,npara) = mod(qxx,qLen)+qLen
                         endif
			 if(qxy.gt.qLen) then 
			                    nyMaxBound = nyMaxBound +1
			!print *, "Out of range y-max",nyMaxBound,rhobd(nxh,1,nzh),Qmd(nxh,1,nzh)
                                q(2,i,npara) = mod(qxy,qLen)-qLen
                         endif
			 if(qxy.lt.-qLen) then 
								nyMinBound = nyMinBound +1
			!print *, "Out of range y-min",nyMinBound,rhobd(nxh,ny,nzh),Qmd(nxh,ny,nzh)
                                q(2,i,npara) = mod(qxy,qLen)+qLen
                         endif
			 if(qxz.gt.qLen) then 
                                q(3,i,npara) = mod(qxz,qLen)-qLen
                        endif
			 if(qxz.lt.-qLen) then 
                         q(3,i,npara) = mod(qxz,qLen)+qLen
                         endif
			 
		  
		  endif
		  
          if(ievolution.eq.1)then
          write(31,2031) int(q(9,i,npara)),q(5,i,npara),q(6,i,npara),q(7,i,npara), &
		  q(8,i,npara),q(1,i,npara),q(2,i,npara),q(3,i,npara),t,ihd(i,npara),icoll(i,npara),EdnPar         
          endif
          
! for p, not for p*

!	      q0(1,i,npara)= Ux(mx,my,mz)
!          q0(2,i,npara)= Uy(mx,my,mz)
!          q0(3,i,npara)= Uz(mx,my,mz)
!------------------------------------ Freeze-out------------------------------------
!<============Determine the hadronization mass by flavor============
          if (abs(q(9,i,npara)).lt.2.5) then
            xmhd=0.3 !0.3   !0.5*xmuv
          else
            xmhd= 0.45  !0.5*xmsv
          end if
!<=============Once a parton is hadronized, record its information in array qq=====
			if(iglobal.eq.2) then 
          if((q(8,i,npara).gt.xmhd).and.(ihd(i,npara).eq.0).and.t.gt.thd) then
            iparthd(npara)=iparthd(npara)+1 !<====# of hadronized partons +1
            NumFt=NumFt+1
            !xNtot = (numbd(mx,my,mz)+numbu(mx,my,mz)+numbs(mx,my,mz))/(dx*dy*dz)/nevent
            !print *,"EdnTL", EdnTL(mx,my,mz),xNtot
            do idx=1, 9
                qq(idx,i,npara)=q(idx,i,npara)
            end do
            qq(4,i,npara)=t-dt
            ihd(i,npara)=1
          end if
			endif
!-----------------------------------Freeze-out--------------------------------------

	   enddo
    enddo          
     
    
	        do i1=1, nx
	        do j1=1, ny
	        do k1=1, nz
               pAxd(i1,j1,k1)=Axd(i1,j1,k1)  !*(1.-1./numbd(i1,j1,k1)) ! to prevent self-interactions
               pAyd(i1,j1,k1)=Ayd(i1,j1,k1)  !*(1.-1./numbd(i1,j1,k1))
               pAzd(i1,j1,k1)=Azd(i1,j1,k1)  !*(1.-1./numbd(i1,j1,k1))
               
               pAxu(i1,j1,k1)=Axu(i1,j1,k1)  !*(1.-1./numbu(i1,j1,k1)) ! to prevent self-interactions
               pAyu(i1,j1,k1)=Ayu(i1,j1,k1)  !*(1.-1./numbu(i1,j1,k1))
               pAzu(i1,j1,k1)=Azu(i1,j1,k1)  !*(1.-1./numbu(i1,j1,k1))
               
               pAxs(i1,j1,k1)=Axs(i1,j1,k1)  !*(1.-1./numbs(i1,j1,k1)) ! to prevent self-interactions
               pAys(i1,j1,k1)=Ays(i1,j1,k1)  !*(1.-1./numbs(i1,j1,k1))
               pAzs(i1,j1,k1)=Azs(i1,j1,k1)  !*(1.-1./numbs(i1,j1,k1))              
               
               
	        enddo
	        enddo
	        enddo      
   
   
   
      ken = 1
	  energykin=0.
	  do npara=1, nevent                 
	   do i=1, itotal(npara) 
	   
	      energykin=energykin+sqrt(q(8,i,npara)**2. +q(5,i,npara)**2. +q(6,i,npara)**2. +q(7,i,npara)**2.)-q(8,i,npara)
		  ken=ken+1
	   enddo
	   enddo
	   energykin = energykin/(ken-1.)
	   write(115,2034)t,energykin
	   
    !print *, "***************--recording----******************" 
    !print *,"recording: t=", t, "nstep = ", Nstep
    t=t+dt
    Nstep=Nstep+1
!******************************** record ************************************************************

    TotCenNum=(2*NCenX+1.)*(2*NCenY+1.0)*(2*NCenZ+1.)
    centerEd=(energyden/(TotCenNum*dx*dy*dz)/npart/nevent +(centerS+centerV)/TotCenNum*5.07**3.) ! unit is GeV/fm^3
    
    totalE=(totalE/npart/nevent +(totalS+totalV)*dV)                         ! unit = GeV
	centerdentotal = dentotal/(TotCenNum*dx*dy*dz)/npart/nevent              ! unit = fm^-3; total quark number density in central cell
	centerdennet = dennet/(TotCenNum*dx*dy*dz)/npart/nevent                  ! unit = fm^-3; net quark number density in central cell
    centerscalU = scalU/(TotCenNum)/npart/nevent                             ! unit = GeV^3; average condensate of u quark
    centerscalD = scalD/(TotCenNum)/npart/nevent                             ! unit = GeV^3; average condensate of d quark
    centerscalS = scalS/(TotCenNum)/npart/nevent                             ! unit = GeV^3; average condensate of d quark
    
    centervectUt = vecUt/(TotCenNum*dx*dy*dz)/npart/nevent                   ! unit = fm^-3; average rhob of u quark
    centervectUx = vecUx/(TotCenNum*dx*dy*dz)/npart/nevent
    centervectUy = vecUy/(TotCenNum*dx*dy*dz)/npart/nevent
    centervectUz = vecUz/(TotCenNum*dx*dy*dz)/npart/nevent

    centervectSt = vecSt/(TotCenNum*dx*dy*dz)/npart/nevent                   ! unit = fm^-3; average rhob of u quark
    centervectSx = vecSx/(TotCenNum*dx*dy*dz)/npart/nevent
    centervectSy = vecSy/(TotCenNum*dx*dy*dz)/npart/nevent
    centervectSz = vecSz/(TotCenNum*dx*dy*dz)/npart/nevent    
    
    centervectDt = vecDt/(TotCenNum*dx*dy*dz)/npart/nevent
    centervectDx = vecDx/(TotCenNum*dx*dy*dz)/npart/nevent
    centervectDy = vecDy/(TotCenNum*dx*dy*dz)/npart/nevent
    centervectDz = vecDz/(TotCenNum*dx*dy*dz)/npart/nevent
    write(1001,*) t,itotal2,NumFt,NumOut
    write(1002,2033)t,Exf(nxh,nyh,nzh),Eyf(nxh,nyh,nzh),Ezf(nxh,nyh,nzh), &
                    Bxf(nxh,nyh,nzh),Byf(nxh,nyh,nzh),Bzf(nxh,nyh,nzh)
    print *,"*************--------Evolution------------**************************"
	print *, "t=           ", t-dt, dt    
    !print *,"Nstep,Ncell   ", Nstep,xNcell,xEzcell/(nx*ny*nz)
    print *,"NumTotal   =           ", itotal2
    print *,"Number of freeze=      ", NumFt
    print *,"Number outside the grid", NumOut
    print *, "Number in central cell", dentotal
    print *, "dx,dy,dz=",dx,dy,dz,"Nevent=",Nevent
    print *, "ntotal = ", ntotal(1),ntotal(2)
    print *, "itotal = ", itotal(1),itotal(2)
    print *, "last paritcle = ", q(1,4642,99)
    print *,"average Num of colls   ", ncoll/nevent
    print *,"centerEd =             ",centerEd
    print *,"centerEdL =            ", EdnTL(nxh,nyh,nzh)
    if(iem.eq.1)then
    print *,"center Exf =            ", Exf(nxh,nyh,nzh)
    print *,"center Eyf =            ", Eyf(nxh,nyh,nzh)
    print *,"center Ezf =            ", Ezf(nxh,nyh,nzh)
    print *,"center Bxf =            ", Bxf(nxh,nyh,nzh)
    print *,"center Byf =            ", Byf(nxh,nyh,nzh)
    print *,"center Bzf =            ", Bzf(nxh,nyh,nzh)
    print *,"center EM =            ", EdnEM(nxh,nyh,nzh)
    print *,"center EMP =          ", EdnEMP(nxh,nyh,nzh)   
    endif
    print *,"number per cell=       ", dentotal/TotCenNum/npart/nevent
    !print *,"number per cell=", (numbd(22,25,20)+numbu(22,25,20)+numbs(22,25,20))/npart/float(nevent)
    print *,"EdnK =                 ", EdnK(nxh,nyh,nzh)      
    !print *,"gam= ", vxs(22,25,20),vys(22,25,20),vzs(22,25,20),gamLV(22,25,20)
    print *,"EdnKL=                 ", EdnKL(nxh,nyh,nzh)
    print *,"*************--------Evolution------------**************************"    
!    print *, "energyden", centerEd,centerdentotal,centerdennet,xmudcen,xmuucen,xmuscen 
    print *,"***************-----properties of central cell------**************** "
    print *,"Energy density in central cell      =",centerEd,"GeV/fm^3"
    print *,"total quark number density          =",centerdentotal,"fm^-3"
    print *,"net quark number density            =",centerdennet,"fm^-3"
    print *,"mass of d quark in central cell     =",xmudcen,"GeV"
    print *,"mass of u quark in central cell     =",xmuucen,"GeV"
    print *,"mass of s quark in central cell     =",xmuscen,"GeV"
    print *,"average condensate of u quark       =",centerscalU,"GeV^3"
    print *,"average condensate of d quark       =",centerscalD,"GeV^3"
    print *,"average condensate of s quark       =",centerscalS,"GeV^3"    
    print *,"average number density of u quark   =",centervectUt,"fm^-3"
    print *,"average x current density of u quark=",centervectUx,"fm^-3"
    print *,"average y current density of u quark=",centervectUy,"fm^-3"
    print *,"average z current density of u quark=",centervectUz,"fm^-3"
    print *,"average number density of d quark   =",centervectDt,"fm^-3"
    print *,"average x current density of d quark=",centervectDx,"fm^-3"
    print *,"average y current density of d quark=",centervectDy,"fm^-3"
    print *,"average z current density of d quark=",centervectDz,"fm^-3"
    print *,"average number density of s quark   =",centervectSt,"fm^-3"
    print *,"average x current density of s quark=",centervectSx,"fm^-3"
    print *,"average y current density of s quark=",centervectSy,"fm^-3"
    print *,"average z current density of s quark=",centervectSz,"fm^-3"    
    print *,"***************------------------------------------**************** "    
	write(101,2014) t, centerEd,centerdentotal,centerdennet,xmudcen,xmuucen,xmuscen, &
    centerscalU, centerscalD, centervectUt, centervectUx, centervectUy, centervectUz, & !xj
    centervectDt, centervectDx, centervectDy, centervectDz
    write(111,2016) t,xmudcen,xmuucen,xmuscen 
    write(110,2015)t,ncoll
    
	xmucen = xmudcen/2. + xmuucen/2.
2014 format(1X,17F15.8)
2016 format(1X,4F10.5)
2015 format(1X,1F15.8,1I8,1F15.8)     
2031 FORMAT(1I6,8(1X,F10.3),1X,1I6,1X,1I6,1X,F10.3)     
2033 FORMAT(7(1X,F13.6))   
2034 FORMAT(2(1X,F13.6))      
     if (iglobal.eq.1)then
         if ((t .gt. thd) .and. (xmucen .gt. 0.3)) goto 99
     end if
     if (iglobal.eq.-1)then
         if ((t .gt. thd) .and. (centerEd .lt. 0.6)) goto 99
     end if     
    
	if (Iglobal.eq.4) then
	if(t.gt.tend) goto 99   
	end if
	if (t .lt. tend) goto 30

99  continue
  
	do npara=1, nevent
       write(100,49) 1,1, ipart(npara) , 0.0, 1, 1, 1, 1
	   do i=1, ipart(npara)

          write(100,2013) q(5,i,npara),q(6,i,npara),q(7,i,npara),&
          int(q(9,i,npara)),q(8,i,npara),q(1,i,npara),q(2,i,npara),q(3,i,npara), t 

	   
	   enddo
    enddo

        !print *, "rhobd :", rhobd(11,11,:)
        ntotald = 0
        rhototald = 0.
        do k=1,nx
           do k1=1,ny
              do k2=1,nz

                ntotald = ntotald + numbd(k,k1,k2)
                rhototald = rhototald + rhobd(k,k1,k2)
                enddo
            enddo
        enddo
        print *,"final total number of d: ",ntotald
        print *, "final total density of d: ",rhototald*5.07**3.

	do i=1, nx
	   do j=1, ny
	      do k=1, nz
             write(103,17) scalarD(i,j,k),scalarU(i,j,k)  !(scalarD(i,j,k)-scalarU(i,j,k))/(scalarD(i,j,k)+scalarU(i,j,k)) 
             write(104,17) rhobd(i,j,k),rhobu(i,j,k)      !(vectDt(i,j,k)-vectUt(i,j,k))/(vectDt(i,j,k)+vectUt(i,j,k))
             write(105,17) rhobxd(i,j,k),rhobxu(i,j,k)    !(vectDx(i,j,k)-vectUx(i,j,k))/(vectDx(i,j,k)+vectUx(i,j,k))
             write(106,17) rhobyd(i,j,k),rhobyu(i,j,k)    !(vectDy(i,j,k)-vectUy(i,j,k))/(vectDy(i,j,k)+vectUy(i,j,k))
             write(107,17) rhobzd(i,j,k),rhobzu(i,j,k)    !(vectDz(i,j,k)-vectUz(i,j,k))/(vectDz(i,j,k)+vectUz(i,j,k))             
          enddo
       enddo
       enddo
    

 15 format(1x,2(i7,1x),1(F10.3,1x),4(I7,1x))
 16 format(1x,3(i5,1x),8(f9.6,1x))
 17 format(1x,2(1x,F14.6))
 19 format(1x,6(1x,F14.6))
 49 format(3(1x,i6),1x,f10.4,4(1x,i6))
2013	FORMAT(3(1X,F10.3),1X,1I6,5(1X,F10.3))    
    
    
    
!    do i=1,4!
!	do j=1,4
!       ETens(i,j)=0.d0!
!	enddo
!    enddo
    
!	ETens(1,1)=1.0
!	ETens(2,2)=0.8
!	ETens(3,3)=0.5
!	ETens(4,4)=0.5
!	ETens(1,2)=0.8
!	ETens(1,3)=0.
!	ETens(1,4)=0.
!	ETens(2,1)=0.8
!	ETens(3,1)=0.
!	ETens(4,1)=0.
 !   vx=0.89    !ETens(1,2)/ETens(1,1)
!	vy=ETens(1,3)/ETens(1,1)
!	vz=ETens(1,4)/ETens(1,1)
!    call xLorMatrix(vx,vy,vz,xLorMat)
!    call xLorMatrix(-vx,-vy,-vz,xLorMatinv)
!	print *, xLorMatinv(1,1)
!    call Cross(xLorMat,ETens,ETens1)
!	print *, ETens1(2,2)
!    call Cross(ETens1,xLorMatinv,ETens)
!    print *, ETens(1,1) 
    
    
    enddo
    close(101)
    end program

		
		
	subroutine xBound(x,nx,ny,nz,Boundx)
	implicit real*8 (a-h,o-z)
	integer,parameter :: Ngx = 31, Ngy = 31, Ngz = 31  
	dimension x(Ngx,Ngy,Ngz),Boundx(Ngx+2,Ngy+2,Ngz+2)
	do ii=1, nx
	   do jj=1, ny
	      do kk=1, nz

			 Boundx(ii+1,jj+1,kk+1) = x(ii,jj,kk)
             
             !iif (numbs(ii,jj,kk) .lt. 1) numbs(ii,jj,kk)=1  ! added by sunkkjj
          enddo
       enddo
    enddo
	
	do ii=1, nx
	   do jj=1, ny
	      !do kk=1, nz
			 kk=nz
			 Boundx(ii+1,jj+1,1) = x(ii,jj,kk)
			 
			 kk=1
			 Boundx(ii+1,jj+1,nz+2) = x(ii,jj,kk)
          !enddo
       enddo
    enddo
	
	do ii=1, nx
	   !do jj=1, ny
	      do kk=1, nz
			 jj=ny
			 Boundx(ii+1,1,kk+1) = x(ii,jj,kk)
			 
			 jj=1
			 Boundx(ii+1,ny+2,kk+1) = x(ii,jj,kk)
          enddo
       !enddo
    enddo
	
	!do ii=1, nx
	   do jj=1, ny
	      do kk=1, nz
			 ii=nx
			 Boundx(1,jj+1,kk+1) = x(ii,jj,kk)
			 
			 ii=1
			 Boundx(nx+2,jj+1,kk+1) = x(ii,jj,kk)
          enddo
       enddo
	   
	   return
    end
		
	  subroutine xLorentz(t1,x1,y1,z1, vx,vy,vz, t2,x2,y2,z2)
	  implicit real*8 (a-h,o-z)

      beta2=vx*vx +vy*vy +vz*vz
      if (beta2.gt.1.d-9) then  
      gam=1./sqrt(1.-beta2)

      t2=gam*(t1 -vx*x1 -vy*y1 -vz*z1)
      x2=x1 -vx*gam*t1 +(gam-1.)*vx*vx/beta2*x1 +(gam-1.)*vx*vy/beta2*y1 +(gam-1.)*vx*vz/beta2*z1
      y2=y1 -vy*gam*t1 +(gam-1.)*vy*vx/beta2*x1 +(gam-1.)*vy*vy/beta2*y1 +(gam-1.)*vy*vz/beta2*z1
      z2=z1 -vz*gam*t1 +(gam-1.)*vz*vx/beta2*x1 +(gam-1.)*vz*vy/beta2*y1 +(gam-1.)*vz*vz/beta2*z1
      else
      t2=t1
      x2=x1
      y2=y1
      z2=z1
      end if

      return
      end


	  subroutine xLorMatrix(vx,vy,vz,xLorMat)
	  implicit real*8 (a-h,o-z)
      dimension xLorMat(4,4)
	  do i=1,4
	  do j=1,4
	  xLorMat(i,j)=0.d0
	  enddo
	  enddo
	  xLorMat(1,1)=1.d0
	  xLorMat(2,2)=1.d0
	  xLorMat(3,3)=1.d0
	  xLorMat(4,4)=1.d0

      beta2=vx*vx +vy*vy +vz*vz
      if (beta2.gt.1.d-9) then  
      gam=1./sqrt(1.-beta2)
	  xLorMat(1,1)=gam
	  xLorMat(1,2)=-gam*vx
	  xLorMat(1,3)=-gam*vy
	  xLorMat(1,4)=-gam*vz
	  xLorMat(2,1)=-gam*vx
	  xLorMat(2,2)= (gam-1.)*vx*vx/beta2+1.
	  xLorMat(2,3)= (gam-1.)*vx*vy/beta2
	  xLorMat(2,4)= (gam-1.)*vx*vz/beta2
	  xLorMat(3,1)=-gam*vy
	  xLorMat(3,2)=(gam-1.)*vy*vx/beta2
	  xLorMat(3,3)=(gam-1.)*vy*vy/beta2+1.
	  xLorMat(3,4)=(gam-1.)*vy*vz/beta2
	  xLorMat(4,1)=-gam*vz
	  xLorMat(4,2)=(gam-1.)*vz*vx/beta2
	  xLorMat(4,3)=(gam-1.)*vz*vy/beta2
	  xLorMat(4,4)=(gam-1.)*vz*vz/beta2+1.
      end if

      return
    end

      subroutine EMForce(id1,x1,y1,z1,vx1,vy1,vz1,id2,x2,y2,z2,vx2,vy2,vz2,Fx,Fy,Fz)
	  implicit real*8 (a-h,o-z)
      charge=1/3.
	  sign1=id1/abs(id1)
	  sign2=id2/abs(id2)
	  ealpha=1/137.          ! electro-magnetic structure constant
	  const1 = sign1*sign2*charge**2.*ealpha
	  hbarc = 0.19733        ! fm*GeV
	  Rx = x1-x2
	  Ry = y1-y2
	  Rz = z1-z2
	  R=sqrt(Rx**2.+Ry**2.+Rz**2.)
	  v22=vx2**2.+vy2**2.+vz2**2.
	  Rv2 = Rx*vx2+Ry*vy2+Rz*vz2
	  
	  eEx=(1-v22)/(R-Rv2)**3. *(Rx-R*vx2)
      eEx =eEx*const1
	  eEy=(1-v22)/(R-Rv2)**3. *(Ry-R*vy2)
      eEy =eEy*const1
	  eEz=(1-v22)/(R-Rv2)**3. *(Rz-R*vz2)
      eEz =eEz*const1

	  eBx=(1-v22)/(R-Rv2)**3. *(vy2*Rz-Ry*vz2)
	  eBy=(1-v22)/(R-Rv2)**3. *(vz2*Rx-Rz*vx2)
	  eBz=(1-v22)/(R-Rv2)**3. *(vx2*Ry-Rx*vy2)
	  eBx=eBx*const1
	  eBy=eBy*const1
	  eBz=eBz*const1
      
	  Fx=eEx + (vy1*eBz-vz1*eBy)            ! unit= GeV/fm
	  Fy=eEy + (vz1*eBx-vx1*eBz)  
	  Fz=eEz + (vx1*eBy-vy1*eBx)
      Fx=Fx*hbarc                           ! unit=GeV/fm
	  Fy=Fy*hbarc
	  Fz=Fz*hbarc
	  return 
    end    

      subroutine EMField(id1,x1,y1,z1,vx1,vy1,vz1,id2,x2,y2,z2,vx2,vy2,vz2,Ex,Ey,Ez,Bx,By,Bz,Eneg) 
	  implicit real*8 (a-h,o-z)
      charge=1/3.
	  sign1=id1/abs(id1)
	  sign2=id2/abs(id2)
	  ealpha=1/137.          ! electro-magnetic structure constant
      charge1 = sign1*charge
      charge2 = sign2*charge
	  const1 = charge1*charge2*ealpha

      v1xv2 = vx1*vx2+vy1*vy2+vz1*vz2
	  hbarc = 0.19733        ! fm*GeV
	  Rx = x1-x2
	  Ry = y1-y2
	  Rz = z1-z2
	  R=sqrt(Rx**2.+Ry**2.+Rz**2.)
	  v22=vx2**2.+vy2**2.+vz2**2.
	  Rv2 = Rx*vx2+Ry*vy2+Rz*vz2
	  
	  Ex=(1-v22)/(R-Rv2)**3. *(Rx-R*vx2)
	  Ey=(1-v22)/(R-Rv2)**3. *(Ry-R*vy2)
	  Ez=(1-v22)/(R-Rv2)**3. *(Rz-R*vz2)
      Ex=Ex*sign2
      Ey=Ey*sign2
      Ez=Ez*sign2

	  Bx=(1-v22)/(R-Rv2)**3. *(vy2*Rz-Ry*vz2)
	  By=(1-v22)/(R-Rv2)**3. *(vz2*Rx-Rz*vx2)
	  Bz=(1-v22)/(R-Rv2)**3. *(vx2*Ry-Rx*vy2)
      Bx=Bx*sign2
      By=By*sign2
      Bz=Bz*sign2
      Eneg = const1*(1-v22)/(R-Rv2)- const1*(1-v22)/(R-Rv2)*(v1xv2)
      Eneg = Eneg * hbarc                ! GeV 
      !if(Ex.gt.10.0)then
         ! print*,R,Ex,Ey,Ez,Bx,By,Bz
          !endif
	  return 
      end      
    
      subroutine EMFieldcell(rx,ry,rz,V,EMJt,EMJx,EMJy,EMJz,Ex,Ey,Ez,Bx,By,Bz) 
	  implicit real*8 (a-h,o-z)
	  ealpha=1/137.          ! electro-magnetic structure constant
      r=(rx**2.+ry**2.+rz**2.)**0.5
	  EMQ=EMJt*V
      
      Ex=EMQ/r**3*rx
      Ey=EMQ/r**3*ry
      Ez=EMQ/r**3*rz
      if(abs(EMQ).le.1E-3)then
          Bx=0.
          By=0.
          Bz=0.
      else
          EMQx=EMJx*V
          EMQy=EMJy*V
          EMQz=EMJz*V
          vx=EMQx/EMQ
          vy=EMQy/EMQ
          vz=EMQz/EMQ
          Bx=vy*Ez-vz*Ey
          By=vz*Ex-vx*Ez
          Bz=vx*Ey-vy*Ex
      endif
	  return 
    end 

      subroutine EMFieldspec(Rx,Ry,Rz,vx2,vy2,vz2,sign,Ex,Ey,Ez,Bx,By,Bz) 
	  implicit real*8 (a-h,o-z)
      R=sqrt(Rx**2.+Ry**2.+Rz**2.)
      
	  v22=vx2**2.+vy2**2.+vz2**2.
	  Rv2 = Rx*vx2+Ry*vy2+Rz*vz2
      
	  Ex=(1-v22)/(R-Rv2)**3. *(Rx-R*vx2)
	  Ey=(1-v22)/(R-Rv2)**3. *(Ry-R*vy2)
	  Ez=(1-v22)/(R-Rv2)**3. *(Rz-R*vz2)
      Ex=(1-v22)/(Rv2**2.+R**2.*(1-v22))**1.5*Rx
      Ey=(1-v22)/(Rv2**2.+R**2.*(1-v22))**1.5*Ry
      Ez=(1-v22)/(Rv2**2.+R**2.*(1-v22))**1.5*Rz 
      Ex=Ex*sign
      Ey=Ey*sign
      Ez=Ez*sign

	  Bx=(1-v22)/(R-Rv2)**3. *(vy2*Rz-Ry*vz2)
	  By=(1-v22)/(R-Rv2)**3. *(vz2*Rx-Rz*vx2)
	  Bz=(1-v22)/(R-Rv2)**3. *(vx2*Ry-Rx*vy2)
      Bx=(1-v22)/(Rv2**2.+R**2.*(1-v22))**1.5*(vy2*Rz-Ry*vz2)
      By=(1-v22)/(Rv2**2.+R**2.*(1-v22))**1.5*(vz2*Rx-Rz*vx2)
      Bz=(1-v22)/(Rv2**2.+R**2.*(1-v22))**1.5*(vx2*Ry-Rx*vy2)     
      Bx=Bx*sign
      By=By*sign
      Bz=Bz*sign
	  return 
    end 
    
	  subroutine Cross(A,B,C)
	  implicit real*8 (a-h,o-z)
      dimension A(4,4),B(4,4),C(4,4)
	  do i=1,4
	  do j=1,4
         C(i,j)=0.
	  enddo
	  enddo
	  do i=1,4
	  do j=1,4
	  do k=1,4
         C(i,j)=C(i,j)+A(i,k)*B(k,j)
	  enddo
	  enddo
	  enddo

      return 
	  end


      SUBROUTINE SRAND(ISEED)
!C
!C  This subroutine sets the integer seed to be used with the
!C  companion RAND function to the value of ISEED.  A flag is
!C  set to indicate that the sequence of pseudo-random numbers
!C  for the specified seed should start from the beginning.
!C
      COMMON /SEED/JSEED,IFRST
!C
      JSEED = ISEED
      IFRST = 0
!C
      RETURN
      END

      REAL*8 FUNCTION RAND()
!C
!C  This function returns a pseudo-random number for each invocation.
!C  It is a FORTRAN 77 adaptation of the "Integer Version 2" minimal
!C  standard number generator whose Pascal code appears in the article:
!C
!C     Park, Steven K. and Miller, Keith W., "Random Number Generators:
!C     Good Ones are Hard to Find", Communications of the ACM,
!C     October, 1988.
!C
      PARAMETER (MPLIER=16807,MODLUS=2147483647,MOBYMP=127773, &
                MOMDMP=2836)
!C
      COMMON  /SEED/JSEED,IFRST
      INTEGER HVLUE, LVLUE, TESTV, NEXTN
      SAVE    NEXTN
!C
      IF (IFRST .EQ. 0) THEN
        NEXTN = JSEED
        IFRST = 1
      ENDIF
!C
      HVLUE = NEXTN / MOBYMP
      LVLUE = MOD(NEXTN, MOBYMP)
      TESTV = MPLIER*LVLUE - MOMDMP*HVLUE
      IF (TESTV .GT. 0) THEN
        NEXTN = TESTV
      ELSE
        NEXTN = TESTV + MODLUS
      ENDIF
      RAND = REAL(NEXTN)/REAL(MODLUS)
!C
      RETURN
      END
      BLOCKDATA RANDBD
      COMMON /SEED/JSEED,IFRST
!C
      DATA JSEED,IFRST/123456789,0/
!C
      END
