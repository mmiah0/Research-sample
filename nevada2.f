      program nevada

*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
* comment out LINES 2844 and 2845 for fewer number of equilibrium iterations   *                                                           *
*                                                                              *
*                                  "NEVADA"                                    *
*                                                                              *
*                                                                              *
*              A FINITE ELEMENT PROGRAM FOR THE THREE DIMENSIONAL              *
*                    NONLINEAR ANALYSIS OF BUILDING STRUCTURES                 *
*                                                                              *
*                                                                              *
*                               DAVID MCCALLEN                                 *
*                    CENTER FOR COMPLEX DISTRIBUTED SYSTEMS                    *
*                           ENGINEERING DIRECTORATE                            *
*                    LAWRENCE LIVERMORE NATIONAL LABORATORY                    *
*                            UNIVERSITY OF CALIFORNIA                          *
*                                                                              *
*                               LARRY SANFORD                                  *
*                COMPUTING APPLICATIONS AND RESEARCH DEPARTMENT                *
*                          COMPUTATIONS DIRECTORATE                            *
*                    LAWRENCE LIVERMORE NATIONAL LABORATORY                    *
*                            UNIVERSITY OF CALIFORNIA                          *
*                                                                              *
*                                August 2002                                   *
*                                                                              *
*                                                                              *
*                                                                              *
*   This program solves for the linear and nonlinear transient response of a   *
*          building subjected to earthquake loading using implicit             * 
*                            Newmark-Beta integration                          *
*                                                                              *
*     This version currently has the following static memory allocation:       *
*                                                                              *
*     30,000 DOF                                                               *
*     10,000 nodes                                                             *
*     5000 beams, 50 beam properties, 35 beam fibers points                    *
*                                                                              *
*  SGI Compile: f77 -r8 -g -o nevada -old_rl -static -O -G 0 nevada.f          *
*  GNU Fortran 4.9.2:  gfortran -o xnevada -g nevada2.f -fdefault-double-8 -fdefault-real-8 
*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

*                             ARRAY DECLARATIONS    

*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
                        
*-----Perform array declarations

*-----Establishment of global parameters

      integer
     :        max_a_array
     :       ,max_beams
     :       ,max_beam_fibers_points
     :       ,max_beam_properties
     :       ,max_beams_times_three
     :       ,max_dof
     :       ,max_dof_plus_1
     :       ,max_iplot
     :       ,max_nodes

      parameter  (
     :           max_a_array            = 1500000
     :          ,max_beams              =    5000
     :          ,max_beam_fibers_points =      20
     :          ,max_beam_properties    =      50
     :          ,max_dof                =   30000
     :          ,max_iplot              =  100000
     :          ,max_nodes              =   10000
     :           )

      parameter  (
     :           max_beams_times_three = max_beams * 3
     :          ,max_dof_plus_1        = max_dof + 1 
     :           )
 
*-----Geometry and displacement block

      integer
     :        id( max_nodes, 6 )

      real
     :      accvctr( max_dof ), dispinc( max_dof ),dispvctr( max_dof )
     :      ,velvctr( max_dof ), xyzloc( max_nodes, 3 ), 
     :      oldvel( max_dof ), oldacc( max_dof ), prevdisp( max_dof )


      common /geometry/
     :        id, xyzloc, dispvctr, dispinc, velvctr, accvctr

*-----Matrix and equation block

      integer
     :        maxa( max_dof_plus_1 ),   maxasave( max_dof_plus_1 ),
     :        ncolht( max_dof_plus_1 ), ncolhtsave( max_dof_plus_1 )

      real
     :     a( max_a_array ),         damp( max_a_array ),
     :     residual( max_dof ),      xkeff( max_a_array ),
     :     xkeffsave( max_a_array ), xmass( max_dof )


      common /matrix/
     :        a,         damp,  residual,  xmass,    maxa,
     :        ncolht,    xkeff, xkeffsave, maxasave,
     :        ncolhtsave

*-----Beam

      integer
     :        nbeamn( max_beams, 3 ), nbeamp( max_beams, 2 ),
     :        numfibrs( max_beam_fibers_points ),
     :        nyldink2( max_beams, 9, max_beam_properties ),
     :        nyldstk2( max_beams, 9, max_beam_properties )

      real
     :      area( max_beam_properties )
     :     ,beamdef( max_beams, 7 )
     :     ,bmdefmst( max_beams, 7 )
     :     ,delepsk2( max_beams, 9, max_beam_properties )
     :     ,e( max_beam_properties )
     :     ,elin( max_beam_properties )
     :     ,etab( max_beam_properties )
     :     ,fibrdefn( max_beam_fibers_points, max_beam_properties, 3 )
     :     ,fmstrk2( max_beams, 12 )
     :     ,forcek2( max_beams, 12 )
     :     ,gj( max_beam_properties )
     :     ,resk2( 12 )
     :     ,sheary( max_beam_properties )
     :     ,shearz( max_beam_properties )
     :     ,sigyld( max_beam_properties )
     :     ,straink2( max_beams, 9, max_beam_properties )
     :     ,stressk2( max_beams, 9, max_beam_properties )
     :     ,strmstk2( max_beams, 9, max_beam_properties )

      real
     :      strsmstk2( max_beams, 9, max_beam_properties )
     :     ,tik2( max_beams_times_three, 3 )
     :     ,tjk2( max_beams_times_three, 3 )
     :     ,tk2( 3, 3 )
     :     ,xj( max_beam_properties )
     :     ,xkk2( 12, 12 )
     :     ,xlamdab( max_beam_properties )
     :     ,xlk2( max_beams )
     :     ,xlnmsk2( max_beam_properties )
     :     ,xlobat( 3 )
     :     ,xlobfres( max_beams, 9, 4 )
     :     ,xlobwt( 3 )
     :     ,xnlnmsk2( max_beam_properties )
     :     ,xnu( max_beam_properties )
     :     ,yi( max_beam_properties )
     :     ,yldcntrk2( max_beams, 9, max_beam_properties, 2 )
     :     ,zi( max_beam_properties )


      common /beam/
     :        nbeamn,   nbeamp,   tk2,      tik2,      tjk2,
     :        xlk2,     xkk2,     area,     e, 
     :        yi,       zi,       xj,       xnu,       sheary,
     :        shearz,   fmstrk2,  gj,       resk2,     forcek2,
     :        xlnmsk2,  beamdef,  xnlnmsk2, elin,      sigyld,
     :        xlamdab,  etab,     straink2, xlobat,    fibrdefn,
     :        numfibrs, xlobfres, nyldink2, nyldstk2,  yldcntrk2,
     :        stressk2, delepsk2, strmstk2, strsmstk2, bmdefmst,
     :        xlobwt

*-----Load vectors

      integer
     :        loads( 10, 3 ), ncurvloc( max_beam_properties )

      real
     :     ar( 50000 ),  blo( 80 ), br( 50000 ),  bup( 80 ),
     :     bupc( 80 ),   dd( 80 ),  eigval( 80 ), eigvec( max_dof, 80 ),
     :     rtolv( 80 ),  ttt( max_dof ),          vec( 80, 80 ),    
     :     w( max_dof ), xload( 200000, 2 ), xnstif( max_a_array )


      common /loads/
     :        loads, xload, eigvec, eigval,   ttt,
     :        w,     ar,    br,     vec,      dd,    rtolv,
     :        bup,   blo,   bupc,   ncurvloc, xnstif


      real dispold( max_dof )

*-----Initialization of variables

      real
     :     time

      data
     :     time /0.0/


*-----Griz/Taurus files

      integer
     :        iob12( 1080 )
     :       ,iplot( max_iplot )

      real * 4
     :         plot( max_iplot )


      common /taurus/
     :        iob12, plot

c     equivalence( plot, iplot )


      data mxfsiz / 2621440 /
      character*8 fname, ckeep
      fname='plots'
      ckeep='keep'

*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

*                               ASSIGN I/O FILES   

*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     
*-----Assign input and output files

      open(unit=5,file='input',status='old',access='sequential',
     $    form='formatted')
      open(unit=6,file='output',status='unknown',access='sequential',
     $    form='formatted')
      open(unit=8,file='igeninfo', status='unknown',access='sequential',
     $     form='formatted')

*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

*                        ASSIGN PROGRAM CONTROL PARAMETERS   

*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

*-----Problem solution type control information

*-----itype=0 - linear dynamic time-history analysis with implicit 
*               time integration
*-----itype=1 - nonlinear dynamic time-history analysis with implicit
*               time integration
*-----itype=2 - eigensolution
*-----itype=3 - nonlinear static analysis
*-----ngravity=1 - first load step is gravity initialization of the model
*-----ngravity=0 - no gravity initialiation of model

*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

*                                READ INPUT DATA   

*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

*-----Read input data

*-----Master control cards

      read(5,5000)numnp,nbeam,itype,tolrance
      read(5,5002)nsteps,stepsize,numfreq,ngravity,ngrndacc
      read(5,5003)npntload,nloadcrv,accgrav,dampstiff,dampmass
      read(5,5004)nprint,nplotcntrl,gamma,beta

*-----Nodal point data

      do 10 i=1,numnp
      read(5,5010)(id(i,j),j=1,6),(xyzloc(i,j),j=1,3)
   10 continue

      if(itype .eq. 2)then
      nmasdof=0
      do 12 i=1,numnp
      do 11 j=1,3
      if(id(i,j) .eq. 0)nmasdof=nmasdof+1
   11 continue
   12 continue
      endif


*-----Beam Element data

      numbmp=0
      nfbrbmp=0
      nfbrsctn=0

      if(nbeam .ne. 0)then
      read(5,5045)numbmp,nfbrbmp,nfbrsctn
      read(5,5050)(nbeamn(i,1),nbeamn(i,2),nbeamn(i,3),
     $             nbeamp(i,1),nbeamp(i,2),i=1,nbeam)  
      if(numbmp .ne. 0)then
      read(5,5055)(e(i),area(i),yi(i),zi(i),xj(i),xnu(i),
     $            sheary(i),shearz(i),xlnmsk2(i),i=1,numbmp)
      endif
      if(nfbrbmp .ne. 0)then
      read(5,5059)(elin(i),sigyld(i),xlamdab(i),etab(i),gj(i),
     $            xnlnmsk2(i),i=1,nfbrbmp) 
      do 16 i=1,nfbrsctn
      read(5,5060)numfibrs(i)
      do 15 j=1,numfibrs(i)
      read(5,5061)(fibrdefn(i,j,k),k=1,3)
   15 continue
   16 continue     
      endif
      endif

*-----Read load locations
       
      if(npntload .ne. 0)then
      read(5,5070)((loads(i,j),j=1,3),i=1,npntload)
      endif

*-----NOTE...
*-----First column of "loads" gives node number
*-----Second column of "loads" gives direction (1=x,2=y,3=z)
*-----Third column of "loads" gives load curve number

*-----Read in load curves - store curves stacked sequentially in a vector

      if(nloadcrv .ne. 0)then
      lshift=0
      do 23 k=1,nloadcrv
      read(5,5080)npts,curvscal
      nlow=lshift+1
      nhigh=lshift+npts
      do 22 kk=nlow,nhigh
      read(5,*)(xload(kk,j),j=1,2)
      xload(kk,2)=xload(kk,2)*curvscal
   22 continue
      ncurvloc(k)=nlow
      lshift=lshift+npts
   23 continue
      endif

*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

*                               ECHO INPUT DATA   

*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

*-----Echo input data to output file

      write(6,6000)
      write(6,6010)
      write(6,6020)itype
      write(6,6030)numnp,nbeam
      write(6,6040)nsteps,stepsize
      write(6,6005)
      write(6,6050)
      write(6,6060)
      write(6,6070)(id(i,1),id(i,2),id(i,3),id(i,4),id(i,5),id(i,6),
     $              xyzloc(i,1),xyzloc(i,2),xyzloc(i,3),i=1,numnp)
      write(6,6005)

      write(6,6096)
      write(6,6097)(i,nbeamn(i,1),nbeamn(i,2),nbeamn(i,3),nbeamp(i,1),
     $i=1,nbeam)
      write(6,6098)
      write(6,7000)(i,e(i),xnu(i),area(i),yi(i),zi(i),xj(i),sheary(i),
     $shearz(i),i=1,numbmp)

      write(6,6005)
      write(6,6099)

*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

*                          INITIALIZE TAURUS DATABASE   

*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

*-----Set up Taurus database and initialize model geometry 

*-----Initialize Taurus database

c     call rwabsf(iob12,fname,12,512,mxfsiz)

*-----Determine state length for plot file adjustments

      if(itype .eq. 2)then
      lenst=3*numnp+1
      elseif(nplotcntrl .eq. 0)then
      lenst=3*numnp+1
      elseif(nplotcntrl .eq. 1)then
      lenst=(3*numnp+1)+6*nbeam
      endif

c      if((itype .ne. 1) .or. (itype .eq. 1. .and. 
c     $     ngravity .eq. 0))then
c     call bullinit(numnp,nbeam,
c    $              numbmp,nfbrbmp,iaddr,itype,ngravity,
c    $              nplotcntrl)
c      endif

*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

*                              INITIALIZE ARRAYS   

*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

*-----Initialize parameters and arrays


*-----Switch ID array from boundary condition codes to equation numbers

      nactdof=0
      numeqn=0

      do 32 i=1,numnp
      do 31 j=1,6
      if(id(i,j) .eq. 1)then
      id(i,j)=0
      elseif(id(i,j) .lt. 0)then
      id(i,j)=id(abs(id(i,j)),j)
      elseif(id(i,j) .eq. 0)then
      numeqn=numeqn+1
      id(i,j)=numeqn
      endif
   31 continue
   32 continue
      nactdof=numeqn
      write(6,6100)nactdof

*-----Perform intializations for beam elements

      if(nbeam .ne. 0)then
      call beaminit(nbeam)
      endif

*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

*                          ENTER LOAD (TIME) STEPPING LOOP  

*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

*-----Load (time) increment stepping loop
     
*-----For linear eigensolution without gravity initialization skip load set do loop 

      if(itype .eq. 2)nsteps=1

*-----For gravity initialization increment steps by 1

      if(ngravity .eq. 1 .and. itype .ne. 2)nsteps=nsteps+1
      if(ngravity .eq. 1 .and. itype .eq. 2)nsteps=2

      do 1000 kload=1,nsteps

      iter=0

*-----Set continue statement to jump back to when solution is unconverged
*-----during equilibrium iterations

   45 continue

*-----Initialize matrices

      do 50 i=1,max_dof,1
      residual(i)=0.
   50 continue

      If(itype .eq. 0 .and. kload .gt. 1)go to 70

      nelembs=nactdof

      do 51 i=1,max_dof_plus_1,1
      maxa(i)=i
      ncolht(i)=1
   51 continue

      do 55 i=1,max_a_array,1
      a(i)=0.
      xnstif(i)=0.
      xkeff(i)=0.
   55 continue

   70 continue
       
*-----Enter beam element loop

      if(nbeam .ne. 0)then

      do 200 k2=1,nbeam

*-----Form beam translation matrix based on displacements from last 
*-----equilibrium state

      if(kload .eq. 1 .and. iter .eq. 0)then
      call trninit(k2)
      else
      call trnbeam(k2,iter,itype,kload)
      endif

*-----Form beam instantaneous stiffness and residual vector 
*-----in updated coordinate system (order of call depends on
*-----whether using standard or fiber beam element)

      if(nbeamp(k2,2) .eq. 0)then
      call stfbeam(k2)
      if((itype .eq. 3 .and. iter .gt. 0) .or.
     $(itype .eq. 1 .and. iter .gt. 0) .or.
     $(itype .eq. 0 .and. kload .gt. 1) .or.
     $(itype .eq. 2 .and. ngravity .eq. 1 .and. iter .gt. 0))then
      call forcbeam(k2)
      endif
      else
      call frcbeamf(k2,itype,iter)
      call stfresk2(k2)
      endif

*-----Translate beam stiffness and residual to global coordinates

      call trntbeam(k2,iter,itype,kload)

*-----Form beam mass matrix

      if(itype .lt. 3)then
      if(kload .le. 1 .and. iter .eq. 0)then
      call massbeam(k2,beammass,beamrmass,beamtmass)
      endif
      endif    
      if(itype .eq. 3 .and. ngravity .eq. 1)then
      if(kload .le. 1 .and. iter .eq. 0)then
      call massbeam(k2,beammass,beamrmass,beamtmass)
      endif
      endif 

*-----Load beam element contribution into global stiffness matrix and 
*-----residual vector

      if(itype .eq. 0 .and. kload .gt. 1)go to 200

      call loadbeam(kload,k2,itype,beammass,beamrmass,beamtmass,
     $             nelembs,nactdof,iter,ngravity)

  200 continue
      endif

*-----End beam element loop


*-----Form system damping matrix

      if(dampstiff .ne. 0 .or. dampmass .ne. 0)then
      if((itype .eq. 0 .and. kload .eq. 1 .and. iter .eq. 0) .or.
     $(itype .eq. 1 .and. kload .le. 50))then
*     The above line allows the system stiffness matrix to become fully 
*     populated after finite deformation of the structure (i.e. the 
*     matrix to fill-out with off-diagonal coupling terms, before settling 
*     on the final damping matrix
      call dampmatrix(nelembs,numnp,dampstiff,dampmass)
      endif
      endif

*-----Jump to subspace for eigensolution

      if(itype .eq. 2 .and. ngravity .ne. 1)go to 1010
      if(itype .eq. 2  .and. ngravity .eq. 1 .and. kload .eq. 2)
     $go to 1010

*-----Apply point loads to model 

      if((npntload .ne. 0 .and. ngravity .eq. 1 .and.
     $kload .gt. 1) .or. (npntload .ne. 0 .and. ngravity .eq. 0))
     $then
      call pointloads(npntload,ngravity,iter,kload)
      endif

*-----Add ground accelerations to residual vector

      if((ngrndacc .gt. 0 .and. ngravity .eq. 1 .and. 
     $kload .gt. 1) .or. (ngrndacc .gt. 0 .and. ngravity .eq. 0))
     $then 
      call groundmotion(numnp,kload,iter)
      endif
     
      if(ngravity .eq. 1)then
      if(iter .gt. 0 .or. kload .eq. 1)then
      do 306 i=1,numnp
      if(id(i,3) .ne. 0)then
      residual(id(i,3))=residual(id(i,3))-
     $xmass(id(i,3))*accgrav
      endif
  306 continue
      endif
      endif

*-----For dynamic analysis form effective stiffness and residual 
*-----vector
      
      if(ngravity .eq. 1 .and. kload .eq. 1)go to 307
      if(itype .le. 1)then
      call timehistory(gamma,beta,stepsize,iter,numnp,nelembs,
     $nactdof,kload,itype)
      endif
  307 continue

*-----For linear time history solution, save matrices

      if(itype .eq. 0 .and. kload .eq. 1)then
      do 308 i=1,nelembs
      xkeffsave(i)=xkeff(i)
  308 continue
      ntopn=nactdof+1
      do 309 i=1,ntopn
      maxasave(i)=maxa(i)
      ncolhtsave(i)=ncolht(i)
  309 continue
      endif

*-----For linear time history solution, skip equilibrium check

      if(itype .eq. 0)go to 400

*-----Check convergence by comparison of the Eucledian norms of the 
*-----residual and incremental displacement vectors

        
      sum=0.
      do 311 i=1,nactdof
      sum=sum+(residual(i)**2)
  311 continue
      sum=sqrt(sum)
      if(iter .eq. 0)then
      suminit=sum 
      write(6,6105)suminit
      endif
      conv=sum/suminit
      write(6,6106)conv

*-----If solution converged, update resultants and move on

      if(conv .lt. tolrance)go to 500
 
  400 continue

*-----Solve equations

      if(itype .eq. 0)then
      do 320 i=1,nelembs
      a(i)=xkeffsave(i)
  320 continue
      do 322 i=1,15001
      maxa(i)=maxasave(i)
      ncolht(i)=ncolhtsave(i)
  322 continue
      endif

      if((itype .eq. 1 .and. ngravity .eq. 0) .or.
     $(itype .eq. 1 .and. ngravity .eq. 1 .and. kload .gt. 1))
     $then
      do 323 i=1,nelembs
      a(i)=xkeff(i)
  323 continue
      endif

      nnm=nactdof+1
      maxa(nactdof+1)=nelembs+1
      call colsol(a,residual,maxa,nactdof,nelembs,1,6)
      call colsol(a,residual,maxa,nactdof,nelembs,2,6)
      
*-----Add incremental displacements to displacement vector

      do 350 i=1,nactdof
      dispinc(i)=residual(i)
      dispold(i)=dispvctr(i)
      dispvctr(i)=dispvctr(i)+residual(i)
  350 continue

*-----For transient analysis update velocity and acceleration
*-----vectors

      if(ngravity .eq. 1 .and. kload .eq. 1)go to 354

      if((itype .eq. 0) .or.
     $(itype .eq. 1 .and. iter .eq. 0))then
      do 352 i=1,nactdof
      oldvel(i)=velvctr(i)
      oldacc(i)=accvctr(i)
      prevdisp(i)=dispold(i)
      vel=(gamma/(stepsize*beta))*dispvctr(i)+velvctr(i)
     $+(stepsize*(1.-gamma))*accvctr(i)-(gamma/(stepsize*beta))
     $*dispold(i)-(gamma/beta)*velvctr(i)
     $-(stepsize*gamma)*((1./(2.*beta))-1.)*accvctr(i)
      acc=(1./((stepsize**2)*beta))*dispvctr(i)-
     $(1./((stepsize**2)*beta))*dispold(i)-(1./(stepsize*beta))
     $*velvctr(i)-((1./(2.*beta))*(1.-(2.*beta)))*accvctr(i)
      velvctr(i)=vel
      accvctr(i)=acc
  352 continue
      endif

      if(itype .eq. 1 .and. iter .gt. 0)then
      do 353 i=1,nactdof
      vel=(gamma/(stepsize*beta))*dispvctr(i)+oldvel(i)
     $+(stepsize*(1.-gamma))*oldacc(i)-(gamma/(stepsize*beta))
     $*prevdisp(i)-(gamma/beta)*oldvel(i)
     $-(stepsize*gamma)*((1./(2.*beta))-1.)*oldacc(i)
      acc=(1./((stepsize**2)*beta))*dispvctr(i)-
     $(1./((stepsize**2)*beta))*prevdisp(i)-(1./(stepsize*beta))
     $*oldvel(i)-((1./(2.*beta))*(1.-(2.*beta)))*oldacc(i)
      velvctr(i)=vel
      accvctr(i)=acc
  353 continue
      endif

  354 continue

      if(itype .eq. 0)go to 500

*-----Increment iteration counter and check max iterations

      iter=iter+1
      if(iter .gt. 50)then
      write(6,6200)kload
      stop
      endif

*-----Solution was unconverged, so go back to top of loop

      go to 45

  500 continue

*-----Update beam element information

      if(nbeam .gt. 0)then
      call beamupdt(nbeam)
      endif

*-----Print converged solution

      if(itype .le. 1)then
      time=time+stepsize
      endif

      write(6,6110)kload,iter
      if(nprint .eq. 1)then
      write(6,6125)
      do 600 i=1,numnp
      write(6,6150)i,dispvctr(id(i,1)),dispvctr(id(i,2)),
     $dispvctr(id(i,3))
  600 continue
      endif
      if(ngravity .eq. 1 .and. itype .eq. 2)nprnt=1
c     call statedump(numnp,numfreq,eigvec,eigval,iaddr,
c    $              itype,kload,ngravity,time,nbeam,
c    $              mxfsiz,lenst,nprnt)

 1010 continue

*-----For ITYPE=2 call eigensolver

      if(itype .eq. 2)then
      if(itype .eq. 2 .and. ngravity .eq. 1 .and. kload .eq. 1)
     $go to 1000

      if(ngravity .eq. 1)then
      nprnt=2
      endif

      actmasx=0.
      actmasy=0.
      actmasz=0.
      do 1014 i=1,numnp
      if(id(i,1) .ne. 0)then
      nskip=0
      do 1011 j=1,i-1
      if(id(i,1) .eq. id(j,1))then
      nskip=1
      endif
 1011 continue
      if (nskip .eq. 0)then
      actmasx=actmasx+xmass(id(i,1))
      endif
      endif

      if(id(i,2) .ne. 0)then
      nskip=0
      do 1012 j=1,i-1
      if(id(i,2) .eq. id(j,2))then
      nskip=1
      endif
 1012 continue
      if (nskip .eq. 0)then
      actmasy=actmasy+xmass(id(i,2))
      endif
      endif

      if(id(i,3) .ne. 0)then
      nskip=0
      do 1013 j=1,i-1
      if(id(i,3) .eq. id(j,3))then
      nskip=1
      endif
 1013 continue
      if (nskip .eq. 0)then
      actmasz=actmasz+xmass(id(i,3))
      endif
      endif

 1014 continue

      write(6,6520)actmasx,actmasy,actmasz
      nc=min(2*numfreq,8+numfreq)
      if(nc .gt. nmasdof)nc=nmasdof
      rtol=1.e-12
      nnc=nc*(nc+1)/2
      nnm=nactdof+1
      iiout=8
      maxa(nactdof+1)=nelembs+1
      nrowdim=15000
      ncoldim=80
      nitem=16
      ifss=1
      ifpr=1
      call sspace(a,xmass,maxa,eigvec,eigval,ttt,w,ar,br,
     $vec,dd,rtolv,bup,blo,bupc,nactdof,nnm,nelembs,
     $nactdof,numfreq,rtol,nc,nnc,nitem,ifss,ifpr,
     $xnstif,iiout,nrowdim,ncoldim)
      write(6,6500) 
      do 1015 i=1,numfreq
      eigval(i)=(1./(8.*atan(1.0)))*sqrt(eigval(i))
      write(6,6510)i,eigval(i)
 1015 continue
c     call statedump(numnp,numfreq,eigvec,eigval,iaddr,
c    $itype,kload,ngravity,time,nbeam,
c    $mxfsiz,lenst,nprnt)

      endif  
 
*-----end time stepping do loop

 1000 continue

*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

*                           END OF TIME STEPPING LOOP   

*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

*-----End of main program

*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

*                              END OF MAIN PROGRAM   

*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

*-----Format statements, subroutine etc.

*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

*                               FORMAT STATEMENTS   

*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

*-----Read format statements

 5000 format(3i5,f10.0)
 5002 format(i5,f10.0,3i5)
 5003 format(2i5,3f10.0)
 5004 format(2i5,2f10.0)
 5010 format(5x,6i5,3f10.0)
 5045 format(3i5)
 5050 format(5x,5i5)
 5055 format(8f10.0,/,f10.0)
 5059 format(6f10.0)
 5060 format(i5)
 5061 format(3f10.0)
 5070 format(3i5)
 5080 format(i5,f10.0)

*-----Write format statements

 6000 format(//,'*****************************************************',
     $'**********************',///,
     $11x,'NONLINBLDG VERSION 1.0',18x,'DEC 2001',
     $////,29x,'DAVID MCCALLEN',/,
     $17x,'CENTER FOR COMPLEX DISTRIBUTED SYSTEMS',/,
     $17x,'LAWRENCE LIVERMORE NATIONAL LABORATORY',/,
     $25x,'UNIVERSITY OF CALIFORNIA',/,
     $26x,'LIVERMORE CALIFORNIA',////,'****************************',
     $'***********************************************',//)
 6005 format(//,'------------------------------------------------',
     $'-------------------------',//)
 6010 format(5x,'PROBLEM CONTROL DATA...',//)
 6020 format(5x,'Analysis type =',i3,/,
     $5x,'(itype=0 yields linear dynamic analysis)',/,
     $5x,'(itype=1 yields nonlinear dynamic analysis)',/,
     $5x,'(itype=2 yields eigensolution)',/,
     $5x,'(itype=3 yields nonlinear static analysis)',//)
 6030 format(5x,'Number of nodes =',i5,/,
     $5x,'Number of beam elements =',i5,/,
     $5x,'Number of truss elements =',i5,/,
     $5x,'Number of cable elements =',i5,/)
 6040 format(5x,'Number of load steps = ',i5,/,
     $5x,'Step size =',f10.3)
 6050 format(5x,'NODAL POINT BOUNDARY CONDITION ',
     $'CODES AND LOCATIONS...',//)
 6060 format(10x,'Boundary condition codes',11x,'X',11x,'Y',11x,'Z',/)
 6070 format(5x,6i5,5x,3e12.5)
 6096 format(5x,'BEAM ELEMENT INPUT...',//,5x,'Beam number',
     $5x,'I-node',5x,'J-node',5x,'K-node',5x,'Material Property',/)
 6097 format(9x,i5,6x,i5,6x,i5,6x,i5,13x,i5)
 6098 format(//,5x,'Beam material properties...',//,1x,'Mat#',
     $5x,'E',8x,'NU',10x,'A',10x,'YI',9x,'ZI',9x,'J',9x,'YSH',8x,
     $'ZSH',/)
 7000 format(1x,i2,1x,e10.3,1x,e10.3,1x,e10.3,1x,e10.3,1x,e10.3,
     $1x,e10.3,1x,e10.3,1x,e10.3/)
 6099 format(5x,'SOLUTION...',//)
 6100 format(5x,'Number of active degrees of freedom =',i6)
 6105 format(/,5x,'First Eucledian norm of residual vector =',e12.5,/)
 6106 format(5x,'Residual norm ratio =',e12.5)
 6110 format(//,5x,'Load step = ',i6,5x,
     $'Number of equilibrium iterations =',i5,//)
 6125 format(5x,'Node number',5x,'X-displacement',5x,
     $'Y-displacement',5x,'Z-displacement',//)
 6150 format(5x,i5,8x,e15.5,4x,e15.5,4x,e15.5)
 6200 format(1x,//,'No convergence after 50 iterations at load step',i5)
 6500 format(/,5x,'NATURAL FREQUENCIES OBTAINED FROM EIGSOLUTION',//,
     $5x,'Mode Number',12x,'Frequency',/)
 6510 format(5x,i6,14x,e12.4)
 6520 format(//,5x,'Active x mass = ',e15.5,5x,
     $ 'Active y mass = ',e15.5,5x,'Active z mass = ',e15.5)
c6535 format(e11.4,1x,e11.4)

*-----Close out Taurus data base

c     call rwabsf(iob12,ckeep,12,512,mxfsiz)

      stop
      end

*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

*                           SHAPE FUNCTIONS FOR BEAMS  

*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

*______________________________________________________________
*--------------------------------------------------------------
*
*                           SHAPE FUNCTIONS
*
*______________________________________________________________
*--------------------------------------------------------------

*-----These function subroutines provide the shape function and shape 
*-----function derivatives of the beam elements

      function xm1pp(x)
      xm1pp=1.5*x
      return
      end

      function xm2pp(x,xl)
      xm2pp=((3.*xl)/4.)*x-(xl/4.)
      return
      end

      function xm3pp(x)
      xm3pp=-1.5*x
      return
      end

      function xm4pp(x,xl)
      xm4pp=((3.*xl)/4.)*x+(xl/4.)
      return
      end

      function xl1pp(x)
      xl1pp=1.5*x
      return   
      end

      function xl2pp(x,xl)
      xl2pp=-((3.*xl)/4.)*x+(xl/4.)
      return
      end

      function xl3pp(x)
      xl3pp=-1.5*x
      return
      end

      function xl4pp(x,xl)
      xl4pp=-((3.*xl)/4.)*x-(xl/4.)
      return
      end

*_______________________________________________________________________________
*-------------------------------------------------------------------------------
*
*                                    BULLINIT
*
*_______________________________________________________________________________
*-------------------------------------------------------------------------------

      subroutine bullinit(numnp,nbeam,
     $                    numbmp,nfbrbmp,iaddr,itype,
     $                    ngravity,nplotcntrl)

*-----This subroutine initalizes the GRIZ plot file database


*-----Variable list...
 
*-----LIST OF VARIABLES AND ARRAYS PASSED INTO ROUTINE AS INPUT...


 
*-----LIST OF VARIABLES PASSED OUT OF ROUTINE AS OUTPUT...



*-----Perform common and array declarations

*-----Establishment of global parameters

      integer
     :        max_beams
     :       ,max_beam_fibers_points
     :       ,max_beam_properties
     :       ,max_beams_times_three
     :       ,max_dof
     :       ,max_iplot
     :       ,max_nodes


      parameter (
     :           max_beams              =    5000
     :          ,max_beam_fibers_points =      20
     :          ,max_beam_properties    =      50
     :          ,max_dof                =   30000
     :          ,max_iplot              =  100000
     :          ,max_nodes              =   10000
     :          )


      parameter (
     :           max_beams_times_three = max_beams * 3
     :          )

*-----Geometry and displacement block

      integer
     :        id( max_nodes, 6 )

      real
     :      accvctr( max_dof ), dispinc( max_dof ), dispvctr( max_dof )
     :     ,velvctr( max_dof ), xyzloc( max_nodes, 3 )


      common /geometry/
     :        id, xyzloc, dispvctr, dispinc, velvctr, accvctr

*-----Beam

      integer
     :        nbeamn( max_beams, 3 ), nbeamp( max_beams, 2 ),
     :        numfibrs( max_beam_fibers_points ),
     :        nyldink2( max_beams, 9, max_beam_properties ),
     :        nyldstk2( max_beams, 9, max_beam_properties )

      real
     :      area( max_beam_properties )
     :     ,beamdef( max_beams, 7 )
     :     ,bmdefmst( max_beams, 7 )
     :     ,delepsk2( max_beams, 9, max_beam_properties )
     :     ,e( max_beam_properties )
     :     ,elin( max_beam_properties )
     :     ,etab( max_beam_properties )
     :     ,fibrdefn( max_beam_fibers_points, max_beam_properties, 3 )
     :     ,fmstrk2( max_beams, 12 )
     :     ,forcek2( max_beams, 12 )
     :     ,gj( max_beam_properties )
     :     ,resk2( 12 )
     :     ,sheary( max_beam_properties )
     :     ,shearz( max_beam_properties )
     :     ,sigyld( max_beam_properties )
     :     ,straink2( max_beams, 9, max_beam_properties )
     :     ,stressk2( max_beams, 9, max_beam_properties )
     :     ,strmstk2( max_beams, 9, max_beam_properties )

      real
     :      strsmstk2( max_beams, 9, max_beam_properties )
     :     ,tik2( max_beams_times_three, 3 )
     :     ,tjk2( max_beams_times_three, 3 )
     :     ,tk2( 3, 3 )
     :     ,xj( max_beam_properties )
     :     ,xkk2( 12, 12 )
     :     ,xlamdab( max_beam_properties )
     :     ,xlk2( max_beams )
     :     ,xlnmsk2( max_beam_properties )
     :     ,xlobat( 3 )
     :     ,xlobfres( max_beams, 9, 4 )
     :     ,xlobwt( 3 )
     :     ,xnlnmsk2( max_beam_properties )
     :     ,xnu( max_beam_properties )
     :     ,yi( max_beam_properties )
     :     ,yldcntrk2( max_beams, 9, max_beam_properties, 2 )
     :     ,zi( max_beam_properties )


      common /beam/
     :        nbeamn,   nbeamp,   tk2,      tik2,      tjk2,
     :        xlk2,     xkk2,     area,     e, 
     :        yi,       zi,       xj,       xnu,       sheary,
     :        shearz,   fmstrk2,  gj,       resk2,     forcek2,
     :        xlnmsk2,  beamdef,  xnlnmsk2, elin,      sigyld,
     :        xlamdab,  etab,     straink2, xlobat,    fibrdefn,
     :        numfibrs, xlobfres, nyldink2, nyldstk2,  yldcntrk2,
     :        stressk2, delepsk2, strmstk2, strsmstk2, bmdefmst,
     :        xlobwt


      integer
     :        iob12( 1080 )
     :       ,iplot( max_iplot )

      real * 4
     :         plot( max_iplot )


      common /taurus/
     :        iob12, plot

      equivalence( plot, iplot )


c
c     Initialization transferred to "blkdat"
c
c     data iplot(16),iplot(18),iplot(21),iplot(31),iplot(34) 
c    $     / 4,2,1,0,0 /
c     data iplot(19) / 0 /


*-----Write problem control data

      numprop=numbmp+nfbrbmp
      numstick=nbeam
     
      iplot(17)=numnp
      iplot(29)=numstick
      iplot(30)=numprop
      if(nplotcntrl .eq. 1)then
      iplot(31)=6
      else
      iplot(31)=0.
      endif
      iaddr=0
c     call wrabsf(iob12,iplot,64,iaddr)
      iaddr=64
      nshift=0

*-----Write initial model geometry

      do 100 i=1,numnp
      plot(3*i)=xyzloc(i,3)
      plot(3*i-1)=xyzloc(i,2)
      plot(3*i-2)=xyzloc(i,1)
  100 continue
      length=3*numnp
c     call wrabsf(iob12,plot,length,iaddr)
      iaddr=iaddr+length

*-----Write element connectivities and material numbers

      if(nbeam .gt. 0)then
      do 200 i=1,nbeam
      iplot(6*i+nshift)=nbeamp(i,1)
      iplot(6*i-1+nshift)=0
      iplot(6*i-2+nshift)=0
      iplot(6*i-3+nshift)=nbeamn(i,3)
      iplot(6*i-4+nshift)=nbeamn(i,2)
      iplot(6*i-5+nshift)=nbeamn(i,1)
  200 continue
      nshift=nshift+6*nbeam
      endif

      length=nshift
c     call wrabsf(iob12,iplot,length,iaddr)
      iaddr=iaddr+length

      return
      end

*_______________________________________________________________________________
*-------------------------------------------------------------------------------
*
*                                   STATEDUMP
*
*_______________________________________________________________________________
*-------------------------------------------------------------------------------

      subroutine statedump(numnp,numfreq,eigvec,eigval,iaddr,
     $                     itype,kload,ngravity,time,
     $                     nbeam,mxfsiz,lenst,nprnt)

*-----This subroutine dumps a state to the GRIZ database


*-----Variable list...
 
*-----LIST OF VARIABLES AND ARRAYS PASSED INTO ROUTINE AS INPUT...


 
*-----LIST OF VARIABLES PASSED OUT OF ROUTINE AS OUTPUT...



*-----Perform common and array declarations

*-----Establishment of global parameters

      integer
     :        max_beams
     :       ,max_beam_fibers_points
     :       ,max_beam_properties
     :       ,max_beams_times_three
     :       ,max_dof
     :       ,max_iplot
     :       ,max_nodes


      parameter (
     :           max_beams              =    5000
     :          ,max_beam_fibers_points =      20
     :          ,max_beam_properties    =      50
     :          ,max_dof                =   30000
     :          ,max_iplot              =  100000
     :          ,max_nodes              =   10000
     :          )


      parameter (
     :           max_beams_times_three = max_beams * 3
     :          )

*-----Geometry and displacement block

      integer
     :        id( max_nodes, 6 )

      real
     :      accvctr( max_dof ), dispinc( max_dof ), dispvctr( max_dof )
     :     ,velvctr( max_dof ), xyzloc( max_nodes, 3 )


      common /geometry/
     :        id, xyzloc, dispvctr, dispinc, velvctr, accvctr

*-----Beam

      integer
     :        nbeamn( max_beams, 3 ), nbeamp( max_beams, 2 ),
     :        numfibrs( max_beam_fibers_points ),
     :        nyldink2( max_beams, 9, max_beam_properties ),
     :        nyldstk2( max_beams, 9, max_beam_properties )

      real
     :      area( max_beam_properties )
     :     ,beamdef( max_beams, 7 )
     :     ,bmdefmst( max_beams, 7 )
     :     ,delepsk2( max_beams, 9, max_beam_properties )
     :     ,e( max_beam_properties )
     :     ,elin( max_beam_properties )
     :     ,etab( max_beam_properties )
     :     ,fibrdefn( max_beam_fibers_points, max_beam_properties, 3 )
     :     ,fmstrk2( max_beams, 12 )
     :     ,forcek2( max_beams, 12 )
     :     ,gj( max_beam_properties )
     :     ,resk2( 12 )
     :     ,sheary( max_beam_properties )
     :     ,shearz( max_beam_properties )
     :     ,sigyld( max_beam_properties )
     :     ,straink2( max_beams, 9, max_beam_properties )
     :     ,stressk2( max_beams, 9, max_beam_properties )
     :     ,strmstk2( max_beams, 9, max_beam_properties )

      real
     :      strsmstk2( max_beams, 9, max_beam_properties )
     :     ,tik2( max_beams_times_three, 3 )
     :     ,tjk2( max_beams_times_three, 3 )
     :     ,tk2( 3, 3 )
     :     ,xj( max_beam_properties )
     :     ,xkk2( 12, 12 )
     :     ,xlamdab( max_beam_properties )
     :     ,xlk2( max_beams )
     :     ,xlnmsk2( max_beam_properties )
     :     ,xlobat( 3 )
     :     ,xlobfres( max_beams, 9, 4 )
     :     ,xlobwt( 3 )
     :     ,xnlnmsk2( max_beam_properties )
     :     ,xnu( max_beam_properties )
     :     ,yi( max_beam_properties )
     :     ,yldcntrk2( max_beams, 9, max_beam_properties, 2 )
     :     ,zi( max_beam_properties )


      common /beam/
     :        nbeamn,   nbeamp,   tk2,      tik2,      tjk2,
     :        xlk2,     xkk2,     area,     e, 
     :        yi,       zi,       xj,       xnu,       sheary,
     :        shearz,   fmstrk2,  gj,       resk2,     forcek2,
     :        xlnmsk2,  beamdef,  xnlnmsk2, elin,      sigyld,
     :        xlamdab,  etab,     straink2, xlobat,    fibrdefn,
     :        numfibrs, xlobfres, nyldink2, nyldstk2,  yldcntrk2,
     :        stressk2, delepsk2, strmstk2, strsmstk2, bmdefmst,
     :        xlobwt


      integer
     :        iob12( 1080 )
     :       ,iplot( max_iplot )

      real * 4
     :         plot( max_iplot )


      common /taurus/
     :        iob12, plot

      equivalence( plot, iplot )


      dimension eigvec(15000,80), eigval(80)


*-----If state will not fit in current plot file, jump to next file

      nxtfil=(iaddr/mxfsiz+1)*mxfsiz
      if((iaddr+lenst) .gt. nxtfil)iaddr=nxtfil

      if(itype .ne. 2)go to 110
      if(ngravity .eq. 1 .and. nprnt .eq. 1)then
      plot(1)=0
      do 25 j=1,numnp
      plot(3*j+1)=xyzloc(j,3)+dispvctr(id(j,3))
      plot(3*j)=xyzloc(j,2)+dispvctr(id(j,2))
      plot(3*j-1)=xyzloc(j,1)+dispvctr(id(j,1))
   25 continue
      length=3*numnp+1
c     call wrabsf(iob12,plot,length,iaddr)
      iaddr=iaddr+length
      else
      do 100 i=1,numfreq
      plot(1)=eigval(i)
      do 50 j=1,numnp
      plot(3*j+1)=xyzloc(j,3)+eigvec(id(j,3),i)
      plot(3*j)=xyzloc(j,2)+eigvec(id(j,2),i)
      plot(3*j-1)=xyzloc(j,1)+eigvec(id(j,1),i)
   50 continue
      length=3*numnp+1
c     call wrabsf(iob12,plot,length,iaddr)
      iaddr=iaddr+length
  100 continue
      endif
      go to 260

  110 continue    

      if(itype .eq. 3)then
      plot(1)=kload
      elseif((itype .le. 1 .and. ngravity .eq. 0) .or.
     $(itype .le. 1 .and. ngravity .eq. 1 .and. kload .gt. 1))
     $then
      plot(1)=time
      endif
      do 150 j=1,numnp
      plot(3*j+1)=xyzloc(j,3)+dispvctr(id(j,3))
      plot(3*j)=xyzloc(j,2)+dispvctr(id(j,2))
      plot(3*j-1)=xyzloc(j,1)+dispvctr(id(j,1))
  150 continue
      length=3*numnp+1
c     call wrabsf(iob12,plot,length,iaddr)
      iaddr=iaddr+length

*-----Write element force resultants

      if(nbeam .gt. 0)then
      do 250 i=1,nbeam
      plot(6*i-5)=fmstrk2(i,1)
      plot(6*i-4)=fmstrk2(i,2)
      plot(6*i-3)=fmstrk2(i,3)
      plot(6*i-2)=fmstrk2(i,5)
      plot(6*i-1)=fmstrk2(i,6)
      plot(6*i)=fmstrk2(i,4)
  250 continue
      length=6*nbeam
c     call wrabsf(iob12,plot,length,iaddr)
      iaddr=iaddr+length
      endif

  260 continue
           
      return
      end

*_______________________________________________________________________________
*-------------------------------------------------------------------------------
*
*                                    BEAMINIT
*
*_______________________________________________________________________________
*-------------------------------------------------------------------------------

      subroutine beaminit(nbeam)

*-----This subroutine initializes the beam elements


*-----Variable list...
 
*-----LIST OF VARIABLES AND ARRAYS PASSED INTO ROUTINE AS INPUT...

 
*-----LIST OF VARIABLES PASSED OUT OF ROUTINE AS OUTPUT...


*-----Perform common and array declarations

*-----Establishment of global parameters

      integer
     :        max_beams
     :       ,max_beam_fibers_points
     :       ,max_beam_properties
     :       ,max_beams_times_three
     :       ,max_dof
     :       ,max_nodes


      parameter (
     :           max_beams              =    5000
     :          ,max_beam_fibers_points =      20
     :          ,max_beam_properties    =      50
     :          ,max_dof                =   30000
     :          ,max_nodes              =   10000
     :          )


      parameter (
     :           max_beams_times_three = max_beams * 3
     :          )

*-----Geometry and displacement block

      integer
     :        id( max_nodes, 6 )

      real
     :      accvctr( max_dof ), dispinc( max_dof ), dispvctr( max_dof )
     :     ,velvctr( max_dof ), xyzloc( max_nodes, 3 )


      common /geometry/
     :        id, xyzloc, dispvctr, dispinc, velvctr, accvctr

*-----Beam

      integer
     :        nbeamn( max_beams, 3 ), nbeamp( max_beams, 2 ),
     :        numfibrs( max_beam_fibers_points ),
     :        nyldink2( max_beams, 9, max_beam_properties ),
     :        nyldstk2( max_beams, 9, max_beam_properties )

      real
     :      area( max_beam_properties )
     :     ,beamdef( max_beams, 7 )
     :     ,bmdefmst( max_beams, 7 )
     :     ,delepsk2( max_beams, 9, max_beam_properties )
     :     ,e( max_beam_properties )
     :     ,elin( max_beam_properties )
     :     ,etab( max_beam_properties )
     :     ,fibrdefn( max_beam_fibers_points, max_beam_properties, 3 )
     :     ,fmstrk2( max_beams, 12 )
     :     ,forcek2( max_beams, 12 )
     :     ,gj( max_beam_properties )
     :     ,resk2( 12 )
     :     ,sheary( max_beam_properties )
     :     ,shearz( max_beam_properties )
     :     ,sigyld( max_beam_properties )
     :     ,straink2( max_beams, 9, max_beam_properties )
     :     ,stressk2( max_beams, 9, max_beam_properties )
     :     ,strmstk2( max_beams, 9, max_beam_properties )

      real
     :      strsmstk2( max_beams, 9, max_beam_properties )
     :     ,tik2( max_beams_times_three, 3 )
     :     ,tjk2( max_beams_times_three, 3 )
     :     ,tk2( 3, 3 )
     :     ,xj( max_beam_properties )
     :     ,xkk2( 12, 12 )
     :     ,xlamdab( max_beam_properties )
     :     ,xlk2( max_beams )
     :     ,xlnmsk2( max_beam_properties )
     :     ,xlobat( 3 )
     :     ,xlobfres( max_beams, 9, 4 )
     :     ,xlobwt( 3 )
     :     ,xnlnmsk2( max_beam_properties )
     :     ,xnu( max_beam_properties )
     :     ,yi( max_beam_properties )
     :     ,yldcntrk2( max_beams, 9, max_beam_properties, 2 )
     :     ,zi( max_beam_properties )


      common /beam/
     :        nbeamn,   nbeamp,   tk2,      tik2,      tjk2,
     :        xlk2,     xkk2,     area,     e, 
     :        yi,       zi,       xj,       xnu,       sheary,
     :        shearz,   fmstrk2,  gj,       resk2,     forcek2,
     :        xlnmsk2,  beamdef,  xnlnmsk2, elin,      sigyld,
     :        xlamdab,  etab,     straink2, xlobat,    fibrdefn,
     :        numfibrs, xlobfres, nyldink2, nyldstk2,  yldcntrk2,
     :        stressk2, delepsk2, strmstk2, strsmstk2, bmdefmst,
     :        xlobwt


c     call bounds_check( 'beaminit'
c    :                  ,'k2'
c    :                  ,1
c    :                  ,nbeam
c    :                  ,max_beams
c    :                 )

      do 200 k2=1,nbeam

*-----Get locations of element i,j and k nodes

      nodei=nbeamn(k2,1)
      nodej=nbeamn(k2,2)
      nodek=nbeamn(k2,3)

c     call bounds_check( 'beaminit'
c    :                  ,'nodei'
c    :                  ,1
c    :                  ,nodei
c    :                  ,max_nodes
c    :                 )

c     call bounds_check( 'beaminit'
c    :                  ,'nodej'
c    :                  ,1
c    :                  ,nodej
c    :                  ,max_nodes
c    :                 )

c     call bounds_check( 'beaminit'
c    :                  ,'nodek'
c    :                  ,1
c    :                  ,nodek
c    :                  ,max_nodes
c    :                 )

      xloci=xyzloc(nodei,1)
      yloci=xyzloc(nodei,2)
      zloci=xyzloc(nodei,3)
      xlocj=xyzloc(nodej,1)
      ylocj=xyzloc(nodej,2)
      zlocj=xyzloc(nodej,3)
      xlock=xyzloc(nodek,1)
      ylock=xyzloc(nodek,2)
      zlock=xyzloc(nodek,3)
      delx=xlocj-xloci
      dely=ylocj-yloci
      delz=zlocj-zloci
      
*-----Determine length of beam element

      xlk2(k2)=sqrt((delx)**2+(dely)**2
     $+(delz)**2)

*-----Initialize orientation of nodal coordinate axes

      c1=((ylocj-yloci)*(zlock-zloci))-
     $   ((zlocj-zloci)*(ylock-yloci))
      c2=(-(xlocj-xloci)*(zlock-zloci))+
     $   ((zlocj-zloci)*(xlock-xloci))
      c3=((xlocj-xloci)*(ylock-yloci))-
     $   ((ylocj-yloci)*(xlock-xloci))
   
      denom=sqrt(c1**2+c2**2+c3**2)

      loc=3*k2-2   

      tik2(loc,1)=delx/xlk2(k2)
      tik2(loc,2)=dely/xlk2(k2)
      tik2(loc,3)=delz/xlk2(k2)
      tjk2(loc,1)=tik2(loc,1)
      tjk2(loc,2)=tik2(loc,2)
      tjk2(loc,3)=tik2(loc,3)
      tik2(loc+2,1)=c1/denom
      tik2(loc+2,2)=c2/denom           
      tik2(loc+2,3)=c3/denom
      tik2(loc+1,1)=tik2(loc+2,2)*tik2(loc,3)-tik2(loc+2,3)*tik2(loc,2)
      tik2(loc+1,2)=-tik2(loc+2,1)*tik2(loc,3)+tik2(loc+2,3)*tik2(loc,1)
      tik2(loc+1,3)=tik2(loc+2,1)*tik2(loc,2)-tik2(loc+2,2)*tik2(loc,1)
      tjk2(loc+2,1)=c1/denom
      tjk2(loc+2,2)=c2/denom           
      tjk2(loc+2,3)=c3/denom
      tjk2(loc+1,1)=tik2(loc+1,1)
      tjk2(loc+1,2)=tik2(loc+1,2)
      tjk2(loc+1,3)=tik2(loc+1,3)

*-----For the fiber elements, initialize the yield center information 

      if(nbeamp(k2,2) .gt. 0)then
      do 180 i=1,3 
      do 175 k=1,4
      xlobfres(k2,i,k)=0.
  175 continue
      if(sigyld(nbeamp(k2,1)) .gt. 0)then
      do 177 j=1,numfibrs(nbeamp(k2,2))
      nyldink2(k2,i,j)=0
      nyldstk2(k2,i,j)=0
      yldcntrk2(k2,i,j,1)=0.
      yldcntrk2(k2,i,j,2)=0.
  177 continue 
      endif        
  180 continue
      endif

  200 continue

      return
      end

*_______________________________________________________________________________
*-------------------------------------------------------------------------------
*                                 
*                                     TRNINIT
*
*_______________________________________________________________________________
*-------------------------------------------------------------------------------

      subroutine trninit(k2) 

*-----This subroutine initializes the beam element direction cosines


*-----Variable list...
 
*-----LIST OF VARIABLES AND ARRAYS PASSED INTO ROUTINE AS INPUT...

 
*-----LIST OF VARIABLES PASSED OUT OF ROUTINE AS OUTPUT...


*-----Perform common and array declarations

*----Establishment of global parameters

      integer
     :        max_beams
     :       ,max_beam_fibers_points
     :       ,max_beam_properties
     :       ,max_beams_times_three
     :       ,max_dof
     :       ,max_nodes


      parameter (
     :           max_beams              =    5000
     :          ,max_beam_fibers_points =      20
     :          ,max_beam_properties    =      50
     :          ,max_dof                =   30000
     :          ,max_nodes              =   10000
     :          )


      parameter (
     :           max_beams_times_three = max_beams * 3
     :          )

*----Geometry and displacement block

      integer
     :        id( max_nodes, 6 )

      real
     :      accvctr( max_dof ), dispinc( max_dof ), dispvctr( max_dof )
     :     ,velvctr( max_dof ), xyzloc( max_nodes, 3 )


      common /geometry/
     :        id, xyzloc, dispvctr, dispinc, velvctr, accvctr

*-----Beam

      integer
     :        nbeamn( max_beams, 3 ), nbeamp( max_beams, 2 ),
     :        numfibrs( max_beam_fibers_points ),
     :        nyldink2( max_beams, 9, max_beam_properties ),
     :        nyldstk2( max_beams, 9, max_beam_properties )

      real
     :      area( max_beam_properties )
     :     ,beamdef( max_beams, 7 )
     :     ,bmdefmst( max_beams, 7 )
     :     ,delepsk2( max_beams, 9, max_beam_properties )
     :     ,e( max_beam_properties )
     :     ,elin( max_beam_properties )
     :     ,etab( max_beam_properties )
     :     ,fibrdefn( max_beam_fibers_points, max_beam_properties, 3 )
     :     ,fmstrk2( max_beams, 12 )
     :     ,forcek2( max_beams, 12 )
     :     ,gj( max_beam_properties )
     :     ,resk2( 12 )
     :     ,sheary( max_beam_properties )
     :     ,shearz( max_beam_properties )
     :     ,sigyld( max_beam_properties )
     :     ,straink2( max_beams, 9, max_beam_properties )
     :     ,stressk2( max_beams, 9, max_beam_properties )
     :     ,strmstk2( max_beams, 9, max_beam_properties )

      real
     :      strsmstk2( max_beams, 9, max_beam_properties )
     :     ,tik2( max_beams_times_three, 3 )
     :     ,tjk2( max_beams_times_three, 3 )
     :     ,tk2( 3, 3 )
     :     ,xj( max_beam_properties )
     :     ,xkk2( 12, 12 )
     :     ,xlamdab( max_beam_properties )
     :     ,xlk2( max_beams )
     :     ,xlnmsk2( max_beam_properties )
     :     ,xlobat( 3 )
     :     ,xlobfres( max_beams, 9, 4 )
     :     ,xlobwt( 3 )
     :     ,xnlnmsk2( max_beam_properties )
     :     ,xnu( max_beam_properties )
     :     ,yi( max_beam_properties )
     :     ,yldcntrk2( max_beams, 9, max_beam_properties, 2 )
     :     ,zi( max_beam_properties )


      common /beam/
     :        nbeamn,   nbeamp,   tk2,      tik2,      tjk2,
     :        xlk2,     xkk2,     area,     e, 
     :        yi,       zi,       xj,       xnu,       sheary,
     :        shearz,   fmstrk2,  gj,       resk2,     forcek2,
     :        xlnmsk2,  beamdef,  xnlnmsk2, elin,      sigyld,
     :        xlamdab,  etab,     straink2, xlobat,    fibrdefn,
     :        numfibrs, xlobfres, nyldink2, nyldstk2,  yldcntrk2,
     :        stressk2, delepsk2, strmstk2, strsmstk2, bmdefmst,
     :        xlobwt


      nodei=nbeamn(k2,1)
      nodej=nbeamn(k2,2)
      nodek=nbeamn(k2,3)
      xloci=xyzloc(nodei,1)
      yloci=xyzloc(nodei,2)
      zloci=xyzloc(nodei,3)
      xlocj=xyzloc(nodej,1)
      ylocj=xyzloc(nodej,2)
      zlocj=xyzloc(nodej,3)
      xlock=xyzloc(nodek,1)
      ylock=xyzloc(nodek,2)
      zlock=xyzloc(nodek,3)
      delx=xlocj-xloci
      dely=ylocj-yloci
      delz=zlocj-zloci
      
*-----Determine length of beam element

      xlk2(k2)=sqrt((delx)**2+(dely)**2
     $+(delz)**2)

*-----Initialize orientation of nodal coordinate axes

      c1=((ylocj-yloci)*(zlock-zloci))-
     $   ((zlocj-zloci)*(ylock-yloci))
      c2=(-(xlocj-xloci)*(zlock-zloci))+
     $   ((zlocj-zloci)*(xlock-xloci))
      c3=((xlocj-xloci)*(ylock-yloci))-
     $   ((ylocj-yloci)*(xlock-xloci))
   
      denom=sqrt(c1**2+c2**2+c3**2)

      tk2(1,1)=delx/xlk2(k2)
      tk2(1,2)=dely/xlk2(k2)
      tk2(1,3)=delz/xlk2(k2)
      tk2(3,1)=c1/denom
      tk2(3,2)=c2/denom
      tk2(3,3)=c3/denom
      tk2(2,1)=tk2(3,2)*tk2(1,3)-tk2(3,3)*tk2(1,2)
      tk2(2,2)=-tk2(3,1)*tk2(1,3)+tk2(3,3)*tk2(1,1)
      tk2(2,3)=tk2(3,1)*tk2(1,2)-tk2(3,2)*tk2(1,1)

      return
      end

*_______________________________________________________________________________
*-------------------------------------------------------------------------------
*                                 
*                                     TRNBEAM
*
*_______________________________________________________________________________
*-------------------------------------------------------------------------------

      subroutine trnbeam(k2,iter,itype,kload) 

*-----This subroutine computes the updated element direction cosines and
*-----element angular deformations

*-----Variable list...
 
*-----LIST OF VARIABLES AND ARRAYS PASSED INTO ROUTINE AS INPUT...


*-----LIST OF VARIABLES PASSED OUT OF ROUTINE AS OUTPUT...


*-----Perform common and array declarations

*-----Establishment of global parameters

      integer
     :        max_beams
     :       ,max_beam_fibers_points
     :       ,max_beam_properties
     :       ,max_beams_times_three
     :       ,max_dof
     :       ,max_nodes


      parameter (
     :           max_beams              =    5000
     :          ,max_beam_fibers_points =      20
     :          ,max_beam_properties    =      50
     :          ,max_dof                =   30000
     :          ,max_nodes              =   10000
     :          )


      parameter (
     :           max_beams_times_three = max_beams * 3
     :          )

*-----Geometry and displacement block

      integer
     :        id( max_nodes, 6 )

      real
     :      accvctr( max_dof ), dispinc( max_dof ), dispvctr( max_dof )
     :     ,velvctr( max_dof ), xyzloc( max_nodes, 3 )


      common /geometry/
     :        id, xyzloc, dispvctr, dispinc, velvctr, accvctr

*----Bbeam

      integer
     :        nbeamn( max_beams, 3 ), nbeamp( max_beams, 2 ),
     :        numfibrs( max_beam_fibers_points ),
     :        nyldink2( max_beams, 9, max_beam_properties ),
     :        nyldstk2( max_beams, 9, max_beam_properties )

      real
     :      area( max_beam_properties )
     :     ,beamdef( max_beams, 7 )
     :     ,bmdefmst( max_beams, 7 )
     :     ,delepsk2( max_beams, 9, max_beam_properties )
     :     ,e( max_beam_properties )
     :     ,elin( max_beam_properties )
     :     ,etab( max_beam_properties )
     :     ,fibrdefn( max_beam_fibers_points, max_beam_properties, 3 )
     :     ,fmstrk2( max_beams, 12 )
     :     ,forcek2( max_beams, 12 )
     :     ,gj( max_beam_properties )
     :     ,resk2( 12 )
     :     ,sheary( max_beam_properties )
     :     ,shearz( max_beam_properties )
     :     ,sigyld( max_beam_properties )
     :     ,straink2( max_beams, 9, max_beam_properties )
     :     ,stressk2( max_beams, 9, max_beam_properties )
     :     ,strmstk2( max_beams, 9, max_beam_properties )

      real
     :      strsmstk2( max_beams, 9, max_beam_properties )
     :     ,tik2( max_beams_times_three, 3 )
     :     ,tjk2( max_beams_times_three, 3 )
     :     ,tk2( 3, 3 )
     :     ,xj( max_beam_properties )
     :     ,xkk2( 12, 12 )
     :     ,xlamdab( max_beam_properties )
     :     ,xlk2( max_beams )
     :     ,xlnmsk2( max_beam_properties )
     :     ,xlobat( 3 )
     :     ,xlobfres( max_beams, 9, 4 )
     :     ,xlobwt( 3 )
     :     ,xnlnmsk2( max_beam_properties )
     :     ,xnu( max_beam_properties )
     :     ,yi( max_beam_properties )
     :     ,yldcntrk2( max_beams, 9, max_beam_properties, 2 )
     :     ,zi( max_beam_properties )


      common /beam/
     :        nbeamn,   nbeamp,   tk2,      tik2,      tjk2,
     :        xlk2,     xkk2,     area,     e, 
     :        yi,       zi,       xj,       xnu,       sheary,
     :        shearz,   fmstrk2,  gj,       resk2,     forcek2,
     :        xlnmsk2,  beamdef,  xnlnmsk2, elin,      sigyld,
     :        xlamdab,  etab,     straink2, xlobat,    fibrdefn,
     :        numfibrs, xlobfres, nyldink2, nyldstk2,  yldcntrk2,
     :        stressk2, delepsk2, strmstk2, strsmstk2, bmdefmst,
     :        xlobwt


      dimension dthxipp(3), dthxjpp(3), xnewtit(3,3), xnewtjt(3,3),
     $          dthxi(3), dthxj(3), anglesi(3,3), anglesj(3,3), 
     $          tit(3,3), tjt(3,3), gammai(2), gammaj(2)


*-----Zero out working matrices

      do 5 i=1,2
      gammai(i)=0.
      gammaj(i)=0.
    5 continue
    
      do 20 i=1,3
      do 10 j=1,3
      dthxipp(i)=0.
      dthxjpp(j)=0.
      xnewtit(i,j)=0.
      xnewtjt(i,j)=0.
   10 continue
   20 continue

      tau3i=0.
      tau3j=0.

*-----Get current location of element nodes
      
      nodei=nbeamn(k2,1)
      nodej=nbeamn(k2,2)
      if(id(nodei,1) .gt. 0)then
      xloci=xyzloc(nodei,1)+dispvctr(id(nodei,1))
      else
      xloci=xyzloc(nodei,1)
      endif
      if(id(nodei,2) .gt. 0)then
      yloci=xyzloc(nodei,2)+dispvctr(id(nodei,2))
      else
      yloci=xyzloc(nodei,2)
      endif
      if(id(nodei,3) .gt. 0)then
      zloci=xyzloc(nodei,3)+dispvctr(id(nodei,3))
      else
      zloci=xyzloc(nodei,3)
      endif
      if(id(nodej,1) .gt. 0)then
      xlocj=xyzloc(nodej,1)+dispvctr(id(nodej,1))
      else
      xlocj=xyzloc(nodej,1)
      endif
      if(id(nodej,2) .gt. 0)then
      ylocj=xyzloc(nodej,2)+dispvctr(id(nodej,2))
      else 
      ylocj=xyzloc(nodej,2)
      endif
      if(id(nodej,3) .gt. 0)then
      zlocj=xyzloc(nodej,3)+dispvctr(id(nodej,3))
      else
      zlocj=xyzloc(nodej,3)
      endif  

*-----Get element coordinate changes and change of length

      delx=xlocj-xloci
      dely=ylocj-yloci
      delz=zlocj-zloci
      xlength=sqrt(delx**2+dely**2+delz**2)
      loc=3*k2-2

      if((itype .eq. 3 .and. iter .eq. 0) .or.
     $(itype .eq. 1 .and. iter .eq. 0) .or.
     $(itype .eq. 0 .and. kload .eq. 1) .or.
     $(itype .eq. 2 .and. iter .eq. 0))go to 95

      dthxi(1)=dispinc(id(nodei,4))
      dthxi(2)=dispinc(id(nodei,5))
      dthxi(3)=dispinc(id(nodei,6))
      dthxj(1)=dispinc(id(nodej,4))
      dthxj(2)=dispinc(id(nodej,5))
      dthxj(3)=dispinc(id(nodej,6))
   
      do 40 i=1,3
      do 30 j=1,3
      tit(j,i)=tik2(loc+i-1,j)
      tjt(j,i)=tjk2(loc+i-1,j)
   30 continue
   40 continue

*-----Translate rotation increments to nodal coordinates

      do 48 i=1,3
      do 45 j=1,3
      dthxipp(i)=dthxipp(i)+tik2(loc+i-1,j)*dthxi(j)
      dthxjpp(i)=dthxjpp(i)+tjk2(loc+i-1,j)*dthxj(j)
   45 continue
   48 continue

      anglesi(1,1)=sqrt(1.-(dthxipp(3)**2)-(dthxipp(2)**2))
      anglesi(2,1)=dthxipp(3)
      anglesi(3,1)=-dthxipp(2)
      anglesi(1,2)=-dthxipp(3)
      anglesi(2,2)=sqrt(1.-(dthxipp(3)**2)-(dthxipp(1)**2)) 
      anglesi(3,2)=dthxipp(1)
      anglesi(1,3)=dthxipp(2)
      anglesi(2,3)=-dthxipp(1)
      anglesi(3,3)=sqrt(1.-(dthxipp(2)**2)-(dthxipp(1)**2))

      anglesj(1,1)=sqrt(1.-(dthxjpp(3)**2)-(dthxjpp(2)**2))
      anglesj(2,1)=dthxjpp(3)
      anglesj(3,1)=-dthxjpp(2)
      anglesj(1,2)=-dthxjpp(3)
      anglesj(2,2)=sqrt(1.-(dthxjpp(3)**2)-(dthxjpp(1)**2)) 
      anglesj(3,2)=dthxjpp(1)
      anglesj(1,3)=dthxjpp(2)
      anglesj(2,3)=-dthxjpp(1)
      anglesj(3,3)=sqrt(1.-(dthxjpp(2)**2)-(dthxjpp(1)**2))

*-----Generate new transformation matrix at node I and node J

      do 70 i=1,3
      do 60 j=1,3
      do 50 k=1,3
      xnewtit(i,j)=xnewtit(i,j)+tit(i,k)*anglesi(k,j)
      xnewtjt(i,j)=xnewtjt(i,j)+tjt(i,k)*anglesj(k,j)
   50 continue
   60 continue
   70 continue

      do 90 i=1,3
      do 80 j=1,3
      tik2(loc+i-1,j)=xnewtit(j,i)
      tjk2(loc+i-1,j)=xnewtjt(j,i)
   80 continue
   90 continue

*-----Get element updated coordinate system

*-----To get direction cosines of element y axis, establish a new "k" node
*-----which is determined from the sum of the nodal "j" direction vectors

   95 continue

      xmag=sqrt(((tik2(loc+1,1)+tjk2(loc+1,1))**2)+
     $((tik2(loc+1,2)+tjk2(loc+1,2))**2)+
     $((tik2(loc+1,3)+tjk2(loc+1,3))**2))

      etax=xlength*((tik2(loc+1,1)+tjk2(loc+1,1))/xmag)
      etay=xlength*((tik2(loc+1,2)+tjk2(loc+1,2))/xmag)
      etaz=xlength*((tik2(loc+1,3)+tjk2(loc+1,3))/xmag)

      xlock=xloci+etax
      ylock=yloci+etay
      zlock=zloci+etaz

      c1=((ylocj-yloci)*(zlock-zloci))-
     $   ((zlocj-zloci)*(ylock-yloci))
      c2=(-(xlocj-xloci)*(zlock-zloci))+
     $   ((zlocj-zloci)*(xlock-xloci))
      c3=((xlocj-xloci)*(ylock-yloci))-
     $   ((ylocj-yloci)*(xlock-xloci))
   
      denom=sqrt(c1**2+c2**2+c3**2)

      tk2(1,1)=delx/xlength
      tk2(1,2)=dely/xlength
      tk2(1,3)=delz/xlength
      tk2(3,1)=c1/denom
      tk2(3,2)=c2/denom           
      tk2(3,3)=c3/denom
      tk2(2,1)=tk2(3,2)*tk2(1,3)-tk2(3,3)*tk2(1,2)
      tk2(2,2)=-tk2(3,1)*tk2(1,3)+tk2(3,3)*tk2(1,1)
      tk2(2,3)=tk2(3,1)*tk2(1,2)-tk2(3,2)*tk2(1,1)

      if((itype .eq. 3 .and. iter .eq. 0) .or.
     $(itype .eq. 1 .and. iter .eq. 0) .or.
     $(itype .eq. 0 .and. kload .eq. 1) .or.
     $(itype .eq. 2 .and. iter .eq. 0))go to 250

*-----Determine element deformations based on the updated element
*-----and nodal coordinate systems

*-----Translate nodal "i" and "j" vectors into the element convected
*-----system

      do 110 i=1,2
      do 100 j=1,3
      gammai(i)=gammai(i)+tk2(i+1,j)*tik2(loc,j)
      gammaj(i)=gammaj(i)+tk2(i+1,j)*tjk2(loc,j)
  100 continue
  110 continue

      do 120 i=1,3
      tau3i=tau3i+tk2(3,i)*tik2(loc+1,i)
      tau3j=tau3j+tk2(3,i)*tjk2(loc+1,i)
  120 continue   


*-----Note: beam element strain measures are in the following order
*           beamdef(k2,1)=element rotation about x' axis at node I
*           beamdef(k2,2)=element rotation about y' axis at node I
*           beamdef(k2,3)=element rotation about z' axis at node I
*           beamdef(k2,4)=element rotation about x' axis at node J
*           beamdef(k2,5)=element rotation about y' axis at node J
*           beamdef(k2,6)=element rotation about z' axis at node J
*           beamdef(k2,7)=element elongation as measured by longitudinal 
*                         displacement at node J

      beamdef(k2,1)=tau3i
      beamdef(k2,2)=-gammai(2)
      beamdef(k2,3)=gammai(1)
      beamdef(k2,4)=tau3j
      beamdef(k2,5)=-gammaj(2)
      beamdef(k2,6)=gammaj(1)
      beamdef(k2,7)=xlength-xlk2(k2)

*-----Determine fiber strains and yield indicator for stiffness
*-----for fiber elements

      if(nbeamp(k2,2) .gt. 0)then
      xji=2./xlk2(k2)
      do 200 i=1,3
      do 150 j=1,numfibrs(nbeamp(k2,2))
      dudx=0.5*beamdef(k2,7)*xji
*-----(Note...for the simple approximation employed for dvdx and dwdx
*-----these terms vanish because of the particular convected coordinate
*-----system employed)
      dv2dx2=(xm2pp(xlobat(i),xlk2(k2))*beamdef(k2,3)+
     $xm4pp(xlobat(i),xlk2(k2))*beamdef(k2,6))*(xji**2)
      dw2dx2=(xl2pp(xlobat(i),xlk2(k2))*beamdef(k2,2)+
     $xl4pp(xlobat(i),xlk2(k2))*beamdef(k2,5))*(xji**2)
      y=fibrdefn(nbeamp(k2,2),j,1)
      z=fibrdefn(nbeamp(k2,2),j,2)
      straink2(k2,i,j)=dudx-dv2dx2*y-dw2dx2*z
      delepsk2(k2,i,j)=straink2(k2,i,j)-strmstk2(k2,i,j)
      if(iter .gt. 0 .and. sigyld(nbeamp(k2,1)) .gt. 0.)then
      upyldstr=yldcntrk2(k2,i,j,1)+xlamdab(nbeamp(k2,1))
      btyldstr=yldcntrk2(k2,i,j,1)-xlamdab(nbeamp(k2,1))
      if(straink2(k2,i,j) .gt. btyldstr .and. 
     $straink2(k2,i,j) .lt. upyldstr)then
      nyldstk2(k2,i,j)=0
      else
      nyldstk2(k2,i,j)=1
      endif
      endif
  150 continue
  200 continue
      endif
  250 continue

      return
      end
     
*_______________________________________________________________________________
*-------------------------------------------------------------------------------
*                                 
*                                   FRCBEAMF
*
*_______________________________________________________________________________
*-------------------------------------------------------------------------------

      subroutine frcbeamf(k2,itype,iter)      

*-----This subroutine computes the beam element forces


*-----Variable list...
 
*-----LIST OF VARIABLES AND ARRAYS PASSED INTO ROUTINE AS INPUT...

 
*-----LIST OF VARIABLES PASSED OUT OF ROUTINE AS OUTPUT...


*-----Establishment of global parameters

      integer
     :        max_beams
     :       ,max_beam_fibers_points
     :       ,max_beam_properties
     :       ,max_beams_times_three


      parameter (
     :           max_beams              =    5000
     :          ,max_beam_fibers_points =      20
     :          ,max_beam_properties    =      50
     :          )


      parameter (
     :           max_beams_times_three = max_beams * 3
     :          )

*-----Beam

      integer
     :        nbeamn( max_beams, 3 ), nbeamp( max_beams, 2 ),
     :        numfibrs( max_beam_fibers_points ),
     :        nyldink2( max_beams, 9, max_beam_properties ),
     :        nyldstk2( max_beams, 9, max_beam_properties )

      real
     :      area( max_beam_properties )
     :     ,beamdef( max_beams, 7 )
     :     ,bmdefmst( max_beams, 7 )
     :     ,delepsk2( max_beams, 9, max_beam_properties )
     :     ,e( max_beam_properties )
     :     ,elin( max_beam_properties )
     :     ,etab( max_beam_properties )
     :     ,fibrdefn( max_beam_fibers_points, max_beam_properties, 3 )
     :     ,fmstrk2( max_beams, 12 )
     :     ,forcek2( max_beams, 12 )
     :     ,gj( max_beam_properties )
     :     ,resk2( 12 )
     :     ,sheary( max_beam_properties )
     :     ,shearz( max_beam_properties )
     :     ,sigyld( max_beam_properties )
     :     ,straink2( max_beams, 9, max_beam_properties )
     :     ,stressk2( max_beams, 9, max_beam_properties )
     :     ,strmstk2( max_beams, 9, max_beam_properties )

      real
     :      strsmstk2( max_beams, 9, max_beam_properties )
     :     ,tik2( max_beams_times_three, 3 )
     :     ,tjk2( max_beams_times_three, 3 )
     :     ,tk2( 3, 3 )
     :     ,xj( max_beam_properties )
     :     ,xkk2( 12, 12 )
     :     ,xlamdab( max_beam_properties )
     :     ,xlk2( max_beams )
     :     ,xlnmsk2( max_beam_properties )
     :     ,xlobat( 3 )
     :     ,xlobfres( max_beams, 9, 4 )
     :     ,xlobwt( 3 )
     :     ,xnlnmsk2( max_beam_properties )
     :     ,xnu( max_beam_properties )
     :     ,yi( max_beam_properties )
     :     ,yldcntrk2( max_beams, 9, max_beam_properties, 2 )
     :     ,zi( max_beam_properties )


      common /beam/
     :        nbeamn,   nbeamp,   tk2,      tik2,      tjk2,
     :        xlk2,     xkk2,     area,     e, 
     :        yi,       zi,       xj,       xnu,       sheary,
     :        shearz,   fmstrk2,  gj,       resk2,     forcek2,
     :        xlnmsk2,  beamdef,  xnlnmsk2, elin,      sigyld,
     :        xlamdab,  etab,     straink2, xlobat,    fibrdefn,
     :        numfibrs, xlobfres, nyldink2, nyldstk2,  yldcntrk2,
     :        stressk2, delepsk2, strmstk2, strsmstk2, bmdefmst,
     :        xlobwt


      if(sigyld(nbeamp(k2,1)) .gt. 0.)then
      eyld=elin(nbeamp(k2,1))/etab(nbeamp(k2,1))
      endif
     
      do 50 i=1,3
      xlobfres(k2,i,1)=0.
      xlobfres(k2,i,2)=0.
      xlobfres(k2,i,3)=0.
      xlobfres(k2,i,4)=0.
   50 continue

       
      do 300 i=1,3
      do 275 j=1,numfibrs(nbeamp(k2,2))

*.. If linear stress-strain relationship, or linear analysis, skip integration 
*.. of nonlinear stress-strain relationship and determine stress from incremental
*.. relationship after line 250

      if(iter .eq. 0)go to 260
      if(sigyld(nbeamp(k2,1)) .eq. 0.)go to 250
      if(itype .eq. 0)go to 250
      upyldstr=yldcntrk2(k2,i,j,1)+xlamdab(nbeamp(k2,1))
      btyldstr=yldcntrk2(k2,i,j,1)-xlamdab(nbeamp(k2,1))
      if(nyldink2(k2,i,j) .eq. 0)then
      if(straink2(k2,i,j) .ge. btyldstr .and. straink2(k2,i,j) 
     $.le. upyldstr)then
      delsig=elin(nbeamp(k2,1))*delepsk2(k2,i,j)
      elseif(straink2(k2,i,j) .lt. btyldstr)then
      delepsy=yldcntrk2(k2,i,j,1)-xlamdab(nbeamp(k2,1))
     $-strmstk2(k2,i,j)
      rho=delepsy/delepsk2(k2,i,j)
      delsig=eyld*(1.+rho*(etab(nbeamp(k2,1))-1.))*delepsk2(k2,i,j)
      elseif(straink2(k2,i,j) .gt. upyldstr)then
      delepsy=yldcntrk2(k2,i,j,1)+xlamdab(nbeamp(k2,1))
     $-strmstk2(k2,i,j)
      rho=delepsy/delepsk2(k2,i,j)
      delsig=eyld*(1.+rho*(etab(nbeamp(k2,1))-1.))*delepsk2(k2,i,j)
      endif
      endif
      if(nyldink2(k2,i,j) .eq. 1)then
      if(straink2(k2,i,j) .gt. upyldstr)then
      delsig=eyld*delepsk2(k2,i,j)
      elseif(straink2(k2,i,j) .ge. btyldstr .and.
     $straink2(k2,i,j) .le. upyldstr)then
      delsig=elin(nbeamp(k2,1))*delepsk2(k2,i,j)
      elseif(straink2(k2,i,j) .lt. btyldstr)then
      delepsy=yldcntrk2(k2,i,j,1)-xlamdab(nbeamp(k2,1))
     $-strmstk2(k2,i,j)
      rho=delepsy/delepsk2(k2,i,j)
      delsig=eyld*(1.+rho*(etab(nbeamp(k2,1))-1.))*delepsk2(k2,i,j)
      endif
      endif
      if(nyldink2(k2,i,j) .eq. -1)then
      if(straink2(k2,i,j) .lt. btyldstr)then
      delsig=eyld*delepsk2(k2,i,j)
      elseif(straink2(k2,i,j) .ge. btyldstr .and. straink2(k2,i,j)
     $.le. upyldstr)then
      delsig=elin(nbeamp(k2,1))*delepsk2(k2,i,j)
      elseif(straink2(k2,i,j) .gt. upyldstr)then
      delepsy=yldcntrk2(k2,i,j,1)+xlamdab(nbeamp(k2,1))
     $-strmstk2(k2,i,j)
      rho=delepsy/delepsk2(k2,i,j)
      delsig=eyld*(1.+rho*(etab(nbeamp(k2,1))-1.))*delepsk2(k2,i,j)
      endif
      endif
      go to 260
  250 continue
      delsig=elin(nbeamp(k2,1))*delepsk2(k2,i,j)
  260 continue
      if(iter .gt. 0)then
      stressk2(k2,i,j)=strsmstk2(k2,i,j)+delsig
      else 
      stressk2(k2,i,j)=strsmstk2(k2,i,j)
      endif
      y=fibrdefn(nbeamp(k2,2),j,1)
      z=fibrdefn(nbeamp(k2,2),j,2)
      a=fibrdefn(nbeamp(k2,2),j,3)
      xlobfres(k2,i,1)=xlobfres(k2,i,1)+stressk2(k2,i,j)*a
      xlobfres(k2,i,2)=xlobfres(k2,i,2)-stressk2(k2,i,j)*y*a
      xlobfres(k2,i,3)=xlobfres(k2,i,3)-stressk2(k2,i,j)*z*a
  275 continue
      xlobfres(k2,i,4)=(gj(nbeamp(k2,1))/xlk2(k2))*
     $(beamdef(k2,4)-beamdef(k2,1))
  300 continue
  
      return
      end

*_______________________________________________________________________________
*-------------------------------------------------------------------------------
*                                 
*                                     STFBEAM
*
*_______________________________________________________________________________
*-------------------------------------------------------------------------------

      subroutine stfbeam(k2) 

*-----This subroutine computes the element stiffness matrix for simple linear
*-----beam elements

*-----Variable list...
 
*-----LIST OF VARIABLES AND ARRAYS PASSED INTO ROUTINE AS INPUT...

 
*-----LIST OF VARIABLES PASSED OUT OF ROUTINE AS OUTPUT...


*-----Perform common and array declarations

*-----Establishment of global parameters

      integer
     :        max_beams
     :       ,max_beam_fibers_points
     :       ,max_beam_properties
     :       ,max_beams_times_three
     :       ,max_dof
     :       ,max_nodes


      parameter (
     :           max_beams              =    5000
     :          ,max_beam_fibers_points =      20
     :          ,max_beam_properties    =      50
     :          ,max_dof                =   30000
     :          ,max_nodes              =   10000
     :          )


      parameter (
     :           max_beams_times_three = max_beams * 3
     :          )

*-----Geometry and displacement block

      integer
     :        id( max_nodes, 6 )

      real
     :      accvctr( max_dof ), dispinc( max_dof ), dispvctr( max_dof )
     :     ,velvctr( max_dof ), xyzloc( max_nodes, 3 )


      common /geometry/
     :        id, xyzloc, dispvctr, dispinc, velvctr, accvctr

*-----Beam

      integer
     :        nbeamn( max_beams, 3 ), nbeamp( max_beams, 2 ),
     :        numfibrs( max_beam_fibers_points ),
     :        nyldink2( max_beams, 9, max_beam_properties ),
     :        nyldstk2( max_beams, 9, max_beam_properties )

      real
     :      area( max_beam_properties )
     :     ,beamdef( max_beams, 7 )
     :     ,bmdefmst( max_beams, 7 )
     :     ,delepsk2( max_beams, 9, max_beam_properties )
     :     ,e( max_beam_properties )
     :     ,elin( max_beam_properties )
     :     ,etab( max_beam_properties )
     :     ,fibrdefn( max_beam_fibers_points, max_beam_properties, 3 )
     :     ,fmstrk2( max_beams, 12 )
     :     ,forcek2( max_beams, 12 )
     :     ,gj( max_beam_properties )
     :     ,resk2( 12 )
     :     ,sheary( max_beam_properties )
     :     ,shearz( max_beam_properties )
     :     ,sigyld( max_beam_properties )
     :     ,straink2( max_beams, 9, max_beam_properties )
     :     ,stressk2( max_beams, 9, max_beam_properties )
     :     ,strmstk2( max_beams, 9, max_beam_properties )

      real
     :      strsmstk2( max_beams, 9, max_beam_properties )
     :     ,tik2( max_beams_times_three, 3 )
     :     ,tjk2( max_beams_times_three, 3 )
     :     ,tk2( 3, 3 )
     :     ,xj( max_beam_properties )
     :     ,xkk2( 12, 12 )
     :     ,xlamdab( max_beam_properties )
     :     ,xlk2( max_beams )
     :     ,xlnmsk2( max_beam_properties )
     :     ,xlobat( 3 )
     :     ,xlobfres( max_beams, 9, 4 )
     :     ,xlobwt( 3 )
     :     ,xnlnmsk2( max_beam_properties )
     :     ,xnu( max_beam_properties )
     :     ,yi( max_beam_properties )
     :     ,yldcntrk2( max_beams, 9, max_beam_properties, 2 )
     :     ,zi( max_beam_properties )


      common /beam/
     :        nbeamn,   nbeamp,   tk2,      tik2,      tjk2,
     :        xlk2,     xkk2,     area,     e, 
     :        yi,       zi,       xj,       xnu,       sheary,
     :        shearz,   fmstrk2,  gj,       resk2,     forcek2,
     :        xlnmsk2,  beamdef,  xnlnmsk2, elin,      sigyld,
     :        xlamdab,  etab,     straink2, xlobat,    fibrdefn,
     :        numfibrs, xlobfres, nyldink2, nyldstk2,  yldcntrk2,
     :        stressk2, delepsk2, strmstk2, strsmstk2, bmdefmst,
     :        xlobwt

*-----Linear beam element stiffness coefficients are given by Przemieniecki
*-----in "Theory of Matrix Structural Analysis" (McGraw-Hill 1968) pg. 70-79
 
      n=nbeamp(k2,1)
      s1=(e(n)*area(n))/xlk2(k2)
      s2=(12.*e(n)*zi(n))/((xlk2(k2)**3)*(1.+sheary(n)))
      s3=(12.*e(n)*yi(n))/((xlk2(k2)**3)*(1.+shearz(n)))
      s4=((e(n)/(2.*(1.+xnu(n))))*xj(n))/xlk2(k2)
      s5=(-6.*e(n)*yi(n))/((xlk2(k2)**2)*(1.+shearz(n)))
      s6=(-6.*e(n)*zi(n))/((xlk2(k2)**2)*(1.+sheary(n)))
      s7=((4.+shearz(n))*(e(n)*yi(n)))/(xlk2(k2)*(1.+shearz(n)))
      s8=((4.+sheary(n))*(e(n)*zi(n)))/(xlk2(k2)*(1.+sheary(n)))
      s9=((2.-shearz(n))*(e(n)*yi(n)))/(xlk2(k2)*(1.+shearz(n)))
      s10=((2.-sheary(n))*(e(n)*zi(n)))/(xlk2(k2)*(1.+sheary(n)))

      do 5 i=1,12
      do 4 j=1,12
      xkk2(i,j)=0.
    4 continue
    5 continue

      xkk2(1,1)=s1
      xkk2(7,1)=-s1
      xkk2(2,2)=s2
      xkk2(6,2)=-s6
      xkk2(8,2)=-s2
      xkk2(12,2)=-s6
      xkk2(3,3)=s3
      xkk2(5,3)=s5
      xkk2(9,3)=-s3
      xkk2(11,3)=s5
      xkk2(4,4)=s4
      xkk2(10,4)=-s4
      xkk2(5,5)=s7
      xkk2(9,5)=-s5
      xkk2(11,5)=s9
      xkk2(6,6)=s8
      xkk2(8,6)=s6
      xkk2(12,6)=s10
      xkk2(7,7)=s1
      xkk2(8,8)=s2
      xkk2(12,8)=s6
      xkk2(9,9)=s3
      xkk2(11,9)=-s5
      xkk2(10,10)=s4
      xkk2(11,11)=s7
      xkk2(12,12)=s8

      do 20 i=1,12
      do 10 j=i,12
      xkk2(i,j)=xkk2(j,i)
   10 continue
   20 continue

*-----Add geometric stiffness components

      nodei=nbeamn(k2,1)
      nodej=nbeamn(k2,2)
      xloci=xyzloc(nodei,1)+dispvctr(id(nodei,1))
      yloci=xyzloc(nodei,2)+dispvctr(id(nodei,2))
      zloci=xyzloc(nodei,3)+dispvctr(id(nodei,3))
      xlocj=xyzloc(nodej,1)+dispvctr(id(nodej,1))
      ylocj=xyzloc(nodej,2)+dispvctr(id(nodej,2))
      zlocj=xyzloc(nodej,3)+dispvctr(id(nodej,3))
      delx=xlocj-xloci
      dely=ylocj-yloci
      delz=zlocj-zloci
      xlength=sqrt(delx**2+dely**2+delz**2)
      axlforc=((e(n)*area(n))/xlk2(k2))*(xlength-xlk2(k2))
      geostiff=axlforc/xlk2(k2)
      xkk2(2,2)=xkk2(2,2)+geostiff
      xkk2(8,2)=xkk2(8,2)-geostiff
      xkk2(2,8)=xkk2(2,8)-geostiff
      xkk2(8,8)=xkk2(8,8)+geostiff
      xkk2(3,3)=xkk2(3,3)+geostiff
      xkk2(3,9)=xkk2(3,9)-geostiff
      xkk2(9,3)=xkk2(9,3)-geostiff
      xkk2(9,9)=xkk2(9,9)+geostiff

      return
      end

*_______________________________________________________________________________
*-------------------------------------------------------------------------------
*                                 
*                                   FORCBEAM
*
*_______________________________________________________________________________
*-------------------------------------------------------------------------------

      subroutine forcbeam(k2)

*-----This subroutine computes the element residual contributions for
*-----the simple linear beam elements


*-----Variable list...
 
*-----LIST OF VARIABLES AND ARRAYS PASSED INTO ROUTINE AS INPUT...

 
*-----LIST OF VARIABLES PASSED OUT OF ROUTINE AS OUTPUT...


*-----Perform common and array declarations     

*----Establishment of global parameters

      integer
     :        max_beams
     :       ,max_beam_fibers_points
     :       ,max_beam_properties
     :       ,max_beams_times_three


      parameter (
     :           max_beams              =    5000
     :          ,max_beam_fibers_points =      20
     :          ,max_beam_properties    =      50
     :          )


      parameter (
     :           max_beams_times_three = max_beams * 3
     :          )

*-----Beam

      integer
     :        nbeamn( max_beams, 3 ), nbeamp( max_beams, 2 ),
     :        numfibrs( max_beam_fibers_points ),
     :        nyldink2( max_beams, 9, max_beam_properties ),
     :        nyldstk2( max_beams, 9, max_beam_properties )

      real
     :      area( max_beam_properties )
     :     ,beamdef( max_beams, 7 )
     :     ,bmdefmst( max_beams, 7 )
     :     ,delepsk2( max_beams, 9, max_beam_properties )
     :     ,e( max_beam_properties )
     :     ,elin( max_beam_properties )
     :     ,etab( max_beam_properties )
     :     ,fibrdefn( max_beam_fibers_points, max_beam_properties, 3 )
     :     ,fmstrk2( max_beams, 12 )
     :     ,forcek2( max_beams, 12 )
     :     ,gj( max_beam_properties )
     :     ,resk2( 12 )
     :     ,sheary( max_beam_properties )
     :     ,shearz( max_beam_properties )
     :     ,sigyld( max_beam_properties )
     :     ,straink2( max_beams, 9, max_beam_properties )
     :     ,stressk2( max_beams, 9, max_beam_properties )
     :     ,strmstk2( max_beams, 9, max_beam_properties )

      real
     :      strsmstk2( max_beams, 9, max_beam_properties )
     :     ,tik2( max_beams_times_three, 3 )
     :     ,tjk2( max_beams_times_three, 3 )
     :     ,tk2( 3, 3 )
     :     ,xj( max_beam_properties )
     :     ,xkk2( 12, 12 )
     :     ,xlamdab( max_beam_properties )
     :     ,xlk2( max_beams )
     :     ,xlnmsk2( max_beam_properties )
     :     ,xlobat( 3 )
     :     ,xlobfres( max_beams, 9, 4 )
     :     ,xlobwt( 3 )
     :     ,xnlnmsk2( max_beam_properties )
     :     ,xnu( max_beam_properties )
     :     ,yi( max_beam_properties )
     :     ,yldcntrk2( max_beams, 9, max_beam_properties, 2 )
     :     ,zi( max_beam_properties )


      common /beam/
     :        nbeamn,   nbeamp,   tk2,      tik2,      tjk2,
     :        xlk2,     xkk2,     area,     e, 
     :        yi,       zi,       xj,       xnu,       sheary,
     :        shearz,   fmstrk2,  gj,       resk2,     forcek2,
     :        xlnmsk2,  beamdef,  xnlnmsk2, elin,      sigyld,
     :        xlamdab,  etab,     straink2, xlobat,    fibrdefn,
     :        numfibrs, xlobfres, nyldink2, nyldstk2,  yldcntrk2,
     :        stressk2, delepsk2, strmstk2, strsmstk2, bmdefmst,
     :        xlobwt


      dimension deldisp(12),delfk2(12)

      do 10 i=1,12
      delfk2(i)=0.
   10 continue

*-----Get beam element deformation increments from last converged state

      deldisp(1)=0.
      deldisp(2)=0.
      deldisp(3)=0.
      deldisp(4)=beamdef(k2,1)-bmdefmst(k2,1)
      deldisp(5)=beamdef(k2,2)-bmdefmst(k2,2)
      deldisp(6)=beamdef(k2,3)-bmdefmst(k2,3)
      deldisp(7)=beamdef(k2,7)-bmdefmst(k2,7)
      deldisp(8)=0.
      deldisp(9)=0.
      deldisp(10)=beamdef(k2,4)-bmdefmst(k2,4)
      deldisp(11)=beamdef(k2,5)-bmdefmst(k2,5)
      deldisp(12)=beamdef(k2,6)-bmdefmst(k2,6)

*-----Get beam element force increments

      do 100 i=1,12
      do 90 j=1,12
      delfk2(i)=delfk2(i)+xkk2(i,j)*deldisp(j)
   90 continue
  100 continue

*-----Add force increment to total force

      do 110 i=1,12
      forcek2(k2,i)=fmstrk2(k2,i)+delfk2(i)
  110 continue

      return
      end     

*_______________________________________________________________________________
*-------------------------------------------------------------------------------
*                                 
*                                   STFRESK2
*
*_______________________________________________________________________________
*-------------------------------------------------------------------------------

      subroutine stfresk2(k2)

*-----This subroutine forms the element stiffness and residual contributions
*-----for nonlinear fiber beam elements

*-----Variable list...
 
*-----LIST OF VARIABLES AND ARRAYS PASSED INTO ROUTINE AS INPUT...

 
*-----LIST OF VARIABLES PASSED OUT OF ROUTINE AS OUTPUT...


*-----Perform common and array declarations

*-----Establishment of global parameters

      integer
     :        max_beams
     :       ,max_beam_fibers_points
     :       ,max_beam_properties
     :       ,max_beams_times_three


      parameter (
     :           max_beams              =    5000
     :          ,max_beam_fibers_points =      20
     :          ,max_beam_properties    =      50
     :          )


      parameter (
     :           max_beams_times_three = max_beams * 3
     :          )

*-----Beam

      integer
     :        nbeamn( max_beams, 3 ), nbeamp( max_beams, 2 ),
     :        numfibrs( max_beam_fibers_points ),
     :        nyldink2( max_beams, 9, max_beam_properties ),
     :        nyldstk2( max_beams, 9, max_beam_properties )

      real
     :      area( max_beam_properties )
     :     ,beamdef( max_beams, 7 )
     :     ,bmdefmst( max_beams, 7 )
     :     ,delepsk2( max_beams, 9, max_beam_properties )
     :     ,e( max_beam_properties )
     :     ,elin( max_beam_properties )
     :     ,etab( max_beam_properties )
     :     ,fibrdefn( max_beam_fibers_points, max_beam_properties, 3 )
     :     ,fmstrk2( max_beams, 12 )
     :     ,forcek2( max_beams, 12 )
     :     ,gj( max_beam_properties )
     :     ,resk2( 12 )
     :     ,sheary( max_beam_properties )
     :     ,shearz( max_beam_properties )
     :     ,sigyld( max_beam_properties )
     :     ,straink2( max_beams, 9, max_beam_properties )
     :     ,stressk2( max_beams, 9, max_beam_properties )
     :     ,strmstk2( max_beams, 9, max_beam_properties )

      real
     :      strsmstk2( max_beams, 9, max_beam_properties )
     :     ,tik2( max_beams_times_three, 3 )
     :     ,tjk2( max_beams_times_three, 3 )
     :     ,tk2( 3, 3 )
     :     ,xj( max_beam_properties )
     :     ,xkk2( 12, 12 )
     :     ,xlamdab( max_beam_properties )
     :     ,xlk2( max_beams )
     :     ,xlnmsk2( max_beam_properties )
     :     ,xlobat( 3 )
     :     ,xlobfres( max_beams, 9, 4 )
     :     ,xlobwt( 3 )
     :     ,xnlnmsk2( max_beam_properties )
     :     ,xnu( max_beam_properties )
     :     ,yi( max_beam_properties )
     :     ,yldcntrk2( max_beams, 9, max_beam_properties, 2 )
     :     ,zi( max_beam_properties )


      common /beam/
     :        nbeamn,   nbeamp,   tk2,      tik2,      tjk2,
     :        xlk2,     xkk2,     area,     e, 
     :        yi,       zi,       xj,       xnu,       sheary,
     :        shearz,   fmstrk2,  gj,       resk2,     forcek2,
     :        xlnmsk2,  beamdef,  xnlnmsk2, elin,      sigyld,
     :        xlamdab,  etab,     straink2, xlobat,    fibrdefn,
     :        numfibrs, xlobfres, nyldink2, nyldstk2,  yldcntrk2,
     :        stressk2, delepsk2, strmstk2, strsmstk2, bmdefmst,
     :        xlobwt


      dimension b(4,12), pfreps(4,4), pbsubgrdf(12,12), bt(12,4), 
     $          scratch(12,4), scratch1(12,12), resultant(4)

      do 50 i=1,4
      do 45 j=1,12
      b(i,j)=0.
   45 continue
   50 continue

      do 60 i=1,12
      do 55 j=1,12
      xkk2(i,j)=0.
      pbsubgrdf(i,j)=0.
   55 continue
   60 continue

      do 63 i=1,12
      forcek2(k2,i)=0.
   63 continue

      xjcobn=xlk2(k2)/2.
      xji=2./xlk2(k2)
 
      do 200 i=1,3
     
*-----Form constitutive matrix

      do 70 ii=1,4
      do 65 jj=1,4
      pfreps(ii,jj)=0.
   65 continue
   70 continue

      do 82 j=1,numfibrs(nbeamp(k2,2))
      if(sigyld(nbeamp(k2,1)) .eq. 0.)then
      ee=elin(nbeamp(k2,1))
      elseif(sigyld(nbeamp(k2,1)) .gt. 0.)then
      if(nyldstk2(k2,i,j) .eq. 0)then
      ee=elin(nbeamp(k2,1))
      elseif(nyldstk2(k2,i,j) .eq. 1)then
      xk1=elin(nbeamp(k2,1))/etab(nbeamp(k2,1))
      xk2=0.7*elin(nbeamp(k2,1))
      ee=(xk1+xk2)/2.
      ee=xk1
      endif
      endif
      y=fibrdefn(nbeamp(k2,2),j,1)
      z=fibrdefn(nbeamp(k2,2),j,2)
      a=fibrdefn(nbeamp(k2,2),j,3)
      pfreps(1,1)=pfreps(1,1)+ee*a
      pfreps(1,2)=pfreps(1,2)-ee*y*a
      pfreps(1,3)=pfreps(1,3)-ee*z*a
      pfreps(2,1)=pfreps(2,1)-ee*y*a
      pfreps(2,2)=pfreps(2,2)+ee*(y**2)*a
      pfreps(2,3)=pfreps(2,3)+ee*y*z*a
      pfreps(3,1)=pfreps(3,1)-ee*z*a
      pfreps(3,2)=pfreps(3,2)+ee*y*z*a
      pfreps(3,3)=pfreps(3,3)+ee*(z**2)*a   
   82 continue
      pfreps(4,4)=gj(nbeamp(k2,1))
 
      do 90 ii=1,12
      do 85 jj=1,4
      scratch(ii,jj)=0.
   85 continue
   90 continue

      do 95 ii=1,12
      do 93 jj=1,12
      scratch1(ii,jj)=0.
   93 continue
   95 continue

*-----Form geometric stiffness contribution

      pbsubgrdf(2,2)=0.25*(xji**2)*xlobfres(k2,i,1)
      pbsubgrdf(8,2)=-0.25*(xji**2)*xlobfres(k2,i,1)
      pbsubgrdf(3,3)=0.25*(xji**2)*xlobfres(k2,i,1)
      pbsubgrdf(9,3)=-0.25*(xji**2)*xlobfres(k2,i,1)
      pbsubgrdf(2,8)=-0.25*(xji**2)*xlobfres(k2,i,1)
      pbsubgrdf(8,8)=0.25*(xji**2)*xlobfres(k2,i,1)
      pbsubgrdf(3,9)=-0.25*(xji**2)*xlobfres(k2,i,1)
      pbsubgrdf(9,9)=0.25*(xji**2)*xlobfres(k2,i,1)
    

*-----Form B transpose
*-----Form B matrix

      b(1,1)=-0.5*xji
      b(1,7)=0.5*xji
      b(2,2)=xm1pp(xlobat(i))*(xji**2)
      b(2,6)=xm2pp(xlobat(i),xlk2(k2))*(xji**2)
      b(2,8)=xm3pp(xlobat(i))*(xji**2)
      b(2,12)=xm4pp(xlobat(i),xlk2(k2))*(xji**2)
      b(3,3)=xl1pp(xlobat(i))*(xji**2)
      b(3,5)=xl2pp(xlobat(i),xlk2(k2))*(xji**2)
      b(3,9)=xl3pp(xlobat(i))*(xji**2)
      b(3,11)=xl4pp(xlobat(i),xlk2(k2))*(xji**2)
      b(4,4)=-0.5*xji
      b(4,10)=0.5*xji

      do 115 ii=1,4
      do 110 jj=1,12
      bt(jj,ii)=b(ii,jj)
  110 continue
  115 continue


*-----Form product B transpose x constitutive matrix

      do 130 ii=1,12
      do 125 jj=1,4
      do 120 k=1,4
      scratch(ii,jj)=scratch(ii,jj)+bt(ii,k)*pfreps(k,jj)
  120 continue
  125 continue
  130 continue

*-----Form product B transpose x constitutive matrix x B

      do 145 ii=1,12
      do 140 jj=1,12
      do 135 k=1,4
      scratch1(ii,jj)=scratch1(ii,jj)+scratch(ii,k)*b(k,jj)
  135 continue
  140 continue
  145 continue

*-----Add into the stiffness matrix with Lobatto weights 

      do 155 ii=1,12
      do 150 jj=1,12
      xkk2(ii,jj)=xkk2(ii,jj)+
     $((scratch1(ii,jj)+pbsubgrdf(ii,jj))*xlobwt(i)*xjcobn)
  150 continue
  155 continue

*-----Form element end forces vector

      resultant(1)=xlobfres(k2,i,1)
      resultant(2)=xlobfres(k2,i,2)
      resultant(3)=xlobfres(k2,i,3)
      resultant(4)=xlobfres(k2,i,4)

      do 165 ii=1,12
      do 160 k=1,4
      forcek2(k2,ii)=forcek2(k2,ii)+
     $(bt(ii,k)*resultant(k)*xlobwt(i)*xjcobn)
  160 continue
  165 continue
  
  200 continue
  
      return
      end

*_______________________________________________________________________________
*-------------------------------------------------------------------------------
*                                 
*                                   TRNTBEAM
*
*_______________________________________________________________________________
*-------------------------------------------------------------------------------

      subroutine trntbeam(k2,iter,itype,kload)

*-----This subroutine translates the beam element stiffness and residual
*-----matrices into global coordinates


*-----Variable list...
 
*-----LIST OF VARIABLES AND ARRAYS PASSED INTO ROUTINE AS INPUT...

 
*-----LIST OF VARIABLES PASSED OUT OF ROUTINE AS OUTPUT...


*-----Perform common and array declarations

*-----Establishment of global parameters

      integer
     :        max_beams
     :       ,max_beam_fibers_points
     :       ,max_beam_properties
     :       ,max_beams_times_three


      parameter (
     :           max_beams              =    5000
     :          ,max_beam_fibers_points =      20
     :          ,max_beam_properties    =      50
     :          )


      parameter (
     :           max_beams_times_three = max_beams * 3
     :          )

*-----Beam

      integer
     :        nbeamn( max_beams, 3 ), nbeamp( max_beams, 2 ),
     :        numfibrs( max_beam_fibers_points ),
     :        nyldink2( max_beams, 9, max_beam_properties ),
     :        nyldstk2( max_beams, 9, max_beam_properties )

      real
     :      area( max_beam_properties )
     :     ,beamdef( max_beams, 7 )
     :     ,bmdefmst( max_beams, 7 )
     :     ,delepsk2( max_beams, 9, max_beam_properties )
     :     ,e( max_beam_properties )
     :     ,elin( max_beam_properties )
     :     ,etab( max_beam_properties )
     :     ,fibrdefn( max_beam_fibers_points, max_beam_properties, 3 )
     :     ,fmstrk2( max_beams, 12 )
     :     ,forcek2( max_beams, 12 )
     :     ,gj( max_beam_properties )
     :     ,resk2( 12 )
     :     ,sheary( max_beam_properties )
     :     ,shearz( max_beam_properties )
     :     ,sigyld( max_beam_properties )
     :     ,straink2( max_beams, 9, max_beam_properties )
     :     ,stressk2( max_beams, 9, max_beam_properties )
     :     ,strmstk2( max_beams, 9, max_beam_properties )

      real
     :      strsmstk2( max_beams, 9, max_beam_properties )
     :     ,tik2( max_beams_times_three, 3 )
     :     ,tjk2( max_beams_times_three, 3 )
     :     ,tk2( 3, 3 )
     :     ,xj( max_beam_properties )
     :     ,xkk2( 12, 12 )
     :     ,xlamdab( max_beam_properties )
     :     ,xlk2( max_beams )
     :     ,xlnmsk2( max_beam_properties )
     :     ,xlobat( 3 )
     :     ,xlobfres( max_beams, 9, 4 )
     :     ,xlobwt( 3 )
     :     ,xnlnmsk2( max_beam_properties )
     :     ,xnu( max_beam_properties )
     :     ,yi( max_beam_properties )
     :     ,yldcntrk2( max_beams, 9, max_beam_properties, 2 )
     :     ,zi( max_beam_properties )


      common /beam/
     :        nbeamn,   nbeamp,   tk2,      tik2,      tjk2,
     :        xlk2,     xkk2,     area,     e, 
     :        yi,       zi,       xj,       xnu,       sheary,
     :        shearz,   fmstrk2,  gj,       resk2,     forcek2,
     :        xlnmsk2,  beamdef,  xnlnmsk2, elin,      sigyld,
     :        xlamdab,  etab,     straink2, xlobat,    fibrdefn,
     :        numfibrs, xlobfres, nyldink2, nyldstk2,  yldcntrk2,
     :        stressk2, delepsk2, strmstk2, strsmstk2, bmdefmst,
     :        xlobwt


      dimension tt(12,12), scratch(12,12), t(12,12), scratch1(12,12)

      do 2 i=1,12
      do 1 j=1,12
      scratch(i,j)=0.
      scratch1(i,j)=0.
      resk2(i)=0.
      tt(i,j)=0.
    1 continue
    2 continue

c     UNUSED:  loc=3*k2-2

      do 20 i=1,3
      do 10 j=1,3
      tt(j,i)=tk2(i,j)
      tt(3+j,3+i)=tk2(i,j)
      tt(6+j,6+i)=tk2(i,j)
      tt(9+j,9+i)=tk2(i,j)
   10 continue
   20 continue

      do 25 i=1,12
      do 24 j=1,12
      t(j,i)=tt(i,j)
   24 continue
   25 continue

      do 50 i=1,12
      do 40 j=1,12
      do 30 k=1,12
      scratch(i,j)=scratch(i,j)+tt(i,k)*xkk2(k,j)
   30 continue
   40 continue
   50 continue

      do 53 i=1,12
      do 52 j=1,12
      do 51 k=1,12
      scratch1(i,j)=scratch1(i,j)+scratch(i,k)*t(k,j)
   51 continue
   52 continue
   53 continue

      do 60 i=1,12
      do 55 j=1,12
      xkk2(i,j)=scratch1(i,j)
   55 continue
   60 continue

      if((itype .eq. 3 .and. iter .gt. 0) .or.
     $(itype .eq. 1 .and. iter .gt. 0) .or.
     $(itype .eq. 0 .and. kload .gt. 1) .or. 
     $(itype .eq. 2 .and. iter .gt. 0))then  
      do 70 i=1,12 
      do 65 j=1,12
      resk2(i)=resk2(i)+tt(i,j)*forcek2(k2,j)
   65 continue
   70 continue
      endif


      return
      end
      
*_______________________________________________________________________________
*-------------------------------------------------------------------------------
*                                 
*                                   MASSBEAM
*
*_______________________________________________________________________________
*-------------------------------------------------------------------------------

      subroutine massbeam(k2,beammass,beamrmass,beamtmass)
 
*-----This subroutine computes the beam element mass matrix

*-----Variable list...
 
*-----LIST OF VARIABLES AND ARRAYS PASSED INTO ROUTINE AS INPUT...

 
*-----LIST OF VARIABLES PASSED OUT OF ROUTINE AS OUTPUT...


*-----Perform common and array declarations

*-----Establishment of global parameters

      integer
     :        max_beams
     :       ,max_beam_fibers_points
     :       ,max_beam_properties
     :       ,max_beams_times_three


      parameter (
     :           max_beams              =    5000
     :          ,max_beam_fibers_points =      20
     :          ,max_beam_properties    =      50
     :          )


      parameter (
     :           max_beams_times_three = max_beams * 3
     :          )

*-----Beam

      integer
     :        nbeamn( max_beams, 3 ), nbeamp( max_beams, 2 ),
     :        numfibrs( max_beam_fibers_points ),
     :        nyldink2( max_beams, 9, max_beam_properties ),
     :        nyldstk2( max_beams, 9, max_beam_properties )

      real
     :      area( max_beam_properties )
     :     ,beamdef( max_beams, 7 )
     :     ,bmdefmst( max_beams, 7 )
     :     ,delepsk2( max_beams, 9, max_beam_properties )
     :     ,e( max_beam_properties )
     :     ,elin( max_beam_properties )
     :     ,etab( max_beam_properties )
     :     ,fibrdefn( max_beam_fibers_points, max_beam_properties, 3 )
     :     ,fmstrk2( max_beams, 12 )
     :     ,forcek2( max_beams, 12 )
     :     ,gj( max_beam_properties )
     :     ,resk2( 12 )
     :     ,sheary( max_beam_properties )
     :     ,shearz( max_beam_properties )
     :     ,sigyld( max_beam_properties )
     :     ,straink2( max_beams, 9, max_beam_properties )
     :     ,stressk2( max_beams, 9, max_beam_properties )
     :     ,strmstk2( max_beams, 9, max_beam_properties )

      real
     :      strsmstk2( max_beams, 9, max_beam_properties )
     :     ,tik2( max_beams_times_three, 3 )
     :     ,tjk2( max_beams_times_three, 3 )
     :     ,tk2( 3, 3 )
     :     ,xj( max_beam_properties )
     :     ,xkk2( 12, 12 )
     :     ,xlamdab( max_beam_properties )
     :     ,xlk2( max_beams )
     :     ,xlnmsk2( max_beam_properties )
     :     ,xlobat( 3 )
     :     ,xlobfres( max_beams, 9, 4 )
     :     ,xlobwt( 3 )
     :     ,xnlnmsk2( max_beam_properties )
     :     ,xnu( max_beam_properties )
     :     ,yi( max_beam_properties )
     :     ,yldcntrk2( max_beams, 9, max_beam_properties, 2 )
     :     ,zi( max_beam_properties )


      common /beam/
     :        nbeamn,   nbeamp,   tk2,      tik2,      tjk2,
     :        xlk2,     xkk2,     area,     e, 
     :        yi,       zi,       xj,       xnu,       sheary,
     :        shearz,   fmstrk2,  gj,       resk2,     forcek2,
     :        xlnmsk2,  beamdef,  xnlnmsk2, elin,      sigyld,
     :        xlamdab,  etab,     straink2, xlobat,    fibrdefn,
     :        numfibrs, xlobfres, nyldink2, nyldstk2,  yldcntrk2,
     :        stressk2, delepsk2, strmstk2, strsmstk2, bmdefmst,
     :        xlobwt


      if(nbeamp(k2,2) .eq. 0)then
      beammass=(xlnmsk2(nbeamp(k2,1))*xlk2(k2))/2.
*      beamrmass=(xlnmsk2(nbeamp(k2,1))*(xlk2(k2)**3))/24.
*      xxj=yi(nbeamp(k2,1))+zi(nbeamp(k2,1))
      beamrmass=0.
*      beamtmass=((xlnmsk2(nbeamp(k2,1)))/area(nbeamp(k2,1))
*     $*xxj)*(xlk2(k2))/2.
      beamtmass=0.
      else
      beammass=(xnlnmsk2(nbeamp(k2,1))*xlk2(k2))/2.
      beamrmass=0.
      endif

      return
      end

*_______________________________________________________________________________
*-------------------------------------------------------------------------------
*                                 
*                                   LOADBEAM
*
*_______________________________________________________________________________
*-------------------------------------------------------------------------------

      subroutine loadbeam(kload,k2,itype,beammass,beamrmass,beamtmass,
     $nelembs,nactdof,iter,ngravity)
 
*-----This subroutine loads the beam stiffness, residual and mass matrices into
*-----the global matrices

*-----Variable list...
 
*-----LIST OF VARIABLES AND ARRAYS PASSED INTO ROUTINE AS INPUT...

 
*-----LIST OF VARIABLES PASSED OUT OF ROUTINE AS OUTPUT...


*-----Perform common and array declarations

*-----Establishment of global parameters

      integer
     :        max_a_array
     :       ,max_beams
     :       ,max_beam_fibers_points
     :       ,max_beam_properties
     :       ,max_beams_times_three
     :       ,max_dof
     :       ,max_dof_plus_1
     :       ,max_nodes


      parameter (
     :           max_a_array            = 1500000
     :          ,max_beams              =    5000
     :          ,max_beam_fibers_points =      20
     :          ,max_beam_properties    =      50
     :          ,max_dof                =   30000
     :          ,max_nodes              =   10000
     :          )


      parameter (
     :           max_beams_times_three = max_beams * 3
     :          ,max_dof_plus_1        = max_dof + 1
     :          )

*-----Geometry and displacement block

      integer
     :        id( max_nodes, 6 )

      real
     :      accvctr( max_dof ), dispinc( max_dof ), dispvctr( max_dof )
     :     ,velvctr( max_dof ), xyzloc( max_nodes, 3 )


      common /geometry/
     :        id, xyzloc, dispvctr, dispinc, velvctr, accvctr

*-----Matrix and equation block

      integer
     :        maxa( max_dof_plus_1 ),   maxasave( max_dof_plus_1 ),
     :        ncolht( max_dof_plus_1 ), ncolhtsave( max_dof_plus_1 )

      real
     :     a( max_a_array ),         damp( max_a_array ),
     :     residual( max_dof ),      xkeff( max_a_array ),
     :     xkeffsave( max_a_array ), xmass( max_dof )


      common /matrix/
     :        a,         damp,  residual,  xmass,    maxa,
     :        ncolht,    xkeff, xkeffsave, maxasave,
     :        ncolhtsave

*-----Beam

      integer
     :        nbeamn( max_beams, 3 ), nbeamp( max_beams, 2 ),
     :        numfibrs( max_beam_fibers_points ),
     :        nyldink2( max_beams, 9, max_beam_properties ),
     :        nyldstk2( max_beams, 9, max_beam_properties )

      real
     :      area( max_beam_properties )
     :     ,beamdef( max_beams, 7 )
     :     ,bmdefmst( max_beams, 7 )
     :     ,delepsk2( max_beams, 9, max_beam_properties )
     :     ,e( max_beam_properties )
     :     ,elin( max_beam_properties )
     :     ,etab( max_beam_properties )
     :     ,fibrdefn( max_beam_fibers_points, max_beam_properties, 3 )
     :     ,fmstrk2( max_beams, 12 )
     :     ,forcek2( max_beams, 12 )
     :     ,gj( max_beam_properties )
     :     ,resk2( 12 )
     :     ,sheary( max_beam_properties )
     :     ,shearz( max_beam_properties )
     :     ,sigyld( max_beam_properties )
     :     ,straink2( max_beams, 9, max_beam_properties )
     :     ,stressk2( max_beams, 9, max_beam_properties )
     :     ,strmstk2( max_beams, 9, max_beam_properties )

      real
     :      strsmstk2( max_beams, 9, max_beam_properties )
     :     ,tik2( max_beams_times_three, 3 )
     :     ,tjk2( max_beams_times_three, 3 )
     :     ,tk2( 3, 3 )
     :     ,xj( max_beam_properties )
     :     ,xkk2( 12, 12 )
     :     ,xlamdab( max_beam_properties )
     :     ,xlk2( max_beams )
     :     ,xlnmsk2( max_beam_properties )
     :     ,xlobat( 3 )
     :     ,xlobfres( max_beams, 9, 4 )
     :     ,xlobwt( 3 )
     :     ,xnlnmsk2( max_beam_properties )
     :     ,xnu( max_beam_properties )
     :     ,yi( max_beam_properties )
     :     ,yldcntrk2( max_beams, 9, max_beam_properties, 2 )
     :     ,zi( max_beam_properties )


      common /beam/
     :        nbeamn,   nbeamp,   tk2,      tik2,      tjk2,
     :        xlk2,     xkk2,     area,     e, 
     :        yi,       zi,       xj,       xnu,       sheary,
     :        shearz,   fmstrk2,  gj,       resk2,     forcek2,
     :        xlnmsk2,  beamdef,  xnlnmsk2, elin,      sigyld,
     :        xlamdab,  etab,     straink2, xlobat,    fibrdefn,
     :        numfibrs, xlobfres, nyldink2, nyldstk2,  yldcntrk2,
     :        stressk2, delepsk2, strmstk2, strsmstk2, bmdefmst,
     :        xlobwt


      dimension ncupl(12),bmmass(12)

      ncupl(1)=id(nbeamn(k2,1),1)
      ncupl(2)=id(nbeamn(k2,1),2)
      ncupl(3)=id(nbeamn(k2,1),3)
      ncupl(4)=id(nbeamn(k2,1),4)
      ncupl(5)=id(nbeamn(k2,1),5)
      ncupl(6)=id(nbeamn(k2,1),6)
      ncupl(7)=id(nbeamn(k2,2),1)
      ncupl(8)=id(nbeamn(k2,2),2)
      ncupl(9)=id(nbeamn(k2,2),3)
      ncupl(10)=id(nbeamn(k2,2),4)
      ncupl(11)=id(nbeamn(k2,2),5)
      ncupl(12)=id(nbeamn(k2,2),6)

      bmmass(1)=beammass
      bmmass(2)=beammass
      bmmass(3)=beammass
      bmmass(4)=beamtmass
      bmmass(5)=beamrmass
      bmmass(6)=beamrmass
      bmmass(7)=beammass
      bmmass(8)=beammass
      bmmass(9)=beammass
      bmmass(10)=beamtmass
      bmmass(11)=beamrmass
      bmmass(12)=beamrmass

      do 100 i=1,12
      do 90 j=1,12

*-----If DOF are null skip to next index

      if(ncupl(i) .eq. 0)go to 100
      if(ncupl(j) .eq. 0)go to 90

*-----If stiffness coefficient falls left of diagonal
*-----don't load the term (because of symmetry only using upper triangle)

      if(ncupl(j) .lt. ncupl(i))go to 90

*-----If stiffness coefficient falls below existing profile load it directly

      if(((ncupl(j)-ncupl(i))+1) .gt. ncolht(ncupl(j)))go to 20

*-----If stiffness coefficient falls below existing profile load it directly

      loc=maxa(ncupl(j))+ncupl(j)-ncupl(i)
      a(loc)=a(loc)+xkk2(i,j)
      go to 90

   20 continue

*-----If stiffness coefficient is zero do not put it at top of column
*-----(sspace will puke!)

*-----I don't think its needed (check it some time)
     
      if(xkk2(i,j) .eq. 0.)go to 90

*-----Determine the new column height

      noldclht=ncolht(ncupl(j))
      ncolht(ncupl(j))=(ncupl(j)-ncupl(i))+1
      nshift=ncolht(ncupl(j))-noldclht

*-----Increase the number of elements below skyline 

      nelembs=nelembs+nshift
      nbottom=nelembs-(maxa(ncupl(j))+ncolht(ncupl(j))-1)

*-----Shift the elements below the new profile terms to make
*-----room for the new column height

      do 30 kk=1,nbottom
      loc=nelembs+1-kk
      a(loc)=a(loc-nshift)
   30 continue

*-----Add new skyline term into the correct position in the stiffness
*-----matrix

      a(maxa(ncupl(j))+ncolht(ncupl(j))-1)=xkk2(i,j)

*-----Fill in zeros between new profile term and old profile term

      if(ncolht(ncupl(j))-noldclht .eq. 1)go to 41
      lower=maxa(ncupl(j))+(noldclht-1)+1
      nupper=maxa(ncupl(j))+(ncolht(ncupl(j))-1)-1

      do 40 kk=lower,nupper
      a(kk)=0.
   40 continue

   41 continue

*-----Adjust pointers to stiffness matrix diagonal elements for the "columns"
*-----to the right of column j

      loc=ncupl(j)+1
      do 12 kk=loc,nactdof
      maxa(kk)=maxa(kk)+nshift
   12 continue
   
   90 continue
  100 continue

*-----Load residual vector into global residual vector
      
      if((itype .eq. 3 .and. iter .gt. 0) .or.
     $(itype .eq. 1 .and. iter .gt. 0) .or.
     $(itype .eq. 2 .and. iter .gt. 0))then 
      do 120 i=1,12
      if(ncupl(i) .ne. 0)then
      residual(ncupl(i))=residual(ncupl(i))-resk2(i)
      endif
  120 continue
      endif

*-----Load mass vector into global mass vector

      if(itype .lt. 3)then
      if(kload .le. 1 .and. iter .eq. 0)then
      do 125 i=1,12
      if(ncupl(i) .ne. 0)then
      xmass(ncupl(i))=xmass(ncupl(i))+bmmass(i)
      endif
  125 continue
      endif
      endif

      if(itype .eq. 3 .and. ngravity .eq. 1)then
      if(kload .le. 1 .and. iter .eq. 0)then
      do 126 i=1,12
      if(ncupl(i) .ne. 0)then
      xmass(ncupl(i))=xmass(ncupl(i))+bmmass(i)
      endif
  126 continue
      endif
      endif

      return
      end

*_______________________________________________________________________________
*-------------------------------------------------------------------------------
*                                 
*                                   BEAMUPDT
*
*_______________________________________________________________________________
*-------------------------------------------------------------------------------

      subroutine beamupdt(nbeam)
 
*-----This subroutine updates the nonlinear beam element plasticity parameters

*-----Variable list...
 
*-----LIST OF VARIABLES AND ARRAYS PASSED INTO ROUTINE AS INPUT...

 
*-----LIST OF VARIABLES PASSED OUT OF ROUTINE AS OUTPUT...


*-----Perform common and array declarations

*-----Establishment of global parameters

      integer
     :        max_beams
     :       ,max_beam_fibers_points
     :       ,max_beam_properties
     :       ,max_beams_times_three


      parameter (
     :           max_beams              =    5000
     :          ,max_beam_fibers_points =      20
     :          ,max_beam_properties    =      50
     :          )


      parameter (
     :           max_beams_times_three = max_beams * 3
     :          )

*-----Beam

      integer
     :        nbeamn( max_beams, 3 ), nbeamp( max_beams, 2 ),
     :        numfibrs( max_beam_fibers_points ),
     :        nyldink2( max_beams, 9, max_beam_properties ),
     :        nyldstk2( max_beams, 9, max_beam_properties )

      real
     :      area( max_beam_properties )
     :     ,beamdef( max_beams, 7 )
     :     ,bmdefmst( max_beams, 7 )
     :     ,delepsk2( max_beams, 9, max_beam_properties )
     :     ,e( max_beam_properties )
     :     ,elin( max_beam_properties )
     :     ,etab( max_beam_properties )
     :     ,fibrdefn( max_beam_fibers_points, max_beam_properties, 3 )
     :     ,fmstrk2( max_beams, 12 )
     :     ,forcek2( max_beams, 12 )
     :     ,gj( max_beam_properties )
     :     ,resk2( 12 )
     :     ,sheary( max_beam_properties )
     :     ,shearz( max_beam_properties )
     :     ,sigyld( max_beam_properties )
     :     ,straink2( max_beams, 9, max_beam_properties )
     :     ,stressk2( max_beams, 9, max_beam_properties )
     :     ,strmstk2( max_beams, 9, max_beam_properties )

      real
     :      strsmstk2( max_beams, 9, max_beam_properties )
     :     ,tik2( max_beams_times_three, 3 )
     :     ,tjk2( max_beams_times_three, 3 )
     :     ,tk2( 3, 3 )
     :     ,xj( max_beam_properties )
     :     ,xkk2( 12, 12 )
     :     ,xlamdab( max_beam_properties )
     :     ,xlk2( max_beams )
     :     ,xlnmsk2( max_beam_properties )
     :     ,xlobat( 3 )
     :     ,xlobfres( max_beams, 9, 4 )
     :     ,xlobwt( 3 )
     :     ,xnlnmsk2( max_beam_properties )
     :     ,xnu( max_beam_properties )
     :     ,yi( max_beam_properties )
     :     ,yldcntrk2( max_beams, 9, max_beam_properties, 2 )
     :     ,zi( max_beam_properties )


      common /beam/
     :        nbeamn,   nbeamp,   tk2,      tik2,      tjk2,
     :        xlk2,     xkk2,     area,     e, 
     :        yi,       zi,       xj,       xnu,       sheary,
     :        shearz,   fmstrk2,  gj,       resk2,     forcek2,
     :        xlnmsk2,  beamdef,  xnlnmsk2, elin,      sigyld,
     :        xlamdab,  etab,     straink2, xlobat,    fibrdefn,
     :        numfibrs, xlobfres, nyldink2, nyldstk2,  yldcntrk2,
     :        stressk2, delepsk2, strmstk2, strsmstk2, bmdefmst,
     :        xlobwt


      do 300 k2=1,nbeam
      if(nbeamp(k2,2) .eq. 0)then
      do 50 i=1,12
      fmstrk2(k2,i)=forcek2(k2,i)
   50 continue
      do 60 i=1,7
      bmdefmst(k2,i)=beamdef(k2,i)
   60 continue
      else
      do 250 i=1,3
      do 200 j=1,numfibrs(nbeamp(k2,1))
      if(sigyld(nbeamp(k2,1)) .gt. 0.)then
      upyldstr=yldcntrk2(k2,i,j,1)+xlamdab(nbeamp(k2,1))
      btyldstr=yldcntrk2(k2,i,j,1)-xlamdab(nbeamp(k2,1))
      if(straink2(k2,i,j) .ge. upyldstr)then
      nyldink2(k2,i,j)=1
      delstress=stressk2(k2,i,j)-(yldcntrk2(k2,i,j,2)
     $+sigyld(nbeamp(k2,1)))
      yldcntrk2(k2,i,j,2)=yldcntrk2(k2,i,j,2)+delstress
      yldcntrk2(k2,i,j,1)=straink2(k2,i,j)-xlamdab(nbeamp(k2,1))
      endif
      if(straink2(k2,i,j) .le. btyldstr)then
      nyldink2(k2,i,j)=-1
      delstress=stressk2(k2,i,j)-(yldcntrk2(k2,i,j,2)
     $-sigyld(nbeamp(k2,1)))
      yldcntrk2(k2,i,j,2)=yldcntrk2(k2,i,j,2)+delstress
      yldcntrk2(k2,i,j,1)=straink2(k2,i,j)+xlamdab(nbeamp(k2,1))
      endif
      if(straink2(k2,i,j) .gt. btyldstr .and. 
     $straink2(k2,i,j) .lt. upyldstr)nyldink2(k2,i,j)=0
      endif
      strmstk2(k2,i,j)=straink2(k2,i,j)
      strsmstk2(k2,i,j)=stressk2(k2,i,j)
  200 continue
  250 continue
      endif
  300 continue
   
      return
      end

*_______________________________________________________________________________
*-------------------------------------------------------------------------------
*                                 
*                                  DAMPMATRIX
*
*_______________________________________________________________________________
*-------------------------------------------------------------------------------

      subroutine dampmatrix(nelembs,numnp,dampstiff,dampmass)
 
*-----This subroutine forms the system damping matrix 


*-----Variable list...
 
*-----LIST OF VARIABLES AND ARRAYS PASSED INTO ROUTINE AS INPUT...

 
*-----LIST OF VARIABLES PASSED OUT OF ROUTINE AS OUTPUT...


*-----Perform common and array declarations
 
*-----Establishment of global parameters

      integer
     :        max_a_array
     :       ,max_dof
     :       ,max_dof_plus_1
     :       ,max_nodes


      parameter (
     :           max_a_array            = 1500000
     :          ,max_dof                =   30000
     :          ,max_nodes              =   10000
     :          )


      parameter (
     :           max_dof_plus_1        = max_dof + 1
     :          )

*-----Geometry and displacement block

      integer
     :        id( max_nodes, 6 )

      real
     :      accvctr( max_dof ), dispinc( max_dof ), dispvctr( max_dof )
     :     ,velvctr( max_dof ), xyzloc( max_nodes, 3 )


      common /geometry/
     :        id, xyzloc, dispvctr, dispinc, velvctr, accvctr

*-----Matrix and equation block

      integer
     :        maxa( max_dof_plus_1 ),   maxasave( max_dof_plus_1 ),
     :        ncolht( max_dof_plus_1 ), ncolhtsave( max_dof_plus_1 )

      real
     :     a( max_a_array ),         damp( max_a_array ),
     :     residual( max_dof ),      xkeff( max_a_array ),
     :     xkeffsave( max_a_array ), xmass( max_dof )


      common /matrix/
     :        a,         damp,  residual,  xmass,    maxa,
     :        ncolht,    xkeff, xkeffsave, maxasave,
     :        ncolhtsave


      do 10 i=1,nelembs
      damp(i)=dampstiff*a(i)
   10 continue

      do 30 i=1,numnp
      do 20 j=1,6
      if(id(i,j) .ne. 0)then
      loc=maxa(id(i,j))
      damp(loc)=damp(loc)+dampmass*xmass(id(i,j))
      endif
   20 continue
   30 continue

      return
      end

*_______________________________________________________________________________
*-------------------------------------------------------------------------------
*                                 
*                                POINTLOADS
*
*_______________________________________________________________________________
*-------------------------------------------------------------------------------

      subroutine pointloads(npntload,ngravity,iter,kload)
 
*-----This subroutine computes the beam element mass matrix


*-----Variable list...
 
*-----LIST OF VARIABLES AND ARRAYS PASSED INTO ROUTINE AS INPUT...

 
*-----LIST OF VARIABLES PASSED OUT OF ROUTINE AS OUTPUT...


*-----Perform common and array declarations

*-----Establishment of global parameters

      integer
     :        max_a_array
     :       ,max_beam_properties
     :       ,max_dof
     :       ,max_dof_plus_1
     :       ,max_nodes


      parameter (
     :           max_a_array            = 1500000
     :          ,max_beam_properties    =      50
     :          ,max_dof                =   30000
     :          ,max_nodes              =   10000
     :          )


      parameter (
     :           max_dof_plus_1        = max_dof + 1
     :          )

*-----Geometry and displacement block

      integer
     :        id( max_nodes, 6 )

      real
     :      accvctr( max_dof ), dispinc( max_dof ), dispvctr( max_dof )
     :     ,velvctr( max_dof ), xyzloc( max_nodes, 3 )


      common /geometry/
     :        id, xyzloc, dispvctr, dispinc, velvctr, accvctr

*-----Matrix and equation block

      integer
     :        maxa( max_dof_plus_1 ),   maxasave( max_dof_plus_1 ),
     :        ncolht( max_dof_plus_1 ), ncolhtsave( max_dof_plus_1 )

      real
     :     a( max_a_array ),         damp( max_a_array ),
     :     residual( max_dof ),      xkeff( max_a_array ),
     :     xkeffsave( max_a_array ), xmass( max_dof )


      common /matrix/
     :        a,         damp,  residual,  xmass,    maxa,
     :        ncolht,    xkeff, xkeffsave, maxasave,
     :        ncolhtsave

*-----Load vectors

      integer
     :        loads( 10, 3 ), ncurvloc( max_beam_properties )

      real
     :     ar( 50000 ),  blo( 80 ), br( 50000 ),  bup( 80 ),
     :     bupc( 80 ),   dd( 80 ),  eigval( 80 ), eigvec( max_dof, 80 ),
     :     rtolv( 80 ),  ttt( max_dof ),          vec( 80, 80 ),    
     :     w( max_dof ), xload( 200000, 2 ), xnstif( max_a_array )


      common /loads/
     :        loads, xload, eigvec, eigval,   ttt,
     :        w,     ar,    br,     vec,      dd,    rtolv,
     :        bup,   blo,   bupc,   ncurvloc, xnstif


      do 10 i=1,npntload
      loc=id(loads(i,1),loads(i,2))
      if(ngravity .eq. 1)then
      loadloc=(ncurvloc(loads(i,3))+(kload-2))
      else
      loadloc=(ncurvloc(loads(i,3))+(kload-1))
      endif
      if(iter .eq. 0 .and. kload .gt. 1)then
      residual(loc)=residual(loc)+xload(loadloc,2)
     $-xload(loadloc-1,2)
      else
      residual(loc)=residual(loc)+xload(loadloc,2)
      endif
   10 continue

      return
      end

*_______________________________________________________________________________
*-------------------------------------------------------------------------------
*                                 
*                                GROUNDMOTION
*
*_______________________________________________________________________________
*-------------------------------------------------------------------------------

      subroutine groundmotion(numnp,kload,iter)
 
*-----This subroutine modifies the residual vector to include ground accelerations


*-----Variable list...
 
*-----LIST OF VARIABLES AND ARRAYS PASSED INTO ROUTINE AS INPUT...

 
*-----LIST OF VARIABLES PASSED OUT OF ROUTINE AS OUTPUT...


*-----Perform common and array declarations

*-----Establishment of global parameters

      integer
     :        max_a_array
     :       ,max_beam_properties
     :       ,max_dof
     :       ,max_dof_plus_1
     :       ,max_nodes


      parameter (
     :           max_a_array            = 1500000
     :          ,max_beam_properties    =      50
     :          ,max_dof                =   30000
     :          ,max_nodes              =   10000
     :          )


      parameter (
     :           max_dof_plus_1        = max_dof + 1
     :          )

*-----Geometry and displacement block

      integer
     :        id( max_nodes, 6 )

      real
     :      accvctr( max_dof ), dispinc( max_dof ), dispvctr( max_dof )
     :     ,velvctr( max_dof ), xyzloc( max_nodes, 3 )


      common /geometry/
     :        id, xyzloc, dispvctr, dispinc, velvctr, accvctr

*-----Matrix and equation block

      integer
     :        maxa( max_dof_plus_1 ),   maxasave( max_dof_plus_1 ),
     :        ncolht( max_dof_plus_1 ), ncolhtsave( max_dof_plus_1 )

      real
     :     a( max_a_array ),         damp( max_a_array ),
     :     residual( max_dof ),      xkeff( max_a_array ),
     :     xkeffsave( max_a_array ), xmass( max_dof )


      common /matrix/
     :        a,         damp,  residual,  xmass,    maxa,
     :        ncolht,    xkeff, xkeffsave, maxasave,
     :        ncolhtsave

*-----Load vectors

      integer
     :        loads( 10, 3 ), ncurvloc( max_beam_properties )

      real
     :     ar( 50000 ),  blo( 80 ), br( 50000 ),  bup( 80 ),
     :     bupc( 80 ),   dd( 80 ),  eigval( 80 ), eigvec( max_dof, 80 ),
     :     rtolv( 80 ),  ttt( max_dof ),          vec( 80, 80 ),    
     :     w( max_dof ), xload( 200000, 2 ), xnstif( max_a_array )


      common /loads/
     :        loads, xload, eigvec, eigval,   ttt,
     :        w,     ar,    br,     vec,      dd,    rtolv,
     :        bup,   blo,   bupc,   ncurvloc, xnstif


      do 10 i=1,numnp
      if(id(i,1) .ne. 0)then
      loc=id(i,1)
      if(iter .eq. 0 .and. kload .gt. 1)then
      residual(loc)=residual(loc)-xmass(loc)*(xload(kload,2)
     $-xload(kload-1,2))
      else
      residual(loc)=residual(loc)-xmass(loc)*xload(kload,2)
      endif
      endif
   10 continue

      return
      end

*-------------------------------------------------------------------------------
*                                 
*                                MATRIXMULT
*
*_______________________________________________________________________________
*-------------------------------------------------------------------------------

      subroutine matrixmult(nactdof,coefmatrix,ncolht,maxa,vector,prod)
 
*-----This subroutine computes the product of a symmetric banded matrix 
*-----matrix stored in skyline form with a vector


*-----Variable list...
 
*-----LIST OF VARIABLES AND ARRAYS PASSED INTO ROUTINE AS INPUT...

*-----nactdof=number of active degrees of freedom

*-----coefmatrix=left hand side matrix (skyline storage)

*-----ncolht=vector containing column heights of matrix

*-----maxa=vector of locations of diagonal terms of matrix

*-----rhs=right hand side vector

 
*-----LIST OF VARIABLES PASSED OUT OF ROUTINE AS OUTPUT...


*-----prod=vector containing product of coefficient matrix and rhs vector
 

*-----Perform common and array declarations

*-----Establishment of global parameters

      integer
     :        max_a_array
     :       ,max_beam_properties
     :       ,max_dof
     :       ,max_dof_plus_1
     :       ,max_nodes


      parameter (
     :           max_a_array            = 1500000
     :          ,max_beam_properties    =      50
     :          ,max_dof                =   30000
     :          ,max_nodes              =   10000
     :          )


      parameter (
     :           max_dof_plus_1        = max_dof + 1
     :          )
 
      real coefmatrix( max_a_array ), vector( max_dof ), 
     :      prod( max_dof )

      integer ncolht( max_dof_plus_1 ), maxa( max_dof_plus_1 )

      do 10 i=1,max_dof,1
      prod(i)=0.
   10 continue


*-----Perform row x vector calculation for each row
         
      do 50 i=1,nactdof  
      
*-----Take products for row elements left of (and including) diagonal

      do 20 j=1,ncolht(i)
      prod(i)=prod(i)+coefmatrix(maxa(i)+ncolht(i)-j)
     $*vector(i-ncolht(i)+j)
   20 continue

*-----Take products for row elements right of diagonal

      nstart=i+1

      do 30 j=nstart,nactdof
      nrow=j-i+1
      if(ncolht(j) .ge. nrow)then
      prod(i)=prod(i)+coefmatrix(maxa(j)+(j-i))*vector(j)
      endif
   30 continue

   50 continue

      return
      end

*-------------------------------------------------------------------------------
*                                 
*                                TIMEHISTORY
*
*_______________________________________________________________________________
*-------------------------------------------------------------------------------

      subroutine timehistory(gamma,beta,stepsize,iter,numnp,nelembs,
     $                       nactdof,kload,itype)
 
*-----This subroutine computes the effective stiffness matrix and effective 
*-----residual vector for the direct integration time history solution


*-----Variable list...
 
*-----LIST OF VARIABLES AND ARRAYS PASSED INTO ROUTINE AS INPUT...

*-----gamma=Newmark integration constant

*-----beta=Newmark integration constant

*-----stepsize=integration time step size

*-----iter=equilibrium iteration counter

 
*-----LIST OF VARIABLES PASSED OUT OF ROUTINE AS OUTPUT...


*-----Perform common and array declarations

*-----Establishment of global parameters

      integer
     :        max_a_array
     :       ,max_beams
     :       ,max_dof
     :       ,max_dof_plus_1
     :       ,max_nodes


      parameter (
     :           max_a_array            = 1500000
     :          ,max_beams              =    5000
     :          ,max_dof                =   30000
     :          ,max_nodes              =   10000
     :          )


      parameter (
     :           max_dof_plus_1        = max_dof + 1
     :          )

*-----Geometry and displacement block

      integer
     :        id( max_nodes, 6 )

      real
     :      accvctr( max_dof ), dispinc( max_dof ), dispvctr( max_dof )
     :     ,velvctr( max_dof ), xyzloc( max_nodes, 3 )


      common /geometry/
     :        id, xyzloc, dispvctr, dispinc, velvctr, accvctr

*-----Matrix and equation block

      integer
     :        maxa( max_dof_plus_1 ),   maxasave( max_dof_plus_1 ),
     :        ncolht( max_dof_plus_1 ), ncolhtsave( max_dof_plus_1 )

      real
     :     a( max_a_array ),         damp( max_a_array ),
     :     residual( max_dof ),      xkeff( max_a_array ),
     :     xkeffsave( max_a_array ), xmass( max_dof )


      common /matrix/
     :        a,         damp,  residual,  xmass,    maxa,
     :        ncolht,    xkeff, xkeffsave, maxasave,
     :        ncolhtsave


      real xmatrix1( max_a_array ), xmatrix2( max_a_array ), 
     :     residual1( max_dof ), residual2( max_dof )


*-----Form required constants
         
      c1=1./((stepsize**2)*beta)
      c2=gamma/(stepsize*beta)
      c3=1./(stepsize*beta)
      c4=gamma/beta
      c5=1./(2.*beta)
      c6=(stepsize*gamma)/(2.*beta)-stepsize 

      do 5 i=1,max_a_array,1
      xmatrix1(i)=0.
      xmatrix2(i)=0.
    5 continue
      do 6 i=1,max_dof,1
      residual1(i)=0.
      residual2(i)=0.
    6 continue

      if((itype .eq. 0 .and. kload .eq. 1) .or.
     $(itype .eq. 1))then 

*-----Form effective stiffness matrix 

      do 10 i=1,nelembs
      xkeff(i)=c2*damp(i)+a(i)
   10 continue
  
      do 20 i=1,numnp     
      do 15 j=1,6
      if(id(i,j) .ne. 0)then
      loc=maxa(id(i,j))
      xkeff(loc)=xkeff(loc)+c1*xmass(id(i,j))
      endif
   15 continue
   20 continue

      endif

*-----Form effective residual vector

      if(iter .eq. 0)then

      do 35 i=1,numnp
      do 30 j=1,6
      if(id(i,j) .ne. 0)then
      loc=maxa(id(i,j))
      xmatrix1(loc)=xmatrix1(loc)+c3*xmass(id(i,j))
      xmatrix2(loc)=xmatrix2(loc)+c5*xmass(id(i,j))
      endif
   30 continue
   35 continue

      do 40 i=1,nelembs
      xmatrix1(i)=xmatrix1(i)+c4*damp(i)
      xmatrix2(i)=xmatrix2(i)+c6*damp(i)
   40 continue

      call matrixmult(nactdof,xmatrix1,ncolht,maxa,velvctr,residual1)
      call matrixmult(nactdof,xmatrix2,ncolht,maxa,accvctr,residual2)

      do 45 i=1,nactdof
      residual(i)=residual(i)+residual1(i)+residual2(i)
   45 continue
 
      else

      do 50 i=1,nactdof
      residual1(i)=xmass(i)*accvctr(i)
   50 continue

      call matrixmult(nactdof,damp,ncolht,maxa,velvctr,residual2)

      do 55 i=1,nactdof
      residual(i)=residual(i)-residual1(i)-residual2(i)
   55 continue

      endif
      
      return
      end

*_______________________________________________________________________________
*-------------------------------------------------------------------------------
*                                 
*                                   COLSOL
*
*_______________________________________________________________________________
*-------------------------------------------------------------------------------

      subroutine colsol(a,v,maxa,nn,nwk,kkk,iout)
 
*-----This subroutine solves the finite element equilibrium equations
*-----in-core using compacted storage and column reduction scheme
*-----(From "Numerical Methods in Finite Element Analysis"
*-----by Bathe and Wilson, Prentice-Hall 1976


*-----Variable list...
 
*-----LIST OF VARIABLES AND ARRAYS PASSED INTO ROUTINE AS INPUT...

*-----a(nwk)=stiffness matrix stored in compacted form

*-----v(nn)=right-hand-side load vector
*-----maxa(nnm)=vector containing addresses of diagonal
*               elements of stiffness matrix in A
*-----nn=number of equations
*-----nwk=number of elements below skyline of matrix
*-----nnm=nn+1
*-----kkk=input flag
*         EQ. 1 Triangularization of stiffness matrix
*         EQ. 2 Reduction and back-substitution of load
*               vector
*-----iout=number of output device
 
*-----LIST OF VARIABLES PASSED OUT OF ROUTINE AS OUTPUT...
      
*-----a(nwk)=D and L factors of stiffness matrix
*-----v(nn)=displacement vector

*-----Perform common and array declarations

      dimension a(nwk),v(*),maxa(*)

c     perform L*D*L(T) factorization of stiffness matrix

      if(kkk-2)40,150,150
   40 do 140 n=1,nn
      kn=maxa(n)
      kl=kn+1
      ku=maxa(n+1)-1
      kh=ku-kl
      if(kh)110,90,50
   50 k=n-kh
      ic=0
      klt=ku
      do 80 j=1,kh
      ic=ic+1
      klt=klt-1
      ki=maxa(k)
      nd=maxa(k+1)-ki-1
      if(nd)80,80,60
   60 kk=min0(ic,nd)
      c=0.
      do 70 l=1,kk
   70 c=c+a(ki+l)*a(klt+l)
      a(klt)=a(klt)-c
   80 k=k+1
   90 k=n
      b=0.
      do 100 kk=kl,ku
      k=k-1
      ki=maxa(k)
      c=a(kk)/a(ki)
      b=b+c*a(kk)
  100 a(kk)=c
      a(kn)=a(kn)-b
  110 if(a(kn))120,120,140
  120 write(iout,2000)n,a(kn)
      stop
  140 continue
      return
 
c     reduce right-hand-side load vector

  150 do 180 n=1,nn
      kl=maxa(n)+1
      ku=maxa(n+1)-1
      if(ku-kl)180,160,160
  160 k=n
      c=0.
      do 170 kk=kl,ku
      k=k-1
  170 c=c+a(kk)*v(k)
      v(n)=v(n)-c
  180 continue
 
c     back-substitute

      do 200 n=1,nn
      k=maxa(n)
  200 v(n)=v(n)/a(k)
      if(nn.eq.1)return
      n=nn
      do 230 l=2,nn
      kl=maxa(n)+1
      ku=maxa(n+1)-1
      if(ku-kl)230,210,210
  210 k=n
      do 220 kk=kl,ku
      k=k-1
  220 v(k)=v(k)-a(kk)*v(n)
  230 n=n-1
      return
 2000 format(//,'Stop - stiffness matrix not positive definite',//,
     $'Nonpositive pivot for equation',i4,//,
     $'Pivot=',e20.12)
      end
 

c-------------------------------------------------------------------------------

      subroutine sspace(a,b,maxa,r,eigv,tt,w,ar,br,vec,d,rtolv,bup,blo,
     $bupc,nn,nnm,nwk,nwm,nroot,rtol,nc,nnc,nitem,ifss,ifpr,xnstif,
     $iout,nrowdim,ncoldim)

c-------------------------------------------------------------------------------
c
c Program
c This program solves for the smallest eigenvalues and the corresponding
c eigenvectors in the generalized eigenproblem using the subspace
c iteration method
c
c From "Numerical Methods in Finite Element Analysis"
c by Bathe and Wilson, Prentice-Hall 1976
c
c-----Input Variables...
c
c     a(nwk)      = stiffness matrix in compacted form (assumed
c                   positive definite)
c     b(nwm)      = mass matrix in compacted form
c     maxa(nnm)   = vector containing  addresses of diagonal
c                   elements of stiffness matrix a
c     r(n1,n2)    = eigenvectors on solution exit
c      where n1=nrodim, n2=ncoldim (so same as main dimensioning)
c     eigv(nc)    = eigenvalues on solution exit
c     tt(nn)      = working vector
c     w(nn)       = working vector
c     ar(nnc)     = working matrix storing projection of k
c     br(nnc)     = working matrix storing projection of m
c     vec(nc,nc)  = working matrix
c     d(nc)       = working vector
c     rtolv(nc)   = working vector
c     bup(nc)     = working vector
c     blo(nc)     = working vector
c     bupc(nc)    = working vector
c     nn          = order of stiffness and mass matrices
c     nnm         = nn+1
c     nwk         = number of elements below skyline of
c                   stiffness matrix
c     nwm         = number of elements below skyline of
c                   mass matrix
c                   i.e. nwm=nwk for consistent mass matrix
c                        nwm=nn  for lumped mass matrix
c     nroot       = number of required eigenvalues and eigenvectors
c     rtol        = convergence tolerance on eigenvalues
c                   (1.e-6 or smaller)
c     nc          = number of iteration vectors used
c                   (usually set to min(2*nroot, nroot+8), but nc
c                   cannot be larger than the number of mass
c                   degrees of freedom)
c     nnc         = nc*(nc+1)/2 dimension of storage vectors ar,br
c     nitem       = maximum number of subspace iterations permitted
c                   (usually set to 16)
c                   the parameters nc and/or nitem must be
c                   increased if a solution has not converged
c     ifss        = flag for Sturm sequence check
c                   eq.0 no check
c                   eq.1 check
c     ifpr        = flag for printing during iteration
c                   eq.0 no printing
c                   eq.1 print
c     nstif       = scratch file to store stiffness matrix
c     iout        = outout printing file
c
c-----Output...
c
c     eigv(nroot) = eigenvalues
c     r(nn,nroot) = eigenvectors
c
c-------------------------------------------------------------------------------
 
      dimension a(nwk),b(nwm),r(nrowdim,ncoldim),tt(nn),w(nn),eigv(nc),
     $d(nc),vec(nc,nc),ar(nnc),br(nnc),rtolv(nc),bup(nc),
     $blo(nc),bupc(nc),xnstif(nwk)
      integer maxa(nnm)
 
c
c     Set tolerance for Jacobi iteration
c
 
      tolj=1.e-12
 
c
c     initialization
c
 
      iconv=0
      nsch=0
      nsmax=12
      n1=nc+1
      nc1=nc-1
      do 9999 iiii=1,nwk
      xnstif(iiii)=a(iiii)
 9999 continue
      do 60 i=1,nc
   60 d(i)=0.
 
c
c      Establish starting iteration vectors
c
 
      nd=nn/nc
      if(nwm.gt.nn)go to 4
      j=0
      do 2 i=1,nn
      ii=maxa(i)
      r(i,1)=b(i)
      if(b(i).gt.0)j=j+1
    2 w(i)=b(i)/a(ii)
      if(nc.le.j)go to 16
      write(iout,1007)
      stop
    4 do 10 i=1,nn
      ii=maxa(i)
      r(i,1)=b(ii)
   10 w(i)=b(ii)/a(ii)
   16 do 20 j=2,nc
      do 20 i=1,nn
   20 r(i,j)=0.
c
      l=nn-nd
      do 30 j=2,nc
      rt=0.
      do 40 i=1,l
      if (w(i).lt.rt)go to 40
      rt=w(i)
      ij=i
   40 continue
      do 50 i=l,nn
      if (w(i).le.rt)go to 50
      rt=w(i)
      ij=i
   50 continue
      tt(j)=float(ij)
      w(ij)=0.
      l=l-nd
   30 r(ij,j)=1.
c
      write(iout,1008)
      write(iout,1002)(tt(j),j=2,nc)
c
c      Factorize matrix A into (L)*(D)*(L(T))
c
      ish=0
      call decomp(a,maxa,nn,ish,iout)
c
c-----START OF ITERATION LOOP
c
      nite=0
  100 nite=nite+1
      if(ifpr.eq.0)go to 90
      write(iout,1010)nite
c
c     Calculate the projections of A and B
c
   90 ij=0
      do 110 j=1,nc
      do 120 k=1,nn
  120 tt(k)=r(k,j)
      call redbak(a,tt,maxa,nn)
      do 130 i=j,nc
      art=0.
      do 140 k=1,nn
  140 art=art+r(k,i)*tt(k)
      ij=ij+1
  130 ar(ij)=art
      do 150 k=1,nn
  150 r(k,j)=tt(k)
  110 continue
      ij=0
      do 160 j=1,nc
      call mult(tt,b,r(1,j),maxa,nn,nwm)
      do 180 i=j,nc
      brt=0.
      do 190 k=1,nn
  190 brt=brt+r(k,i)*tt(k)
      ij=ij+1
  180 br(ij)=brt
      if(iconv.gt.0)go to 160
      do 200 k=1,nn
  200 r(k,j)=tt(k)
  160 continue
c
c      Solve for eigensystem of subspace operators
c
      if(ifpr.eq.0)go to 320
      ind=1
  210 write(iout,1020)
      ii=1
      do 300 i=1,nc
      itemp=ii+nc-i
      write(iout,1005)(ar(j),j=ii,itemp)
  300 ii=ii+n1-i
      write(iout,1030)
      ii=1
      do 310 i=1,nc
      itemp=ii+nc-i
      write(iout,1005)(br(j),j=ii,itemp)
  310 ii=ii+n1-i
      if(ind.eq.2)go to 350
c
  320 call jacobi(ar,br,vec,eigv,w,nc,nnc,tolj,nsmax,ifpr,iout)
c
c     e1=eigv(1)
c     e2=eigv(2)
c     e3=eigv(3)

      if(ifpr.eq.0)go to 350
      write(iout,1040)
      ind=2
      go to 210
c
c     Arrange eigenvalues in ascending order
c
  350 is=0
      ii=1
      do 360 i=1,nc1
      itemp=ii+n1-i
      if(eigv(i+1).ge.eigv(i))go to 360
      is=is+1
      eigvt=eigv(i+1)
      eigv(i+1)=eigv(i)
      eigv(i)=eigvt
      bt=br(itemp)
      br(itemp)=br(ii)
      br(ii)=bt
      do 370 k=1,nc
      rt=vec(k,i+1)
      vec(k,i+1)=vec(k,i)
  370 vec(k,i)=rt
  360 ii=itemp
      if(is.gt.0)go to 350
      if(ifpr.eq.0)go to 375
      write(iout,1035)
      write(iout,1006)(eigv(i),i=1,nc)
c
c     Calculate B times approximate eigenvectors (iconv.eq.0)
c     or final eigenvector approximations (iconv.gt.0)
c
  375 do 420 i=1,nn
      do 422 j=1,nc
  422 tt(j)=r(i,j)
      do 424 k=1,nc
      rt=0.
      do 430 l=1,nc
  430 rt=rt+tt(l)*vec(l,k)
  424 r(i,k)=rt
  420 continue
      if(iconv.gt.0)go to 500
c
c     Check for convergence of eigenvalues
c
      do 380 i=1,nc
      dif=abs(eigv(i)-d(i))
  380 rtolv(i)=dif/eigv(i)
      if(ifpr.eq.0)go to 385
      write(iout,1050)
      write(iout,1005)(rtolv(i),i=1,nc)
c
  385 do 390 i=1,nroot
      if(rtolv(i).gt.rtol)go to 400
  390 continue
      write(iout,1060)rtol
      iconv=1
      go to 100
  400 if(nite.lt.nitem)go to 410
      write(iout,1070)
      iconv=2
      ifss=0
      go to 100
c
  410 do 440 i=1,nc
  440 d(i)=eigv(i)
      go to 100
c
c---END OF ITERATION LOOP
c
  500 write(iout,1100)
      write(iout,1006)(eigv(i),i=1,nroot)
      write(iout,1110)
 
c-----****** comment out writing of eigenvectors ******
 
      do 530 j=1,nroot
  530 write(iout,1005)(r(k,j),k=1,nn)
c
c      Calculate and print error norms
c
      do 8888 iiii=1,nwk
      a(iiii)=xnstif(iiii)
 8888 continue
c
      do 580 l=1,nroot
      rt=eigv(l)
      call mult(tt,a,r(1,l),maxa,nn,nwk)
      vnorm=0.
c
      do 590 i=1,nn
  590 vnorm=vnorm+tt(i)*tt(i)
      call mult(w,b,r(1,l),maxa,nn,nwm)
      wnorm=0.
      do 600 i=1,nn
      tt(i)=tt(i)-rt*w(i)
  600 wnorm=wnorm+tt(i)*tt(i)
      vnorm=sqrt(vnorm)
      wnorm=sqrt(wnorm)
      d(l)=wnorm/vnorm
  580 continue
      write(iout,1115)
      write(iout,1006)(d(i),i=1,nroot)
c
c     Apply Sturm sequence check
c
      if(ifss.eq.0)go to 700
      call scheck(eigv,rtolv,bup,blo,bupc,d,nc,nei,rtol,shift,iout)
c
      write(iout,1120)shift
c
c     Shift matrix A
c
      do 7777 iiii=1,nwk
      a(iiii)=xnstif(iiii)
 7777 continue
      if(nwm.gt.nn)go to 645
      do 640 i=1,nn
      ii=maxa(i)
  640 a(ii)=a(ii)-b(i)*shift
      go to 660
  645 do 650 i=1,nwk
  650 a(i)=a(i)-b(i)*shift
c
c     Factorize shifted matrix
c
  660 ish=1
      call decomp(a,maxa,nn,ish,iout)
c
c     count number of negative diagonal elements
c
      nsch=0
      do 664 i=1,nn
      ii=maxa(i)
      if(a(ii).lt.0)nsch=nsch+1
  664 continue
      if(nsch.eq.nei)go to 670
      nmis=nsch-nei
      write(iout,1130)nmis
      go to 700
  670 write(iout,1140)nsch
  700 return
c
 1002 format(/,5f10.0)
 1005 format(/,5e11.4)
 1006 format(/,3e22.14)
 1007 format(///,'Stop, NC is larger than the number of mass ',
     $'degrees of freedom')
 1008 format(///,'degrees of freedom excited by unit starting',
     $' iteration vectors')
 1010 format(/,'iteration number',i4)
 1020 format(/,'projection of A (/,matrix AR)')
 1030 format(/,'projection of B (/,matrix BR)')
 1035 format(/,'eigenvalues of AR-lambda*BR')
 1040 format(/,'AR and BR after Jacobi diagonalization')
 1050 format(/,'relative tolerance reached on eigenvalues')
 1060 format(///,'convergence reached for rtol',1x,e11.4)
 1070 format(1h1,'***No convergence in maximum number of iterations',
     $' permitted',/,'we accept the current iteration values',/,
     $'the Sturm sequence check is not performed')
 1100 format(///,'the calculated eigenvalues are')
 1115 format(//,'print error norms on the eigenvalues')
 1110 format(//,'the calculated eigenvectors are',//)
 1120 format(///,'check applied at shift ',e22.14)
 1130 format(//,'there are ',i4,'eigenvalues mising ')
 1140 format(//,'we found the lowest ',i4,' eigenvalues')
c
      end

*2 4 6 8(1)2 4 6 8(2)2 4 6 8(3)2 4 6 8(4)2 4 6 8(5)2 4 6 8(6)2 4 6 8(7)2 4 6 8(8
*-------------------------------------------------------------------------------

      subroutine decomp(a,maxa,nn,ish,iout)

*-------------------------------------------------------------------------------
c
c  Program
c      To calculate (L)*(D)*(L)(T) factorization of
c      stiffness matrix
c
*-------------------------------------------------------------------------------
      dimension a(*),maxa(*)
      if(nn.eq.1)return
      do 200 n=1,nn
      kn=maxa(n)
      kl=kn+1
      ku=maxa(n+1)-1
      kh=ku-kl
      if(kh)304,240,210
  210 k=n-kh
      ic=0
      klt=ku
      do 260 j=1,kh
      ic=ic+1
      klt=klt-1
      ki=maxa(k)
      nd=maxa(k+1)-ki-1
      if(nd)260,260,270
  270 kk=min0(ic,nd)
      c=0.
      do 280 l=1,kk
  280 c=c+a(ki+l)*a(klt+l)
      a(klt)=a(klt)-c
  260 k=k+1
  240 k=n
      b=0.
      do 300 kk=kl,ku
      k=k-1
      ki=maxa(k)
      c=a(kk)/a(ki)
      if(abs(c).lt.1.e7)go to 290
      write(iout,2010)n,c
      stop
  290 b=b+c*a(kk)
  300 a(kk)=c
      a(kn)=a(kn)-b
  304 if(a(kn))310,310,200
  310 if(ish.eq.0)go to 320
      if(a(kn).eq.0)a(kn)=-1.e-16
      go to 200
  320 write(iout,2000)n,a(kn)
      stop
  200 continue
c
      return
 2000 format(//,'Stop-stiffness matrix not positive definite',//,
     $'nonpositive pivot for equation',i4,
     $'pivot = ',e20.12)
 2010 format(//,'Stop Sturm sequence check failed because of '
     $'multiplier growth for column',i4,//,' multiplier =',e20.8)
      end

 
*2 4 6 8(1)2 4 6 8(2)2 4 6 8(3)2 4 6 8(4)2 4 6 8(5)2 4 6 8(6)2 4 6 8(7)2 4 6 8(8
*-------------------------------------------------------------------------------

      subroutine redbak(a,v,maxa,nn)

*-------------------------------------------------------------------------------
c
c   Program
c     To reduce and back-substitute iteration vectors
c
*-------------------------------------------------------------------------------
 
       dimension a(*),v(*),maxa(*)
c
      do 400 n=1,nn
      kl=maxa(n)+1
      ku=maxa(n+1)-1
      if(ku-kl)400,410,410
  410 k=n
      c=0.
      do 420 kk=kl,ku
      k=k-1
  420 c=c+a(kk)*v(k)
      v(n)=v(n)-c
  400 continue
c
      do 480 n=1,nn
      k=maxa(n)
  480 v(n)=v(n)/a(k)
      if(nn.eq.1)return
      n=nn
      do 500 l=2,nn
      kl=maxa(n)+1
      ku=maxa(n+1)-1
      if(ku-kl)500,510,510
  510 k=n
      do 520 kk=kl,ku
      k=k-1
  520 v(k)=v(k)-a(kk)*v(n)
  500 n=n-1
c
      return
      end
 
*2 4 6 8(1)2 4 6 8(2)2 4 6 8(3)2 4 6 8(4)2 4 6 8(5)2 4 6 8(6)2 4 6 8(7)2 4 6 8(8
*-------------------------------------------------------------------------------

      subroutine mult(tt,b,rr,maxa,nn,nwm)

*-------------------------------------------------------------------------------
c
c   Program
c     To evaluate product of B times RR and store result in TT
c
*-------------------------------------------------------------------------------
 
      dimension tt(*),b(*),rr(*),maxa(*)
c
      if(nwm.gt.nn)go to 20
      do 10 i=1,nn
   10 tt(i)=b(i)*rr(i)
      return
c
   20 do 40 i=1,nn
   40 tt(i)=0.
      do 100 i=1,nn
      kl=maxa(i)
      ku=maxa(i+1)-1
      ii=i+1
      cc=rr(i)
      do 100 kk=kl,ku
      ii=ii-1
  100 tt(ii)=tt(ii)+b(kk)*cc
      if(nn.eq.1)return
      do 200 i=2,nn
      kl=maxa(i)+1
      ku=maxa(i+1)-1
      if(ku-kl)200,210,210
  210 ii=i
      aa=0.
      do 220 kk=kl,ku
      ii=ii-1
  220 aa=aa+b(kk)*rr(ii)
      tt(i)=tt(i)+aa
  200 continue
c
      return
      end
 
*2 4 6 8(1)2 4 6 8(2)2 4 6 8(3)2 4 6 8(4)2 4 6 8(5)2 4 6 8(6)2 4 6 8(7)2 4 6 8(8
*-------------------------------------------------------------------------------

      subroutine scheck(eigv,rtolv,bup,blo,bupc,neiv,nc,nei,rtol,shift,
     $iout)

*-------------------------------------------------------------------------------
c
c   Program
c      To evaluate shift for Sturm sequence check
c
*-------------------------------------------------------------------------------
      dimension eigv(nc),rtolv(nc),bup(nc),blo(nc),bupc(nc),neiv(nc)
c
      ftol=0.01
c
      do 100 i=1,nc
      bup(i)=eigv(i)*(1.+ftol)
  100 blo(i)=eigv(i)*(1.-ftol)
      nroot=0
      do 120 i=1,nc
  120 if(rtolv(i).lt.rtol)nroot=nroot+1
      if(nroot.ge.1)go to 200
      write(iout,1010)
      stop
c
c     Find upper bounds on eigenvalue clusters
c
  200 do 240 i=1,nroot
  240 neiv(i)=1
      if(nroot.ne.1)go to 260
      bupc(1)=bup(1)
      lm=1
      l=1
      i=2
      go to 295
  260 l=1
      i=2
  270 if(bup(i-1).le.blo(i))go to 280
      neiv(l)=neiv(l)+1
      i=i+1
      if(i.le.nroot)go to 270
  280 bupc(l)=bup(i-1)
      if(i.gt.nroot)go to 290
      l=l+1
      i=i+1
      if(i.le.nroot)go to 270
      bupc(l)=bup(i-1)
  290 lm=l
      if(nroot.eq.nc)go to 300
  295 if(bup(i-1).le.blo(i))go to 300
      if(rtolv(i).gt.rtol)go to 300
      bupc(l)=bup(i)
      neiv(l)=neiv(l)+1
      nroot=nroot+1
      if(nroot.eq.nc)go to 300
      i=i+1
      go to 295
c
c     Find shift
c
  300 write(iout,1020)
      write(iout,1005)(bupc(i),i=1,lm)
      write(iout,1030)
      write(iout,1006)(neiv(i),i=1,lm)
      ll=lm-1
      if(lm.eq.1)go to 310
  330 do 320 i=1,ll
  320 neiv(l)=neiv(l)+neiv(i)
      l=l-1
      ll=ll-1
      if(l.ne.1)go to 330
  310 write(iout,1040)
      write(iout,1006)(neiv(i),i=1,lm)
      l=0
      do 340 i=1,lm
      l=l+1
      if(neiv(i).ge.nroot)go to 350
  340 continue
  350 shift=bupc(l)
      nei=neiv(l)
c
      return
c
 1005 format(/,6e22.14)
 1006 format(/,6i22)
 1010 format(/,'****Error solution stop in SCHECK',/,10x,
     $'no eigenvalues found',/)
 1020 format(///,'upper bounds on eigenvalue clusters')
 1030 format(/,'No of eigenvalues in each cluster')
 1040 format(/,'No of eigenvalues less than upper bounds')
      end
 
*2 4 6 8(1)2 4 6 8(2)2 4 6 8(3)2 4 6 8(4)2 4 6 8(5)2 4 6 8(6)2 4 6 8(7)2 4 6 8(8
*-------------------------------------------------------------------------------

      subroutine jacobi(a,b,x,eigv,d,n,nwa,rtol,nsmax,ifpr,iout)

*-------------------------------------------------------------------------------
c
c  Program
c     To solve the generalized eigenproblem using the generalized
c     Jacobi iteration
c
*-------------------------------------------------------------------------------
 
      dimension a(nwa),b(nwa),x(n,n),eigv(n),d(n)
c
c     Initialize eigenvalue and eigenvector matrices
c
      n1=n+1
      ii=1
      do 10 i=1,n
      if(a(ii).gt.0. .and. b(ii).gt.0.)go to 4
      write(iout,2020)ii,a(ii),b(ii)
      stop
    4 d(i)=a(ii)/b(ii)
      eigv(i)=d(i)
   10 ii=ii+n1-i
      do 30 i=1,n
      do 20 j=1,n
   20 x(i,j)=0.
   30 x(i,i)=1.
      if(n.eq.1)return
c
c     Initialize sweep counter and begin iteration
c
      nsweep=0
      nr=n-1
   40 nsweep=nsweep+1
      if(ifpr.eq.1)write(iout,2000)nsweep
c
c     Check if present off-diagonal element is large enough to
c     require zeroing
c
      eps=(0.01**nsweep)**2
      do 210 j=1,nr
      jp1=j+1
      jm1=j-1
      ljk=jm1*n-jm1*j/2
      jj=ljk+j
      do 210 k=jp1,n
      kp1=k+1
      km1=k-1
      jk=ljk+k
      kk=km1*n-km1*k/2+k
      eptola=(a(jk)*a(jk))/(a(jj)*a(kk))
      eptolb=(b(jk)*b(jk))/(b(jj)*b(kk))
      if((eptola.lt.eps).and.(eptolb.lt.eps))go to 210
c
c     If zeroing is required, calculate the rotation matrix elements
c     CA and CG
c
      akk=a(kk)*b(jk)-b(kk)*a(jk)
      ajj=a(jj)*b(jk)-b(jj)*a(jk)
      ab=a(jj)*b(kk)-a(kk)*b(jj)
      check=(ab*ab+4.*akk*ajj)/4.
      if(check)50,60,60
   50 write(iout,2020)
      stop
   60 sqch=sqrt(check)
      d1=ab/2.+sqch
      d2=ab/2.-sqch
      den=d1
      if(abs(d2).gt.abs(d1))den=d2
      if(den)80,70,80
   70 ca=0.
      cg=-a(jk)/a(kk)
      go to 90
   80 ca=akk/den
      cg=-ajj/den
c
c     Perform the generalized rotation to zero the present
c     off-diagonal element
c
   90 if(n-2)100,190,100
  100 if(jm1-1)130,110,110
  110 do 120 i=1,jm1
      im1=i-1
      ij=im1*n-im1*i/2+j
      ik=im1*n-im1*i/2+k
      aj=a(ij)
      bj=b(ij)
      ak=a(ik)
      bk=b(ik)
      a(ij)=aj+cg*ak
      b(ij)=bj+cg*bk
      a(ik)=ak+ca*aj
  120 b(ik)=bk+ca*bj
  130 if(kp1-n)140,140,160
  140 lji=jm1*n-jm1*j/2
      lki=km1*n-km1*k/2
      do 150 i=kp1,n
      ji=lji+i
      ki=lki+i
      aj=a(ji)
      bj=b(ji)
      ak=a(ki)
      bk=b(ki)
      a(ji)=aj+cg*ak
      b(ji)=bj+cg*bk
      a(ki)=ak+ca*aj
  150 b(ki)=bk+ca*bj
  160 if(jp1-km1)170,170,190
  170 lji=jm1*n-jm1*j/2
      do 180 i=jp1,km1
      ji=lji+i
      im1=i-1
      ik=im1*n-im1*i/2+k
      aj=a(ji)
      bj=b(ji)
      ak=a(ik)
      bk=b(ik)
      a(ji)=aj+cg*ak
      b(ji)=bj+cg*bk
      a(ik)=ak+ca*aj
  180 b(ik)=bk+ca*bj
  190 ak=a(kk)
      bk=b(kk)
      a(kk)=ak+2.*ca*a(jk)+ca*ca*a(jj)
      b(kk)=bk+2.*ca*b(jk)+ca*ca*b(jj)
      a(jj)=a(jj)+2.*cg*a(jk)+cg*cg*ak
      b(jj)=b(jj)+2.*cg*b(jk)+cg*cg*bk
      a(jk)=0.
      b(jk)=0.
c
c     Update the eigenvector matrix after each rotation
c
      do 200 i=1,n
      xj=x(i,j)
      xk=x(i,k)
      x(i,j)=xj+cg*xk
  200 x(i,k)=xk+ca*xj
  210 continue
c
c     Update the eigenvalues after each sweep
c
      ii=1
      do 220 i=1,n
c     test1=a(ii)
c     test2=b(ii)
      if(a(ii).gt.0 .and. b(ii).gt.0.)go to 215
      write(iout,2020)ii,a(ii),b(ii)
      stop
  215 eigv(i)=a(ii)/b(ii)
      ii=ii+n1-i
  220 continue
      if(ifpr.eq.0)go to 230
      write(iout,2030)
      write(iout,2010)(eigv(i),i=1,n)
c
c     Check for convergence
c
  230 do 240 i=1,n
      tol=rtol*d(i)
      dif=abs(eigv(i)-d(i))
      if(dif.gt.tol)go to 280
  240 continue
c
c     Check all off-diagonal elements to see if another
c     sweep is required
c
      eps=rtol**2
      do 250 j=1,nr
      jm1=j-1
      jp1=j+1
      ljk=jm1*n-jm1*j/2
      jj=ljk+j
      do 250 k=jp1,n
      km1=k-1
      jk=ljk+k
      kk=km1*n-km1*k/2+k
      epsa=(a(jk)*a(jk))/(a(jj)*a(kk))
      epsb=(b(jk)*b(jk))/(b(jj)*b(kk))
      if((epsa.lt.eps).and.(epsb.lt.eps))go to 250
      go to 280
  250 continue
c
c      Fill out bottom triangle of resultant matrices and
c      scale eigenvectors
c
  255 ii=1
      do 275 i=1,n
      bb=sqrt(b(ii))
      do 270 k=1,n
  270 x(k,i)=x(k,i)/bb
  275 ii=ii+n1-i
      return
c
c     Update D matrix and start new sweep, if allowed
c
  280 do 290 i=1,n
  290 d(i)=eigv(i)
      if(nsweep.lt.nsmax)go to 40
      go to 255
 2000 format(/,'sweep number in JACOBI =',i4)
 2010 format(/,6e20.12)
 2020 format(/,'***Error solution stop***',/,
     $'matrix not positive definite',/,
     $'ii=',i4,'a(ii)=',e20.12,'b(ii)=',e20.12)
 2030 format(/,'current eigenvalues in JACOBI are',/)
      end


      subroutine rwabsf (fit,lfn,us,bfs,fcs)                                  
c
c     ansi fortran equivalents to llnl familied file random i/o subrs.
c
c     steven j. sackett
c
c     lawrence livermore laboratory
c     livermore, california 94550
c
c     january 3, 1984
c
c-----------------------------------------------------------------------
c     this package consists of seven subrs. for handling word
c     addressable random i/o on familied files:
c
c               asgrfm  (not called directly)
c               nrfnam  (not called directly)
c               rdabsf  (fit,w,nw,da)
c               rdiska  (not called directly)
c                wdiska
c               riosta  (fit)
c               rwabsf  (fit,lfn,us,bfs,fcs)
c               wrabsf  (fit,w,nw,da)
c
c
c     if the disk address given in the argument list to rdabsf/wrabsf
c     is greater than or equal to the family (file) size, it is biased
c     to access the correct family member. here the family size is
c     defined to be the size at which family members are created.
c
c     the root name for a family, which is the name of the first family
c     member, is taken to be the name associated with the given unit
c     specifier (logical unit). names for succeeding family members are
c     generated by appending a two digit integer to the root name (or
c     to its first six characters). the disk address bias for a member,
c     which is also the first word address for the member, is equal to
c     this integer times the family (file) size. assuming, for example,
c     the root name 'diska' and a family size of 1000000b words, a
c     family with five members would appear as follows:
c             --------------------------------------
c             member     name     first word address
c             --------------------------------------
c               1        diska         0
c               2        diska01       1000000b
c               3        diska02       2000000b
c               4        diska03       3000000b
c               5        diska04       4000000b
c             --------------------------------------
c     note that with this naming scheme errors can occur if an eight
c     character root name of the form axxxxxnn is used, where 'a'
c     denotes an alphabetic character, 'x' denotes any character, and
c     'n' denotes a numeric character.
c
c     if a family member exists, it is opened and used. if it does
c     not exist, a new file is created at the family (file) size.
c     note that this will result in an error exit if a write is
c     attempted on an existing file that is read-only.
c
c-----------------------------------------------------------------------
c
c.... open/close a file for word addressable random i/o
c
c     calling sequence: call rwabsf(fit,lfn,us,bfs,fcs)
c
c     input arguments
c            fit      array to use for the file information table
c            lfn      the file name for an open call; the file
c                     disposition status ('keep' or 'delete') for
c                     a close call
c            us       the file unit specifier (logical unit no.)
c                     for an open call; zero (0) or omitted for a
c                     close call
c            bfs      the buffer size in words for an open call;
c                     zero (0) or omitted for a close call
c            fcs      file creation size (family size) for an open
c                     call; zero (0) or omitted for a close call
c
c     fit must be dimensioned as an array of at least bfs+7 words in
c     the user's program and must not be changed while the file is open
c     bfs must be a multiple of 512 (1000b). fcs must be a multiple of
c     bfs.
c
      implicit integer (a-z)                                                  
      dimension fit(8)                                                        
      character*8 lfn                                                         
      common/frfcm1/mxfrf,ifrf,buflen,fcsize,diskloc,curlen,kop,ier           
c            mxfrf     dimension of the familied random file name
c                      table (currently set to 16) - this is the
c                      maximum number of familied random files
c                      allowed to be open at the same time
c            ifrf      index in familied random file name table for
c                      the file accessed last
c
      character*8 frfn,frn,kfn                                                
      common/frfcm2/frfn(2,16),frn,kfn                                        
c            frfn      familied random file name table
c
      logical fxist                                                           
      logical
     :        opened_status

c
c     parameter giving number of record units per integer word
c     for most systems a single character is used as a record unit

c     parameter (ncpw=4)                                                sun
c     parameter (ncpw=1)                                                sgi
      parameter (ncpw=1)                                                dec
c


c
c     Initialization transferred to "blkdat"
c
c     data mxfrf/16/                                                          
c     data frfn/32*'        '/                                                


c
c     fit(1) = us
c     fit(2) = index in random file root name table for this file
c     fit(3) = bfs
c     fit(4) = fcs
c     fit(5) = disk address of first word in the buffer
c     fit(6) = number of words of data currently in the buffer
c     fit(7) = disk address of last word in file + 1
c


      if ( 0 .ge. bfs ) then

          write (*,1111) bfs
1111      format( 'rwabsf:  bfs size:  ', i10, ' no action performed.' )

          return

      endif


      if ((lfn.eq.'keep').or.(lfn.eq.'delete')) go to 50                      
c.... initialize fit and put name in familied random file name table
      do 10 n=1,mxfrf                                                         
      if (frfn(1,n).ne.'        ') go to 10                                   
      frfn(1,n)=lfn                                                           
      fit(2)=n                                                                
      go to 12                                                                
   10 continue                                                                
      stop        ' rwabsf open error- too many random files '                
   12 continue                                                                
      buflen=bfs-mod(bfs,512)                                                 

c
c     ORIGINAL:
c
c     inquire (file=lfn,recl=rcl)                                             
c     rcl=rcl/ncpw                                                            
c
      inquire( file = lfn, opened = opened_status, recl = rcl )

      if ( opened_status .eqv. .false. ) then

          rcl = 0

      endif

      rcl = rcl / ncpw

c     END REPLACEMENT


      if (buflen.lt.rcl) then                                                 
      stop        ' rwabsf open error- buffer too small '                     
      endif                                                                   
      if (rcl.ne.0) buflen=rcl                                                
      fit(1)=us                                                               
      fit(3)=buflen                                                           
      fit(4)=fcs-mod(fcs,buflen)                                              
      fit(5)=fit(4)                                                           
      fit(6)=0                                                                
      fit(7)=0                                                                
      return                                                                  
c.... flush the buffer if data is present which is not on disk
   50 continue                                                                
      if (fit(3).lt.0) then                                                   
      fit(3)=-fit(3)                                                          
      call wdiska (fit(1),fit(8),fit(3),fit(5))                               
      fit(7)=max0(fit(7),fit(5)+fit(3))                                       
      endif                                                                   
c.... close the file
c     commented out 7/19/90 Tom Spelce
c     if (fit(7).eq.0) return
      close (fit(1),status='keep')                                            
      ifrf=fit(2)                                                             
      frn=frfn(1,ifrf)                                                        
      frfn(1,ifrf)='        '                                                 
      frfn(2,ifrf)='        '                                                 
      if (lfn.ne.'delete') return                                             
c.... destroy the family if requested
      n=0                                                                     
      kfn=frn                                                                 
   60 continue                                                                
      ifrf=fit(2)                                                             
      inquire (file=kfn,exist=fxist)                                          
      if (.not.fxist) return                                                  
      buflen=fit(3)                                                           

c
c     ORIGINAL:
c
c     open (fit(1),file=kfn,access='direct',form='unformatted',               
c    1recl=buflen,status='old')                                               

      open( fit(1), file = kfn, access = 'direct', form = 'formatted',
     :      recl = buflen * ncpw, status = 'old' )

c
c     END REPLACEMENT

      close (fit(1),status='delete')                                          
      n=n+1                                                                   
      call nrfnam (frn,n,kfn)                                                 
      go to 60                                                                
      end                                                                     

*2 4 6 8(1)2 4 6 8(2)2 4 6 8(3)2 4 6 8(4)2 4 6 8(5)2 4 6 8(6)2 4 6 8(7)2 4 6 8(8
*------------------------------- SUBROUTINE nrfnam -----------------------------

      subroutine nrfnam(frn,i,nfn)                                            
c***********************************************************************c
c                                                                       c
c     get name of i+1st member of random access file family -frn-       c
c                                                                       c
c-----------------------------------------------------------------------c
c                                                                       c
c     input arguments:                                                  c
c           frn    family root name (name of the first family member)   c
c            i     family member index for member i+1                   c
c                                                                       c
c     output arguments:                                                 c
c           nfn    name of family member i+1                            c
c                                                                       c
c-----------------------------------------------------------------------c
c                                                                       c
c      Note: this version of 8/22/88 was created for Robert Whirley     c
c            to handle naming problems in DYNA3D                        c
c                                                                       c
c***********************************************************************c
      character*8 frn,nfn                                                     
      character*1 ni(10)                                                      
      data ni/'0','1','2','3','4','5','6','7','8','9'/                        
c
      if (i.ne.0) go to 11                                                    
c.....i = 0
      nfn=frn                                                                 
      return                                                                  
   11 if (i.lt.100) go to 21                                                  
c.....i >= 100
   20 format(' family index exceeds 99 ')                                     
      write(* ,20)                                                            
      call exit                                                               
c.....0 < i < 100
   21 do 30 k=1,6                                                             
         if (frn(k:k).eq.' ') go to 40                                        
   30    continue                                                             
   40 k = k-1                                                                 
      j=i/10                                                                  
      nfn=frn(1:k)//ni(j+1)//ni(i-10*j+1)                                     
      return                                                                  
      end 
     
*2 4 6 8(1)2 4 6 8(2)2 4 6 8(3)2 4 6 8(4)2 4 6 8(5)2 4 6 8(6)2 4 6 8(7)2 4 6 8(8
*------------------------------- SUBROUTINE rdiska -----------------------------

      subroutine rdiska (lus,w,nw,da)                                         
c
c     interface to direct access i/o for rwabsf random i/o subrs.
c
c.... entry to transfer one record from disk to a buffer
c
c     input arguments
c            lus      the file unit specifier (logical unit no.)
c            nw       number of words to read from disk
c                     (must be a multiple of 512)
c            da       zero base disk word.address
c                     (must be on a sector boundary)
c
c     output arguments
c            w        data read from disk
c
      implicit integer(a-z)                                                   
      dimension w(nw)                                                         
c
      common/frfcm1/mxfrf,ifrf,buflen,fcsize,diskloc,curlen,kop,ier           
c
      lda=da/buflen+1                                                         
      read (lus,rec=lda,iostat=ios) w                                         
      return                                                                  
c
c.... entry to transfer one record from a buffer to disk
c
c     input arguments
c            lus      the file unit specifier (logical unit no.)
c            w        data to be written to disk
c            nw       number of words to write to disk
c                     (must be a multiple of 512)
c            da       zero base disk word.address
c                     (must be on a sector boundary)
c
      entry wdiska(lus,w,nw,da)                                               
      lda=da/buflen+1                                                         
      write (lus,rec=lda) w                                                   
      return                                                                  
c
      end                                                                     

*2 4 6 8(1)2 4 6 8(2)2 4 6 8(3)2 4 6 8(4)2 4 6 8(5)2 4 6 8(6)2 4 6 8(7)2 4 6 8(8
*------------------------------- SUBROUTINE wrabsf -----------------------------

      subroutine wrabsf (fit,w,nw,da)                                         
c
c.... entry for random write
c
c     calling sequence: call wrabsf(fit,w,nw,da)
c
c     input arguments
c            fit      the file information table
c            w        data to be written to disk
c            nw       number of words to write to disk
c            da       zero base disk address
c
      implicit integer (a-z)                                                  
      dimension fit(8),w(nw)                                                  
c
      common/frfcm1/mxfrf,ifrf,buflen,fcsize,diskloc,curlen,kop,ier           
      character*8 frfn,frn,kfn                                                
      common/frfcm2/frfn(2,16),frn,kfn                                        
c
c.... get family size, family root name, and name of open family member
      ifrf=fit(2)                                                             
      fcsize=fit(4)                                                           
      frn=frfn(1,ifrf)                                                        
      kfn=frfn(2,ifrf)                                                        
c.... get buffer pointers
      buflen=iabs(fit(3))                                                     
      diskloc=fit(5)                                                          
      curlen=fit(6)                                                           
      kop=1                                                                   
      l=nw                                                                    
      m=0                                                                     
c.... set up access to correct family member
   30 kd=da+m                                                                 
c     call asgrfm (kd,fit)                                                    
      ll=min0(l,fit(4)-kd)                                                    
c.... move data into buffer
   40 i=kd-diskloc                                                            
      if (i.lt.0) go to 50                                                    
      bloc=i                                                                  
      blen=min0(ll,buflen-i)                                                  
      i=0                                                                     
      go to 60                                                                
   50 bloc=0                                                                  
      blen=ll+i                                                               
      if (blen.gt.buflen) go to 80                                            
c.... branch if no overlap
   60 if (blen.le.0) go to 80                                                 
      ll=blen                                                                 
      do 70 k=1,ll                                                            
      fit(k+bloc+7)=w(k+m-i)                                                  
   70 continue                                                                
      fit(3)=-buflen                                                          
      if (i.lt.0) m=m-ll                                                      
      l=l-ll                                                                  
      m=m+ll                                                                  
c.... loop if all requested data has not been transferred
      if (l.ne.0) go to 30                                                    
      return                                                                  
c.... flush the buffer if data is present which is not on disk
   80 if (fit(3).lt.0) then                                                   
      fit(3)=-fit(3)                                                          
      call wdiska (fit(1),fit(8),buflen,diskloc)                              
      fit(7)=max0(fit(7),diskloc+buflen)                                      
      endif                                                                   
c.... for blocks larger than the buffer, write the data directly
   90 if (ll.lt.buflen) go to 110                                             
      if (mod(kd,buflen).ne.0) go to 110                                      
      nr=ll/buflen                                                            
      do 100 n=1,nr                                                           
      call wdiska (fit(1),w(m+1),buflen,kd)                                   
      m=m+buflen                                                              
      kd=kd+buflen                                                            
  100 continue                                                                
      fit(7)=max0(fit(7),kd)                                                  
      l=l-nr*buflen                                                           
      if (l.eq.0) return                                                      
      ll=ll-nr*buflen                                                         
      if (ll.eq.0) go to 30                                                   
c.... initialize the buffer if required
  110 continue                                                                
      diskloc=kd-mod(kd,buflen)                                               
      fit(5)=diskloc                                                          
      curlen=min0(idim(fit(7),diskloc),buflen)                                
      if (curlen.lt.1) go to 120                                              
      if ((kd.ne.diskloc).or.(ll.lt.curlen)) then                             
      call rdiska (fit(1),fit(8),curlen,diskloc)                              
      endif                                                                   
  120 curlen=buflen                                                           
      fit(6)=curlen                                                           
      go to 40                                                                
c
      end                                                                     

*2 4 6 8(1)2 4 6 8(2)2 4 6 8(3)2 4 6 8(4)2 4 6 8(5)2 4 6 8(6)2 4 6 8(7)2 4 6 8(8
*------------------------------- SUBROUTINE asgrfm -----------------------------

      subroutine asgrfm(da,fit)                                               
c
c     assign next family member for random i/o
c
      implicit integer(a-z)                                                   
      dimension fit(8)                                                        
c
      common/frfcm1/mxfrf,ifrf,buflen,fcsize,diskloc,curlen,kop,ier           
      character*8 frfn,frn,kfn                                                
      common/frfcm2/frfn(2,16),frn,kfn                                        
c
      logical fxist                                                           
      character nfn*8,msg*49                                                  
c
c     parameter giving number of record units per integer word
c     for most systems a single character is used as a record unit

c     parameter (ncpw=4)                                                sun
c     parameter (ncpw=1)                                                sgi
      parameter (ncpw=1)                                                dec
c
      data msg/' read requested from nonexistent family member - '/           
c
c.... compute family member index & bias disk address for correct access
      i=da/fcsize                                                             
      da=da-i*fcsize                                                          
c.... get the name of the requested family member
      call nrfnam(frn,i,nfn)                                                  
c.... return if current family member is the desired one
      if(kfn.eq.nfn) return                                                   
c.... flush the buffer if data is present which is not on disk
      if (fit(3).lt.0) then                                                   
      fit(3)=-fit(3)                                                          
      call wdiska (fit(1),fit(8),buflen,diskloc)                              
      fit(7)=max0(fit(7),diskloc+buflen)                                      
      endif                                                                   
c.... determine if requested family member exists
      if (kop.ne.0) go to 20                                                  
      inquire (file=nfn,exist=fxist)                                          
      if (fxist) go to 20                                                     
      ier=-ier                                                                
      if (ier.lt.0) return                                                    
      stop        ' read attempted from nonexistent file '                    
c.... close the current family member
   20 continue                                                                
      if (kfn.ne.'        ') then                                             
      close (fit(1),status='keep')                                            
      endif                                                                   
c.... open/create the requested family member
c     open (fit(1),file=nfn,access='direct',form='unformatted'          vms
c    1,recl=buflen,status='unknown')                                    vms
      open (fit(1),file=nfn,access='direct',form='unformatted'          unix
     1,recl=ncpw*buflen,status='unknown')                               unix
      kfn=nfn                                                                 
      frfn(2,ifrf)=kfn                                                        
      fit(5)=fit(4)                                                           
      diskloc=fit(5)                                                          
      fit(6)=0                                                                
      curlen=fit(6)                                                           
      fit(7)=0                                                                
      if (fxist) fit(7)=fit(4)                                                
      return                                                                  
      end                         


*_______________________________________________________________________________
*-------------------------------------------------------------------------------
*
*                                     BLKDAT
*
*_______________________________________________________________________________
*-------------------------------------------------------------------------------

      block data blkdat

c******************************  PARABLDG ******************************
c
c ROUTINE:  blkdat
c
c PURPOSE:  initialize global variables
c
c PROGRAMMER:  L. Sanford
c
c NOTES:
c
c MODIFICATION HISTORY:  [L. Sanford]
c
c     05/10/02:  Establishment of block data for initialization of
c                global parameters
c
c***********************************************************************

      implicit none

c
c     establishment of global parameters
c

      integer
     :        max_a_array
     :       ,max_beams
     :       ,max_beam_fibers_points
     :       ,max_beam_properties
     :       ,max_beams_times_three
     :       ,max_dof
     :       ,max_dof_plus_1
     :       ,max_iplot
     :       ,max_nodes


      parameter (
     :           max_a_array            = 1500000
     :          ,max_beams              =    5000
     :          ,max_beam_fibers_points =      20
     :          ,max_beam_properties    =      50
     :          ,max_dof                =   30000
     :          ,max_iplot              =  100000
     :          ,max_nodes              =   10000
     :          )


      parameter (
     :           max_beams_times_three = max_beams * 3
     :          ,max_dof_plus_1        = max_dof + 1
     :          )





c
c     geometry and displacement block
c

      integer
     :        id( max_nodes, 6 )

      real
     :      accvctr( max_dof ), dispinc( max_dof ), dispvctr( max_dof )
     :     ,velvctr( max_dof ), xyzloc( max_nodes, 3 )


      common /geometry/
     :        id, xyzloc, dispvctr, dispinc, velvctr, accvctr


c
c     matrix and equation block
c

      integer
     :        maxa( max_dof_plus_1 ),   maxasave( max_dof_plus_1 ),
     :        ncolht( max_dof_plus_1 ), ncolhtsave( max_dof_plus_1 )

      real
     :     a( max_a_array ),         damp( max_a_array ),
     :     residual( max_dof ),      xkeff( max_a_array ),
     :     xkeffsave( max_a_array ), xmass( max_dof )


      common /matrix/
     :        a,      damp,  residual,  xmass,    maxa,
     :        ncolht, xkeff, xkeffsave, maxasave, ncolhtsave


c
c     beam
c

      integer
     :        nbeamn( max_beams, 3 ), nbeamp( max_beams, 2 ),
     :        numfibrs( max_beam_fibers_points ),
     :        nyldink2( max_beams, 9, max_beam_properties ),
     :        nyldstk2( max_beams, 9, max_beam_properties )

      real
     :      area( max_beam_properties )
     :     ,beamdef( max_beams, 7 )
     :     ,bmdefmst( max_beams, 7 )
     :     ,delepsk2( max_beams, 9, max_beam_properties )
     :     ,e( max_beam_properties )
     :     ,elin( max_beam_properties )
     :     ,etab( max_beam_properties )
     :     ,fibrdefn( max_beam_fibers_points, max_beam_properties, 3 )
     :     ,fmstrk2( max_beams, 12 )
     :     ,forcek2( max_beams, 12 )
     :     ,gj( max_beam_properties )
     :     ,resk2( 12 )
     :     ,sheary( max_beam_properties )
     :     ,shearz( max_beam_properties )
     :     ,sigyld( max_beam_properties )
     :     ,straink2( max_beams, 9, max_beam_properties )
     :     ,stressk2( max_beams, 9, max_beam_properties )
     :     ,strmstk2( max_beams, 9, max_beam_properties )

      real
     :      strsmstk2( max_beams, 9, max_beam_properties )
     :     ,tik2( max_beams_times_three, 3 )
     :     ,tjk2( max_beams_times_three, 3 )
     :     ,tk2( 3, 3 )
     :     ,xj( max_beam_properties )
     :     ,xkk2( 12, 12 )
     :     ,xlamdab( max_beam_properties )
     :     ,xlk2( max_beams )
     :     ,xlnmsk2( max_beam_properties )
     :     ,xlobat( 3 )
     :     ,xlobfres( max_beams, 9, 4 )
     :     ,xlobwt( 3 )
     :     ,xnlnmsk2( max_beam_properties )
     :     ,xnu( max_beam_properties )
     :     ,yi( max_beam_properties )
     :     ,yldcntrk2( max_beams, 9, max_beam_properties, 2 )
     :     ,zi( max_beam_properties )


      common /beam/
     :        nbeamn,   nbeamp,   tk2,      tik2,      tjk2,
     :        xlk2,     xkk2,     area,     e, 
     :        yi,       zi,       xj,       xnu,       sheary,
     :        shearz,   fmstrk2,  gj,       resk2,     forcek2,
     :        xlnmsk2,  beamdef,  xnlnmsk2, elin,      sigyld,
     :        xlamdab,  etab,     straink2, xlobat,    fibrdefn,
     :        numfibrs, xlobfres, nyldink2, nyldstk2,  yldcntrk2,
     :        stressk2, delepsk2, strmstk2, strsmstk2, bmdefmst,
     :        xlobwt


c
c     load vectors
c

      integer
     :        loads( 10, 3 ), ncurvloc( max_beam_properties )

      real
     :     ar( 50000 ),  blo( 80 ), br( 50000 ),  bup( 80 ),
     :     bupc( 80 ),   dd( 80 ),  eigval( 80 ), eigvec( max_dof, 80 ),
     :     rtolv( 80 ),  ttt( max_dof ),          vec( 80, 80 ),    
     :     w( max_dof ), xload( 200000, 2 ), xnstif( max_a_array )


      common /loads/
     :        loads, xload, eigvec, eigval,   ttt,
     :        w,     ar,    br,     vec,      dd,    rtolv,
     :        bup,   blo,   bupc,   ncurvloc, xnstif


c

      integer
     :        buflen, curlen, diskloc, fcsize, ier,
     :        ifrf,   kop,    mxfrf


      common /frfcm1/
     :        mxfrf,  ifrf, buflen, fcsize, diskloc,
     :        curlen, kop,  ier



      character * 8
     :              frfn( 2, 16 ), frn, kfn

      common /frfcm2/
     :        frfn, frn, kfn



      integer
     :        iob12( 1080 )
     :       ,iplot( max_iplot )

      real * 4
     :         plot( max_iplot )


      common /taurus/
     :        iob12, plot

      equivalence( plot, iplot )



c
c     initialization of global parameters:  common blocks
c

c     geometry and displacement block

      integer
     :        bnd_id
     :       ,bnd_xyzloc


      parameter (
     :           bnd_id     = max_nodes * 6
     :          ,bnd_xyzloc = max_nodes * 3
     :          )


      data
     :     accvctr  /max_dof * 0.0/
     :    ,dispinc  /max_dof * 0.0/
     :    ,dispvctr /max_dof * 0.0/
     :    ,id       /bnd_id * 0/
     :    ,velvctr  /max_dof * 0.0/
     :    ,xyzloc   /bnd_xyzloc * 0.0/


c     matrix and equation block

      data
     :     a          /max_a_array * 0.0/
     :    ,damp       /max_a_array * 0.0/
     :    ,maxa       /max_dof_plus_1 * 0/
     :    ,maxasave   /max_dof_plus_1 * 0/
     :    ,ncolht     /max_dof_plus_1 * 0/
     :    ,ncolhtsave /max_dof_plus_1 * 0/
     :    ,residual   /max_dof * 0.0/
     :    ,xkeff      /max_a_array * 0.0/
     :    ,xkeffsave  /max_a_array * 0.0/
     :    ,xmass      /max_dof * 0.0/


c     beam

      integer
     :        bnd_beamdef
     :       ,bnd_bmdefmst
     :       ,bnd_delepsk2
     :       ,bnd_fibrdefn
     :       ,bnd_fmstrk2
     :       ,bnd_forcek2
     :       ,bnd_nbeamn
     :       ,bnd_nbeamp
     :       ,bnd_nyldink2
     :       ,bnd_nyldstk2
     :       ,bnd_straink2
     :       ,bnd_stressk2
     :       ,bnd_strmstk2
     :       ,bnd_strsmstk2
     :       ,bnd_tik2
     :       ,bnd_tjk2
     :       ,bnd_xlobfres
     :       ,bnd_yldcntrk2


      parameter (
     :           bnd_beamdef   = max_beams * 7
     :          ,bnd_bmdefmst  = max_beams * 7
     :          ,bnd_delepsk2  = max_beams * 9 * max_beam_properties
     :          ,bnd_fibrdefn  = max_beam_fibers_points * 
     :                           max_beam_properties *
     :                           3
     :          ,bnd_fmstrk2   = max_beams * 12
     :          ,bnd_forcek2   = max_beams * 12
     :          ,bnd_nbeamn    = max_beams * 3
     :          ,bnd_nbeamp    = max_beams * 2 
     :          )
      parameter (
     :           bnd_nyldink2  = max_beams * 9 * max_beam_properties
     :          ,bnd_nyldstk2  = max_beams * 9 * max_beam_properties
     :          ,bnd_straink2  = max_beams * 9 * max_beam_properties
     :          ,bnd_stressk2  = max_beams * 9 * max_beam_properties
     :          ,bnd_strmstk2  = max_beams * 9 * max_beam_properties
     :          ,bnd_strsmstk2 = max_beams * 9 * max_beam_properties
     :          ,bnd_tik2      = max_beams_times_three * 3
     :          ,bnd_tjk2      = max_beams_times_three * 3
     :          ,bnd_xlobfres  = max_beams * 9 * 4
     :          ,bnd_yldcntrk2 = max_beams * 9 * max_beam_properties * 2
     :          )


c
c     NOTE:  "sigyld":  yield value vector for beam elements
c                       is used for branching decisions and
c                       must be initialized to "0.0"

      data
     :     area     /max_beam_properties * 0.0/
     :    ,beamdef  /bnd_beamdef * 0.0/
     :    ,bmdefmst /bnd_bmdefmst * 0.0/
     :    ,delepsk2 /bnd_delepsk2 * 0.0/
     :    ,e        /max_beam_properties * 0.0/
     :    ,elin     /max_beam_properties * 0.0/
     :    ,etab     /max_beam_properties * 0.0/
     :    ,fibrdefn /bnd_fibrdefn * 0.0/
     :    ,fmstrk2  /bnd_fmstrk2 * 0.0/
     :    ,forcek2  /bnd_forcek2 * 0.0/
     :    ,gj       /max_beam_properties * 0.0/
     :    ,nbeamn   /bnd_nbeamn * 0/
     :    ,nbeamp   /bnd_nbeamp * 0/
     :    ,numfibrs /max_beam_fibers_points * 0/
     :    ,nyldink2 /bnd_nyldink2 * 0/
     :    ,nyldstk2 /bnd_nyldstk2 * 0/
     :    ,resk2    /12 * 0.0/
     :    ,sheary   /max_beam_properties * 0.0/
     :    ,shearz   /max_beam_properties * 0.0/

      data
     :     sigyld   /max_beam_properties * 0.0/
     :    ,straink2 /bnd_straink2 * 0.0/
     :    ,stressk2 /bnd_stressk2 * 0.0/
     :    ,strmstk2 /bnd_strmstk2 * 0.0/
     :    ,strsmstk2 /bnd_strsmstk2 * 0.0/
     :    ,tik2      /bnd_tik2 * 0.0/
     :    ,tjk2      /bnd_tjk2 * 0.0/
     :    ,tk2       /9 * 0.0/
     :    ,xj        /max_beam_properties * 0.0/
     :    ,xkk2      /144 * 0.0/
     :    ,xlamdab   /max_beam_properties * 0.0/
     :    ,xlk2      /max_beams * 0.0/
     :    ,xlobat    /-1.0, 0.0, 1.0/
     :    ,xlobfres  /bnd_xlobfres * 0/
     :    ,xlobwt    /0.33333, 1.33333, 0.33333/
     :    ,xnlnmsk2  /max_beam_properties * 0.0/
     :    ,xnu       /max_beam_properties * 0.0/
     :    ,yldcntrk2 /bnd_yldcntrk2 * 0.0/
     :    ,zi        /max_beam_properties * 0.0/


c     load vectors

      integer
     :        bnd_eigvec


      parameter (
     :           bnd_eigvec = max_dof * 80
     :          )


      data
     :      ar       /50000 * 0.0/
     :     ,blo      /80 * 0.0/
     :     ,br       /50000 * 0.0/
     :     ,bup      /80 * 0.0/
     :     ,bupc     /80 * 0.0/
     :     ,dd       /80 * 0.0/
     :     ,eigval   /80 * 0.0/
     :     ,eigvec   /bnd_eigvec * 0.0/
     :     ,loads    /30 * 0/
     :     ,ncurvloc /max_beam_properties * 0/
     :     ,rtolv    /80 * 0.0/
     :     ,ttt      /max_dof * 0.0/
     :     ,vec      /6400 * 0.0/
     :     ,w        /max_dof * 0.0/
     :     ,xload    /400000 * 0.0/
     :     ,xnstif   /max_a_array * 0.0/





c     miscellaneous global initialization

      data
     :     frfn / 32 * '        '/


c
c     iplot initialization template:
c
c     data
c    :     iplot(16) /4/
c    :    ,iplot(18) /2/
c    :    ,iplot(21) /1/
c    :    ,iplot(19) /0/
c    :    ,iplot(31) /0/
c    :    ,iplot(34) /0/


      data iplot / 15 * 0,
     :                  4,
     :                  0,
     :                  2,
     :                  0,
     :                  0,
     :                  1,
     :              6 * 0,
     :                  7,
     :              5 * 0,
     :                  9,
     :          99966 * 0 /


      data
     :     mxfrf /16/


      end

