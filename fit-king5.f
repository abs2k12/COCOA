c
c Thanks to Douglas Haggie for providing this routine which fits the
c king 1966 model to surface brightness and velocity dispersion profiles.
c
cProgramme to fit a King model to the surface brightness profile of
cTrager et al and the velocity dispersion data of Gebhardt
cThis version separates the creation of the model and the construction
cof the goodness of fit
cfrom the specific case of trager
cThis version allows reading a Monte Carlo profiles.dat file
c     x(1) = sigma (km/s)
c     x(2) = King scale radius r_c in arcsec
c     x(3) = W0
c     x(4) = m/l
       real m2l,aaa,bbb,ccc,ddd,xxx
       character input*35,output*35,input1*35,output1*35
      common /cluster/ m2l,dist,A_V,iline,iline1,input,output,
     &                 input1,output1
      common /options/ iout,option
coption = 0: original code (Trager, Gebhardt), 1: temp* files
      external trager,gof
      parameter (ndim = 4)
      real p(ndim + 1,ndim),y(ndim+1),x(ndim)
c      m2l = 2.d0
      open(10,file='fit-king5.conf')
      read (10,*) option
      read (10,*) dist
      read (10,*) A_V
      read (10,*) iline
      read (10,*) iline1
      read (10,*) input
      read (10,*) output
      read (10,*) input1
      read (10,*) output1
      print*,'option,dist,A_V =',option,dist,A_V,iline,iline1,
     & input,output,input1,output1
      read (10,*) (p(1,i),i=1,4)
      read (10,*) (p(i+1,i),i=1,4)
      print*,'p(4,3) =', (p(1,i), i=1,ndim),(p(i+1,i), i=1,ndim)
      close(10)
      iout = 0
ciout = 2 for diagnostic output
      do i = 2,ndim+1
         do j = 1,4
            if (i.ne.j+1) p(i,j) = p(1,j)
         enddo
      enddo
c$$$         p(i,1) = 10.d0
c$$$         p(i,2) = 20.d0
c$$$         p(i,3) = 8.d0
c$$$         p(i,4) = 2.d0
c$$$      enddo
c$$$      p(2,1) = 11.d0
c$$$      p(3,2) = 40.d0
c$$$      p(4,3) = 9.d0
c$$$      p(5,4) = 1.5d0
      do i = 1,ndim+1
         do j = 1,ndim
            x(j) = p(i,j)
         enddo
         if (iout.eq.1) print*,x
         y(i) = gof(x)
c         stop
      enddo
      mp = ndim + 1
      np = ndim
      ftol = 1.d-6
      call amoeba(p,y,mp,np,ndim,ftol,gof,iter)
      print*,'#',iter
      do i = 1,ndim+1
         print*,'#',(p(i,j),j=1,ndim),y(i)
      enddo
cNow output the final profiles
      iout = 1
      do i = 1,ndim
         x(i) = p(1,i)
      enddo
      result = gof(x)
      stop
      end


      function gof(x)
      parameter (ndim = 4)
      real x(ndim),r(206),sb(206),ystart(2)
      character*6 name
      data ifirst/0/
      PARAMETER (KMAXX=1000)
      REAL xp(KMAXX),yp(2,KMAXX),jking,kking,rking(kmaxx),
     &     rhoking(kmaxx),rhoking2(kmaxx),sd(kmaxx),sd2(kmaxx),
     &     v2king(kmaxx),v2king2(kmaxx),vlos(kmaxx),vlos2(kmaxx)
      COMMON /path/ kmax,kount,dxsav,xp,yp
      common /king/W0,rking,rhoking,rhoking2,pradius,sd,sd2,jking,kking,
     &     v2king,v2king2,vlos,vlos2
      common /options/iout,option
      external derivs,rkqs,sdfunc
cMake the King model
      if (iout.eq.1) print*,'gof calling makemodel ',x
      call makemodel(x)
      if (iout.eq.1) print*,'gof calling trager ',x
      sbresult = trager(x)
      vresult = gebhardt(x)
c      gof = 1.0d0*sbresult + 1.d0*vresult
      gof = sbresult + vresult
      print*,'#gof ',x,sbresult,vresult
      return
      end

      subroutine makemodel(x)
      parameter (ndim = 4)
      real x(ndim),r(206),sb(206),ystart(2),msun
      character*6 name
      data ifirst/0/
      PARAMETER (KMAXX=1000,g=6.67d-8,msun=1.989d33,pc=3.085678d18)
      REAL xp(KMAXX),yp(2,KMAXX),jking,kking,rking(kmaxx),
     &     rhoking(kmaxx),rhoking2(kmaxx),sd(kmaxx),sd2(kmaxx),
     &     vlos(kmaxx),vlos2(kmaxx),v2king(kmaxx),v2king2(kmaxx)
      real m2l
       character input*35,output*35,input1*35,output1*35
      common /cluster/ m2l,dist,A_V,iline,iline1,input,output,
     &                 input1,output1
      COMMON /path/ kmax,kount,dxsav,xp,yp
      common /king/W0,rking,rhoking,rhoking2,pradius,sd,sd2,jking,kking,
     &     v2king,v2king2,vlos,vlos2
      common /options/iout,option
      external derivs,rkqs,sdfunc,vlosfunc
cDensity as in King eq.(11)
      rho(W) = pi * kking / jking ** 3 * exp(W - W0) * (-0.6D1 * sqrt(W
     #) + 0.3D1 * sqrt(pi) * erf(sqrt(W)) * exp(W) - 0.4D1 *
     # W ** (0.3D1 / 0.2D1)) * exp(-W) / 0.3D1
cOne-dimensional mean  value of j^2v^2:
      v2(W) = 0.1D1 / 0.10D2 * (-dble(20 * W ** (0.3D1 / 0.2D1)) -
     &     0.30D2
     #* sqrt(dble(W)) + 0.15D2 * sqrt(pi) * erf(sqrt(dble(W)
     #)) * exp(dble(W)) - dble(8 * W ** (0.5D1 / 0.2D1))) / (-0.6D1 * sq
     #rt(dble(W)) + 0.3D1 * sqrt(pi) * erf(sqrt(dble(W))) *
     #exp(dble(W)) - dble(4 * W ** (0.3D1 / 0.2D1)))
      pi = 4.d0 *atan(1.d0)
cThe following are reset later
      kking = 1.d0
      jking = 1.d0
cMake the King model; sigma in km/s,rc in arcsec,m2l in solar units
      sigma = x(1)
      rc = x(2)
      W0 = x(3)
      m2l = x(4)
      jking = 1.d0/(sigma*sqrt(2.d0))  !units 1/(km/s)
      rho0 = (9.d0*sigma**2/(4.d0*pi*g*rc**2))*(1.d10/msun)*pc*
c     &     (pi*dist/648000)
     &     (pi*dist/648000)**(-2)     !units: sigma km/s, g cgs, rc arcsec,
!msun grams, pc cm, dist pc => rho0
      print*, 'rho0 (msun/pc^3) =',rho0
cSince kking = 1 at present:
      rho0onkking = rho(W0)
      kking = rho0/rho0onkking
c      print*,'makemodel jking,kking ',jking,kking
cVariables in use: see King eqs.(19), (21) with some differences of notation
c     x = W
c     y1 = R^2
c     y2 = dR^2/dW
      nvar = 2
      ystart(1) = 0.d0
      ystart(2) = -2.d0/3.d0
      x1 = W0
      x2 = 0.d0
      eps = 1.d-10
      h1 = 0.1d0
      hmin = 0.d0
      dxsav = -0.1d0
      kmax = kmaxx
      if (iout.eq.2) print*,'Entering odeint with parameters',
     &     ystart,nvar,x1,x2,eps,h1,
     &     hmin
      call odeint(ystart,nvar,x1,x2,eps,h1,hmin,nok,nbad,derivs,rkqs)
c      print*,ystart
c      print*,kmax,kount
c      print*,'results from odeint'
C       do i = 1,kount
C          print*,xp(i),(yp(j,i),j=1,2)
C       enddo
c      stop
      if (kount.ge.kmax) then
         print*,'kmax too small in trager'
         stop
      endif
cForm a grid of density versus radius
c      print*,'density profile in makemodel'
      do i = 1,kount
         rking(i) = sqrt(yp(1,i))*rc !units arcsec
         rhoking(i) = rho(xp(i)) !units msun/pc^3
c         print*, rking(i),rhoking(i)
      enddo
c      stop
      rmax = rking(kount)       !units arcsec
      yp1 = 0.d0
      ypn = 0.d0
      call spline(rking,rhoking,kount,yp1,ypn,rhoking2)
!units: rking arcsec, rhoking msun/pc^3
cNow construct table of surface densities
c      print*,'Surface density profile'
      do i = 1,kount
         if (i.eq.kount) then
            s = 0.d0
         else
            a = 0.d0
            b = sqrt(rmax**2 - rking(i)**2)
            pradius = rking(i)  !units of b,pradius: arcsec
c            print*,'Entering qtrap with params ',a,b
            call qtrap(sdfunc,a,b,s)
         endif
         sd(i) = s              !units msun*arcsec/pc^3
c         print*,rking(i),sd(i)
c         if (i.gt.1) sd(i) = sd(i)/sd(1)
      enddo
c      sd(1) = 1.d0
      yp1 = 0.d0
      ypn = 0.d0
      call spline(rking,sd,kount,yp1,ypn,sd2)
!units: rking arcsec, sd msun*arcsec/pc^3
cNext a table of 1d velocity dispersion
c      print*,'1d velocity dispersion'
      do i = 1,kount
         if (i.eq.kount) then
            v2king(i) = 0.d0
         else
            v2king(i) = v2(xp(i))/jking**2 !units (km/s)^2
         endif
c         print*,rking(i),v2king(i),jking,xp(i)
      enddo
      yp1 = 0.d0
      dwdr = 2.d0*rmax/yp(2,kount) !units arcsec/r_c^2
      dwdr = dwdr/rc**2               !new: units arcsec^(-1)
c      ypn = (1.d0/7.d0)*dwdr
      ypn = (1.d0/7.d0)*dwdr/jking**2
c      ypn = 0.d0
c      print*,'ypn=',ypn
      call spline(rking,v2king,kount,yp1,ypn,v2king2)
c      iout = 1
c      if (iout.eq.1) print*,'#Line of sight velocity profile'
      do i = 1,kount
         if (i.eq.kount) then
            s = 0.d0
            vlos(i) = 0.d0
         else
            a = 0.d0
            b = sqrt(rmax**2 - rking(i)**2)
            pradius = rking(i)  !units arcsec
c            print*,'Entering qtrap with params ',a,b
            call qtrap(vlosfunc,a,b,s)
!units of s: msun*(km/s)^2*arcsec/pc^3
!units of sd: msun*pc^(-3)*arcsec
            vlos(i) = sqrt(s/sd(i))
         endif
c         if (iout.eq.1) print*,rking(i),vlos(i)
      enddo
      if (iout.eq.1) print*
c      stop
      yp1 = 0.d0
cThe following limit is derived by inspection of a graph!
      ypn = 0.d0
      call spline(rking,vlos,kount,yp1,ypn,vlos2)
      return
      end

      function vlosfunc(x)
      parameter (kmaxx=1000)
      real rking(kmaxx),rhoking(kmaxx),rhoking2(kmaxx),sd(kmaxx),
     &     sd2(kmaxx)
      real dxsav,xp(kmaxx),yp(2,kmaxx),jking,kking,v2king(kmaxx),
     &     v2king2(kmaxx),vlos(kmaxx),vlos2(kmaxx)
      common /king/W0,rking,rhoking,rhoking2,pradius,sd,sd2,jking,kking,
     &     v2king,v2king2,vlos,vlos2
      COMMON /path/ kmax,kount,dxsav,xp,yp
      radius = sqrt(x**2 + pradius**2)
c      print*,'sdfunc calling splint ',rking(1),rhoking(1),rhoking2(1),kount,radius
      call splint(rking,v2king,v2king2,kount,radius,y1)
      call splint(rking,rhoking,rhoking2,kount,radius,y2)
!units: rking arcsec, v2king (km/s)^2, rhoking msun/pc^3
c      print*,'In sdfunc ',rking(1),rhoking(1),rhoking2(1),kount,radius,y
      vlosfunc = 2.d0*y1*y2     !units msun*(km/s)^2/pc^3
c      print*,x,vlosfunc,y1,y2,radius,x,pradius
c      stop
      return
      end


      function gebhardt(x)
cComparison with the los velocity dispersion profile
      parameter (ndim = 4)
      real x(ndim),r(206),sb(206),ystart(2),vdp(100)
      save r,sb,vdp,ndata
      character*6 name
      data ifirst/0/
      PARAMETER (KMAXX=1000,Vsun = 4.83d0)
      REAL xp(KMAXX),yp(2,KMAXX),jking,kking,rking(kmaxx),
     &     rhoking(kmaxx),rhoking2(kmaxx),sd(kmaxx),sd2(kmaxx)
      real m2l
       character input*35,output*35,input1*35,output1*35
      common /cluster/ m2l,dist,A_V,iline,iline1,input,output,
     &                 input1,output1
      COMMON /path/ kmax,kount,dxsav,xp,yp
      real v2king(kmaxx),v2king2(kmaxx),vlos(kmaxx),vlos2(kmaxx),xxx
      common /king/W0,rking,rhoking,rhoking2,pradius,sd,sd2,jking,kking,
     &     v2king,v2king2,vlos,vlos2
      common /options/iout,option
      external derivs,rkqs,sdfunc
c      print*,'Entered gebhardt'
      pi = 4.d0*atan(1.d0)
      if (ifirst.eq.0) then
cRead trager et al data
         if (option.eq.0) then
            ndata = iline1
            open (10,file=input1)
            read (10,*)
            do i = 1,iline1
               read (10,*) r(i),vdp(i)
c     Convert to arcsec
               r(i) = 60.d0*r(i)
            enddo
         else if (option.eq.1) then
            ndata = iline1
            open (10,file=input1)
 10         continue
            do i = 1,ndata
              read (10,*,err=10,end=20) r(i),vdp(i)
C                read (10,*,err=10,end=20) aaa,bbb,ccc,r(i),xxx,vdp(i),
C      &               ddd,yyy,zzz
               vdp(i) = sqrt(vdp(i))
c     Convert to arcsec
c               r(i) = 3600.d0*180*(r(i)/dist)/pi
            enddo
         else
            print*,'option',option,'not recognised'
         endif
 20      continue
         ifirst = 1
         close(10)
      endif
      if (iout.eq.1) then
	print*,''
      	print*,'Velocity dispersion profile'
      endif
      gebhardt = 0.d0
      open(10,file=output1)
      do i = 1,ndata
         radius = r(i)
         call splint(rking,vlos,vlos2,kount,radius,y)
         if (radius.gt.rking(kount)) y = 0.0
cHere y is the surface density in solar masses per square arcsec.
cIt is then converted to solar luminosities, then to absolute magnitude then apparent magnitude.
cWe neglect absorption at present
         if (iout.eq.1) then
            write(10,*) radius,vdp(i),y
         endif
         gebhardt = gebhardt + (vdp(i) - y)**2
      enddo
      close(10)
c     print*,'#gebhardt: ',x,gebhardt
C       if (iout.eq.1) then
C          print*,'Normalised surface density'
C          do i = 1,kount
C             print*, rking(i),sd(i)
C          enddo
C       endif
c      stop
      return
      end


      function trager(x)
      parameter (ndim = 4)
      real x(ndim),r(206),sb(206),ystart(2)
      save r,sb,ndata
      character*6 name
      data ifirst/0/
      PARAMETER (KMAXX=1000,Vsun = 4.83d0)
      REAL xp(KMAXX),yp(2,KMAXX),jking,kking,rking(kmaxx),
     &     rhoking(kmaxx),rhoking2(kmaxx),sd(kmaxx),sd2(kmaxx)
      real m2l,halfsd,rhalf,err
       character input*35,output*35,input1*35,output1*35
      common /cluster/ m2l,dist,A_V,iline,iline1,input,output,
     &                 input1,output1
      COMMON /path/ kmax,kount,dxsav,xp,yp
      real v2king(kmaxx),v2king2(kmaxx),vlos(kmaxx),vlos2(kmaxx)
      common /king/W0,rking,rhoking,rhoking2,pradius,sd,sd2,jking,kking,
     &     v2king,v2king2,vlos,vlos2
      common /options/iout,option
      external derivs,rkqs,sdfunc
c      print*,'Entered trager'
      pi = 4.d0*atan(1.d0)
      if (ifirst.eq.0) then
         if(option.eq.0) then
cRead trager et al data
            ndata = iline
            open (10,file=input)
            do i = 1,iline
               read (10,*) name,r(i),sb(i) !units sb: mV/arcsec^2
cConvert to arcsec
               r(i) = 10.d0**r(i) !units arcsec
            enddo
            close(10)
         else if (option.eq.1) then
            ndata = iline
            open (10,file=input)
 10         continue
            do i = 1,ndata
              read (10,*,err=10,end=20) r(i),sb(i)
C               read (10,*,err=10,end=20) aaa,bbb,ccc,r(i),sb(i),ddd,
C      &              xxx,yyy,zzz
cHere sb is assumed to be in
c V mag per arcsec**2
            enddo
         else
            print*,'option',option,'not recognised'
         endif
 20      continue
         close (10)
         ifirst = 1
      endif
C Compare in magnitudes
c      print*,'Surface brightness data'
      trager = 0.d0
      open(10,file=output)
      do i = 1,ndata
         radius = r(i)          !units arcsec
c         print*,'trager calling splint ',rking(i),sd(i),sd2(i),kount,radius
c         call splint(rking,sd,sd2,kount,radius,y)
         if (radius.lt.rking(kount)) then
            call splint(rking,sd,sd2,kount,radius,y)
!units sd is msun*arcsec/pc^3
            pc2arcsec = (1.d0/dist)*(648000/pi)
            y = y/pc2arcsec**3
!units now msun/arcsec**2 (new)
         else
           y = 0.d0
         endif
c         print*,'trager =',y
cHere y is the surface density in solar masses per square arcsec.
cIt is then converted to solar luminosities, then to absolute magnitude then apparent magnitude.
cWe neglect absorption at present
         if (y.gt.0.d0) then
            y = y/m2l
            sbking = Vsun - 2.5d0*log10(y)
            sbking = sbking + 5.d0*log10(dist/10.d0) + 3.d0*A_V
         else
            sbking = 100.d0
         endif
         if (iout.eq.1) then
            write(10,*) radius,sb(i),sbking
         endif
         trager = trager + (sb(i) - sbking)**2
      enddo
      close(10)
      halfsd = sd(1)/2.0
C       if (iout.eq.1) then
c      rhalf = 0.0
c      err = 0.0
c      print*, kount
      call ratint(sd,rking,kount,halfsd,rhalf,err)
      print*,'2D-rc =',rhalf,'error=',err
c      print*,'Normalised surface density'
c      do i = 1,kount
c        print*, rking(i),sd(i)
c      enddo
C       endif
c      stop
      return
      end

      function sdfunc(x)
      parameter (kmaxx=1000)
      real rking(kmaxx),rhoking(kmaxx),rhoking2(kmaxx),sd(kmaxx),
     &     sd2(kmaxx)
      real dxsav,xp(kmaxx),yp(2,kmaxx),jking,kking,v2king(kmaxx),
     &     v2king2(kmaxx),vlos(kmaxx),vlos2(kmaxx)
      common /king/W0,rking,rhoking,rhoking2,pradius,sd,sd2,jking,kking,
     &     v2king,v2king2,vlos,vlos2
      COMMON /path/ kmax,kount,dxsav,xp,yp
      radius = sqrt(x**2 + pradius**2) !units arcsec
c      print*,'sdfunc calling splint ',rking(1),rhoking(1),rhoking2(1),kount,radius
      call splint(rking,rhoking,rhoking2,kount,radius,y)
!units: rking arcsec, rhoking msun/pc^3
c      print*,'In sdfunc ',rking(1),rhoking(1),rhoking2(1),kount,radius,y
      sdfunc = 2.d0*y
c      print*,x,sdfunc
      return
      end

      SUBROUTINE trapzd(func,a,b,s,n)
      INTEGER n
      REAL a,b,s,func
      EXTERNAL func
      INTEGER it,j
      REAL del,sum,tnm,x
      if (n.eq.1) then
        s=0.5*(b-a)*(func(a)+func(b))
      else
        it=2**(n-2)
        tnm=it
        del=(b-a)/tnm
        x=a+0.5*del
        sum=0.
        do 11 j=1,it
          sum=sum+func(x)
          x=x+del
11      continue
        s=0.5*(s+(b-a)*sum/tnm)
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software ]vrt43D04.


      SUBROUTINE qtrap(func,a,b,s)
      INTEGER JMAX
      REAL a,b,func,s,EPS
      EXTERNAL func
c      PARAMETER (EPS=1.e-6, JMAX=20)
      PARAMETER (EPS=1.e-5, JMAX=20)
CU    USES trapzd
      INTEGER j
      REAL olds
      olds=-1.e30
      do 11 j=1,JMAX
        call trapzd(func,a,b,s,j)
        if (abs(s-olds).lt.EPS*abs(olds)) return
        oldolds = olds
        olds=s
11    continue
      print*,'too many steps in qtrap, s, olds',s,oldolds
c      pause 'too many steps in qtrap'
      END
C  (C) Copr. 1986-92 Numerical Recipes Software ]vrt43D04.


      subroutine derivs(x,y,dydx)
      parameter (kmaxx=1000)
      real y(2),dydx(2),j,k,jking,kking,rking(kmaxx),rhoking(kmaxx),
     &     rhoking2(kmaxx),sd(kmaxx),sd2(kmaxx),v2king(kmaxx),
     &     v2king2(kmaxx),vlos(kmaxx),vlos2(kmaxx)
      common /king/W0,rking,rhoking,rhoking2,pradius,sd,sd2,jking,kking,
     &     v2king,v2king2,vlos,vlos2
      common /options/iout,option
      rho(W) = pi * k / j ** 3 * exp(W - W0) * (-0.6D1 * sqrt(W
     #) + 0.3D1 * sqrt(pi) * erf(sqrt(W)) * exp(W) - 0.4D1 *
     # W ** (0.3D1 / 0.2D1)) * exp(-W) / 0.3D1
      rhodiff(W) = pi * k / j ** 3 * exp(W - W0) * (-0.3D1 * W ** (-
     #0.1D1 / 0.2D1) + 0.3D1 * exp(-W) * W ** (-0.1D1 / 0.2D1) * exp(W)
     #+ 0.3D1 * sqrt(pi) * erf(sqrt(W)) * exp(W) - 0.6D1 * s
     #qrt(W)) * exp(-W) / 0.3D1
cThe following two assignments are not important in this code, but would be to get the correct
cscaling of other parameters in the King model
      j = jking
      k = kking
      pi = 4.d0*atan(1.d0)
      nvar = 2
      dydx(1) = y(2)
      if (y(1).eq.0.d0) then
         a = -1.5d0
         b = (27.d0/40.d0)*rhodiff(W0)/rho(W0)
         dydx(2) = -2.d0*b/a**3
      else
         dydx(2) = 1.5d0*y(2)**2/y(1) +
     &        2.25d0*(rho(x)/rho(W0))*y(2)**3/y(1)
      endif
      if (iout.eq.2) print*,'derivs ',x,y,dydx,rho(x),rho(W0),W0,pi,k,j
      return
      end


      SUBROUTINE rkck(y,dydx,n,x,h,yout,yerr,derivs)
      INTEGER n,NMAX
      REAL h,x,dydx(n),y(n),yerr(n),yout(n)
      EXTERNAL derivs
      PARAMETER (NMAX=50)
CU    USES derivs
      INTEGER i
      REAL ak2(NMAX),ak3(NMAX),ak4(NMAX),ak5(NMAX),ak6(NMAX),
     *ytemp(NMAX),A2,A3,A4,A5,A6,B21,B31,B32,B41,B42,B43,B51,B52,B53,
     *B54,B61,B62,B63,B64,B65,C1,C3,C4,C6,DC1,DC3,DC4,DC5,DC6
      PARAMETER (A2=.2,A3=.3,A4=.6,A5=1.,A6=.875,B21=.2,B31=3./40.,
     *B32=9./40.,B41=.3,B42=-.9,B43=1.2,B51=-11./54.,B52=2.5,
     *B53=-70./27.,B54=35./27.,B61=1631./55296.,B62=175./512.,
     *B63=575./13824.,B64=44275./110592.,B65=253./4096.,C1=37./378.,
     *C3=250./621.,C4=125./594.,C6=512./1771.,DC1=C1-2825./27648.,
     *DC3=C3-18575./48384.,DC4=C4-13525./55296.,DC5=-277./14336.,
     *DC6=C6-.25)
      do 11 i=1,n
        ytemp(i)=y(i)+B21*h*dydx(i)
11    continue
      call derivs(x+A2*h,ytemp,ak2)
      do 12 i=1,n
        ytemp(i)=y(i)+h*(B31*dydx(i)+B32*ak2(i))
12    continue
      call derivs(x+A3*h,ytemp,ak3)
      do 13 i=1,n
        ytemp(i)=y(i)+h*(B41*dydx(i)+B42*ak2(i)+B43*ak3(i))
13    continue
      call derivs(x+A4*h,ytemp,ak4)
      do 14 i=1,n
        ytemp(i)=y(i)+h*(B51*dydx(i)+B52*ak2(i)+B53*ak3(i)+B54*ak4(i))
14    continue
      call derivs(x+A5*h,ytemp,ak5)
      do 15 i=1,n
        ytemp(i)=y(i)+h*(B61*dydx(i)+B62*ak2(i)+B63*ak3(i)+B64*ak4(i)+
     *B65*ak5(i))
15    continue
      call derivs(x+A6*h,ytemp,ak6)
      do 16 i=1,n
        yout(i)=y(i)+h*(C1*dydx(i)+C3*ak3(i)+C4*ak4(i)+C6*ak6(i))
16    continue
      do 17 i=1,n
        yerr(i)=h*(DC1*dydx(i)+DC3*ak3(i)+DC4*ak4(i)+DC5*ak5(i)+DC6*
     *ak6(i))
17    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software ]vrt43D04.


      SUBROUTINE rkqs(y,dydx,n,x,htry,eps,yscal,hdid,hnext,derivs)
      INTEGER n,NMAX
      REAL eps,hdid,hnext,htry,x,dydx(n),y(n),yscal(n)
      EXTERNAL derivs
      PARAMETER (NMAX=50)
CU    USES derivs,rkck
      INTEGER i
      REAL errmax,h,xnew,yerr(NMAX),ytemp(NMAX),SAFETY,PGROW,PSHRNK,
     *ERRCON
      PARAMETER (SAFETY=0.9,PGROW=-.2,PSHRNK=-.25,ERRCON=1.89e-4)
      h=htry
1     call rkck(y,dydx,n,x,h,ytemp,yerr,derivs)
      errmax=0.
      do 11 i=1,n
        errmax=max(errmax,abs(yerr(i)/yscal(i)))
11    continue
      errmax=errmax/eps
      if(errmax.gt.1.)then
        h=SAFETY*h*(errmax**PSHRNK)
        if(h.lt.0.1*h)then
          h=.1*h
        endif
        xnew=x+h
        if(xnew.eq.x)pause 'stepsize underflow in rkqs'
        goto 1
      else
        if(errmax.gt.ERRCON)then
          hnext=SAFETY*h*(errmax**PGROW)
        else
          hnext=5.*h
        endif
        hdid=h
        x=x+h
        do 12 i=1,n
          y(i)=ytemp(i)
12      continue
        return
      endif
      END
C  (C) Copr. 1986-92 Numerical Recipes Software ]vrt43D04.


      SUBROUTINE odeint(ystart,nvar,x1,x2,eps,h1,hmin,nok,nbad,derivs,
     *rkqs)
      INTEGER nbad,nok,nvar,KMAXX,MAXSTP,NMAX
      REAL eps,h1,hmin,x1,x2,ystart(nvar),TINY
      EXTERNAL derivs,rkqs
      PARAMETER (MAXSTP=10000,NMAX=2,KMAXX=1000,TINY=1.e-30)
      INTEGER i,kmax,kount,nstp
      REAL dxsav,h,hdid,hnext,x,xsav,dydx(NMAX),xp(KMAXX),y(NMAX),
     *yp(NMAX,KMAXX),yscal(NMAX)
      COMMON /path/ kmax,kount,dxsav,xp,yp
      x=x1
      h=sign(h1,x2-x1)
      nok=0
      nbad=0
      kount=0
      do 11 i=1,nvar
        y(i)=ystart(i)
11    continue
      if (kmax.gt.0) xsav=x-2.*dxsav
      do 16 nstp=1,MAXSTP
        call derivs(x,y,dydx)
        do 12 i=1,nvar
          yscal(i)=abs(y(i))+abs(h*dydx(i))+TINY
12      continue
        if(kmax.gt.0)then
          if(abs(x-xsav).gt.abs(dxsav)) then
            if(kount.lt.kmax-1)then
              kount=kount+1
              xp(kount)=x
              do 13 i=1,nvar
                yp(i,kount)=y(i)
13            continue
              xsav=x
            endif
          endif
        endif
        if((x+h-x2)*(x+h-x1).gt.0.) h=x2-x
c        print*,'Calling
        call rkqs(y,dydx,nvar,x,h,eps,yscal,hdid,hnext,derivs)
        if(hdid.eq.h)then
          nok=nok+1
        else
          nbad=nbad+1
        endif
        if((x-x2)*(x2-x1).ge.0.)then
          do 14 i=1,nvar
            ystart(i)=y(i)
14        continue
          if(kmax.ne.0)then
            kount=kount+1
            xp(kount)=x
            do 15 i=1,nvar
              yp(i,kount)=y(i)
15          continue
          endif
          return
        endif
        if(abs(hnext).lt.hmin) pause
     *'stepsize smaller than minimum in odeint'
        h=hnext
16    continue
      pause 'too many steps in odeint'
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software ]vrt43D04.








      SUBROUTINE amoeba(p,y,mp,np,ndim,ftol,funk,iter)
      INTEGER iter,mp,ndim,np,NMAX,ITMAX
      REAL ftol,p(mp,np),y(mp),funk
      PARAMETER (NMAX=20,ITMAX=5000)
      EXTERNAL funk
CU    USES amotry,funk
      INTEGER i,ihi,ilo,inhi,j,m,n
      REAL rtol,sum,swap,ysave,ytry,psum(NMAX),amotry
      iter=0
1     do 12 n=1,ndim
        sum=0.
        do 11 m=1,ndim+1
          sum=sum+p(m,n)
11      continue
        psum(n)=sum
12    continue
2     ilo=1
      if (y(1).gt.y(2)) then
        ihi=1
        inhi=2
      else
        ihi=2
        inhi=1
      endif
      do 13 i=1,ndim+1
        if(y(i).le.y(ilo)) ilo=i
        if(y(i).gt.y(ihi)) then
          inhi=ihi
          ihi=i
        else if(y(i).gt.y(inhi)) then
          if(i.ne.ihi) inhi=i
        endif
13    continue
      rtol=2.*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo)))
      if (rtol.lt.ftol) then
        swap=y(1)
        y(1)=y(ilo)
        y(ilo)=swap
        do 14 n=1,ndim
          swap=p(1,n)
          p(1,n)=p(ilo,n)
          p(ilo,n)=swap
c To ensure positivity:
          p(1,n)=abs(p(1,n))
          p(ilo,n)=abs(p(ilo,n))

14      continue
        return
      endif
      if (iter.ge.ITMAX) pause 'ITMAX exceeded in amoeba'
      iter=iter+2
      ytry=amotry(p,y,psum,mp,np,ndim,funk,ihi,-1.0)
      if (ytry.le.y(ilo)) then
        ytry=amotry(p,y,psum,mp,np,ndim,funk,ihi,2.0)
      else if (ytry.ge.y(inhi)) then
        ysave=y(ihi)
        ytry=amotry(p,y,psum,mp,np,ndim,funk,ihi,0.5)
        if (ytry.ge.ysave) then
          do 16 i=1,ndim+1
            if(i.ne.ilo)then
              do 15 j=1,ndim
                psum(j)=0.5*(p(i,j)+p(ilo,j))
                p(i,j)=psum(j)
c To ensure positivity:
                p(i,j) = abs(p(i,j))
15            continue
              y(i)=funk(psum)
            endif
16        continue
          iter=iter+ndim
          goto 1
        endif
      else
        iter=iter-1
      endif
      goto 2
      END
C  (C) Copr. 1986-92 Numerical Recipes Software ]vrt43D04.


      FUNCTION amotry(p,y,psum,mp,np,ndim,funk,ihi,fac)
      INTEGER ihi,mp,ndim,np,NMAX
      REAL amotry,fac,p(mp,np),psum(np),y(mp),funk
      PARAMETER (NMAX=20)
      EXTERNAL funk
CU    USES funk
      INTEGER j
      REAL fac1,fac2,ytry,ptry(NMAX)
      fac1=(1.-fac)/ndim
      fac2=fac1-fac
      do 11 j=1,ndim
        ptry(j)=psum(j)*fac1-p(ihi,j)*fac2
c To ensure positivity:
        ptry(j) = abs(ptry(j))
11    continue

      ytry=funk(ptry)
      if (ytry.lt.y(ihi)) then
        y(ihi)=ytry
        do 12 j=1,ndim
          psum(j)=psum(j)-p(ihi,j)+ptry(j)
          p(ihi,j)=ptry(j)
12      continue
      endif
      amotry=ytry
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software ]vrt43D04.

      SUBROUTINE splint(xa,ya,y2a,n,x,y)
      INTEGER n
      REAL x,y,xa(n),y2a(n),ya(n)
      INTEGER k,khi,klo
      REAL a,b,h
      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.) pause 'bad xa input in splint'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*
     &   (h**2)/6.
c      print*,h,a,b,y
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software ]vrt43D04.



      SUBROUTINE spline(x,y,n,yp1,ypn,y2)
      INTEGER n,NMAX
      REAL yp1,ypn,x(n),y(n),y2(n)
      PARAMETER (NMAX=500)
      INTEGER i,k
      REAL p,qn,sig,un,u(NMAX)
      if (yp1.gt..99e30) then
        y2(1)=0.
        u(1)=0.
      else
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+
     *1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*
     *u(i-1))/p
11    continue
      if (ypn.gt..99e30) then
        qn=0.
        un=0.
      else
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software ]vrt43D04.

      SUBROUTINE ratint(xa,ya,n,x,y,dy)
      INTEGER n,NMAX
      REAL dy,x,y,xa(n),ya(n),TINY
      PARAMETER (NMAX=250,TINY=1.e-25)
      INTEGER i,m,ns
      REAL dd,h,hh,t,w,c(NMAX),d(NMAX)
      ns=1
      hh=abs(x-xa(1))
      do 11 i=1,n
        h=abs(x-xa(i))
        if (h.eq.0.)then
          y=ya(i)
          dy=0.0
          return
        else if (h.lt.hh) then
          ns=i
          hh=h
        endif
        c(i)=ya(i)
        d(i)=ya(i)+TINY
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          w=c(i+1)-d(i)
          h=xa(i+m)-x
          t=(xa(i)-x)*d(i)/h
          dd=t-c(i+1)
          if(dd.eq.0.)pause 'failure in ratint'
          dd=w/dd
          d(i)=c(i+1)*dd
          c(i)=t*dd
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software ]vrt43D04.
