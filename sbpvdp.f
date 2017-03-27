****************************************************************************
*
*     description of output file with SBP and VDP
*
*  1   - r         - poition of the geometrical average for the shell (pc)
*  2   - sbp       - SBP/pc^2
*  3   - vdp       - VDP/pc^2
*  4   - r         - poition of the geometrical average for the shell (arcsec)
*  5   - sbp       - SBP/arcsec^2
*  6   - vdp       - VDP/arcsec^2
*
******************************************************************************
*
      integer i,k,im,idd1,idd2,ikb,ik1,ik2,irc,itot,
     &        ngrid,nstar,indb,iinter1,iinter2,j
      parameter (ngrid = 51,nstar=20000000)
      real*8 rr,rp,rpx,rpy,vrp,vtx,vty,sm1,sm2,slum1,slum2,rad1,rad2,
     &     spin1,spin2,ap,ecc,mv,mbv,mi,mv1,mbv1,mi1,mv2,mbv2,mi2,
     &     R_Sun,rmin,rmax,modulus,redenning,scal,rtt,pi,vr,vt,
     &     max_mag,
     &     r(nstar),slum(nstar),vs(nstar),rlog(ngrid),sb(ngrid),
     &     vd(ngrid),rob(ngrid),den(ngrid),mag_v(nstar),vsx(nstar),
     &     vsy(nstar),vsxy(nstar),vdx(ngrid),vdy(ngrid),vdxy(ngrid)
      character input*20,output_sbp_vdp*20
      parameter (V_sun = 4.83d0)
*
      write(6,*) 'input file'
      read(5,*) input
      write(6,*) 'SBP-VDP output file'
      read(5,*) output_sbp_vdp
      write(6,*) 'distance to Sun (pc)'
      read(5,*) R_Sun
      write(6,*) 'redenning'
      read(5,*) redenning
      write(6,*) 'max_mag'
      read(5,*) max_mag
*
      open (7,file=input)
      open (8,file=output_sbp_vdp)
      open (100,file="check.dat")		
*
      rtt = 0.d0
      pi = 4.0d0*atan(1.0d0)
      scal = 180.d0*3600.d0/R_sun/pi
      modulus = 5.d0*log10(R_Sun/10.d0) + 3.d0*redenning + V_sun
*      
*   read projected snapshot
*
* 10   format(1x,i9,1p5e16.8,2i9,3i4,1p8e16.8,2i11,1p14e16.8,i6)
 40   format(1x,1p6e16.8)
*
      read(7,*)
      read(7,*)
      i = 0
 
C  21   read(7,*)
 20   i = i + 1     
*



      read(7,*,end=30) indb,im,rpx,rpy,rpz,vrp,vtx,vty,idd1,idd2,ikb,
     &     ik1,ik2,sm1,
     &     sm2,slum1,slum2,rad1,rad2,spin1,spin2,iinter1,iinter2,
     &     ap,ecc,mv,mbv,mi,mv1,mbv1,mi1,mv2,mbv2,mi2,rr,vr,vt


    
*      write(100,*), indb,im,rpx,rpy,rpz,vrp,vtx,vty,idd1,idd2,
*     &     ikb,ik1,ik2,
*     &     sm1,sm2,slum1,slum2,rad1,rad2,spin1,spin2,iinter1,iinter2,
*     &     ap,ecc,mv,mbv,mi,mv1,mbv1,mi1,mv2,mbv2,mi2,rr,vr,vt

      r(i) = sqrt(rpx**2 + rpy**2)
      if(rtt.lt.r(i)) rtt = r(i)
*      slum(i) = slum1 + slum2
        if (j.eq.im) then
          mag_v(i) = mv2
          slum(i) = 10.0d0**(0.4d0*V_sun)*10.0d0**(-0.4d0*mv2)
        else
          mag_v(i) = mv1
          slum(i) = 10.0d0**(0.4d0*V_sun)*10.0d0**(-0.4d0*mv1)
        endif
C       slum(i) = 10.0d0**(0.4d0*V_sun)*10.0d0**(-0.4d0*mv)
      vs(i) = vrp**2
      vsx(i) = rpx**2
      vsy(i) = rpy**2
      vsxy(i) = rpx**2 + rpy**2

      
      write(100,*), indb,mag_v(i),slum(i)

      if (indb.eq.1) j = im
     
      go to 20		    			
*
 30   continue
*
      itot = i - 1
*
      do i = 1,ngrid
         rlog(i) = (-2.8d0 + 2.8d0*dfloat((i-1))/dfloat((ngrid - 1)))
         rlog(i) = 10.d0**rlog(i)*rtt
         sb(i) = 0.d0
         vd(i) = 0.d0
         den(i) = 0.d0
         vdx(i) = 0.d0
         vdy(i) = 0.d0
         vdxy(i) = 0.d0
      enddo
*      
      do i = 1,itot
         do k = 1,ngrid
            if(k.eq.1) then
               rmin = 0.d0
               rmax = rlog(k)
               rob(k) = 0.5d0*rmax
               sur = pi*(rmax*rmax - rmin*rmin)
            else
               rmin = rlog(k-1)
               rmax = rlog(k)
c               rob(k) = sqrt(rmax*rmin)
c               rob(k) = 0.5d0*(rmax + rmin)
c               rob(k) = sqrt(0.5d0*(rmax**2 - rmin**2)/log(rmax/rmin))
c               rob(k) = rmin
               rob(k) = 1.d0/sqrt(1.d0/(rmin*rmax))
               sur = pi*(rmax*rmax - rmin*rmin)
            endif 
            if(r(i).gt.rmin.and.r(i).le.rmax) then
              sb(k) = sb(k) + slum(i)/sur
              if(mag_v(i).lt.max_mag) then
                vd(k) = vd(k) + vs(i)/sur
                den(k) = den(k) + 1.d0/sur
                vdx(k) = vd(k) + vsx(i)/sur
                vdy(k) = vd(k) + vsy(i)/sur
                vdxy(k) = vd(k) + vsxy(i)/sur
              endif
            endif
         enddo
      enddo 
*
      do i = 1,ngrid  
         write(8,40) rob(i),-2.5d0*log10(sb(i))+modulus,vd(i)/den(i),
     &               rob(i)*scal,-2.5d0*log10(sb(i)/scal**2)+modulus,   
     &               vd(i)/den(i)
      enddo
* 
      close (7)
      close (8)
      close (10)
	
*
      stop
      end
*
