*
*    Program to delete overlapping stars from different fields
*
*****************************************************************************

      integer i,k,im,idd1,idd2,ikb,ik1,ik2,iinter1,iinter2,irc,itot,
     &        ngrid,nstar,indb,j

      parameter (nstar=1000000)

      real*8 x,y,mag,difm,num,rtx,rty,difx,dify,dife,
     &     mag_v(nstar),fwhm,pixs,sizx,sizy,ovlerap,frame1,
     &     overframe,r(nstar),difpix,
     &     rpx(nstar),rpy(nstar),errmag(nstar),starnum(nstar)

      parameter (V_sun = 4.83d0)
*
*
      open (7,file="catalogue.dat")
      open (8,file="overlap.dat")
      open (9,file="all.dat")
      open (10,file="fwhm.dat")
      rtx = 0.d0
      rty = 0.d0

      read(10,*,end=30) fwhm,pixs,sizx,sizy,ovlerap
      frame1=sizx*pixs
      overframe=frame1*overlap
      difpix=fwhm*pixs

C  10   format(1x,i9,1p5e16.8,2i9,3i4,1p8e16.8,2i11,1p14e16.8,i6)
c 40   format(1x,1p7e16.8)
*
      i = 0
      j = -1

 20   i = i + 1
*

      read(7,*,end=30) x,y,mag,difm,num
c      r(i)= sqrt(rpx**2 + rpy**2)
      rpx(i) = x
      if(rtx.lt.rpx(i)) rtx = rpx(i)
      rpy(i) = y
      if(rty.lt.rpy(i)) rty = rpy(i)
      mag_v(i) = mag
      errmag(i) = difm
      starnum(i) = num
      write(9,*) x,y,mag,difm,num

*
      go to 20
*
 30   continue
*
      itot = i - 1
*
      do i = 1,itot
         do k = i+1,itot
c          if(((rpx.ge.overframe).and.(rpx.le.frame1)).or
c     &    ((rpx.ge.overframe*2).and.(rpx.le.frame1*2))
            difx = abs(rpx(i))-abs(rpx(k))
            dify = abs(rpy(i))-abs(rpy(k))
            dife = abs(errmag(i))-abs(errmag(k))

            if(((difx.le.difpix).and.(dify.le.difpix))
     &      .and.((difx.gt.-difpix).and.(dify.gt.-difpix))) then
              j = 1
              dife = abs(errmag(i))-abs(errmag(k))
              if (dife.lt.0) then
               write(8,*) rpx(i),rpy(i),mag_v(i),errmag(i),starnum(i)
              else
               write(8,*) rpx(k),rpy(k),mag_v(k),errmag(k),starnum(k)
              endif
            endif
         enddo
      enddo

*
      close (7)
      close (8)
*
      stop
      end
*
