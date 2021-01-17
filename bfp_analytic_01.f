      program bfpa

      real x(10000),y(10000),z(10000)
      real ma(3,3),mx(3,1),mb(3,1),mm(3,3)
      real gxl,gyl,tzp
      integer n
      character*48 infile

c     PROGRAM bfp_analytic.f
c     Steve Morice, August 2013
c     Copyright Steve Morice.
c     Best-fit plane to a set of arbitary x,y,z points
c     using Gauss-Jordan inversion of solution matrix.

      call getarg(1,infile)

      if(infile.eq.'') then
       write(*,*)'Usage: bfp_analytic.exe [infile]'
       goto 99
      endif

      open(8,file=infile,status='old')
      n=1
10    read(8,*,end=11) x(n),y(n),z(n)
      n=n+1
      goto 10
11    n=n-1
      close(8)

      do 100 il=1,3
       mb(il,1)=0.
       do 110 jl=1,3
        ma(il,jl)=0.
        mm(il,jl)=0.
110    enddo
100   enddo

      do 120 nl=1,n
       ma(1,1)=ma(1,1)+(x(nl)*x(nl))
       ma(1,2)=ma(1,2)+(x(nl)*y(nl))
       ma(1,3)=ma(1,3)+(x(nl))
       ma(2,1)=ma(1,2)
       ma(2,2)=ma(2,2)+(y(nl)**2)
       ma(2,3)=ma(2,3)+(y(nl))
       ma(3,1)=ma(1,3)
       ma(3,2)=ma(2,3)
       ma(3,3)=real(n)
       mb(1,1)=mb(1,1)+(x(nl)*z(nl))
       mb(2,1)=mb(2,1)+(y(nl)*z(nl))
       mb(3,1)=mb(3,1)+(z(nl))
120   enddo

      call gaussj(ma,3,3,mm,3,3)

      mx(1,1)=(ma(1,1)*mb(1,1))+(ma(1,2)*mb(2,1))+(ma(1,3)*mb(3,1))
      mx(2,1)=(ma(2,1)*mb(1,1))+(ma(2,2)*mb(2,1))+(ma(2,3)*mb(3,1))
      mx(3,1)=(ma(3,1)*mb(1,1))+(ma(3,2)*mb(2,1))+(ma(3,3)*mb(3,1))

      open(unit=9,file='bfp_analytic_output.txt',status='unknown')

      write(*,*)'Solution z = Ax + By + C:'
      write(*,*)'A = ',mx(1,1)
      write(*,*)'B = ',mx(2,1)
      write(*,*)'C = ',mx(3,1)
      write(9,*)'Solution z = Ax + By + C:'
      write(9,*)'A = ',mx(1,1)
      write(9,*)'B = ',mx(2,1)
      write(9,*)'C = ',mx(3,1)

      Write(*,*)'Best fit plane to input points:'
      Write(9,*)'Best fit plane to input points:'
      do 200, wl=1,n
       write(*,*) wl,z(wl),(mx(1,1)*x(wl))+(mx(2,1)*y(wl))+mx(3,1)
       write(9,*) wl,z(wl),(mx(1,1)*x(wl))+(mx(2,1)*y(wl))+mx(3,1)
200   enddo

c     At z=0 and x=0...
      py=(-1.*mx(3,1))/mx(2,1)
c     At z=0 and y=0...
      px=(-1.*mx(3,1))/mx(1,1)
c     At x=0 and y=0...
      pz=mx(3,1)
      str=atan2((-1.*px),py)
      dipd=str+1.57079632679490
      r=py*sin(str)
      dip=3.14159265358979-atan2(pz,r)

c     Aaaarrgh! Quadrants!!!

      if(dip.ge.4.712388980384) then
        if(dipd.lt.3.14159265358979) then
         dipd=dipd+3.14159265358979
         dip=6.2831853071296-dip
        else
         dipd=dipd-3.14159265358979
         dip=6.2831853071296-dip
        endif
      elseif((dip.ge.3.14159265358979).and.(dip.lt.4.712388980384)) then
        dip=dip-3.14159265358979
      elseif((dip.ge.1.570796326795).and.(dip.lt.3.14159265358979)) then
       if(dipd.lt.3.14159265358979) then
         dipd=dipd+3.14159265358979
         dip=3.14159265358979-dip
        else
         dipd=dipd-3.14159265358979
         dip=3.14159265358979-dip
        endif
      endif
      if(dipd.lt.0.) then
       dipd=dipd+6.2831853071296
      elseif(dipd.gt.6.2831853071296) then
       dipd=dipd-6.2831853071296
      endif

      write(*,*)'Best fit strike:        ',str 
      write(*,*)'Best fit dip direction: ',dipd,dipd*57.29578
      write(*,*)'Best fit dip:           ',dip,dip*57.29578
      write(9,*)'Best fit strike:        ',str 
      write(9,*)'Best fit dip direction: ',dipd,dipd*57.29578
      write(9,*)'Best fit dip:           ',dip,dip*57.29578

      close(9)

c     Write plane grid

      open(unit=11,file='bfp_analytic_grid.csv',status='unknown')
      do 400 gxl=1697563.,1721464.,25.
       do 410 gyl=5677453.,5704254.,25.
        tzp=(gxl*mx(1,1))+(gyl*mx(2,1))+mx(3,1)
        write(11,*) gxl,',',gyl,',',tzp
410    enddo
400   enddo
      close(11)

99    end

c     **************************************************************
      SUBROUTINE gaussj(a,n,np,b,m,mp)  
 
      INTEGER m,mp,n,np,NMAX  
      REAL a(np,np),b(np,mp)  
      PARAMETER (NMAX=50)  
      INTEGER i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),ipiv(NMAX)  
      REAL big,dum,pivinv  
      do 11 j=1,n  
        ipiv(j)=0  
11    continue  
      do 22 i=1,n  
        big=0.  
        do 13 j=1,n  
          if(ipiv(j).ne.1)then  
            do 12 k=1,n  
              if (ipiv(k).eq.0) then  
                if (abs(a(j,k)).ge.big)then  
                  big=abs(a(j,k))  
                  irow=j  
                  icol=k  
                endif  
              else if (ipiv(k).gt.1) then  
                pause 'singular matrix in gaussj'  
              endif  
12          continue  
          endif  
13      continue  
        ipiv(icol)=ipiv(icol)+1  
        if (irow.ne.icol) then  
          do 14 l=1,n  
            dum=a(irow,l)  
            a(irow,l)=a(icol,l)  
            a(icol,l)=dum  
14        continue  
          do 15 l=1,m  
            dum=b(irow,l)  
            b(irow,l)=b(icol,l)  
            b(icol,l)=dum  
15        continue  
        endif  
        indxr(i)=irow  
        indxc(i)=icol  
        if (a(icol,icol).eq.0.) pause 'singular matrix in gaussj'  
        pivinv=1./a(icol,icol)  
        a(icol,icol)=1.  
        do 16 l=1,n  
          a(icol,l)=a(icol,l)*pivinv  
16      continue  
        do 17 l=1,m  
          b(icol,l)=b(icol,l)*pivinv  
17      continue  
        do 21 ll=1,n  
          if(ll.ne.icol)then  
            dum=a(ll,icol)  
            a(ll,icol)=0.  
            do 18 l=1,n  
              a(ll,l)=a(ll,l)-a(icol,l)*dum  
18          continue  
            do 19 l=1,m  
              b(ll,l)=b(ll,l)-b(icol,l)*dum  
19          continue  
          endif  
21      continue  
22    continue  
      do 24 l=n,1,-1  
        if(indxr(l).ne.indxc(l))then  
          do 23 k=1,n  
            dum=a(k,indxr(l))  
            a(k,indxr(l))=a(k,indxc(l))  
            a(k,indxc(l))=dum  
23        continue  
        endif  
24    continue  
      return  
      END  

  
 
