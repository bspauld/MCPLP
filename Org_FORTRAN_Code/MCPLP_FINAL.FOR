$LARGE
$PACK
      integer*2 ihr,imin,isec,i100,seedval
      real*4 ranval
      dimension ipop(60),icov(60,60),num(60),ib(60),int(60)
      open(unit=1,file='c:\newben\datainp3.txt',status='old')
      open(unit=2,file='c:\newben\output.txt',status='old')
      open(unit=3,file='c:\newben\diagno2.txt',status='old')
      call gettim(ihr,imin,isec,i100)
      seedval=4
      call seed(i100)  
      n=55
      do 1 i=1,n
      read(1,100) ipop(i),(icov(i,j),j=1,n)
  100 format(i3,55i4)
    1 continue
      do 500 kelf=1,100
      lp=6
      lp2=lp+1
      do 2 i=1,n
      num(i)=0
    2 continue   
      do 3 i=1,lp
   41 call random(ranval)
      xp=(ranval*n)+1.0
      mp=xp
      if(mp.eq.0) mp=1
      if(mp.gt.n) mp=n
      ib(i)=mp
      if(i.eq.1) go to 42
      iminus=i-1
      do 40 i4=1,iminus
      if(ib(i4).eq.ib(i)) go to 41
   40 continue
   42 ii=ib(i)
      do 4 j=1,n
      num(j)=num(j)+icov(ii,j)
    4 continue
    3 continue
      isum=0
      do 6 i=1,n
      if(num(i).eq.0) go to 6
      isum=isum+ipop(i)
    6 continue
      i1=lp+1
      do 7 i=i1,n
      ib(i)=i
    7 continue
      do 43 i=1,lp
      ii=ib(i)
      ib(ii)=i
   43 continue
      kp=3
      lp1=lp-2
      ksum=isum
      do 5 i=1,lp1
      int(1)=ib(i)
      i1=i+1
      lp3=lp-1
      do 8 j=i1,lp3
      int(2)=ib(j)
      i2=j+1
      do 44 jk=i2,lp
      int(3)=ib(jk)
C
C Start of coverage calculation after interdiction
C
      do 9 k=1,n
      num(k)=0
    9 continue
      do 12 k1=1,lp
      ii=ib(k1)
      if(ii.eq.int(1)) go to 12
      if(ii.eq.int(2)) go to 12
      if(ii.eq.int(3)) go to 12
      do 11 k2=1,n
      num(k2)=num(k2)+icov(ii,k2)
   11 continue
   12 continue
      jsum=0
      do 10 i3=1,n
      if(num(i3).eq.0) go to 10
      jsum=jsum+ipop(i3)
   10 continue
C
C End of calculating coverage after interdiction
C      
c      write(3,101) (ib(i3),i3=1,lp),int(1),int(2),int(3),jsum,ksum
      if(jsum.ge.ksum) go to 8
      ksum=jsum
      ks1=int(1)
      ks2=int(2)
      ks3=int(3)
   44 continue
    8 continue
    5 continue 
C
C Start the Teitz and Bart 1-OPT Swaps
C 
c      write(3,101) (ib(i),i=1,lp),ks1,ks2,ks3,ksum
   33 ic=0
      do 20 i=1,lp
      ip=ib(i)
      do 21 j=lp2,n
      jp=ib(j)
      ib(i)=jp
      ib(j)=ip
c      if(i.eq.1.and.j.eq.6) write(6,101) (ib(i3),i3=1,lp)
C
C Calculate new coverage value after potential swap
C
      do 22 k=1,n
      num(k)=0
   22 continue
      do 23 k=1,lp
      ii=ib(k)
      do 24 l=1,n
      num(l)=num(l)+icov(ii,l)
   24 continue
   23 continue
      isum=0
      do 25 k=1,n
      if(num(k).eq.0) go to 25
      isum=isum+ipop(k)
   25 continue
C
C Start of coverage calculation after interdiction
C
      kksum=isum
      do 26 k=1,lp1
      int(1)=ib(k)
      i1=k+1
      do 27 l=i1,lp3
      int(2)=ib(l)
      i2=l+1
      do 45 mo=i2,lp
      int(3)=ib(mo)
C
c      write(3,101) int(1),int(2),int(3),kksum
      do 28 m=1,n
      num(m)=0
   28 continue
      do 29 k1=1,lp
      ii=ib(k1)
      if(ii.eq.int(1)) go to 29
      if(ii.eq.int(2)) go to 29
      if(ii.eq.int(3)) go to 29
      do 30 k2=1,n
      num(k2)=num(k2)+icov(ii,k2)
   30 continue
   29 continue
      jsum=0
      do 31 k3=1,n
      if(num(k3).eq.0) go to 31
      jsum=jsum+ipop(k3)
   31 continue
C
c     write(3,101) (ib(i3),i3=1,lp),int(1),int(2),int(3),jsum,kksum,ksum
      if(jsum.ge.kksum) go to 27
c      write(3,101) kksum,jsum
      kksum=jsum
      kks1=int(1)
      kks2=int(2)
      kks3=int(3)
   45 continue
   27 continue
   26 continue
      if(kksum.le.ksum) go to 32
      ic=1
      ip=ib(i)
      jp=ib(j)
      it1=kks1
      it2=kks2
      it3=kks3
      ksum=kksum
c      write(3,101) (ib(i3),i3=1,lp),kks1,kks2,kks3,ksum
  101 format(15i5)
      go to 21
   32 continue
      ib(i)=ip
      ib(j)=jp
   21 continue
   20 continue
      if(ic.eq.1) go to 33
      write(3,101) (ib(i),i=1,lp),it1,it2,it3,ksum
  500 continue
      stop
      end
