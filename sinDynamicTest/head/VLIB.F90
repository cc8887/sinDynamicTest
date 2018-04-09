module vlib
contains                                         
 subroutine analy4(kma,coord,e,v)
!
! analytical stiffness of a 4-node quadrilateral element
! plane strain formulation
!
 implicit none
 real,intent(in):: coord(:,:),e,v
 real,intent(out):: kma(:,:)
 real:: e1,e2,g,x1,x2,x3,x4,y1,y2,y3,y4,a2,a2st3,alph,beta
 real:: c11,c21,c31,c41,c51,c61,c71,c81,c22,c32,c42,c52,c62,c72,c82
 real:: c33,c43,c53,c63,c73,c83,c44,c54,c64,c74,c84,c55,c65,c75,c85
 real:: c66,c76,c86,c77,c87,c88,f1,f2,s1,s2,s3,s4,t1,t2,t3,t4
 integer:: i,j
!
 e1=e*(1-v)/(1+v)/(1-2*v)
 e2=v*e1/(1-v)
 g=e/2/(1+v)      

 x1=coord(1,1)
 x2=coord(2,1)
 x3=coord(3,1)
 x4=coord(4,1)
 y1=coord(1,2)
 y2=coord(2,2)
 y3=coord(3,2)
 y4=coord(4,2)

 a2=(x4-x2)*(y3-y1)-(x3-x1)*(y4-y2)
 a2st3=3*a2*a2

 call groupa(x1,x2,x3,x4,y1,y2,y3,y4,&
             s1,s2,s3,s4,t1,t2,t3,t4,f1,f2)
 alph=a2*(s1*e1+s2*g)+f1*(s3*e1+s4*g)
 beta=a2*(t1*e1+t2*g)+f2*(t3*e1+t4*g)
 c11=(alph/(a2st3-f1**2)+beta/(a2st3-f2**2))*0.5
 kma(1,1)=c11

 call groupa(x2,x3,x4,x1,y2,y3,y4,y1,&
             s1,s2,s3,s4,t1,t2,t3,t4,f1,f2)
 alph=a2*(s1*e1+s2*g)+f1*(s3*e1+s4*g)
 beta=a2*(t1*e1+t2*g)+f2*(t3*e1+t4*g)
 c33=(alph/(a2st3-f1**2)+beta/(a2st3-f2**2))*0.5
 kma(3,3)=c33

 call groupa(x3,x4,x1,x2,y3,y4,y1,y2,&
             s1,s2,s3,s4,t1,t2,t3,t4,f1,f2)
 alph=a2*(s1*e1+s2*g)+f1*(s3*e1+s4*g)
 beta=a2*(t1*e1+t2*g)+f2*(t3*e1+t4*g)
 c55=(alph/(a2st3-f1**2)+beta/(a2st3-f2**2))*0.5
 kma(5,5)=c55

 call groupa(x4,x1,x2,x3,y4,y1,y2,y3,&
             s1,s2,s3,s4,t1,t2,t3,t4,f1,f2)
 alph=a2*(s1*e1+s2*g)+f1*(s3*e1+s4*g)
 beta=a2*(t1*e1+t2*g)+f2*(t3*e1+t4*g)
 c77=(alph/(a2st3-f1**2)+beta/(a2st3-f2**2))*0.5
 kma(7,7)=c77

 call groupa(y4,y3,y2,y1,x4,x3,x2,x1,&
             s1,s2,s3,s4,t1,t2,t3,t4,f1,f2)
 alph=a2*(s1*e1+s2*g)+f1*(s3*e1+s4*g)
 beta=a2*(t1*e1+t2*g)+f2*(t3*e1+t4*g)
 c88=(alph/(a2st3-f1**2)+beta/(a2st3-f2**2))*0.5
 kma(8,8)=c88

 call groupa(y1,y4,y3,y2,x1,x4,x3,x2,&
             s1,s2,s3,s4,t1,t2,t3,t4,f1,f2)
 alph=a2*(s1*e1+s2*g)+f1*(s3*e1+s4*g)
 beta=a2*(t1*e1+t2*g)+f2*(t3*e1+t4*g)
 c22=(alph/(a2st3-f1**2)+beta/(a2st3-f2**2))*0.5
 kma(2,2)=c22

 call groupa(y2,y1,y4,y3,x2,x1,x4,x3,&
             s1,s2,s3,s4,t1,t2,t3,t4,f1,f2)
 alph=a2*(s1*e1+s2*g)+f1*(s3*e1+s4*g)
 beta=a2*(t1*e1+t2*g)+f2*(t3*e1+t4*g)
 c44=(alph/(a2st3-f1**2)+beta/(a2st3-f2**2))*0.5
 kma(4,4)=c44

 call groupa(y3,y2,y1,y4,x3,x2,x1,x4,&
             s1,s2,s3,s4,t1,t2,t3,t4,f1,f2)
 alph=a2*(s1*e1+s2*g)+f1*(s3*e1+s4*g)
 beta=a2*(t1*e1+t2*g)+f2*(t3*e1+t4*g)
 c66=(alph/(a2st3-f1**2)+beta/(a2st3-f2**2))*0.5
 kma(6,6)=c66
 
 call groupb(x1,x2,x3,x4,y1,y2,y3,y4,&
             s1,s2,s3,s4,t1,t2,t3,t4,f1,f2)
 alph=a2*(s1*e2+s2*g)+f1*(s3*e2+s4*g)
 beta=a2*(t1*e2+t2*g)+f2*(t3*e2+t4*g)
 c21=(alph/(a2st3-f1**2)+beta/(a2st3-f2**2))*0.5
 kma(2,1)=c21

 call groupb(x2,x3,x4,x1,y2,y3,y4,y1,&
             s1,s2,s3,s4,t1,t2,t3,t4,f1,f2)
 alph=a2*(s1*e2+s2*g)+f1*(s3*e2+s4*g)
 beta=a2*(t1*e2+t2*g)+f2*(t3*e2+t4*g)
 c43=(alph/(a2st3-f1**2)+beta/(a2st3-f2**2))*0.5
 kma(4,3)=c43

 call groupb(x3,x4,x1,x2,y3,y4,y1,y2,&
             s1,s2,s3,s4,t1,t2,t3,t4,f1,f2)
 alph=a2*(s1*e2+s2*g)+f1*(s3*e2+s4*g)
 beta=a2*(t1*e2+t2*g)+f2*(t3*e2+t4*g)
 c65=(alph/(a2st3-f1**2)+beta/(a2st3-f2**2))*0.5
 kma(6,5)=c65

 call groupb(x4,x1,x2,x3,y4,y1,y2,y3,&
             s1,s2,s3,s4,t1,t2,t3,t4,f1,f2)
 alph=a2*(s1*e2+s2*g)+f1*(s3*e2+s4*g)
 beta=a2*(t1*e2+t2*g)+f2*(t3*e2+t4*g)
 c87=(alph/(a2st3-f1**2)+beta/(a2st3-f2**2))*0.5
 kma(8,7)=c87

 call groupc(x1,x2,x3,x4,y1,y2,y3,y4,&
             s1,s2,s3,s4,t1,t2,t3,t4,f1,f2)
 alph=a2*(s1*e1+s2*g)+f1*(s3*e1+s4*g)
 beta=a2*(t1*e1+t2*g)+f2*(t3*e1+t4*g)
 c31=(alph/(a2st3-f1**2)+beta/(a2st3-f2**2))*0.5
 kma(3,1)=c31

 call groupc(x2,x3,x4,x1,y2,y3,y4,y1,&
             s1,s2,s3,s4,t1,t2,t3,t4,f1,f2)
 alph=a2*(s1*e1+s2*g)+f1*(s3*e1+s4*g)
 beta=a2*(t1*e1+t2*g)+f2*(t3*e1+t4*g)
 c53=(alph/(a2st3-f1**2)+beta/(a2st3-f2**2))*0.5
 kma(5,3)=c53

 call groupc(x3,x4,x1,x2,y3,y4,y1,y2,&
             s1,s2,s3,s4,t1,t2,t3,t4,f1,f2)
 alph=a2*(s1*e1+s2*g)+f1*(s3*e1+s4*g)
 beta=a2*(t1*e1+t2*g)+f2*(t3*e1+t4*g)
 c75=(alph/(a2st3-f1**2)+beta/(a2st3-f2**2))*0.5
 kma(7,5)=c75

 call groupc(x4,x1,x2,x3,y4,y1,y2,y3,&
             s1,s2,s3,s4,t1,t2,t3,t4,f1,f2)
 alph=a2*(s1*e1+s2*g)+f1*(s3*e1+s4*g)
 beta=a2*(t1*e1+t2*g)+f2*(t3*e1+t4*g)
 c71=(alph/(a2st3-f1**2)+beta/(a2st3-f2**2))*0.5
 kma(7,1)=c71
 
 call groupc(y4,y3,y2,y1,x4,x3,x2,x1,&
             s1,s2,s3,s4,t1,t2,t3,t4,f1,f2)
 alph=a2*(s1*e1+s2*g)+f1*(s3*e1+s4*g)
 beta=a2*(t1*e1+t2*g)+f2*(t3*e1+t4*g)
 c86=(alph/(a2st3-f1**2)+beta/(a2st3-f2**2))*0.5
 kma(8,6)=c86

 call groupc(y1,y4,y3,y2,x1,x4,x3,x2,&
             s1,s2,s3,s4,t1,t2,t3,t4,f1,f2)
 alph=a2*(s1*e1+s2*g)+f1*(s3*e1+s4*g)
 beta=a2*(t1*e1+t2*g)+f2*(t3*e1+t4*g)
 c82=(alph/(a2st3-f1**2)+beta/(a2st3-f2**2))*0.5
 kma(8,2)=c82

 call groupc(y2,y1,y4,y3,x2,x1,x4,x3,&
             s1,s2,s3,s4,t1,t2,t3,t4,f1,f2)
 alph=a2*(s1*e1+s2*g)+f1*(s3*e1+s4*g)
 beta=a2*(t1*e1+t2*g)+f2*(t3*e1+t4*g)
 c42=(alph/(a2st3-f1**2)+beta/(a2st3-f2**2))*0.5
 kma(4,2)=c42

 call groupc(y3,y2,y1,y4,x3,x2,x1,x4,&
             s1,s2,s3,s4,t1,t2,t3,t4,f1,f2)
 alph=a2*(s1*e1+s2*g)+f1*(s3*e1+s4*g)
 beta=a2*(t1*e1+t2*g)+f2*(t3*e1+t4*g)
 c64=(alph/(a2st3-f1**2)+beta/(a2st3-f2**2))*0.5
 kma(6,4)=c64

 call groupd(x1,x2,x3,x4,y1,y2,y3,y4,&
             s1,s2,s3,s4,t1,t2,t3,t4)
 call f1f2(x1,x2,x3,x4,y1,y2,y3,y4,f1,f2)
 alph=a2*(s1*e2+s2*g)+f1*(s3*e2+s4*g)
 beta=a2*(t1*e2+t2*g)+f2*(t3*e2+t4*g)
 c41=(alph/(a2st3-f1**2)+beta/(a2st3-f2**2))*0.5
 kma(4,1)=c41

 call groupd(x2,x3,x4,x1,y2,y3,y4,y1,&
             s1,s2,s3,s4,t1,t2,t3,t4)
 call f1f2(x2,x3,x4,x1,y2,y3,y4,y1,f1,f2)
 alph=a2*(s1*e2+s2*g)+f1*(s3*e2+s4*g)
 beta=a2*(t1*e2+t2*g)+f2*(t3*e2+t4*g)
 c63=(alph/(a2st3-f1**2)+beta/(a2st3-f2**2))*0.5
 kma(6,3)=c63

 call groupd(x3,x4,x1,x2,y3,y4,y1,y2,&
             s1,s2,s3,s4,t1,t2,t3,t4)
 call f1f2(x3,x4,x1,x2,y3,y4,y1,y2,f1,f2)
 alph=a2*(s1*e2+s2*g)+f1*(s3*e2+s4*g)
 beta=a2*(t1*e2+t2*g)+f2*(t3*e2+t4*g)
 c85=(alph/(a2st3-f1**2)+beta/(a2st3-f2**2))*0.5
 kma(8,5)=c85

 call groupd(x4,x1,x2,x3,y4,y1,y2,y3,&
             s1,s2,s3,s4,t1,t2,t3,t4)
 call f1f2(x4,x1,x2,x3,y4,y1,y2,y3,f1,f2)
 alph=a2*(s1*e2+s2*g)+f1*(s3*e2+s4*g)
 beta=a2*(t1*e2+t2*g)+f2*(t3*e2+t4*g)
 c72=(alph/(a2st3-f1**2)+beta/(a2st3-f2**2))*0.5
 kma(7,2)=c72

 call groupd(y1,y2,y3,y4,x1,x2,x3,x4,&
             s1,s2,s3,s4,t1,t2,t3,t4)
 call f1f2(x1,x2,x3,x4,y1,y2,y3,y4,f1,f2)
 alph=a2*(s1*e2+s2*g)+f1*(s3*e2+s4*g)
 beta=a2*(t1*e2+t2*g)+f2*(t3*e2+t4*g)
 c32=(alph/(a2st3-f1**2)+beta/(a2st3-f2**2))*0.5
 kma(3,2)=c32

 call groupd(y2,y3,y4,y1,x2,x3,x4,x1,&
             s1,s2,s3,s4,t1,t2,t3,t4)
 call f1f2(x2,x3,x4,x1,y2,y3,y4,y1,f1,f2)
 alph=a2*(s1*e2+s2*g)+f1*(s3*e2+s4*g)
 beta=a2*(t1*e2+t2*g)+f2*(t3*e2+t4*g)
 c54=(alph/(a2st3-f1**2)+beta/(a2st3-f2**2))*0.5
 kma(5,4)=c54

 call groupd(y3,y4,y1,y2,x3,x4,x1,x2,&
             s1,s2,s3,s4,t1,t2,t3,t4)
 call f1f2(x3,x4,x1,x2,y3,y4,y1,y2,f1,f2)
 alph=a2*(s1*e2+s2*g)+f1*(s3*e2+s4*g)
 beta=a2*(t1*e2+t2*g)+f2*(t3*e2+t4*g)
 c76=(alph/(a2st3-f1**2)+beta/(a2st3-f2**2))*0.5
 kma(7,6)=c76

 call groupd(y4,y1,y2,y3,x4,x1,x2,x3,&
             s1,s2,s3,s4,t1,t2,t3,t4)
 call f1f2(x4,x1,x2,x3,y4,y1,y2,y3,f1,f2)
 alph=a2*(s1*e2+s2*g)+f1*(s3*e2+s4*g)
 beta=a2*(t1*e2+t2*g)+f2*(t3*e2+t4*g)
 c81=(alph/(a2st3-f1**2)+beta/(a2st3-f2**2))*0.5
 kma(8,1)=c81

 call groupe(x1,x2,x3,x4,y1,y2,y3,y4,&
             s1,s2,s3,s4,t1,t2,t3,t4,f1,f2)
 alph=a2*(s1*e1+s2*g)+f1*(s3*e1+s4*g)
 beta=a2*(t1*e1+t2*g)+f2*(t3*e1+t4*g)
 c51=(alph/(a2st3-f1**2)+beta/(a2st3-f2**2))*0.5
 kma(5,1)=c51

 call groupe(x2,x3,x4,x1,y2,y3,y4,y1,&
             s1,s2,s3,s4,t1,t2,t3,t4,f1,f2)
 alph=a2*(s1*e1+s2*g)+f1*(s3*e1+s4*g)
 beta=a2*(t1*e1+t2*g)+f2*(t3*e1+t4*g)
 c73=(alph/(a2st3-f1**2)+beta/(a2st3-f2**2))*0.5
 kma(7,3)=c73

 call groupe(y2,y1,y4,y3,x2,x1,x4,x3,&
             s1,s2,s3,s4,t1,t2,t3,t4,f1,f2)
 alph=a2*(s1*e1+s2*g)+f1*(s3*e1+s4*g)
 beta=a2*(t1*e1+t2*g)+f2*(t3*e1+t4*g)
 c84=(alph/(a2st3-f1**2)+beta/(a2st3-f2**2))*0.5
 kma(8,4)=c84

 call groupe(y3,y2,y1,y4,x3,x2,x1,x4,&
             s1,s2,s3,s4,t1,t2,t3,t4,f1,f2)
 alph=a2*(s1*e1+s2*g)+f1*(s3*e1+s4*g)
 beta=a2*(t1*e1+t2*g)+f2*(t3*e1+t4*g)
 c62=(alph/(a2st3-f1**2)+beta/(a2st3-f2**2))*0.5
 kma(6,2)=c62

 call groupf(x1,x2,x3,x4,y1,y2,y3,y4,&
             s1,s2,s3,s4,t1,t2,t3,t4)
 call f1f2(x1,x2,x3,x4,y1,y2,y3,y4,f1,f2)
 alph=a2*(s1*e2+s2*g)+f1*(s3*e2+s4*g)
 beta=a2*(t1*e2+t2*g)+f2*(t3*e2+t4*g)
 c61=(alph/(a2st3-f1**2)+beta/(a2st3-f2**2))*0.5
 kma(6,1)=c61

 call groupf(x2,x3,x4,x1,y2,y3,y4,y1,&
             s1,s2,s3,s4,t1,t2,t3,t4)
 call f1f2(x2,x3,x4,x1,y2,y3,y4,y1,f1,f2)
 alph=a2*(s1*e2+s2*g)+f1*(s3*e2+s4*g)
 beta=a2*(t1*e2+t2*g)+f2*(t3*e2+t4*g)
 c83=(alph/(a2st3-f1**2)+beta/(a2st3-f2**2))*0.5
 kma(8,3)=c83

 call groupf(y1,y2,y3,y4,x1,x2,x3,x4,&
             s1,s2,s3,s4,t1,t2,t3,t4)
 call f1f2(x1,x2,x3,x4,y1,y2,y3,y4,f1,f2)
 alph=a2*(s1*e2+s2*g)+f1*(s3*e2+s4*g)
 beta=a2*(t1*e2+t2*g)+f2*(t3*e2+t4*g)
 c52=(alph/(a2st3-f1**2)+beta/(a2st3-f2**2))*0.5
 kma(5,2)=c52

 call groupf(y2,y3,y4,y1,x2,x3,x4,x1,&
             s1,s2,s3,s4,t1,t2,t3,t4)
 call f1f2(x2,x3,x4,x1,y2,y3,y4,y1,f1,f2)
 alph=a2*(s1*e2+s2*g)+f1*(s3*e2+s4*g)
 beta=a2*(t1*e2+t2*g)+f2*(t3*e2+t4*g)
 c74=(alph/(a2st3-f1**2)+beta/(a2st3-f2**2))*0.5
 kma(7,4)=c74
 do i=1,8; do j=i+1,8; kma(i,j)=kma(j,i); end do; end do

 return 
 end  subroutine analy4

 subroutine groupa(x1,x2,x3,x4,y1,y2,y3,y4,&
                   s1,s2,s3,s4,t1,t2,t3,t4,f1,f2)
 implicit none
 real,intent(in)::x1,x2,x3,x4,y1,y2,y3,y4
 real,intent(out)::s1,s2,s3,s4,t1,t2,t3,t4,f1,f2
 s1=2*(y4-y2)**2
 s2=2*(x4-x2)**2
 s3=-s1/2
 s4=-s2/2
 t1=(y2-y3)**2+(y3-y4)**2+(y4-y2)**2
 t2=(x2-x3)**2+(x3-x4)**2+(x4-x2)**2
 t3=(y4-y3)**2-(y3-y2)**2
 t4=(x4-x3)**2-(x3-x2)**2
 f1=(x1+x3)*(y4-y2)-(y1+y3)*(x4-x2)-2*(x2*y4-x4*y2)
 f2=(y2+y4)*(x3-x1)-(x2+x4)*(y3-y1)-2*(x3*y1-x1*y3)
 return
 end subroutine groupa

 subroutine groupb(x1,x2,x3,x4,y1,y2,y3,y4,&
                   s1,s2,s3,s4,t1,t2,t3,t4,f1,f2)
 implicit none
 real,intent(in)::x1,x2,x3,x4,y1,y2,y3,y4
 real,intent(out)::s1,s2,s3,s4,t1,t2,t3,t4,f1,f2
 s1=2*(x2-x4)*(y4-y2)
 s2=s1
 s3=-s1/2
 s4=s3
 t1=x2*(y4-2*y2+y3)+x3*(y2-2*y3+y4)+x4*(y2-2*y4+y3)
 t2=t1
 t3=x2*(y2-y3)+x3*(y4-y2)+x4*(y3-y4)
 t4=t3
 f1=(x1+x3)*(y4-y2)-(y1+y3)*(x4-x2)-2*(x2*y4-x4*y2)
 f2=(y2+y4)*(x3-x1)-(x2+x4)*(y3-y1)-2*(x3*y1-x1*y3)
 return
 end subroutine groupb

 subroutine groupc(x1,x2,x3,x4,y1,y2,y3,y4,&
                   s1,s2,s3,s4,t1,t2,t3,t4,f1,f2)
 implicit none
 real,intent(in)::x1,x2,x3,x4,y1,y2,y3,y4
 real,intent(out)::s1,s2,s3,s4,t1,t2,t3,t4,f1,f2
 s1=(y4-y2)*(2*y1-y3-y4)
 s2=(x4-x2)*(2*x1-x3-x4)
 s3=(y4-y2)*(y4-y1)
 s4=(x4-x2)*(x4-x1)
 t1=(y3-y1)*(2*y2-y3-y4)
 t2=(x3-x1)*(2*x2-x3-x4)
 t3=(y3-y1)*(y3-y2)
 t4=(x3-x1)*(x3-x2)
 f1=(x1+x3)*(y4-y2)-(y1+y3)*(x4-x2)-2*(x2*y4-x4*y2)
 f2=(y2+y4)*(x3-x1)-(x2+x4)*(y3-y1)-2*(x3*y1-x1*y3)
 return
 end subroutine groupc

 subroutine groupd(x1,x2,x3,x4,y1,y2,y3,y4,&
                    s1,s2,s3,s4,t1,t2,t3,t4)
 implicit none
 real,intent(in)::x1,x2,x3,x4,y1,y2,y3,y4
 real,intent(out)::s1,s2,s3,s4,t1,t2,t3,t4
 s1=(x3-x1)*(y4-y2)+(x4-x1)*(y4-y2)
 s2=(y3-y1)*(x4-x2)+(y4-y1)*(x4-x2)
 s3=(x4-x1)*(y2-y4)
 s4=(y4-y1)*(x2-x4)
 t1=(x3-x1)*(y4-y2)+(x3-x1)*(y3-y2)
 t2=(y3-y1)*(x4-x2)+(y3-y1)*(x3-x2)
 t3=(x3-x1)*(y2-y3)
 t4=(y3-y1)*(x2-x3)
 return
 end subroutine groupd

 subroutine groupe(x1,x2,x3,x4,y1,y2,y3,y4,&
                   s1,s2,s3,s4,t1,t2,t3,t4,f1,f2)
 implicit none
 real,intent(in)::x1,x2,x3,x4,y1,y2,y3,y4
 real,intent(out)::s1,s2,s3,s4,t1,t2,t3,t4,f1,f2
 s1=-(y4-y2)**2
 s2=-(x4-x2)**2
 s3=0
 s4=0
 t1=(y3+y1)*(y4+y2)-2*(y4-y2)**2-2*(y1*y3+y2*y4)
 t2=(x3+x1)*(x4+x2)-2*(x4-x2)**2-2*(x1*x3+x2*x4)
 t3=(y4-y2)*(y1-y2+y3-y4)
 t4=(x4-x2)*(x1-x2+x3-x4)
 f1=(x1+x3)*(y4-y2)-(y1+y3)*(x4-x2)-2*(x2*y4-x4*y2)
 f2=(y2+y4)*(x3-x1)-(x2+x4)*(y3-y1)-2*(x3*y1-x1*y3)
 return
 end subroutine groupe

 subroutine groupf(x1,x2,x3,x4,y1,y2,y3,y4,&
                   s1,s2,s3,s4,t1,t2,t3,t4)
 implicit none
 real,intent(in)::x1,x2,x3,x4,y1,y2,y3,y4
 real,intent(out)::s1,s2,s3,s4,t1,t2,t3,t4
 s1=(x4-x2)*(y4-y2)
 s2=s1
 s3=0
 s4=0
 t1=(x4-x2)*(y4-y2)+(x2-x1)*(y2-y3)+(x4-x1)*(y4-y3)
 t2=(y4-y2)*(x4-x2)+(y2-y1)*(x2-x3)+(y4-y1)*(x4-x3)
 t3=(x2-x1)*(y3-y2)+(x4-x1)*(y4-y3)
 t4=(y2-y1)*(x3-x2)+(y4-y1)*(x4-x3)
 return
 end subroutine groupf

 subroutine f1f2(x1,x2,x3,x4,y1,y2,y3,y4,f1,f2)
 implicit none
 real,intent(in)::x1,x2,x3,x4,y1,y2,y3,y4
 real,intent(out)::f1,f2
 f1=(x1+x3)*(y4-y2)-(y1+y3)*(x4-x2)-2*(x2*y4-x4*y2)
 f2=(y2+y4)*(x3-x1)-(x2+x4)*(y3-y1)-2*(x3*y1-x1*y3)
 return
 end subroutine f1f2
subroutine bee4(coord,points,i,det,bee)

!  analytical version of the bee matrix for a 4-node quad      

 implicit none
 real,intent(in):: coord(:,:),points(:,:)
 real,intent(out):: det,bee(:,:)
 integer,intent(in):: i
 real:: x1,x2,x3,x4,y1,y2,y3,y4,xi,et,&
        x12,x13,x14,x23,x24,x34,y12,y13,y14,y23,y24,y34,&
        xy12,xy13,xy14,xy23,xy24,xy34,xifac,etfac,const,den,&
        dn1dx,dn2dx,dn3dx,dn4dx,dn1dy,dn2dy,dn3dy,dn4dy  

 x1=coord(1,1); x2=coord(2,1); x3=coord(3,1); x4=coord(4,1)
 y1=coord(1,2); y2=coord(2,2); y3=coord(3,2); y4=coord(4,2)
 xi=points(i,1); et=points(i,2)
 x12=x1-x2; x13=x1-x3; x14=x1-x4; x23=x2-x3; x24=x2-x4; x34=x3-x4
 y12=y1-y2; y13=y1-y3; y14=y1-y4; y23=y2-y3; y24=y2-y4; y34=y3-y4
 xy12=x1*y2-y1*x2; xy13=x1*y3-y1*x3; xy14=x1*y4-y1*x4; xy23=x2*y3-y2*x3
 xy24=x2*y4-y2*x4; xy34=x3*y4-y3*x4
 xifac= xy12-xy13+xy24-xy34; etfac= xy13-xy14-xy23+xy24
 const=-xy12+xy14-xy23-xy34; den=xifac*xi+etfac*et+const; det=0.125*den
 dn1dx=( y23*xi+y34*et-y24)/den; dn2dx=(-y14*xi-y34*et+y13)/den
 dn3dx=( y14*xi-y12*et+y24)/den; dn4dx=(-y23*xi+y12*et-y13)/den
 dn1dy=(-x23*xi-x34*et+x24)/den; dn2dy=( x14*xi+x34*et-x13)/den
 dn3dy=(-x14*xi+x12*et-x24)/den; dn4dy=( x23*xi-x12*et+x13)/den
 bee(1,1)=dn1dx; bee(1,2)=0.; bee(1,3)=dn2dx; bee(1,4)=0.
 bee(1,5)=dn3dx; bee(1,6)=0.; bee(1,7)=dn4dx; bee(1,8)=0.
 bee(2,1)=0.; bee(2,2)=dn1dy; bee(2,3)=0.; bee(2,4)=dn2dy
 bee(2,5)=0.; bee(2,6)=dn3dy; bee(2,7)=0.; bee(2,8)=dn4dy
 bee(3,1)=dn1dy; bee(3,2)=dn1dx; bee(3,3)=dn2dy; bee(3,4)=dn2dx
 bee(3,5)=dn3dy; bee(3,6)=dn3dx; bee(3,7)=dn4dy; bee(3,8)=dn4dx
return
end subroutine bee4           
subroutine seep4(coord,perm,kc)

!     this contains the analytical conductivity mstrix
!     for a four node quadrilateral

 implicit none
 real,intent(in):: coord(:,:),perm(:,:)
 real,intent(out):: kc(:,:)
 integer::i,j
 real:: x1,x2,x3,x4,y1,y2,y3,y4,y31,y42,x31,x42,x34,y34
 real:: a2,f1,f2,kx,ky,ht(4),it(4),temp,s1,s2,t1,t2,t3,gdiag,goppo
 
! input data
 
 x1=coord(1,1); x2=coord(2,1); x3=coord(3,1); x4=coord(4,1)
 y1=coord(1,2); y2=coord(2,2); y3=coord(3,2); y4=coord(4,2)
 kx=perm(1,1); ky=perm(2,2)
 
! procedure definitions
 
 y31=(y3-y1); y42=(y4-y2); x31=(x3-x1); x42=(x4-x2)
 x34=x3-x4; y34=y3-y4; ht(1)=x4-x1; ht(2)=x1-x2
 ht(3)=x2-x3; it(1)=y4-y1; it(2)=y1-y2; it(3)=y2-y3
 
 a2=y42*x31-y31*x42
 f1=2*y42*x34-2*x42*y34-a2; f2=2*y31*x34-2*x31*y34-a2
 s1=kx*y42**2+ky*x42**2; s2=kx*y34**2+ky*x34**2
 t1=kx*y42*y34+ky*x42*x34; t2=kx*y31*y34+ky*x31*x34; t3=kx*y42*y31+x31*x42*ky
 
! calculate parent terms
 
 gdiag=(-(2.*a2+f1)*s1)/(2.*(3.*a2**2-f1**2))
 goppo=-(s1*(2.*a2+f2)+2.*t1*(a2+f2)+2.*a2*s2)/(2.*(3.*a2**2-f2**2))
 kc(1,1)=gdiag+goppo
!
 gdiag=((2.*a2+f1)*t3-(a2+f1)*t1)/(2.*(3.*a2**2-f1**2))
 goppo=((2.*a2+f2)*t3+(a2+f2)*t2)/(2.*(3.*a2**2-f2**2))
 kc(2,1)=gdiag+goppo
!
 gdiag=(s1*a2)/(2.*(3.*a2**2-f1**2))
 goppo=((2.*a2+f2)*s1+(a2+f2)*(2.*t1-t3)+2.*a2*(s2-t2))/(2.*(3.*a2**2-f2**2))
 kc(3,1)=gdiag+goppo
 
! calculate elements using parent terms and transforms
 
 temp=y31; y31=y42; y42=-temp; temp=x31; 
 x31=x42; x42=-temp; x34=ht(1); y34=it(1)
!
 f1=2*y42*x34-2*x42*y34-a2; f2=2*y31*x34-2*x31*y34-a2
 s1=kx*y42**2+ky*x42**2; s2=kx*y34**2+ky*x34**2
 t1=kx*y42*y34+ky*x42*x34; t2=kx*y31*y34+ky*x31*x34; t3=kx*y42*y31+x31*x42*ky
!
 gdiag=(-(2.*a2+f1)*s1)/(2.*(3.*a2**2-f1**2))
 goppo=-(s1*(2.*a2+f2)+2.*t1*(a2+f2)+2.*a2*s2)/(2.*(3.*a2**2-f2**2))
 kc(2,2)=gdiag+goppo
!
 gdiag=((2.*a2+f1)*t3-(a2+f1)*t1)/(2.*(3.*a2**2-f1**2))
 goppo=((2.*a2+f2)*t3+(a2+f2)*t2)/(2.*(3.*a2**2-f2**2))
 kc(3,2)=gdiag+goppo
!
 gdiag=(s1*a2)/(2.*(3.*a2**2-f1**2))
 goppo=((2.*a2+f2)*s1+(a2+f2)*(2.*t1-t3)+2.*a2*(s2-t2))/(2.*(3.*a2**2-f2**2))
 kc(4,2)=gdiag+goppo
!
 temp=y31; y31=y42; y42=-temp; temp=x31
 x31=x42; x42=-temp; x34=ht(2); y34=it(2)
!
 f1=2*y42*x34-2*x42*y34-a2; f2=2*y31*x34-2*x31*y34-a2
 s1=kx*y42**2+ky*x42**2; s2=kx*y34**2+ky*x34**2
 t1=kx*y42*y34+ky*x42*x34; t2=kx*y31*y34+ky*x31*x34; t3=kx*y42*y31+x31*x42*ky
!
 gdiag=(-(2.*a2+f1)*s1)/(2.*(3.*a2**2-f1**2))
 goppo=-(s1*(2.*a2+f2)+2.*t1*(a2+f2)+2.*a2*s2)/(2.*(3.*a2**2-f2**2))
 kc(3,3)=gdiag+goppo
!
 gdiag=((2.*a2+f1)*t3-(a2+f1)*t1)/(2.*(3.*a2**2-f1**2))
 goppo=((2.*a2+f2)*t3+(a2+f2)*t2)/(2.*(3.*a2**2-f2**2))
 kc(4,3)=gdiag+goppo
!
 temp=y31; y31=y42; y42=-temp; temp=x31
 x31=x42; x42=-temp; x34=ht(3); y34=it(3)
!
 f1=2*y42*x34-2*x42*y34-a2; f2=2*y31*x34-2*x31*y34-a2
 s1=kx*y42**2+ky*x42**2; s2=kx*y34**2+ky*x34**2
 t1=kx*y42*y34+ky*x42*x34; t2=kx*y31*y34+ky*x31*x34; t3=kx*y42*y31+x31*x42*ky
!
 gdiag=(-(2.*a2+f1)*s1)/(2.*(3.*a2**2-f1**2))
 goppo=-(s1*(2.*a2+f2)+2.*t1*(a2+f2)+2.*a2*s2)/(2.*(3.*a2**2-f2**2))
 kc(4,4)=gdiag+goppo
!
 gdiag=((2.*a2+f1)*t3-(a2+f1)*t1)/(2.*(3.*a2**2-f1**2))
 goppo=((2.*a2+f2)*t3+(a2+f2)*t2)/(2.*(3.*a2**2-f2**2))
 kc(4,1)=gdiag+goppo
 
! mirror matrix about main diagonal
 
 do i=2,4; do j=1,i-1; kc(j,i)=kc(i,j); enddo; enddo

 return
 end subroutine seep4         
subroutine loc_to_glob(local,global,gamma,coord)
!
!      this subroutine transforms the local end reactions and
!      moments into the global system (3-d)
!
implicit none
real,intent(in)::local(:),gamma,coord(:,:)
real,intent(out)::global(:)
real::t(12,12),r0(3,3),x1,x2,y1,y2,z1,z2,xl,yl,zl,&
      pi,gamrad,cg,sg,den,ell,x,sum
integer::i,j,k
     x1=coord(1,1); y1=coord(1,2); z1=coord(1,3)
     x2=coord(2,1); y2=coord(2,2); z2=coord(2,3)
     xl=x2-x1; yl=y2-y1; zl=z2-z1; ell=sqrt(xl*xl+yl*yl+zl*zl)
     t=0.0
     pi=acos(-1.); gamrad=gamma*pi/180.; cg=cos(gamrad); sg=sin(gamrad)
     den=ell*sqrt(xl*xl+zl*zl)
     if(den /= 0.0)then
      r0(1,1)=xl/ell; r0(2,1)=yl/ell; r0(3,1)=zl/ell
      r0(1,2)=(-xl*yl*cg-ell*zl*sg)/den; r0(2,2)=den*cg/(ell*ell)
      r0(3,2)=(-yl*zl*cg+ell*xl*sg)/den; r0(1,3)=(xl*yl*sg-ell*zl*cg)/den
      r0(2,3)=-den*sg/(ell*ell); r0(3,3)=(yl*zl*sg+ell*xl*cg)/den
     else
      r0(1,1)=0.; r0(3,1)=0.; r0(2,2)=0.; r0(2,3)=0.
      r0(2,1)=1.; r0(1,2)=-cg; r0(3,3)=cg; r0(3,2)=sg; r0(1,3)=sg
     end if
     do i=1,3; do j=1,3 
      x=r0(i,j)
      do k=0,9,3; t(i+k,j+k)=x; end do
     end do; end do
     do i=1,12
      sum=0.
       do j=1,12
        sum=sum+t(i,j)*local(j)
       end do
      global(i)=sum
     end do
      return
end subroutine loc_to_glob    
subroutine glob_to_loc(local,global,gamma,coord)
!
!      this subroutine transforms the global end reactions and
!      moments into the local system (2-d, 3-d)
!
implicit none
real,intent(in)::global(:),gamma,coord(:,:)
real,intent(out)::local(:)
real::t(12,12),r0(3,3),x1,x2,y1,y2,z1,z2,xl,yl,zl,&
      pi,gamrad,cg,sg,den,ell,x,sum
integer::i,j,k,ndim
ndim=ubound(coord,2)
select case(ndim)
  case(2)
      x1=coord(1,1); y1=coord(1,2)
      x2=coord(2,1); y2=coord(2,2)
      ell=sqrt((x2-x1)**2+(y2-y1)**2)
      cg=(x2-x1)/ell
      sg=(y2-y1)/ell
      local(1)=cg*global(1)+sg*global(2)
      local(2)=cg*global(2)-sg*global(1)
      local(3)=global(3)
      local(4)=cg*global(4)+sg*global(5)
      local(5)=cg*global(5)-sg*global(4)
      local(6)=global(6)
  case(3)
     x1=coord(1,1); y1=coord(1,2); z1=coord(1,3)
     x2=coord(2,1); y2=coord(2,2); z2=coord(2,3)
     xl=x2-x1; yl=y2-y1; zl=z2-z1; ell=sqrt(xl*xl+yl*yl+zl*zl)
     t=0.0
     pi=acos(-1.); gamrad=gamma*pi/180.; cg=cos(gamrad); sg=sin(gamrad)
     den=ell*sqrt(xl*xl+zl*zl)
     if(den /= 0.0)then
      r0(1,1)=xl/ell; r0(1,2)=yl/ell; r0(1,3)=zl/ell
      r0(2,1)=(-xl*yl*cg-ell*zl*sg)/den; r0(2,2)=den*cg/(ell*ell)
      r0(2,3)=(-yl*zl*cg+ell*xl*sg)/den; r0(3,1)=(xl*yl*sg-ell*zl*cg)/den
      r0(3,2)=-den*sg/(ell*ell); r0(3,3)=(yl*zl*sg+ell*xl*cg)/den
     else
      r0(1,1)=0.; r0(1,3)=0.; r0(2,2)=0.; r0(3,2)=0.
      r0(1,2)=1.; r0(2,1)=-cg; r0(3,3)=cg; r0(2,3)=sg; r0(3,1)=sg
     end if
     do i=1,3; do j=1,3 
      x=r0(i,j)
      do k=0,9,3; t(i+k,j+k)=x; end do
     end do; end do
     do i=1,12
      sum=0.
       do j=1,12
        sum=sum+t(i,j)*global(j)
       end do
       local(i)=sum
     end do
   end select
   return
end subroutine glob_to_loc                
subroutine kvdet(kv,n,iw,det,ksc)
!
!      this subroutine forms the determinant of the global stiffness
!      matrix stored as a vector (upper triangle)
!
 implicit none
 real,intent(inout)::kv(:)
 integer,intent(in)::n,iw
 real,intent(out)::det
 integer,intent(out)::ksc
 real::const,em
 integer::iwp2,l,ll,j,ma,ii,k,i
 iwp2=iw+2
 l=n-1
 do j=1,l
   if((iw-n+j) <= 0.0)ma=iwp2
   ma=ma-1
   const=1./kv(j)
   ii=j
   do k=2,ma
     ii=ii+1
     if(kv((k-1)*n+j) /= 0.0)then
       em=kv((k-1)*n+j)*const
       ll=0
       do i=k,ma
         ll=ll+1
         kv((ll-1)*n+ii)=kv((ll-1)*n+ii)-em*kv((i-1)*n+j)
       end do
     end if
   end do
 end do
 det=1.
 ksc=0
 do j=1,n
   det=det*kv(j)
   if(kv(j) < 0.)ksc=ksc+1
 end do
 return
end subroutine kvdet            
subroutine hinge(coord,holdr,action,react,prop,iel,etype,gamma)
!
!      this subroutine forms the end forces and moments to be
!      applied to a member if a joint has gone plastic
!
implicit none
!interface
!  subroutine glob_to_loc(local,global,gamma,coord)
!    real,intent(in)::global(:),gamma,coord(:,:)
!    real,intent(out)::local(:)
!  end subroutine glob_to_loc
!  subroutine loc_to_glob(local,global,gamma,coord)
!    real,intent(in)::local(:),gamma,coord(:,:)
!    real,intent(out)::global(:)
!  end subroutine loc_to_glob
!end interface
real,intent(in):: holdr(:,:),coord(:,:),action(:),prop(:,:),gamma(:)
real,intent(out):: react(:)
integer,intent(in)::etype(:),iel
real::ell,x1,x2,y1,y2,z1,z2,csch,snch,bm1,bm2,bm3,bm4,bm5,&
      s1,s2,s3,s4,s5,mpy,mpz,mpx
real::global(12),local(12),total(12),gam
integer::ndim
ndim=ubound(coord,2)
bm1=0.; bm2=0.; bm3=0.; bm4=0.; bm5=0.
total(:)=holdr(:,iel)
select case(ndim)
  case(1)
    mpy=prop(2,etype(iel))
    ell=coord(2,1)-coord(1,1)
    s1=total(2)+action(2); s2=total(4)+action(4)
  if(abs(s1) > mpy)then
    if(s1 >  0.)bm1= mpy-s1
    if(s1 <= 0.)bm1=-mpy-s1
  end if
  if(abs(s2) > mpy)then
    if(s2 >  0.)bm2= mpy-s2
    if(s2 <= 0.)bm2=-mpy-s2
  end if
    react(1)= (bm1+bm2)/ell; react(2)=bm1
    react(3)=-react(1);      react(4)=bm2
  case(2)
    mpy=prop(3,etype(iel))  
    x1=coord(1,1); y1=coord(1,2); x2=coord(2,1); y2=coord(2,2)
    ell=sqrt((y2-y1)**2+(x2-x1)**2)
    csch=(x2-x1)/ell; snch=(y2-y1)/ell
    s1=total(3)+action(3); s2=total(6)+action(6)
  if(abs(s1) > mpy)then
    if(s1 >  0.)bm1= mpy-s1
    if(s1 <= 0.)bm1=-mpy-s1
  end if
  if(abs(s2) > mpy)then
    if(s2 >  0.)bm2= mpy-s2
    if(s2 <= 0.)bm2=-mpy-s2
  end if
    react(1)=-(bm1+bm2)*snch/ell; react(2)= (bm1+bm2)*csch/ell; react(3)=bm1
    react(4)=-react(1); react(5)=-react(2); react(6)=bm2
  case(3)
    gam=gamma(iel)
    mpy=prop(5,etype(iel))
    mpz=prop(6,etype(iel))
    mpx=prop(7,etype(iel))
    x1=coord(1,1); y1=coord(1,2); z1=coord(1,3)
    x2=coord(2,1); y2=coord(2,2); z2=coord(2,3)
    ell=sqrt((z2-z1)**2+(y2-y1)**2+(x2-x1)**2)
    global=total+action
    call glob_to_loc(local,global,gam,coord)
    global=0.0
    s1=local(5); s2=local(11)
    if(abs(s1) > mpy)then
      if(s1 >  0.)bm1= mpy-s1
      if(s1 <= 0.)bm1=-mpy-s1
    end if
      if(abs(s2) > mpy)then
      if(s2 >  0.)bm2= mpy-s2
      if(s2 <= 0.)bm2=-mpy-s2
    end if
      local( 3)=-(bm1+bm2)/ell
      local( 9)=-local(3)
      local( 5)= bm1
      local(11)= bm2
    s3=local(6); s4=local(12)
    if(abs(s3) > mpz)then
      if(s3 >  0.)bm1= mpz-s3
      if(s3 <= 0.)bm1=-mpz-s3
    end if
      if(abs(s4) > mpy)then
      if(s4 >  0.)bm2= mpz-s4
      if(s4 <= 0.)bm2=-mpz-s4
    end if
      local( 2)=(bm3+bm4)/ell
      local( 8)=-local(2)
      local( 6)= bm3
      local(12)= bm4
      s5=local(4)
      if(abs(s5) > mpx)then
        if(s5 >  0.)global(4)= mpx-s5
        if(s5 <= 0.)global(4)=-mpx-s5
      end if
      local(10)=-local(4)
      call loc_to_glob(local,react,gam,coord)
  end select
return
end subroutine hinge                           
 subroutine rigid_jointed(km,prop,gamma,etype,iel,coord) 
!
!      this subroutine forms the stiffness matrix of a
!      general beam/column element(1-d, 2-d or 3-d)
!
 implicit none
   real,intent(in)::gamma(:),coord(:,:),prop(:,:) 
   integer,intent(in)::etype(:),iel
   real,intent(out):: km(:,:)
   integer::ndim,i,j,k
   real::ell,x1,x2,y1,y2,z1,z2,c,s,e1,e2,e3,e4,pi,xl,yl,zl,cg,sg,den
   real::ea,ei,eiy,eiz,gj
   real::a1,a2,a3,a4,a5,a6,a7,a8,sum,gamrad,x
   real::t(12,12),tt(12,12),cc(12,12),r0(3,3)
   ndim=ubound(coord,2)
   select case(ndim)
    case(1)
     ei=prop(1,etype(iel)); ell=coord(2,1)-coord(1,1)
     km(1,1)=12.*ei/(ell*ell*ell) ;km(3,3)=km(1,1)
     km(1,2)=6.*ei/(ell*ell) ;  km(2,1)=km(1,2) ; km(1,4)=km(1,2)
     km(4,1)=km(1,4) ;  km(1,3)=-km(1,1) ; km(3,1)=km(1,3) ; km(3,4)=-km(1,2)
     km(4,3)=km(3,4) ; km(2,3)=km(3,4) ;  km(3,2)=km(2,3); km(2,2)=4.*ei/ell
     km(4,4)=km(2,2) ;km(2,4)=2.*ei/ell ; km(4,2)=km(2,4)
    case(2)
     ea=prop(1,etype(iel)); ei=prop(2,etype(iel))
     x1=coord(1,1); y1=coord(1,2); x2=coord(2,1); y2=coord(2,2)
     ell=sqrt((y2-y1)**2+(x2-x1)**2)
     c=(x2-x1)/ell; s=(y2-y1)/ell
     e1=ea/ell; e2=12.*ei/(ell*ell*ell); e3=ei/ell; e4=6.*ei/(ell*ell)
     km(1,1)=c*c*e1+s*s*e2; km(4,4)=km(1,1); km(1,2)=s*c*(e1-e2)
     km(2,1)=km(1,2); km(4,5)=km(1,2); km(5,4)=km(4,5); km(1,3)=-s*e4
     km(3,1)=km(1,3); km(1,6)=km(1,3); km(6,1)=km(1,6); km(3,4)=s*e4 
     km(4,3)=km(3,4); km(4,6)=km(3,4); km(6,4)=km(4,6); km(1,4)=-km(1,1) 
     km(4,1)=km(1,4); km(1,5)=s*c*(-e1+e2); km(5,1)=km(1,5); km(2,4)=km(1,5)
     km(4,2)=km(2,4); km(2,2)=s*s*e1+c*c*e2; km(5,5)=km(2,2); km(2,5)=-km(2,2)
     km(5,2)=km(2,5); km(2,3)=c*e4; km(3,2)=km(2,3); km(2,6)=km(2,3)
     km(6,2)=km(2,6); km(3,3)=4.*e3; km(6,6)=km(3,3); km(3,5)=-c*e4
     km(5,3)=km(3,5); km(5,6)=km(3,5); km(6,5)=km(5,6); km(3,6)=2.*e3
     km(6,3)=km(3,6)
    case(3)
     ea =prop(1,etype(iel)); eiy=prop(2,etype(iel))
     eiz=prop(3,etype(iel)); gj =prop(4,etype(iel))
     x1=coord(1,1); y1=coord(1,2); z1=coord(1,3)
     x2=coord(2,1); y2=coord(2,2); z2=coord(2,3)
     xl=x2-x1; yl=y2-y1; zl=z2-z1; ell=sqrt(xl*xl+yl*yl+zl*zl)
     km=0.0; t=0.0; tt=0.0
     a1=ea/ell; a2=12.*eiz/(ell*ell*ell); a3=12.*eiy/(ell*ell*ell)
     a4=6.*eiz/(ell*ell); a5=6.*eiy/(ell*ell); a6=4.*eiz/ell
     a7=4.*eiy/ell; a8=gj/ell
     km(1,1)=a1; km(7,7)=a1; km(1,7)=-a1; km(7,1)=-a1; km(2,2)=a2; km(8,8)=a2
     km(2,8)=-a2; km(8,2)=-a2; km(3,3)=a3; km(9,9)=a3; km(3,9)=-a3; km(9,3)=-a3
     km(4,4)=a8; km(10,10)=a8; km(4,10)=-a8; km(10,4)=-a8; km(5,5)=a7
     km(11,11)=a7; km(5,11)=.5*a7; km(11,5)=.5*a7; km(6,6)=a6; km(12,12)=a6
     km(6,12)=.5*a6; km(12,6)=.5*a6; km(2,6)=a4; km(6,2)=a4; km(2,12)=a4
     km(12,2)=a4; km(6,8)=-a4; km(8,6)=-a4; km(8,12)=-a4; km(12,8)=-a4
     km(5,9)=a5; km(9,5)=a5; km(9,11)=a5; km(11,9)=a5; km(3,5)=-a5
     km(5,3)=-a5; km(3,11)=-a5; km(11,3)=-a5
     pi=acos(-1.); gamrad=gamma(iel)*pi/180.; cg=cos(gamrad); sg=sin(gamrad)
     den=ell*sqrt(xl*xl+zl*zl)
     if(den /= 0.0)then
      r0(1,1)=xl/ell; r0(1,2)=yl/ell; r0(1,3)=zl/ell
      r0(2,1)=(-xl*yl*cg-ell*zl*sg)/den; r0(2,2)=den*cg/(ell*ell)
      r0(2,3)=(-yl*zl*cg+ell*xl*sg)/den; r0(3,1)=(xl*yl*sg-ell*zl*cg)/den
      r0(3,2)=-den*sg/(ell*ell); r0(3,3)=(yl*zl*sg+ell*xl*cg)/den
     else
      r0(1,1)=0.; r0(1,3)=0.; r0(2,2)=0.; r0(3,2)=0.
      r0(1,2)=1.; r0(2,1)=-cg; r0(3,3)=cg; r0(2,3)=sg; r0(3,1)=sg
     end if
     do i=1,3; do j=1,3 
      x=r0(i,j)
      do k=0,9,3; t(i+k,j+k)=x; tt(j+k,i+k)=x; end do
     end do; end do
     do i=1,12; do j=1,12
      sum=0.0
      do k=1,12; sum=sum+km(i,k)*t(k,j); end do
      cc(i,j)=sum
     end do; end do
     do i=1,12; do j=1,12
      sum=0.
      do k=1,12; sum=sum+tt(i,k)*cc(k,j); end do
      km(i,j)=sum
     end do; end do
  end select
 return
 end subroutine rigid_jointed                
 subroutine pin_jointed(km,ea,coord)
!
!      this subroutine forms the stiffness matrix of a
!      general rod element(1-d, 2-d or 3-d)
!
 implicit none
   real,intent(in)::ea,coord(:,:); real,intent(out):: km(:,:)
   integer::ndim ,i,j
   real::eaol,ell,cs,sn,x1,x2,y1,y2,z1,z2,a,b,c,d,e,f,xl,yl,zl
   ndim=ubound(coord,2)!确定单元节点坐标第2维的大小，从1~3变化
   select case(ndim)
    case(1)!1维情况
      ell=coord(2,1)-coord(1,1)
      eaol=ea/ell
      km(1,1)=eaol;  km(1,2)=-eaol
      km(2,1)=-eaol; km(2,2)=eaol  
    case(2)!二维情况
      x1=coord(1,1); y1=coord(1,2)
      x2=coord(2,1); y2=coord(2,2)
      ell=sqrt((y2-y1)**2+(x2-x1)**2)
      cs=(x2-x1)/ell; sn=(y2-y1)/ell
      a=cs*cs; b=sn*sn; c=cs*sn
      km(1,1)=a; km(3,3)=a; km(1,3)=-a; km(3,1)=-a; km(2,2)=b; km(4,4)=b
      km(2,4)=-b; km(4,2)=-b; km(1,2)=c; km(2,1)=c; km(3,4)=c; km(4,3)=c
      km(1,4)=-c; km(4,1)=-c; km(2,3)=-c; km(3,2)=-c
      do i=1,4; do j=1,4; km(i,j)=km(i,j)*ea/ell; end do; end do
    case(3)!三维情况
      x1=coord(1,1); y1=coord(1,2); z1=coord(1,3)
      x2=coord(2,1); y2=coord(2,2); z2=coord(2,3)
      xl=x2-x1; yl=y2-y1; zl=z2-z1
      ell=sqrt(xl*xl+yl*yl+zl*zl)
      xl=xl/ell; yl=yl/ell; zl=zl/ell
      a=xl*xl; b=yl*yl; c=zl*zl; d=xl*yl; e=yl*zl; f=zl*xl
      km(1,1)=a; km(4,4)=a; km(2,2)=b; km(5,5)=b; km(3,3)=c; km(6,6)=c
      km(1,2)=d; km(2,1)=d; km(4,5)=d; km(5,4)=d; km(2,3)=e; km(3,2)=e
      km(5,6)=e; km(6,5)=e; km(1,3)=f; km(3,1)=f; km(4,6)=f; km(6,4)=f
      km(1,4)=-a; km(4,1)=-a; km(2,5)=-b; km(5,2)=-b; km(3,6)=-c; km(6,3)=-c
      km(1,5)=-d; km(5,1)=-d; km(2,4)=-d; km(4,2)=-d; km(2,6)=-e; km(6,2)=-e
      km(3,5)=-e; km(5,3)=-e; km(1,6)=-f; km(6,1)=-f; km(3,4)=-f; km(4,3)=-f
      do i=1,6; do j=1,6; km(i,j)=km(i,j)*ea/ell; end do; end do
   end select
 return
 end subroutine pin_jointed       
 subroutine glob_to_axial(axial,global,coord)
!
!      this subroutine transforms the global end reactions
!      into an axial force for rod elements (2-d or 3-d)
!
 implicit none
 real,intent(in)::global(:),coord(:,:)
 real,intent(out)::axial
 real::add,ell
 integer::ndim,i
 ndim=ubound(coord,2)
 add = 0.
 do i = 1 , ndim
   add = add + (coord(2,i)-coord(1,i))**2
 end do
 ell=sqrt(add)
 axial = 0.
 do i = 1 , ndim
   axial = axial + (coord(2,i)-coord(1,i))/ell*global(ndim+i)
 end do
 return
 end subroutine glob_to_axial                       
end module vlib               

