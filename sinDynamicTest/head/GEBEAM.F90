module gebeam
contains
!------------计算单元节点坐标和单元节点编号；--------------------------
!------------整体节点坐标和整体节点编号；------------------------------
!------------单元自由度编号和整体自由度编号；半带宽--------------------
subroutine geometry_get(iel,ell,coord,num)
! element node numbers, nodal coordinates and steering vectors for
! a line of (nonuniform) beam elements
!初始结构位于x轴上
implicit none
	  integer,intent(in)::iel; integer,intent(out)::num(:)
	  real,intent(in)::ell; real,intent(out)::coord(:,:)
	  num=(/iel,iel+1/);
	  if(iel==1)then;
		coord(1,1)=0.0; coord(2,1)=ell; 
	  else
		coord(1,1)=(iel-1)*ell; coord(2,1)=iel*ell;
	  end if
	  coord(1,2)=0;coord(1,3)=0;coord(2,2)=0;coord(2,3)=0;
end subroutine geometry_get 
!----------------------------COMPELEMENT ROTATION---------------------------------------
function comrot(va) result(cva)
implicit none
	real,intent(in)::va(:);
	real::cva(3),fi;
	real::pi=3.1415926;
	integer::n;
	
	fi=sqrt(dot_product(va,va));
	n=floor((fi+pi)/(2*pi));
	cva = va-(2*n*pi/fi)*va;
return
end function comrot
!----------------------------CROSS PRODUCT----------------------------------------------	
function cross_product(va,vb) result(vc)
implicit none
	real,intent(in)::va(:),vb(:)
	real::vc(3)

	vc(1)= va(2)*vb(3)-va(3)*vb(2);
	vc(2)= va(3)*vb(1)-va(1)*vb(3);
	vc(3)= va(1)*vb(2)-va(2)*vb(1);

return
end function cross_product
!--------------------------------ASYMMETRY MATRIX--------------------------------------
function asymmtx(v) result(amtx_v)
implicit none
	real::v(:)
	real::amtx_v(3,3)
	amtx_v = 0; 
	amtx_v(1,2) = -v(3);amtx_v(1,3) = v(2);amtx_v(2,3) = -v(1); 
	amtx_v(2,1) = -amtx_v(1,2);amtx_v(3,1) = -amtx_v(1,3);amtx_v(3,2) = -amtx_v(2,3);
return
end function asymmtx
!---------------------------------矢量A,B的并矢----------------------------------------
function biv(a,b) result(bivv)
implicit none
	real,intent(in)::a(:),b(:);
	real::bivv(3,3);
	integer::i,j
	
	do i=1,3
		do j=1,3
			bivv(i,j)=a(i)*b(j)
		end do
	end do
return
end function biv
!-----------------------------------------------------------------------
!subroutine spfn_get(der,der0,fun,qmtx,nespfn,points,i,jacobi,jacobi0)
subroutine spfn_get(der,fun,qmtx,nespfn,points,i,jacobi)
! this subroutine forms the beam shape functions
! and their 1st derivatives in global coordinates
implicit none
	 real,intent(in)::points(:,:),jacobi;
	 integer,intent(in)::i
	 real,intent(out)::der(:),fun(:),qmtx(:,:),nespfn(:,:)
	 real::xi
	 integer::j,k
	 xi=points(i,1); 
	 fun(1) = (1-xi)/2.0;fun(2) = (1+xi)/2.0; !某GUASS点线性形函数
	 !CALL INVERT(jacobi);!求JACOBIAN矩阵的逆
	 !der0(1)=-0.5*(1.0/jacobi0); der0(2)=0.5*(1.0/jacobi0);!初始时刻的值
	 der(1)=-0.5*(1.0/jacobi); der(2)=0.5*(1.0/jacobi);!全局S坐标系下的形函数导数
	 
	 qmtx = 0;
	 qmtx(1,1)=der(1);qmtx(2,2)=der(1);qmtx(3,3)=der(1);qmtx(4,4)=der(1);qmtx(5,5)=der(1);qmtx(6,6)=der(1);
	 qmtx(1,7)=der(2);qmtx(2,8)=der(2);qmtx(3,9)=der(2);qmtx(4,10)=der(2);qmtx(5,11)=der(2);qmtx(6,12)=der(2);
	 qmtx(7,4)=fun(1);qmtx(8,5)=fun(1);qmtx(9,6)=fun(1);
	 qmtx(7,10)=fun(2);qmtx(8,11)=fun(2);qmtx(9,12)=fun(2);
 
	 nespfn=0;
	 forall(j=1:6,k=1:6,j==k) nespfn(j,k)=fun(1);
	 forall(j=1:6,k=7:12,k==j+6) nespfn(j,k)=fun(2);
 
 end subroutine spfn_get
!------------------------------------BTCBK MATRIX---------------------------------------------------------
!subroutine btcbk_get(btcbk,bmtx,rmtx,tmtx,xc_p,nfc,mr,cst,cst_p,phai_p,phai,cnm,fi)
subroutine btcbk_get(btcbk,bmtx,rmtx,tmtx,rtxc_p,nfc,mr,cst,cst_p,phai_p,phai,cnm,fi)
implicit none

	integer::i,j;
	real,intent(in)::rtxc_p(:),rmtx(:,:),tmtx(:,:),nfc(:),mr(:),cst(:),cst_p(:),phai(:),phai_p(:),cnm(:,:);
	real,intent(out)::btcbk(:,:),bmtx(:,:);
	real::kmtx(9,9),kmat(9,9),c1mtx(3,3),c2mtx_nfc(3,3),c2mtx_mr(3,3),c3mtx(3,3),eyee(3,3);
	real::amtx_nfc_rt(3),fi
	!real::rtxc_p(3);
		
	!rtxc_p=matmul(transpose(rmtx),xc_p);!直接求RTXCP的空间分量
	
	eyee=0;forall(i=1:3,j=1:3,i==j)eyee(i,j)=1;
	amtx_nfc_rt=matmul(asymmtx(nfc),rtxc_p);!C2

	if(fi/=0.0)then
		!调用求解C1矩阵的子程序
		CALL c1mtx_get(c1mtx,cst,phai_p,phai); 
		!调用求解C2矩阵的子程序
		CALL c2mtx_get(c2mtx_nfc,cst,amtx_nfc_rt,phai);
		CALL c2mtx_get(c2mtx_mr,cst,mr,phai);
		!调用求解C3矩阵的子程序
		CALL c3mtx_get(c3mtx,cst,cst_p,mr,phai,phai_p,fi);
	else
		c1mtx=0.5*asymmtx(phai_p);
		c2mtx_nfc=-0.5*asymmtx(amtx_nfc_rt);
		c2mtx_mr=-0.5*asymmtx(mr);
		c3mtx=-1.0/3.0*dot_product(phai_p,mr)*eyee+1.0/6.0*(biv(phai_p,mr)+biv(mr,phai_p));
	end if

	bmtx=0;
	bmtx(1:3,1:3)=transpose(rmtx);
	bmtx(1:3,7:9)=matmul(asymmtx(rtxc_p),tmtx);
	bmtx(4:6,4:6)=tmtx;
	bmtx(4:6,7:9)=c1mtx;

	kmat=0;
	kmat = matmul(matmul(transpose(bmtx),cnm),bmtx);

	kmtx=0;
	kmtx(1:3,7:9)=-matmul(matmul(rmtx,asymmtx(nfc)),tmtx);
	kmtx(4:6,7:9)=c2mtx_mr;
	kmtx(7:9,7:9)=c3mtx+c2mtx_nfc+matmul(matmul(matmul(transpose(tmtx),asymmtx(nfc)),asymmtx(rtxc_p)),tmtx);  
	kmtx(7:9,1:3)=transpose(kmtx(1:3,7:9));
	kmtx(7:9,4:6)=transpose(kmtx(4:6,7:9));

	btcbk=0;
	btcbk=kmat+kmtx;
return
end subroutine btcbk_get
!----------------------------------C1 MATRIX-----------------------------------------------
subroutine c1mtx_get(c1mtx,cst,v,phai)
implicit none
	real,intent(in)::cst(:),v(3),phai(:)
	real,intent(out)::c1mtx(:,:)
	real::amtx_phai_v(3),amtx_phai_v_x_phai(3,3),phai_v
	real::eyee(3,3)
	integer::i,j
	eyee=0;forall(i=1:3,j=1:3,i==j)eyee(i,j)=1;

	amtx_phai_v=matmul(asymmtx(phai),v)!after c2
	phai_v=dot_product(phai,v) !after c3
	
	c1mtx = cst(1)*biv(v,phai)-cst(2)*biv(amtx_phai_v,phai)+         &
			cst(3)*phai_v*biv(phai,phai)-cst(4)*asymmtx(v)+          &
			cst(5)*(phai_v*eyee+biv(phai,v));
return
end subroutine c1mtx_get

!----------------------------------C2 MATRIX----------------------------------------------
subroutine c2mtx_get(c2mtx,cst,v,phai)
implicit none
real,intent(in)::cst(:),v(:),phai(:)
real,intent(out)::c2mtx(:,:)
real::amtx_phai_v(3),phai_v,cst2,cst4
	real::eyee(3,3)
	integer::i,j
	eyee=0;forall(i=1:3,j=1:3,i==j)eyee(i,j)=1;
	
	cst2=-cst(2);cst4=-cst(4);
	amtx_phai_v=matmul(asymmtx(phai),v)!after c2
	phai_v=dot_product(phai,v) !after c3
	c2mtx = cst(1)*biv(v,phai)-cst2*biv(amtx_phai_v,phai)+           &
			cst(3)*phai_v*biv(phai,phai)-cst4*asymmtx(v)+            &
			cst(5)*(phai_v*eyee+biv(phai,v))
return
end subroutine c2mtx_get
!---------------------------------C3 MATRIX---------------------------------------------
subroutine c3mtx_get(c3mtx,cst,cst_p,v,phai,phai_p,fi)
implicit none
	real,intent(out)::c3mtx(:,:)
	real,intent(in)::cst(:),cst_p(:),v(:),phai(:),phai_p(3)
	integer::i,j
	real::amtx_v_phai_p(3),&
	sym_phai_x_amtx_v_php(3,3),&      !after c2
	sym_phai_x_v(3,3),&               !after c3
	sym_phai_x_phai_p(3,3),&
	sym_phai_p_x_v(3,3)               !after c5
	real::eyee(3,3),fi
	eyee=0;forall(i=1:3,j=1:3,i==j)eyee(i,j)=1;
	
	amtx_v_phai_p=matmul(asymmtx(v),phai_p)         !after c2
	do i=1,3
		do j=1,3
			sym_phai_x_amtx_v_php(i,j)=phai(i)*amtx_v_phai_p(j)+amtx_v_phai_p(i)*phai(j)
			sym_phai_x_v(i,j)=phai(i)*v(j)+v(i)*phai(j)
			sym_phai_x_phai_p(i,j)=phai(i)*phai_p(j)+phai_p(i)*phai(j)
			sym_phai_p_x_v(i,j)=phai_p(i)*v(j)+v(i)*phai_p(j)		
		end do
	end do

	c3mtx=(cst(1)*dot_product(phai_p,v)+cst(2)*dot_product(phai,amtx_v_phai_p)          &
		  +cst(3)*dot_product(phai,phai_p)*dot_product(v,phai))*eyee                    &
		  +cst(2)*sym_phai_x_amtx_v_php+cst(3)*dot_product(phai,phai_p)*sym_phai_x_v      &
		  +1.0/fi*(cst_p(1)*dot_product(phai_p,v)+cst_p(2)*dot_product(phai,amtx_v_phai_p)&
		  +cst_p(3)*dot_product(phai,phai_p)*dot_product(v,phai))*biv(phai,phai)          &
		  +cst(3)*dot_product(v,phai)*sym_phai_x_phai_p+cst(5)*sym_phai_p_x_v
return
end subroutine c3mtx_get
!---------------------------------C4 MATRIX------------------------------------------------
subroutine c4mtx_get(c4mtx,dphai,phai,cst)
implicit none
	real,intent(in)::dphai(:),phai(:),cst(:);
	real,intent(out)::c4mtx(3,3);
	real::temp1(3,3),temp2(3,3),eyee(3,3);!after c2
	integer::i,j
	eyee=0;forall(i=1:3,j=1:3,i==j)eyee(i,j)=1;
	
	temp1 = biv(matmul(asymmtx(phai),dphai),phai);
	!temp2 = biv(matmul(phai,dphai),phai);
	temp2=0.0;
	
	c4mtx = (cst(1)+cst(5))*dot_product(phai,dphai)*eyee+2*cst(3)*dot_product(phai,dphai)*biv(phai,phai)+ &
			2*cst(5)*biv(phai,dphai)+(cst(1)+cst(5))*biv(dphai,phai)+ &
			-cst(2)*(temp1+temp2);
return
end subroutine c4mtx_get
!---------------------------------C5 MATRIX------------------------------------------------
subroutine c5mtx_get(c5mtx,phai,v,cst,cst_p,fi)
implicit none
	real,intent(in)::phai(:),v(:),cst(:),cst_p(:),fi;
	real,intent(out)::c5mtx(3,3)
	real::temp(3,3) !after c2 prime
	real::eyee(3,3)
	integer::i,j
	eyee=0;forall(i=1:3,j=1:3,i==j)eyee(i,j)=1;
	temp = biv(matmul(asymmtx(phai),v),phai);
	c5mtx = (cst(3)*dot_product(phai,v)**2+cst(5)*dot_product(v,v))*eyee+ &
			(cst_p(3)/fi*dot_product(phai,v)**2+cst(3)*dot_product(v,v))*biv(phai,phai)+ &
			2*cst(3)*dot_product(phai,v)*biv(phai,v)+(cst_p(1)/fi+cst(3))*dot_product(phai,v)*biv(v,phai)- &
			cst_p(2)/fi*dot_product(phai,v)*temp-cst(2)*temp+cst(2)*dot_product(phai,v)*asymmtx(v)+ &
			(cst(1)+cst(5))*biv(v,v);

return
end subroutine c5mtx_get
!---------------------------------GET N,MR------------------------------------------------
!subroutine nfcmr_get(nfc,mr,cnm,rtxc_p,tphai_p,E1_mat)
subroutine nfcmr_get(nfc,mr,cnm,rtxc_p,tphai_p)
implicit none
!real,intent(in)::cnm(:,:),tphai_p(:),rtxc_p(:),E1_mat(:);
real,intent(in)::cnm(:,:),tphai_p(:),rtxc_p(:);
real,intent(out)::nfc(:),mr(:)
real::E1(3)=(/1.0,0.0,0.0/);
	
	!nfc=MATMUL(cnm(1:3,1:3),rtxc_p-E1_mat);
	nfc=MATMUL(cnm(1:3,1:3),rtxc_p-E1);
	mr=matmul(cnm(4:6,4:6),tphai_p);
	
return
end subroutine nfcmr_get
!------------------------求解单位长梁截面合外力与合外力偶---------------------------------
subroutine nmrbar_get(nbar,mrbar,perwet,rhoaw,diacab,elxc,dxc,d2xc,rmtx,points,k,wtvel0)! dynamic analysis
!subroutine nmrbar_get(nbar,mrbar,perwet,rhoaw,diacab,rmtx,wtvel0);!static analysis，常速海流
!subroutine nmrbar_get(nbar,mrbar,perwet,rhoaw,diacab,rmtx,elxc,points,k,wtvel0);!static analysis,递减海流
implicit none
	!以el-开头的表示单元节点值，如elxc,elphai
	real,intent(out)::nbar(:),mrbar(:);
	integer,intent(in)::k;
	real,intent(in)::rhoaw,diacab,perwet,rmtx(:,:),wtvel0(:);
	real,intent(in)::elxc(:,:),dxc(:),d2xc(:),points(:,:);!动态分析时打开
	!real,intent(in)::elxc(:,:),points(:,:);!递减海流时打开，静态
	real::gravity(3),gforce(3),buofc(3),addmass(3),cm,zguass,z1,z2,xi,vcu(3),wtvel(3),&
		 vr(3),fdl(3),fdg(3),cd1,cd2,cd3,E2(3),zall;
	real::pi =3.1415;
	
	cm = 1.0;!附加质量力系数
	cd1 =0.02; !水阻力系数,切向
	cd2 =2.0;  !水阻力系数,法向
	cd3 =2.0;  !水阻力系数,法向
	zall = 5000;!全海深


	xi=points(k,1); 
	z1 = elxc(1,1);!节点1的z值，当前值
	z2 = elxc(2,1);!节点2的z值，当前值
	!gravity=(/9.8,0.0,0.0/);!or -9.8m/s2
	!E2 = (/0.0,diacab*0.5,0.0/);! X轴垂直向下，Y轴向右，Z轴垂直纸面向外

	!单位体积的重力矢量
	!gforce=rhoac*gravity;
	!单位体积的浮力矢量
	!buofc=rhoaw*gravity;
	!单位缆长的附加质量力矢量,当水流相对缆索切向运动时附加质量力为零，见F. R. Driscoll博士论文APPENDIX C
	!addmass=rhoaw*cm*pi/4.0*diacab*diacab*d2xc;
	
	!单位缆长的水阻力矢量
	!!zguass=(xi+1)/2.0*(z2-z1)+z1;!高斯点对应的实际深度值		
	!全深度h下的海流流速函数
	!!wtvel = wtvel0*(1.0-zguass/zall);!单元高斯点对应的海流
	
	!wtvel = wtvel0![全深度等海流时打开]
	
	!vr=matmul(transpose(rmtx),(-wtvel));!static,无升沉有海流
	!!vr=matmul(transpose(rmtx),(dxc-wtvel));!dynamic，高斯点上的相对速度：[缆索速度-水流速度]；
		!write(11,*)"wtvel",wtvel;
	!write(11,*)rmtx
	!局部坐标系下的水阻力【有海流时考虑】
	!!fdl(1)= -0.5*rhoaw*cd1*pi*diacab*vr(1)*abs(vr(1));!切向
	!!fdl(2)= -0.5*rhoaw*cd2*diacab*vr(2)*sqrt(abs(vr(2)**2+vr(3)**2));!法向
	!!fdl(3)= -0.5*rhoaw*cd3*diacab*vr(3)*sqrt(abs(vr(2)**2+vr(3)**2));!法向
		!write(11,*)"fdl",fdl;	
	!整体坐标系下的水阻力【有海流时考虑】
	!!fdg= matmul(rmtx,fdl);
		!write(11,*)"the hydrodynamic forces:-------------------------"
		!write(11,*)fdg;
	!单位缆长合外力与合外力偶
	!!nbar= (/perwet,0.0,0.0/)+fdg;!有海流情形
	
	!静水条件下的单位缆长水阻力【无海流时考虑】
	fdl(1)= -0.5*rhoaw*cd1*pi*diacab*dxc(1)*abs(dxc(1));!切向
	nbar= (/perwet+fdl(1),0.0,0.0/);!无海流情形
	
	!mrbar= cross_product(matmul(rmtx,E2),fdg);
	!mrbar= cross_product(E2,fdg);
	mrbar =0;
end subroutine nmrbar_get
!--------------------------------MASS MATRIIX-----------------------------------------
!subroutine mass(massx,area,rhoac,tmtx,jtensor)
subroutine mass(massx,tmtx,jtensor,permass)
implicit none
	real,intent(in)::tmtx(:,:),jtensor(:,:),permass;
	real,intent(out)::massx(:,:)
	integer::i,j
	real::addmass,pi = 3.1415; !附加质量

	massx=0;
	addmass =0.0;!无海流时打开
	!!addmass = 1025*1.5*pi/4.0*0.03*0.03;!addmass=rhoaw*cm*pi/4.0*diacab*diacab,cm=2时addmass = 1.449，有海流时打开
	forall(i=1:3,j=1:3,i==j)massx(i,j)=permass+addmass;!单位：kg/m
	massx(4:6,4:6)=matmul(matmul(transpose(tmtx),jtensor),tmtx);
	!do i=1,6;write(11,"(6f10.6)")massx(i,:);end do

return
end subroutine mass		
!--------------------------------KLOAD MATRIX----------------------------------------
subroutine kload(loadx,cst,phai,mrbar)
implicit none
	real,intent(in)::cst(:),phai(:),mrbar(:);
	real,intent(out)::loadx(:,:)
	real::c2mtx_mrbar(3,3)

	CALL c2mtx_get(c2mtx_mrbar,cst,mrbar,phai);
	loadx=0;
	loadx(4:6,4:6) = c2mtx_mrbar;
	
return
end subroutine kload
!---------------------惯性力虚功项产生的离心力矩阵KCENT(6,6)--------------------------
subroutine kcent_gyro(centx,gyrox,phai,dphai,d2phai,tmtx,dertf,jtensor,cst,cst_p,fi)
implicit none
	real,intent(in)::phai(:),dphai(:),d2phai(:),tmtx(:,:),dertf(:,:),jtensor(:,:),cst(:),cst_p(:),fi;
	real,intent(out)::centx(:,:),gyrox(:,:);
	real::c2temp(3),wojwo(3,3),angvel(3),angacc(3),c2mtx(3,3),c1mtx1(3,3),c1mtx2(3,3),c5mtx(3,3),c4mtx(3,3);
	real::eyee(3,3)
	integer::i,j
	eyee=0;forall(i=1:3,j=1:3,i==j)eyee(i,j)=1;
	
	angvel = matmul(tmtx,dphai);
	angacc = matmul(tmtx,d2phai)+matmul(dertf,dphai);
	c2temp = matmul(matmul(asymmtx(angvel),jtensor),angvel)+matmul(jtensor,angacc);
	wojwo = matmul(asymmtx(angvel),jtensor)-asymmtx(matmul(jtensor,angvel));
	
	if(fi/=0)then
		call c2mtx_get(c2mtx,cst,c2temp,phai);
		call c1mtx_get(c1mtx1,cst,dphai,phai);
		call c1mtx_get(c1mtx2,cst,d2phai,phai);
		call c5mtx_get(c5mtx,phai,dphai,cst,cst_p,fi);
		call c4mtx_get(c4mtx,dphai,phai,cst);
	else
		c2mtx=-0.5*asymmtx(c2temp);
		c1mtx1=0.5*asymmtx(dphai);
		c1mtx2=0.5*asymmtx(d2phai);
		c5mtx=-1.0/6.0*biv(dphai,dphai)+1.0/6.0*dot_product(dphai,dphai)*eyee;
		c4mtx=0.0;
	end if
	centx = 0;
	centx(4:6,4:6) = c2mtx+matmul(transpose(tmtx),matmul(wojwo,c1mtx1)+matmul(jtensor,c5mtx+c1mtx2));
	
	gyrox = 0.0;
	gyrox(4:6,4:6) = matmul(transpose(tmtx),matmul(wojwo,tmtx)+matmul(jtensor,c4mtx));
	
return
end subroutine
!-------------------------------转换矩阵T的时间导数--------------------------------------------
subroutine der_tf(dertf,phai,dphai,cst,fi)
implicit none
	real,intent(in)::phai(:),dphai(:),cst(:),fi;
	real,intent(out)::dertf(:,:)
	real::eyee(3,3)
	integer::i,j
	eyee=0;forall(i=1:3,j=1:3,i==j)eyee(i,j)=1;
	
	if(fi/=0.0)then
		dertf=dot_product(phai,dphai)*(cst(1)*eyee-cst(2)*asymmtx(phai)+cst(3)*biv(phai,phai))+ &
			cst(4)*asymmtx(dphai)+cst(5)*(biv(phai,dphai)+biv(dphai,phai))
	else
		dertf=-0.5*asymmtx(dphai);
	end if

end subroutine der_tf
!-------------------------------外部节点力内部矢量-------------------------------------------
subroutine extfcmx_get(extfcmx,nbar,mrbar,tmtx)
implicit none
	real,intent(in)::nbar(:),mrbar(:),tmtx(:,:);
	real,intent(out)::extfcmx(:);
	
	extfcmx(1:3)=nbar;
	extfcmx(4:6)=matmul(transpose(tmtx),mrbar);
return
end subroutine extfcmx_get
!-------------------------------材料内部节点力内部矢量-------------------------------------------
subroutine intfcmx_get(intfcmx,nfc,mr)
implicit none
	real,intent(in)::nfc(:),mr(:);
	real,intent(out)::intfcmx(6);
	
	intfcmx(1:3)=nfc;
	intfcmx(4:6)=mr;
return
end subroutine intfcmx_get
!-------------------------------惯性节点力内部矢量--------------------------------------------
subroutine inefcmx_get(inefcmx,tmtx,phai,dphai,cst,jtensor,fi,dertf)
implicit none
	real,intent(in)::dertf(:,:),phai(:),dphai(:),tmtx(:,:),jtensor(:,:),cst(:),fi;
	real,intent(out)::inefcmx(:);
	real::temp1(3),temp2(3),temp3(3);
	
	temp1 = matmul(tmtx,dphai);
	temp2 = matmul(dertf,dphai);
	temp3 = matmul(matmul(asymmtx(temp1),jtensor),temp1)+matmul(jtensor,temp2);
	
	inefcmx(1:3)=0;
	inefcmx(4:6)= matmul(transpose(tmtx),temp3);
	
return
end subroutine inefcmx_get
!-------------------------------CHECK CONVERGEED----------------------------------------------
subroutine checonverg(ttd,delta,tol,converged)
! sets converged to .false. if relative change in loads and
! oldlds is greater than tol and updates oldlds
implicit none
	real,intent(in)::ttd(0:),tol;real,intent(in out)::delta(0:)
	logical,intent(out)::converged
	  converged=.true.
	  converged=(maxval(abs(delta))/maxval(abs(ttd))<=tol)
return
end subroutine checonverg
!-------------------------------transformation among vector------------------------------------
!---------要注意的是：这在转换节点处的值，不能用GUASS点上的量----------------------------------
!SUBROUTINE t2e_vtr(t2ev,rmtx0_nd1,rmtx0_nd2,rmtx_nd1,rmtx_nd2,tmtx_nd1,tmtx_nd2)
!implicit none
!	real,intent(in)::rmtx0_nd1(:,:),rmtx0_nd2(:,:),rmtx_nd1(:,:),rmtx_nd2(:,:),tmtx_nd1(:,:),tmtx_nd2(:,:);
!	real,intent(out)::t2ev(:,:);
!	real::eyee(3,3);
!	integer::i,j;
!	eyee=0;forall(i=1:3,j=1:3,i==j)eyee(i,j)=1;
!	
!	t2ev=0;
!	!NODE 1:
!	t2ev(1:3,1:3)=matmul(rmtx_nd1,rmtx0_nd1);!static
!	!t2ev(4:6,4:6)=tmtx_nd1;
!	t2ev(4:6,4:6)=eyee;
!	
!	!NODE 2:
!	t2ev(7:9,7:9)=matmul(rmtx_nd2,rmtx0_nd2);!static 	
!	!t2ev(10:12,10:12)=tmtx_nd2;
!	t2ev(10:12,10:12)=eyee;
!end subroutine t2e_vtr
!-------------------------------transformation among matrix---------------------------------------
!SUBROUTINE t2e_mtx(t2em1,t2em2,rmtx0_nd1,rmtx0_nd2,rmtx_nd1,rmtx_nd2,tmtx_nd1,tmtx_nd2)
!implicit none
!	real,intent(in)::rmtx0_nd1(:,:),rmtx0_nd2(:,:),rmtx_nd1(:,:),rmtx_nd2(:,:),tmtx_nd1(:,:),tmtx_nd2(:,:);
!	real,intent(out)::t2em1(:,:),t2em2(:,:);
!	real::eyee(3,3);
!	integer::i,j;
!	eyee=0;forall(i=1:3,j=1:3,i==j)eyee(i,j)=1;
	
	!---------------THE MATRIX PUT ON BEFORE
!	t2em1=0;
!!	!NODE 1:
!	t2em1(1:3,1:3)=matmul(rmtx_nd1,rmtx0_nd1);!static	
!	t2em1(4:6,4:6)=transpose(tmtx_nd1);
!	!t2em1(4:6,4:6)=tmtx_nd1;
!	!t2em1(4:6,4:6)=eyee;
!	!NODE 2:
!	t2em1(7:9,7:9)=matmul(rmtx_nd2,rmtx0_nd2);!static
!	t2em1(10:12,10:12)=transpose(tmtx_nd2);
!	!t2em1(10:12,10:12)=tmtx_nd2;
!	!t2em1(10:12,10:12)=eyee;
	
!	!--------------THE MATRIX PUT ON THE BACK	
!	t2em2=0;
!	!NODE 1:
!	t2em2(1:3,1:3)=transpose(matmul(rmtx_nd1,rmtx0_nd1));!static
!	!t2em2(4:6,4:6)=transpose(tmtx_nd1);
!	t2em2(4:6,4:6)=tmtx_nd1;
!	!t2em2(4:6,4:6)=eyee;

!	!NODE 2:
!	t2em2(7:9,7:9)=transpose(matmul(rmtx_nd2,rmtx0_nd2));!static
	!(10:12,10:12)=transpose(tmtx_nd2);
!	t2em2(10:12,10:12)=tmtx_nd2;
	!t2em2(10:12,10:12)=eyee;	
!end subroutine t2e_mtx
!-------------------------------transformation matrix rmtx0 at nodes---------------------------------------
!节点【材料坐标系】与【惯性坐标系】的转换矩阵——rmtx0;
!subroutine rmtx0_node(rmtx0_nd1,rmtx0_nd2,coord,radius)
!implicit none
!	real,intent(in)::coord(:,:),radius;
!	real,intent(out)::rmtx0_nd1(:,:),rmtx0_nd2(:,:);
!	real::xx1,xx2,yy1,yy2,fi1,fi2,phai1(3),phai2(3);
!	real::eyee(3,3),e2(3)=(/0.0,-1.0,0.0/);
!	integer::i,j;
	
!	eyee=0;forall(i=1:3,j=1:3,i==j)eyee(i,j)=1;	
	
!	xx1=coord(1,1);yy1=coord(1,3);
!	xx2=coord(2,1);yy2=coord(2,3);
	
!	fi1=atan(xx1/(radius-yy1));		
!	fi2=atan(xx2/(radius-yy2));
	
!	phai1=fi1*e2;phai2=fi2*e2; 

!	if(fi1/=0.0)then;
!		rmtx0_nd1 = eyee+sin(fi1)/fi1*asymmtx(phai1)+(1-cos(fi1))/(fi1*fi1)*matmul(asymmtx(phai1),asymmtx(phai1));
!	else
!		rmtx0_nd1 = eyee;		
!	end if
	
!	if(fi2/=0.0)then;
!		rmtx0_nd2 = eyee+sin(fi2)/fi2*asymmtx(phai2)+(1-cos(fi2))/(fi2*fi2)*matmul(asymmtx(phai2),asymmtx(phai2));
!	else
!		rmtx0_nd2 = eyee;		
!	end if				
!end subroutine rmtx0_node
!-------------------------------transformation matrix rmtx at nodes---------------------------------------
!subroutine rtmtx_node(rmtx_nd1,rmtx_nd2,tmtx_nd1,tmtx_nd2,elphai)
!implicit none
!	real,intent(in)::elphai(:,:);
!!	real,intent(out)::rmtx_nd1(:,:),rmtx_nd2(:,:),tmtx_nd1(:,:),tmtx_nd2(:,:);
!	real::phai1(3),phai2(3);
!	real::fi1,fi2;
!	real::eyee(3,3);
!	integer::i,j;
	
!	eyee=0;forall(i=1:3,j=1:3,i==j)eyee(i,j)=1;	
	
!	phai1=elphai(1,:);
!	phai2=elphai(2,:);
	
!	fi1 = sqrt(phai1(1)**2+phai1(2)**2+phai1(3)**2);
!	fi2 = sqrt(phai2(1)**2+phai2(2)**2+phai2(3)**2);
	
!	if(fi1/=0.0)then;
!		rmtx_nd1 = eyee+sin(fi1)/fi1*asymmtx(phai1)+(1-cos(fi1))/(fi1*fi1)*matmul(asymmtx(phai1),asymmtx(phai1));
!		tmtx_nd1 = sin(fi1)/fi1*eyee-(1-cos(fi1))/fi1**2*asymmtx(phai1)+(fi1-sin(fi1))/fi1**3*biv(phai1,phai1);		
!	else
!		rmtx_nd1=eyee;
!		tmtx_nd1=eyee;
!	end if
	
!	if(fi2/=0.0)then;
!		rmtx_nd2 = eyee+sin(fi2)/fi2*asymmtx(phai2)+(1-cos(fi2))/(fi2*fi2)*matmul(asymmtx(phai2),asymmtx(phai2));
!		tmtx_nd2 = sin(fi2)/fi2*eyee-(1-cos(fi2))/fi2**2*asymmtx(phai2)+(fi2-sin(fi2))/fi2**3*biv(phai2,phai2);
!	else 
!		rmtx_nd2=eyee;
!		tmtx_nd2=eyee;
!	end if
!return
!end subroutine rtmtx_node
!---------------------------------------单元内部节点力/力偶矢量变换矩阵（高斯点）---------------------------------------
subroutine E2e_vtr(E2ev,rmtx0)
implicit none
	real,intent(in)::rmtx0(:,:);
	real,intent(out)::E2ev(:,:);
	real::eyee(3,3);
	integer::i,j;
	eyee=0;forall(i=1:3,j=1:3,i==j)eyee(i,j)=1;
	
	E2ev=0;
	!NODE 1:
	E2ev(1:3,1:3)=rmtx0;!static,t->e	
	E2ev(4:6,4:6)=rmtx0;!E->e
	
	!NODE 2:
	E2ev(7:9,7:9)=rmtx0;!static,t->e 
	E2ev(10:12,10:12)=rmtx0;!E->e
	
return
end subroutine E2e_vtr
!---------------------------------------单元刚度矩阵变换矩阵（高斯点）---------------------------------------
subroutine E2e_mtx(E2em1,E2em2,rmtx0)
implicit none
	real,intent(in)::rmtx0(:,:);
	real,intent(out)::E2em1(:,:),E2em2(:,:);
	real::eyee(3,3);
	integer::i,j;
	eyee=0;forall(i=1:3,j=1:3,i==j)eyee(i,j)=1;
	
	!---------------E2em1--------------
	E2em1=0;
	!NODE 1:
	E2em1(1:3,1:3)=rmtx0;
	E2em1(4:6,4:6)=rmtx0;

	
	!NODE 2:
	E2em1(7:9,7:9)=rmtx0;
	E2em1(10:12,10:12)=rmtx0;
	
	!---------------E2em2--------------
	E2em2=0;
	!NODE 1:
	E2em2(1:3,1:3)=transpose(rmtx0);	
	E2em2(4:6,4:6)=transpose(rmtx0);

	!NODE 2:
	E2em2(7:9,7:9)=transpose(rmtx0);
	E2em2(10:12,10:12)=transpose(rmtx0);
	
end subroutine E2e_mtx
!---------------------------------------整体至材料坐标系变换矩阵（高斯点）---------------------------------------
subroutine rmtx00_get(rmtx0,rmtx0_nd1,rmtx0_nd2,coord,radius,fun);
implicit none
	real,intent(in)::coord(:,:),radius,fun(:);
	real,intent(out)::rmtx0(:,:),rmtx0_nd1(:,:),rmtx0_nd2(:,:);	
	real::xx1,xx2,yy1,yy2,fi1,fi2,fi,phai1(3),phai2(3),phai(3);
	real::eyee(3,3),e2(3)=(/0.0,0.0,-1.0/);
	integer::i,j;
	real,parameter::pi=3.1415926;
	
	eyee=0;forall(i=1:3,j=1:3,i==j)eyee(i,j)=1;	
	
	xx1=coord(1,1);yy1=coord(1,2);
	xx2=coord(2,1);yy2=coord(2,2);
	!write(11,*)xx1,yy1;
	
	fi1=atan(yy1/(radius-xx1));	!单元第1节点处的转角大小；	
	fi2=atan(yy2/(radius-xx2)); !单元第2节点处的转角大小；
	
	phai1=fi1*e2;phai2=fi2*e2; !单元第1，2节点处的转动矢量； 
	
	phai = fun(1)*phai1+fun(2)*phai2; !GUASS点处的截面转动矢量（与后面的RMTX保持一致）
    !fi = sqrt(phai(1)**2+phai(2)**2+phai(3)**2);!GUASS点处的截面转动矢量大小。
	fi = fun(1)*fi1 + fun(2)*fi2;
		write(11,*)"initial rotation angle,fi：";
		write(11,*)fi*180/pi;
		write(11,*)"initial rotation angle,fi1：";
		write(11,*)fi1*180/pi;
		write(11,*)"initial rotation angle,fi2：";
		write(11,*)fi2*180/pi;
	
	if(fi/=0.0)then;!GUASS点处的转动矩阵
		rmtx0 = eyee+sin(fi)/fi*asymmtx(phai)+(1-cos(fi))/(fi*fi)*matmul(asymmtx(phai),asymmtx(phai))
	else
		rmtx0 = eyee;
	end if
	
	if(fi1/=0.0)then;!单元节点处的转动矩阵
		rmtx0_nd1 = eyee+sin(fi1)/fi1*asymmtx(phai1)+(1-cos(fi1))/(fi1*fi1)*matmul(asymmtx(phai1),asymmtx(phai1));
		rmtx0_nd2 = eyee+sin(fi2)/fi2*asymmtx(phai2)+(1-cos(fi2))/(fi2*fi2)*matmul(asymmtx(phai2),asymmtx(phai2));
	else
		rmtx0_nd1 = eyee;
		rmtx0_nd2 = eyee;
	end if

end subroutine rmtx00_get
!---------------------------------------材料至运动坐标系变换矩阵（高斯点）---------------------------------------
subroutine rtmtx_get(rmtx,tmtx,phai,fi)
implicit none	
	real,intent(in)::phai(:),fi;
	real,intent(out)::rmtx(:,:),tmtx(:,:);
	real::eyee(3,3);
	integer::i,j;
	
	eyee=0;forall(i=1:3,j=1:3,i==j)eyee(i,j)=1;	
	if(fi/=0.0)then
		rmtx = eyee+sin(fi)/fi*asymmtx(phai)+(1-cos(fi))/(fi*fi)*matmul(asymmtx(phai),asymmtx(phai));
		tmtx = sin(fi)/fi*eyee-(1-cos(fi))/(fi**2)*asymmtx(phai)+(fi-sin(fi))/(fi**3)*biv(phai,phai);
		!tmtx = eyee-(1-cos(fi))/fi**2*asymmtx(phai)+(fi-sin(fi))/fi**3*biv(phai,phai);
	else
		rmtx=eyee;			
		tmtx=eyee;
	end if
end subroutine rtmtx_get
!---------------------------------------总体坐标系->材料坐标系1,2的变换矩阵--------------------------------------
subroutine rmtx0_get(rmtx0_1,rmtx0_2,phai0_1,phai0_2)
implicit none
	real,intent(in)::phai0_1(:),phai0_2(:);
	real,intent(out)::rmtx0_1(:,:),rmtx0_2(:,:);
	real::fi0_1,fi0_2;
	real::eyee(3,3);
	integer::i,j;
	
	eyee=0;forall(i=1:3,j=1:3,i==j)eyee(i,j)=1;	
	fi0_1 = sqrt(dot_product(phai0_1,phai0_1));
	fi0_2 = sqrt(dot_product(phai0_2,phai0_2));
	
	!write(11,*)fi0_1,fi0_2;
	
	rmtx0_1 = eyee+sin(fi0_1)/fi0_1*asymmtx(phai0_1)+(1-cos(fi0_1))/(fi0_1*fi0_1)*matmul(asymmtx(phai0_1),asymmtx(phai0_1));
	rmtx0_2 = eyee+sin(fi0_2)/fi0_2*asymmtx(phai0_2)+(1-cos(fi0_2))/(fi0_2*fi0_2)*matmul(asymmtx(phai0_2),asymmtx(phai0_2));
	
end subroutine rmtx0_get
!---------------------------------------系数C1,C2,C3,C4,C4等-----------------------------------------------------
subroutine coef(cst,cst_p,fi)
implicit none
	real,intent(in)::fi;
	real,intent(out)::cst(:),cst_p(:);
	
	if(fi/=0.0)then
		cst(1)=(fi*cos(fi)-sin(fi))/(fi**3);
		cst(2)=(fi*sin(fi)+2*cos(fi)-2)/(fi**4);
		cst(3)=(3*sin(fi)-2*fi-fi*cos(fi))/(fi**5);
		cst(4)=(cos(fi)-1)/(fi**2);
		cst(5)=(fi-sin(fi))/(fi**3);
		cst_p(1)=((3-fi*fi)*sin(fi)-3*fi*cos(fi))/(fi**4);
		cst_p(2)=((fi*fi-8)*cos(fi)-5*fi*sin(fi)+8)/(fi**5);
		cst_p(3)=((fi*fi-15)*sin(fi)+7*fi*cos(fi)+8*fi)/(fi**6);	
	else
		cst(1)=-1.0/3.0;
		cst(2)=-1.0/12.0;
		cst(3)=-1.0/60.0;
		cst(4)=-0.5;
		cst(5)=1.0/6.0;    
		cst_p(1)=0.0;
		cst_p(2)=0.0;
		cst_p(3)=0.0;	
	end if
end subroutine coef

end module gebeam