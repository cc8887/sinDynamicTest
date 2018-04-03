program cable
!------------------------------------------------------------------------------
! program gebeam: transient analysis of underwater cable by using  
!              3D  geometrically exact beam element.
!------------------------------------------------------------------------------
use new_library ;use geometry_lib;use gebeam ;use libks ;use vlib;implicit none 
!----------------------------static variables----------------------------------
integer::nels,neq,nn,nband,nr,nod=2,nodof=6,ndof=12,iel,i,j,k,ndim=3,  &
		nip1,nip2,LIMIT,iters,loaded_nodes,stat,ocean_step,incs,iy;	
real::beta,gamma,dt,dt0,EA,GA2,GA3,GJ,EI2,EI3,jacobi,fi,fi1,fi2,ell,ts,area,rhoaw,&
		rhoac,diacab,alpha1,alpha2,alpha3,alpha4,alpha5,alpha6,cd1,cd2,TOL,factfun,&
		amp,omego,period,jtensor1,jtensor2,jtensor3,perwet,permass,masscage,cfcage,&
		alpha_m,alpha_f,spradius,factor;
character(len=15) :: element = 'line'; 
real,parameter::pi=3.1415926;
logical::converged;
real::cst(5),cst_p(3),cnm(6,6), xc0_p(3),xc0(3),phai0(3),&
	  c1mtx(3,3),c2mtx(3,3),c3mtx(3,3),c4mtx(3,3),c5mtx(3,3),btcbk(9,9),&
	  qmtx(9,12),nespfn(6,12),tmtx(3,3),massx(6,6),loadx(6,6),centx(6,6),gyrox(6,6),&
	  extfcmx(6),intfcmx(6),inefcmx(6),eextfc(12),eintfc(12),eintfc_g(12),einefc(12),rmtx0(3,3),&
	  rmtx(3,3),rtxc_p(3),tphai_p(3),E1_mat(3),nfc(3),mr(3),nbar(3),mrbar(3),jtensor(3,3),&
	  bmtx(6,9),t2ev(12,12),t2em1(12,12),t2em2(12,12),rmtx0_nd1(3,3),rmtx0_nd2(3,3),&
	  rmtx_nd1(3,3),rmtx_nd2(3,3),tmtx_nd1(3,3),tmtx_nd2(3,3),dertf(3,3),&
	  ocean_disp(1516),ocean_vel(1516),ocean_acc(1516),wtvelo(3),wtvel(3);
!--------------------------dynamic arrays--------------------------------------
!real,allocatable::ekm(:,:),ekm_g(:,:),emm(:,:),egykm(:,:),ectkm(:,:),tkm(:),tmm(:),&
!				  tgykm(:),tctkm(:),d0(:),v0(:),a0(:),acc0(:),points1(:,:),points2(:,:),der(:),der0(:),&
!				  fun(:),kcr(:),eyee(:,:),weights1(:),weights2(:),rloads(:),ttd(:),ttv(:),ttd_0(:),ttv_0(:),&
!				  tta(:),tta_0(:),elxc(:,:),eldd(:,:),eldxc(:,:),eld2xc(:,:),loads(:),&
!				  elphai(:,:),c_elphai(:,:),eldphai(:,:),eld2phai(:,:),xc(:),xc_p(:),dd(:),&
!				  dxc(:),d2xc(:),phai(:),dphai(:),d2phai(:),phai_p(:),no(:),val(:,:),&
!				  tkm_bck(:),tmm_bck(:),tldkm_bck(:),rloads_bck(:),coord(:,:),&
!				  kiter(:),cdmtx(:,:),temp1(:),temp2(:),temp3(:),temp4(:),mtemp(:),ctemp(:),delta(:),g_coord(:,:);
				 				 
real,allocatable::ekm(:,:),ekm_g(:,:),emm(:,:),egykm(:,:),ectkm(:,:),tkm(:,:),tmm(:,:),&
				  tgykm(:,:),tctkm(:,:),d0(:),v0(:),a0(:),acc0(:),points1(:,:),points2(:,:),der(:),der0(:),&
				  fun(:),kcr(:,:),eyee(:,:),weights1(:),weights2(:),rloads(:),tt_intfc(:),ttd(:),ttv(:),ttd_0(:),ttv_0(:),&
				  tta(:),tta_0(:),elxc(:,:),eldd(:,:),eldxc(:,:),eld2xc(:,:),loads(:),ttd_a(:),ttv_a(:),tta_a(:),&
				  elphai(:,:),c_elphai(:,:),eldphai(:,:),eld2phai(:,:),xc(:),xc_p(:),dd(:),&
				  dxc(:),d2xc(:),phai(:),dphai(:),d2phai(:),phai_p(:),no(:),val(:,:),&
				  tkm_bck(:,:),tmm_bck(:,:),tldkm_bck(:,:),rloads_bck(:),coord(:,:),&
				  kiter(:,:),cdmtx(:,:),temp1(:),temp2(:),temp3(:),temp4(:),mtemp(:),ctemp(:),delta(:),g_coord(:,:) ,&
				  work(:,:),copy(:,:),factors(:);
integer,allocatable::nf(:,:),num(:),g_num(:,:),g(:), g_g(:,:);
				 
!-----------------------input and initialisation-------------------------------
open (11 , file = 'cable.res' , status = 'replace', action='write')
open (13 , file = 'G:\����\����ѧģ��\YXYBC3\SIACAB\Dynamic_1\General alpha 1\UMBILICAL\data back\ocean_1\ocean_1_disp.txt' , status = 'old' ,    action ='read')!�������Ĵ�����λ�ƣ�Ƶ��Ϊ10Hz.
open (14 , file = 'G:\����\����ѧģ��\YXYBC3\SIACAB\Dynamic_1\General alpha 1\UMBILICAL\data back\ocean_1\ocean_1_vel.txt' , status = 'old' ,    action ='read')!�������Ĵ������ٶȣ�Ƶ��Ϊ10Hz.
open (15 , file = 'G:\����\����ѧģ��\YXYBC3\SIACAB\Dynamic_1\General alpha 1\UMBILICAL\data back\ocean_1\ocean_1_acc.txt' , status = 'old' ,    action ='read')!�������Ĵ��������ٶȣ�Ƶ��Ϊ10Hz.
!open (13 , file = 'E:\YXYBC3\SIACAB\Dynamic_1\General alpha 1\UMBILICAL\ROPOS1730DATA\ROPOS1730SHIP_Dz.DAT' , status = 'old' ,    action ='read')!ROPOS 1730 ����������λ�ƣ�Ƶ��Ϊ10Hz.
!open (14 , file = 'E:\YXYBC3\SIACAB\Dynamic_1\General alpha 1\UMBILICAL\ROPOS1730DATA\ROPOS1730SHIP_Vz.DAT' , status = 'old' ,    action ='read')!ROPOS 1730 �����������ٶȣ�Ƶ��Ϊ10Hz.
!open (15 , file = 'E:\YXYBC3\SIACAB\Dynamic_1\General alpha 1\UMBILICAL\ROPOS1730DATA\ROPOS1730SHIP_Az.DAT' , status = 'old' ,    action ='read')!ROPOS 1730 �������������ٶȣ�Ƶ��Ϊ10Hz.
open (16 , file = 'D:\���Դ���\Fortran\cableStaticTest1\cableStaticTest1\Data\D17.3 4000�޸��� ��ĩ�˽�� �޺��� ��̬.res', status = 'old' , action ='read');!��ȡ��̬�����Ľ����Ϊ��̬�����ĳ�ʼֵ��
open (10 , file = 'data\ROPOS_4000_17.3.dat' , status = 'old' ,action ='read');!��ROPOSϵͳ����
open (12 , file = 'G:\����\����ѧģ��\YXYBC3\postData\UMBILICAL\DY_DATA\D17.3_4000_�޸���__�޺���_��������_��ROV����.csv' , status = 'replace', action='write')

!read(13, FMT = "(f15.5)", IOSTAT = stat, ADVANCE = 'YES')ocean_disp; !ע�⣺��ȡ�ļ��Ĵ�С����Ӧ�Ĵ洢�����СӦһ�¡�
	if(stat<0)then
		write(*,*)"�ѵ��ļ�ocean_1_disp.txtĩβ,END";!ֱ�Ӵ�ӡ����Ļ
	else if(stat==0)then
		write(*,*)"��ȡ�ļ�ocean_1_disp.txt������NORMAL"
	else
		write(*,*)"��ȡ�ļ�ocean_1_disp.txt��������ERROR"
	end if	
	!write(11,"(f15.5)")ocean_disp;
	!write(11,*)"ocean_disp output endline---------------------------------------------------"
!read(14, FMT = "(f15.5)", IOSTAT = stat, ADVANCE = 'YES')ocean_vel; !ע�⣺��ȡ�ļ��Ĵ�С����Ӧ�Ĵ洢�����СӦһ�¡�
	if(stat<0)then
		write(*,*)"�ѵ��ļ�ocean_1_vel.txtĩβ,END";!ֱ�Ӵ�ӡ����Ļ
	else if(stat==0)then
		write(*,*)"��ȡ�ļ�ocean_1_vel.txt������NORMAL"
	else
		write(*,*)"��ȡ�ļ�ocean_1_vel.txt��������ERROR"
	end if	
	!write(11,"(f15.5)")ocean_vel;
	!write(11,*)"ocean_vel output endline====================================================="
!read(15, FMT = "(f15.5)", IOSTAT = stat, ADVANCE = 'YES')ocean_acc; !ע�⣺��ȡ�ļ��Ĵ�С����Ӧ�Ĵ洢�����СӦһ�¡�
	if(stat<0)then
		write(*,*)"�ѵ��ļ�ocean_1_acc.txtĩβ,END";!ֱ�Ӵ�ӡ����Ļ
	else if(stat==0)then
		write(*,*)"��ȡ�ļ�ocean_1_acc.txt������NORMAL";
	else
		write(*,*)"��ȡ�ļ�ocean_1_acc.txt��������ERROR";
	end if	
	!write(11,"(f15.5)")ocean_acc;
	!write(11,*)"ocean_acc output endline******************************************************"

read(10,*)nels,nn,ell,nip1,nip2,LIMIT,TOL;
!write(*,*)nels,ell,nn,nip,LIMIT,TOL;
!----------Materials parameters-------------
read(10,*)EA,GA2,GA3,GJ,EI2,EI3;
!write(*,*)EA,GA2,GA3,GJ,EI2,EI3,area,rhoaw,rhoac,diacab;
read(10,*)rhoaw,diacab,perwet,permass,masscage,cfcage;
!----------Newmark parameters---------------
read(10,*)cd1,cd2,dt0,spradius;
		cnm =0; 
		cnm(1,1)=EA;cnm(2,2)=GA2;cnm(3,3)=GA3;
		cnm(4,4)=GJ;cnm(5,5)=EI2;cnm(6,6)=EI3;
read(10,*)jtensor1,jtensor2,jtensor3;		
		!���ù�������
		jtensor=0;
		jtensor(1,1)=jtensor1;jtensor(2,2)=jtensor2;jtensor(3,3)=jtensor3;

		!alpha1=1.0/(beta*dt*dt);
		!alpha2=1.0/(beta*dt);
		!alpha3=(1-2*beta)/(2*beta);
		!alpha4=gamma/(beta*dt);
		!alpha5=1.0-gamma/beta;
		!alpha6=(1-gamma/(2*beta))*dt;	
		alpha_f = spradius/(spradius+1.0);
		alpha_m = (2*spradius-1.0)/(spradius+1);
		beta = (1-alpha_m+alpha_f)**2/4.0;
		gamma = 0.5 - alpha_m + alpha_f;
		

allocate(ekm(ndof,ndof),emm(ndof,ndof),egykm(ndof,ndof),ectkm(ndof,ndof),nf(nodof,nn),&
         coord(nod,ndim),g_coord(ndim,nn),num(nod),g_num(nod,nels),g(ndof),g_g(ndof,nels),&
         points1(nip1,ndim),points2(nip2,ndim),der(nod),fun(nod),weights1(nip1),weights2(nip2),eyee(ndim,ndim),&
		 elxc(nod,ndim),eldd(nod,ndim),eldxc(nod,ndim),eld2xc(nod,ndim),&
		 elphai(nod,ndim),c_elphai(nod,ndim),eldphai(nod,ndim),eld2phai(nod,ndim),&
		 xc(ndim),xc_p(ndim),dd(ndim),dxc(ndim),d2xc(ndim),&		 
		 phai(ndim),dphai(ndim),d2phai(ndim),phai_p(ndim));
		 
eyee=0;forall(i=1:3,j=1:3,i==j)eyee(i,j)=1 !���õ�λ����

!read(10,*)g_coord;
!read(10,*)g_num; 
read(10,*)nr;
nf = 1; if(nr>0)read(10,*)(k,nf(:,k),i=1,nr); call formnf(nf); neq=maxval(nf);
!---------------loop the elements to find global array sizes-------------------
nband=0;
elements_1: do iel=1,nels    
              CALL geometry_get(iel,ell,coord,num); CALL num_to_g(num,nf,g)  
              g_num(:,iel)=num; g_coord(:,num)=transpose(coord);
              g_g(:,iel)=g; if(nband<bandwidth(g))nband=bandwidth(g);
			  !num=g_num(:,iel);
			  !CALL num_to_g(num,nf,g);
			  !g_g(:,iel)=g; if(nband<bandwidth(g))nband=bandwidth(g);		  
end do elements_1

!allocate(tkm(neq*(nband+1)),tmm(neq*(nband+1)),tgykm(neq*(nband+1)),tctkm(neq*(nband+1)),&
!		a0(0:neq),acc0(0:neq),d0(0:neq),v0(0:neq),rloads(0:neq),kcr(neq*(nband+1)),                       &
!		ttd(0:neq),ttv(0:neq),tta(0:neq),ttd_0(0:neq),ttv_0(0:neq),tta_0(0:neq),tkm_bck(neq*(nband+1)),tmm_bck(neq*(nband+1)),&
!		tldkm_bck(neq*(nband+1)),rloads_bck(0:neq),kiter(neq*(nband+1)),cdmtx(neq*(nband+1)),&
!		temp1(0:neq),temp2(0:neq),temp3(0:neq),temp4(0:neq),mtemp(0:neq),ctemp(0:neq),delta(0:neq),loads(0:neq),&
!		work(nband+1,neq),copy(nband+1,neq));
allocate(tkm(neq,(2*nband+1)),tmm(neq,(2*nband+1)),tgykm(neq,(2*nband+1)),tctkm(neq,(2*nband+1)),&
		a0(0:neq),acc0(0:neq),d0(0:neq),v0(0:neq),rloads(0:neq),tt_intfc(0:neq),kcr(neq,(2*nband+1)),                       &
		ttd(0:neq),ttv(0:neq),tta(0:neq),ttd_a(0:neq),ttv_a(0:neq),tta_a(0:neq),ttd_0(0:neq),ttv_0(0:neq),tta_0(0:neq),tkm_bck(neq,(2*nband+1)),tmm_bck(neq,(2*nband+1)),&
		tldkm_bck(neq,(2*nband+1)),rloads_bck(0:neq),kiter(neq,(2*nband+1)),cdmtx(neq,(2*nband+1)),&
		temp1(0:neq),temp2(0:neq),temp3(0:neq),temp4(0:neq),mtemp(0:neq),ctemp(0:neq),delta(0:neq),loads(0:neq),&
		work(nband+1,neq),copy(nband+1,neq))
			
	write(11,'(a)')"INITIAL Global coordinates:"
	do k=1,nn 
		write(11,'(a,i5,a,3f12.3)') "Node    ",k,"    ",g_coord(:,k);       			
	end do
	write(11,'(a)')"Global node numbers:"
	do k=1,nels 
		write(11,'(a,i5,a,3i4)') "Element ",k,"      ",g_num(:,k);                    			 
	end do;
	write(11,'(2(a,i5),/)')"There are ",neq,"  equations and the half-bandwidth is ",nband  
!------------------------------GET SEA STATE----------------------------------------------
	!���������˶�
	amp = 0.67; period = 7;! the 4 sea state;
	!amp = 1.0; period = 7;! the 5 sea state;
	omego = 2*pi/period;	
!------------------------------GET LOADS--------------------------------------------------
	!loads=0.0; read(10,*)loaded_nodes
	!if(loaded_nodes/=0)read(10,*)(k,loads(nf(:,k)),i=1,loaded_nodes);
	!write(11,*)loads;
	
	!����ļ��������������غɷ�����������ʱ�䲽������ָ���غ�
	read(10,*)loaded_nodes; allocate(no(loaded_nodes),val(loaded_nodes,nodof)); 
	read(10,*)(no(i),val(i,:),i=1,loaded_nodes);!�ڵ��غ�ֵ
	!�غ���������incs
	!read(10,*)incs; allocate(factors(incs)); 
	!do i=1,incs; factors(i)=i*1.0/incs; end do

	points1=0.0;weights1=0.0;
	points2=0.0;weights2=0.0;
	CALL sample(element,points1,weights1); !��ø�GUASS���ֵ㼰��Ӧ�ļ�Ȩϵ��
	CALL sample(element,points2,weights2); !��ø�GUASS���ֵ㼰��Ӧ�ļ�Ȩϵ��
	!write(11,*)"points1:",points1;
	!write(11,*)"weights1:",weights1;
	!write(11,*)"points2:",points2(1,1),points2(2,1);
	!write(11,*)"weights2:",weights2;
!--------------------------NEWMARK METHOD AND NEWTON-RAPHSON ITERATION---------------------------
	!���ļ��ж�ȡ��ʼֵ
	d0 = 0.0;
	read(16, FMT = "(f15.6)",IOSTAT = stat, ADVANCE = 'YES' )d0(nf(1,1:nn));!��ȡ���ڵ��1���ɶȵ�λ�ơ�
	!read(16, FMT = "(6f15.6)",IOSTAT = stat, ADVANCE = 'YES' )(d0(nf(1:6,I)),I=1,nn);!�к���ʱ������
	if(stat<0)then
		write(*,*)"�ѵ��ļ�initial.resĩβ,END";!ֱ�Ӵ�ӡ����Ļ
	else if(stat==0)then;
		write(*,*)"��ȡ�ļ�initial.res������NORMAL"
	else
		write(*,*)"��ȡ�ļ�initial.res��������ERROR"
	end if	
		!write(11,"(f15.6)")d0(nf(1,1:nn));
		!write(11,"(6f15.6)")(d0(nf(1:6,I)),I=1,nn);
		
	v0=0.0; a0=0.0;!��ʱ�䲽����ֵ��λ�ơ��ٶȡ����ٶȣ�
	!wtvelo = (/0.0,0.5,0.0/);
	wtvel = (/0.0,0.0,0.8/);
	!wtvel = 0.0;
	!ocean_step =1;
	
	write(*,*)"��ʼ����ʱ�䲽ѭ����"
	ts =0.0;
!timeloop:do ts=0,30,dt	
timeloop:do while(ts<=50.0) !����Ӧ���� dt
	if(dt0<=0.0)then;write(*,*)"ʱ�䲽��<0,�˳�!!!";exit timeloop;end if 
	!if(dt0<=0.0001)then;dt0 = 0.01;end if
	dt =dt0;
		alpha1=1.0/(beta*dt*dt);
		alpha2=1.0/(beta*dt);
		alpha3=(1-2*beta)/(2*beta);
		alpha4=gamma/(beta*dt);
		alpha5=1.0-gamma/beta;
		alpha6=(1-gamma/(2*beta))*dt;	
	!ts = ts + dt ;	
	!write(*,*)"The current time is:",ts;
	!ocean_step = ocean_step+1;	
	!n+1ʱ�䲽������ֵԤ�⣨����nʱ�䲽������ֵd0,v0,a0�������������и���У����ֱ������
	!--------------Newmark �㷨��ʱ�䲽������ֵ-------------------
		ttd = d0;
		ttv = alpha5*v0+alpha6*a0;
		tta = -alpha2*v0-alpha3*a0;
	!--------------Newmark �㷨��ʱ�䲽������ֵ-------------------
	!--------------����alpha �㷨��ʱ�䲽������ֵ-----------------
		ttd_a = (1-alpha_f)*ttd+alpha_f*d0;
		ttv_a = (1-alpha_f)*ttv+alpha_f*v0;
		tta_a = (1-alpha_m)*tta+alpha_m*a0;
	!--------------����alpha �㷨��ʱ�䲽������ֵ-----------------		
		iters=0;
		iterations:do
			iters = iters+1;
			!if(iters>3)then;dt0 =dt0-0.001;cycle timeloop;end if !������������ʱ��Ҫ���ٲ���;
			if(iters>5)then;dt0 =dt0/1.01;cycle timeloop;end if !������������ʱ��Ҫ���ٲ���;
				write(11,"(a,i5)")"the main iters:-----",iters;	
			loads=0.0;
			do i = 1,loaded_nodes; loads(nf(:,no(i)))=val(i,:); end do										
			tkm=0.0;tmm=0.0;tgykm=0.0;tctkm=0.0;kcr=0;rloads=0;
			!�����������նȾ�����ÿ�ε����и���
			do iel = 1,nels
				ekm=0.0; emm=0.0; egykm=0;ectkm=0;
				eextfc=0;eintfc=0;einefc=0;
				g=g_g(:,iel);num=g_num(:,iel);
				!if(iel==1)then;write(11,*)g;end if
								
				!���ֽڵ�������λ�ƣ�key point��
				coord=transpose(g_coord(:,num));!��Ԫ�ڵ��ʼ����
				
				!������ȡ�õ�Ԫ�ڵ�ֵ;Լ���ı߽�ֵ
				!�����׶˽ڵ�Լ������
				if(iel==1)then;				
					!�����˶���ע��ĸ�������˶��ǵ����ɶ��˶�����ʼ�˶�����ָ����գ�����X������Ϊ������˳�ʼλ��Ϊ����
					ttd(g(2:3))=0.0;ttd(g(4:6))=0.0;
					ttv(g(2:3))=0.0;ttv(g(4:6))=0.0;
					tta(g(2:3))=0.0;tta(g(4:6))=0.0;
					
					ttd_a(g(2:3))=0.0;ttd_a(g(4:6))=0.0;
					ttv_a(g(2:3))=0.0;ttv_a(g(4:6))=0.0;
					tta_a(g(2:3))=0.0;tta_a(g(4:6))=0.0;
					
					ttd(g(1)) = -amp*sin(omego*(ts+dt));
					ttv(g(1)) = -amp*omego*cos(omego*(ts+dt));
					tta(g(1)) = amp*omego*omego*sin(omego*(ts+dt));	
					
					ttd_a(g(1)) = -((1-alpha_f)*amp*sin(omego*(ts+dt))+alpha_f*amp*sin(omego*(ts)));
					ttv_a(g(1)) = -((1-alpha_f)*amp*omego*cos(omego*(ts+dt))+alpha_f*amp*omego*cos(omego*(ts)));
					tta_a(g(1)) = (1-alpha_m)*amp*omego*omego*sin(omego*(ts+dt))+alpha_m*amp*omego*omego*sin(omego*(ts));	
					
					!��ʵ�ʲ���ĸ���������ݣ�Ƶ��10Hz.
					!ttd(g(1)) = -ocean_disp(ocean_step);
					!ttv(g(1)) = -ocean_vel(ocean_step);
					!tta(g(1)) = -ocean_acc(ocean_step);					
				end if
				
				! extract the displacement of an element from the total
				eldd(1,:)=ttd_a(g(1:3));elphai(1,:)=ttd_a(g(4:6));
				eldd(2,:)=ttd_a(g(7:9));elphai(2,:)=ttd_a(g(10:12));
				
				! extract the velocity of an element from the total
				eldxc(1,:)=ttv_a(g(1:3));eldphai(1,:)=ttv_a(g(4:6));
				eldxc(2,:)=ttv_a(g(7:9));eldphai(2,:)=ttv_a(g(10:12));
				
				! extract the acceleration of an element from the total
				eld2xc(1,:)=tta_a(g(1:3));eld2phai(1,:)=tta_a(g(4:6));
				eld2xc(2,:)=tta_a(g(7:9));eld2phai(2,:)=tta_a(g(10:12));
				
				elxc=coord+eldd;!��Ԫ�ڵ㵱ǰ����											
				jacobi=	0.5*sqrt((coord(2,1)-coord(1,1))**2+(coord(2,2)-coord(1,2))**2+(coord(2,3)-coord(1,3))**2);				
				!write(11,"(a,f6.3)")"the main jacobi:---------",jacobi;
				
				do k=1,nip1
					!���ص�k��GUASS���ֵ��Ӧ�ĵ�Ԫ�ڵ��κ���������
					call spfn_get(der,fun,qmtx,nespfn,points1,k,jacobi);						
					
					!���в�ֵ�õ���Ԫ�ڲ�ĳ�����
					!xc=fun(1)*elxc(1,:)+fun(2)*elxc(2,:);!�ռ������ʾ
					xc_p=der(1)*elxc(1,:)+der(2)*elxc(2,:);!�ռ������ʾ
					dxc=fun(1)*eldxc(1,:)+fun(2)*eldxc(2,:);
					d2xc=fun(1)*eld2xc(1,:)+fun(2)*eld2xc(2,:);
					
					phai=fun(1)*elphai(1,:)+fun(2)*elphai(2,:);
					phai_p=der(1)*elphai(1,:)+der(2)*elphai(2,:);
					!dphai=fun(1)*eldphai(1,:)+fun(2)*eldphai(2,:);		
					!d2phai=fun(1)*eld2phai(1,:)+fun(2)*eld2phai(2,:);
					
					!!!!!�ڵ�Ԫ��Χ�ڶ�����xc,xc_p,phai,phai_p��ص�ʸ�������
					!phaiʸ���Ĵ�С									
					fi = sqrt(dot_product(phai,phai));					
					CALL coef(cst,cst_p,fi);!c1,c2,c3,c4,c5,c6......	
					
					CALL rtmtx_get(rmtx,tmtx,phai,fi);!��Ӳ��ϻ����˶���ת������rmtx,tmtx;																																			
					rtxc_p=matmul(transpose(rmtx),xc_p);																		
					tphai_p = matmul(tmtx,phai_p);
					
					!��ⵥλ�ο����������Ӧ�����Ӧ��ż
					CALL nfcmr_get(nfc,mr,cnm,rtxc_p,tphai_p);!��E1=(1,0,0)ʱ
					!��ⵥλ�ο���������������������ż									
					CALL nmrbar_get(nbar,mrbar,perwet,rhoaw,diacab,elxc,dxc,d2xc,rmtx,points1,k,wtvel);!ˮ������
					!nbar=(/29.7,0.0,0.0/);mrbar=0.0;

					!��ⵥԪ�նȾ����е�BTCBK����
					CALL btcbk_get(btcbk,bmtx,rmtx,tmtx,rtxc_p,nfc,mr,cst,cst_p,phai_p,phai,cnm,fi);					
					!if(iters==1 .and. iel==nels)write(11,"(a,/,81e12.4)")"btcbk:",btcbk
					!��΢�غɾ��� loadx
					!CALL kload(loadx,cst,phai,mrbar);
					
					!��Ԫ�ⲿ�ڵ���(12,1)			
					CALL extfcmx_get(extfcmx,nbar,mrbar,tmtx);		
						!write(11,"(a,/,6e12.4)")"extfcmx",extfcmx;
					eextfc = eextfc + matmul(transpose(nespfn),extfcmx)*jacobi*weights1(k);
						!write(11,"(a,/,12e12.4)")"eextfc",eextfc;
						
					!��Ԫ�ڲ��ڵ���(12,1)
					CALL intfcmx_get(intfcmx,nfc,mr);
					eintfc = eintfc + matmul(transpose(matmul(bmtx,qmtx)),intfcmx)*jacobi*weights1(k);															
									
					!��Ԫ�նȾ���(12,12)
					ekm = ekm+ matmul(matmul(transpose(qmtx),btcbk),qmtx)*jacobi*weights1(k);					
					!write(11,"(a,i5)")"this is ekm of iteration:",iters
					!write(11,*)ekm;
					
					!��Ԫ�غɸնȾ���(12,12),��������������ͨ��Ϊ0;
					!eldkm = eldkm+ matmul(matmul(transpose(nespfn),loadx),nespfn)*jacobi*weights1(k);
					!eldkm=0.0;
					!write(11,"(a,i5)")"this is eldkm of iteration:",iters
					!write(11,*)eldkm;
					
				end do !nip1
				
				do k=1,nip2
					!���ص�k��GUASS���ֵ��Ӧ�ĵ�Ԫ�ڵ��κ���������
					call spfn_get(der,fun,qmtx,nespfn,points2,k,jacobi);	
					!���в�ֵ�õ���Ԫ�ڲ�ĳ�����
					!xc=fun(1)*elxc(1,:)+fun(2)*elxc(2,:);!�ռ������ʾ
					!xc_p=der(1)*elxc(1,:)+der(2)*elxc(2,:);!�ռ������ʾ
					!dxc=fun(1)*eldxc(1,:)+fun(2)*eldxc(2,:);
					!d2xc=fun(1)*eld2xc(1,:)+fun(2)*eld2xc(2,:);
					
					phai=fun(1)*elphai(1,:)+fun(2)*elphai(2,:);
					!phai_p=der(1)*elphai(1,:)+der(2)*elphai(2,:);
					dphai=fun(1)*eldphai(1,:)+fun(2)*eldphai(2,:);
					d2phai=fun(1)*eld2phai(1,:)+fun(2)*eld2phai(2,:);
					
					!!!!!�ڵ�Ԫ��Χ�ڶ�����xc,xc_p,phai,phai_p��ص�ʸ�������
					!phaiʸ���Ĵ�С									
					fi = sqrt(dot_product(phai,phai));
						!write(11,*)"inp2 fi:",fi;
					CALL coef(cst,cst_p,fi);!c1,c2,c3,c4,c5,c6......					
					CALL rtmtx_get(rmtx,tmtx,phai,fi);!��Ӳ��ϻ����˶���ת������rmtx,tmtx;
					
					!��΢����������
					!CALL der_tf(dertf,phai,dphai,cst,fi);
					!CALL kcent_gyro(centx,gyrox,phai,dphai,d2phai,tmtx,dertf,jtensor,cst,cst_p,fi);
					
					!��Ԫ���Խڵ���(12,1)-B����		
					CALL inefcmx_get(inefcmx,tmtx,phai,dphai,cst,jtensor,fi,dertf);	
						!write(11,"(a,/,6e12.4)")"inefcmx",inefcmx;
					einefc = einefc + matmul(transpose(nespfn),inefcmx)*jacobi*weights2(k);
					
					!����phai��ص�΢�������� massx
					!CALL mass(massx,area,rhoac,tmtx,jtensor);
					CALL mass(massx,tmtx,jtensor,permass);					
					!��Ԫ��������(12,12)
					emm = emm+ matmul(matmul(transpose(nespfn),massx),nespfn)*jacobi*weights2(k);						
														
					!��Ԫ�������նȾ���(�ǶԳƣ�����ת��ʱ�ɺ��ԣ�����ʱ����)
					!ectkm = ectkm+ matmul(matmul(transpose(nespfn),centx),nespfn)*jacobi*weights2(k);
					
					!��Ԫ�������նȾ���(�ǶԳƣ�����ת��ʱ�ɺ��ԣ�����ʱ����)
					!egykm = egykm+ matmul(matmul(transpose(nespfn),gyrox),nespfn)*jacobi*weights2(k);					
				end do !nip2				
					g=g_g(:,iel);
					!CALL formkv(tkm,ekm,g,neq) ;!��װ����նȾ���tkm����neq��(nband+1)�еĶ�ά�������洢���Խ���Ԫ�ش洢�ڵ�(nband+1)��
					!CALL formkv(tmm,emm,g,neq) ;!��װ������������tmm
					!CALL formkv(tctkm,ectkm,g,neq);!��װ�����������նȾ���ectkm
					!CALL formkv(tgykm,egykm,g,neq);!��װ�����ת���նȾ���egykm
					if(iel==nels)then!����ˮ��[ĩ�˽����ROV������]����[��������=��������ϵ��*�û�ˮ������(����/g)]����λKG
						!emm(7,7)=emm(7,7)+3.0/3.0*(2790+3000+1.5*4320);
						emm(7,7)=emm(7,7)+ masscage;!masscage��ʾĩ�˽�������������ˮ�������ܺͣ���F. R. Driscoll�в�����ROV����
						do i=1,12;write(11,"(12f12.3)")emm;end do
					end if 
								
					!��װ�ǶԳƾ���
					CALL FORMTB(tkm,ekm,g);
					CALL FORMTB(tmm,emm,g);
					!CALL FORMTB(tctkm,ectkm,g);
					!CALL FORMTB(tgykm,egykm,g);
					
					!ĩ�˽��(�ʵ�)���ظ�����ˮ����
					if(iel==nels)then;
						loads(nf(1,iel+1)) =loads(nf(1,iel+1))-cfcage*eldxc(2,1)*abs(eldxc(2,1));
						!write(11,*)loads;
					end if
					
					!��װ����ڵ��غ�
					rloads(g)=rloads(g)+eextfc - eintfc - einefc;
					!rloads(g)=rloads(g) - eintfc- einefc;
			end do !nels
			!------------------------------equation solution---------------------------------------------
			!------------------------------Newmark Algorithm---------------------------------------------
			!��װ��ɺ��������󱸷ݡ�����ڵ��غɱ���
			tkm_bck=tkm;
			tmm_bck=tmm;
			!tldkm_bck=tldkm;
			rloads_bck=rloads;				
			kcr = 0.0-tkm;!����ϵ������ Kcr
			
			!write(11,"(a,/,31e12.4)")"total rloads:--------------------",rloads
			
			!write(11,"(a,i5)")"this is tkm of iteration:",iters
			!write(11,*)tkm;
			
			!write(11,"(a,i5)")"this is tmm of iteration:",iters
			!write(11,*)tmm;
			
			!write(11,"(a,i5)")"this is tldkm of iteration:",iters
			!write(11,*)tldkm;
					
			!���ʼ���ٶ�
			!call cholin(tmm);!�������������CHOLESKI�ֽ�
			!a=0;!���ٶȳ�ֵΪ��
			!a=rloads;!�������ɶȳ�ʼ�غɸ������ٶȱ���
			!call chobac(tmm,a)!���ٶȽ��������a��		
			
			!NEWTON�������������Է�����
			cdmtx= cd1*tmm+cd2*tkm;!�������
			kiter=alpha1*(1-alpha_m)*tmm +alpha4*(1-alpha_f)*(cdmtx) -(1-alpha_f)*kcr;	!����ϵ������,Ted Belytschko���ϳ�ΪJACOBIAN����
			!����λ�Ʊ߽�����������JACOBIAN����,�Ⱥ���ˣ�
			kiter(1,:)=0.0;kiter(1,nband+1)=10000;
			kiter(:,nband+1)=kiter(:,nband+1)+1.0;!Ϊ��ֹ���Է�����Ax=b�е�A���ֲ�̬ʱ������ֵ���졣
			!write(11,"(a)")"the system matrix :";
			!write(11,*)kiter;
			!�󷽳̵Ⱥš�=���ұ߸���
			temp1=(1-alpha_m)*(alpha1*(ttd-d0)-alpha2*v0-alpha3*a0)+alpha_m*a0; 
			temp2=(1-alpha_f)*(alpha4*(ttd-d0)+alpha5*v0+alpha6*a0)+alpha_f*v0;
			CALL BANTMUL(tmm,temp1,mtemp);
			CALL BANTMUL(cdmtx,temp2,ctemp);
			delta=loads+rloads-mtemp-ctemp;
			delta(0)=0.0;
			delta(1) = 0.0;!����JACOBIAN���󣬵Ⱥ��Ҷ�
		
			!CALL banred(kiter,neq); CALL bacsub(kiter,delta);
			CALL GAUSS_BAND(kiter,WORK);
			COPY = WORK;
			CALL SOLVE_BAND(kiter,COPY,delta);			
			!write(11,"(a)")"the incremental displacement delta:";
			!write(11,*)delta;
			
			!-----------UPDATE PROCESS----------------	
			!delta = delta*2.0;
			ttd = ttd+delta; !��λ�Ƹ���
			ttv = ttv+alpha4*delta;
			tta = tta+alpha1*delta;
			
			ttd_a = ttd_a+(1-alpha_f)*delta;
			ttv_a = ttv_a+alpha4*(1-alpha_f)*delta;
			tta_a = tta_a+alpha1*(1-alpha_m)*delta

			!ttv = alpha4*(ttd-d0)+alpha5*v0+alpha6*a0;
			!tta = alpha1*(ttd-d0)-alpha2*v0-alpha3*a0;
			!ttv = ttv+gamma/(beta*dt)*delta;!���ٶȸ���
			!tta = tta+1.0/(beta*dt*dt)*delta;!�ܼ��ٶȸ���	
			
			!ttd = ttd +delta;
			!tta = (ttd-ttd_0)/(beta*dt*dt);
			!ttv = ttv_0+gamma*dt*tta;

			ttd(0)=0.0;
					
			!-----------check converged---------------
			CALL checonverg(ttd,delta,tol,converged)		
			if(iters==limit .or. converged)exit iterations
			!if(iters==converged)exit iterations
			!write(11,*)"END DO once ITERATIONS!"	
		end do iterations
			!-----------����������ĸ���--------------
			d0 = ttd;
			v0 = ttv;
			a0 = tta;
			
			!����������Ԫ�ڲ��ڵ���,Ҳ����������
			tt_intfc = 0.0;
			do iel = 1,nels
				eintfc=0;
				g=g_g(:,iel);num=g_num(:,iel);
				coord=transpose(g_coord(:,num));
				
				eldd(1,:)=d0(g(1:3));elphai(1,:)=d0(g(4:6));
				eldd(2,:)=d0(g(7:9));elphai(2,:)=d0(g(10:12));
				
				elxc=coord+eldd;!��Ԫ�ڵ㵱ǰ����											
				jacobi=	0.5*sqrt((coord(2,1)-coord(1,1))**2+(coord(2,2)-coord(1,2))**2+(coord(2,3)-coord(1,3))**2);			
				
				do k = 1,nip1
					call spfn_get(der,fun,qmtx,nespfn,points1,k,jacobi);
					xc_p=der(1)*elxc(1,:)+der(2)*elxc(2,:);
					phai=fun(1)*elphai(1,:)+fun(2)*elphai(2,:);
					phai_p=der(1)*elphai(1,:)+der(2)*elphai(2,:);
					
					fi = sqrt(dot_product(phai,phai));					
					CALL coef(cst,cst_p,fi);!c1,c2,c3,c4,c5,c6......	
					
					CALL rtmtx_get(rmtx,tmtx,phai,fi);!��Ӳ��ϻ����˶���ת������rmtx,tmtx;																																			
					rtxc_p=matmul(transpose(rmtx),xc_p);																		
					tphai_p = matmul(tmtx,phai_p);
					
					!��ⵥλ���������Ӧ�����Ӧ��ż
					CALL nfcmr_get(nfc,mr,cnm,rtxc_p,tphai_p);!��E1=(1,0,0)ʱ
					!��ⵥԪ�նȾ����е�BTCBK����
					CALL btcbk_get(btcbk,bmtx,rmtx,tmtx,rtxc_p,nfc,mr,cst,cst_p,phai_p,phai,cnm,fi);
					!��Ԫ�ڲ��ڵ���(12,1)
					CALL intfcmx_get(intfcmx,nfc,mr);
					eintfc = eintfc + matmul(transpose(matmul(bmtx,qmtx)),intfcmx)*jacobi*weights1(k);							
				end do
				
				g=g_g(:,iel);
				tt_intfc(g) = tt_intfc(g)+0.0-eintfc;
				
			end do !nels when converged
			
			!�����к����
			write(11,'(a,F10.5)')"--------------------Computed time histories at time:---------------------------  ",ts;
			write(11,"(a,i5,a,/)")"Oh,my God,Converged need ",iters," iterations!"				
			write(11,'(a)')"NODE:                     DISPLACEMENT:        "
			do i=1,nn
				write(11,"(f10.3,i5,6e12.4),//") ts,i,d0(nf(1:3,i));
				write(12,"(f10.3,',',i5,',',3(e15.6,','),3(e15.6,','),3(e15.6,','),3(e15.6,','),e15.6)") ts,i,d0(nf(1:3,i))+g_coord(:,i),d0(nf(1:3,i)),v0(nf(1:3,i)),a0(nf(1:3,i)),tt_intfc(nf(1,i));!�����ļ�Ϊ.csv��ʽ
			end do
		if(iters == limit)then
			write(11,"(a)")"The iterations has overpassed the LIMIT number!";
			exit timeloop;
		end if 
	!ocean_step = ocean_step+1;	
	ts = ts + dt ;	
		write(*,*)"The current time is:",ts;
	if(iters<2)then;dt0 = dt0*2.0;write(11,*)"The dt is added at last!";end if !������������ʱ��Ҫ���Ӳ���;
		!write(11,"(a,i5)")"the main iters:-----",iters;		
end do timeloop
end program cable