class DQ:
 # Initialize Double Queue with joining & serving probabilities, methods find decay rates of large queue probabilities
 def __init__(self,np=0.1,nn=0.1,no=0.2,up=0.3,un=0.3):
   if (np+nn+no+up+un!=1)|(np<=0)|(nn<=0)|(no<=0)|(up<=0)|(un<=0)|(np>=up)|(nn>=un)|(up+un<=0.5): print("bad input")
   else: self.__Np,self.__Nn,self.__No,self.__Up,self.__Un=np,nn,no,up,un
 #decay rates of either & both queues by Kernal, TailsD private & tailsD public for true Double Queues
 def __TailsK(np,nn,no,up,un):
   from pylab import poly1d
   po,pp,pn=(np+nn+no)/(up+un),np/up,nn/un
   axp,axn,bx,cxp,cxn=poly1d([nn,up]),poly1d([np,un]),poly1d([1.,0.]),poly1d([np+no,un,0.]),poly1d([nn+no,up,0.])
   ayp,ayn,byp,byn,cyp,cyn=np+no,nn+no,poly1d([-nn,1.,-un]),poly1d([-np,1.,-up]),poly1d([up,0.,0.]),poly1d([un,0.,0.])
   dxp,dxn,dyp,dyn=bx**2-4*axp*cxp,bx**2-4*axn*cxn,byp**2-4*ayp*cyp,byn**2-4*ayn*cyn
   kp,kn,ko=min(1/pn,dyp.r[1]),min(1/pp,dyn.r[1]),min(dxp.r[0],po**-2,dxn.r[0])
   kyp,kyn=(byp(1/pn)+dyp(1/pn)**0.5)/2/ayp,(byn(1/pp)+dyn(1/pp)**0.5)/2/ayn
   if kyp>1: ko=min(ko,kyp)
   if kyn>1: ko=min(ko,kyn)
   kxp,kxn=(po**-2+dxp(po**-2)**0.5)/2/axp(po**-2),(po**-2+dxn(po**-2)**0.5)/2/axn(po**-2)
   if kxp>1: kp=min(kp,kxp)
   if kxn>1: kn=min(kn,kxn)
   return [1/kn,1/ko,1/kp,"Kernel"]
 def tailsK(self): return DQ.__TailsK(self.__Np,self.__Nn,self.__No,self.__Up,self.__Un)
 #decay rates of either & both queues by Doublesided Quasi Birth Death, TailsD private & tailsD public for true Double Queues
 def __TailsD(np,nn,no,up,un):
   from math import log,e
   po,pp,pn=(np+nn+no)/(up+un),np/up,nn/un
   yp,yn,qp,qn=up*po**2+nn,un*po**2+np,pp*(nn+no)/un,pn*(np+no)/up
   mn1=log(1-2*(un*(nn+no))**0.5+(1+4*(un*(nn+no))-(un*(nn+no))**0.5-up*np)**0.5)-log(2)-log(np)
   mp1=log(1-2*(up*(np+no))**0.5+(1+4*(up*(np+no))-(up*(np+no))**0.5-un*nn)**0.5)-log(2)-log(nn)
   mn2=log(1-4*(np*up+un*(nn+no))+((1-4*(np*up+un*(nn+no)))**2-64*(nn+no)*np*un*up)**0.5)-log(8*np)-log(nn+no)
   mp2=log(1-4*(nn*un+up*(np+no))+((1-4*(nn*un+up*(np+no)))**2-64*(np+no)*nn*up*un)**0.5)-log(8*nn)-log(np+no)
   sp1,sp2,sn1,sn2,co,cp,cn=mp1**1,mp2**1,mn1**1,mn2**1,-2*log(po),-log(po),-log(po)
   if un>nn+no: sn1,sn2=-log(pp),-log(pp)
   if up>np+no: sp1,sp2=-log(pn),-log(pn)
   rn,ro,rp=e**-sn1,po**2,e**-sp1
   if (yn+no>yp)&(yp+no>yn):
     if (pp<po)&(pn<po):
       if sn2>-2*log(po): rn=po*yn/(yp+no)
       if sp2>-2*log(po): rp=po*yp/(yn+no)
     elif pn<po:
       if sp2>-log(qp): rp=(qp*(nn+up*qp)/(un*qp+np+no))**0.5
       ro,rn=qp,pp
     else:
       if sn2>-log(qn): rn=(qn*(np+un*qn)/(up*qn+nn+no))**0.5
       ro,rp=qn,pn
   elif yp+no>yn:
     co,cn=mn2**1,mn2-log(2)-log(np*e**mn2+un)
     cp=log(1-(np+no/2+un/e**co)/e**cn)-log(up/e**co+nn+no/2)
     if (pp<po)&(pn<po):
       if sp1<cp: ro,rp=qn,pn
       else: ro,rp=e**-mn2,min(rp,((nn*e**co+up)/(un*e**co+np+no))**0.5)
     elif pp<po:
       if sn2>=-log(qn): rn=(qn*(np+un*qn)/(up*qn+nn+no))**0.5
       ro,rp=qn,pn
     else:
       if (cn<-log(pp))&(cp<sp1): ro,rn=e**-mn2,pp
       elif (cp<sp1)|(cn>=-log(pp))&(qp<qn): ro,rn,rp=qp,pp,min(rp,(qp*(nn+up*qp)/(un*qp+np+no))**0.5)
       else: ro,rp,rn=qn,pn,min(rn,(qn*(np+un*qn)/(up*qn+nn+no))**0.5)
   else:
     co,cp=mp2**1,mp2-log(2)-log(nn*e**mp2+up)
     cn=log(1-(nn+no/2+up/e**co)/e**cp)-log(un/e**co+np+no/2)
     if (pn<po)&(pp<po):
       if sn1<cn: ro,rn=qp,pp
       else: ro,rn=e**-mp2,min(rn,((np*e**co+un)/(up*e**co+nn+no))**0.5)
     elif pn<po:
       if sp2>=-log(qp): rp=(qp*(nn+up*qp)/(un*qp+np+no))**0.5
       ro,rn=qp,pp
     else:
       if (cp<-log(pn))&(cn<sn1): ro,rp=e**-mp2,pn
       elif (cn<sn1)|(cp>=-log(pn))&(qn<qp): ro,rp,rn=qn,pn,min(rn,(qn*(np+un*qn)/(up*qn+nn+no))**0.5)
       else: ro,rn,rp=qp,pp,min(rp,(qp*(nn+up*qp)/(un*qp+np+no))**0.5)
   return [rn,ro,rp,"DQBD"]
 def tailsD(self): return DQ.__TailsD(self.__Np,self.__Nn,self.__No,self.__Up,self.__Un)
 #decay rate & error of either queue by Iterative Power, Tailpn private & tailp & tailn public for true Double Queues
 def __Tailpn(m=200,np=0.1,nn=0.1,no=0.2,up=0.3,un=0.3,k=2,mo=[[0.,0.],[1.,0.],[0.,0.]],ms=[[0.3,0.3],[0.3,0.],[0.1,0.]]):
   if k>m: return [ms[k][0]/ms[k-1][0],ms[k][0]/ms[k-1][0]-mo[k][0]/mo[k-1][0],(ms[k][0]/ms[k-1][0])**2]
   mo,ms=[mo[0]+[0.]]+mo+[mo[0]+[0.]],[mo[0]+[0.]]+ms+[mo[0]+[0.]]
   for i in range(1,2*k):
     mo[i],ms[i]=mo[i]+[0.],ms[i]+[0.]
     ms[i][0]=nn*mo[i-1][0]+un*mo[i+1][0]+up*(mo[i][0]+mo[i-1][1])
     for j in range(1,k): ms[i][j]=up*mo[i-1][j+1]+nn*mo[i-1][j]+un*mo[i+1][j]+(np+no)*mo[i+1][j-1]
   return DQ.__Tailpn(m,np,nn,no,up,un,k+1,ms,mo)
 def tailp(self,m=200): return DQ.__Tailpn(m,self.__Np,self.__Nn,self.__No,self.__Up,self.__Un)+["+Power"]
 def tailn(self,m=200): return DQ.__Tailpn(m,self.__Nn,self.__Np,self.__No,self.__Un,self.__Up)+["-Power"]
 #decay rate & error of both queues by Iterative Power, Tailo private & tailo public for true Double Queues
 def __Tailo(m,np=0.1,nn=0.1,no=0.2,up=0.3,un=0.3,k=1,o=1.,mo=[[0.,0.,0.],[0.,0.,0.]],ms=[[1.,0.,0.],[0.,0.,0.]]):
   mo,ms=[mo[0]+[0.,0.,0.,0.]]+mo+[mo[0]]+[mo[0]+[0.,0.,0.,0.]],[mo[0]+[0.,0.,0.,0.]]+ms+[mo[0]]+[mo[0]+[0.,0.,0.,0.]]
   if k>m: return [ms[k-1][0]/ms[k][0],ms[k-1][0]/ms[k][0]-o,(ms[k-1][0]/ms[k][0])**2,"oPower"]
   mo[1],ms[1],o=mo[1]+[0.,0.,0.,0.],ms[1]+[0.,0.,0.,0.],ms[k-1][0]/ms[k][0]
   for i in range(1,3*k):
     mo[i+1],ms[i+1]=mo[i+1]+[0.,0.,0.,0.],ms[i+1]+[0.,0.,0.,0.]
     mo[i][1]=up*ms[i-1][0]+(nn+no/2)*ms[i][0]+un*ms[i][3]+(np+no)*ms[i+1][3]
     mo[i][2]=un*ms[i-1][0]+(np+no/2)*ms[i][0]+up*ms[i][4]+(nn+no)*ms[i+1][4]
     for j in range(1,k):
       mo[i][4*j+1]=up*ms[i-1][4*j-1]+nn*ms[i][4*j-1]+un*ms[i][4*j+3]+(np+no)*ms[i+1][4*j+3]
       mo[i][4*j+2]=un*ms[i-1][4*j]+np*ms[i][4*j]+up*ms[i][4*(j+1)]+(nn+no)*ms[i+1][4*(j+1)]
   ms[0][0]=(np+no)*mo[1][1]+(nn+no)*mo[1][2]
   for i in range(1,3*k+1):
     ms[i][0]=un*mo[i][1]+up*mo[i][2]+(np+no)*mo[i+1][1]+(nn+no)*mo[i+1][2]
     for j in range(1,k+1):
       ms[i][4*j-1]=up*mo[i-1][4*j-3]+nn*mo[i][4*j-3]+un*mo[i][4*j+1]+(np+no)*mo[i+1][4*j+1]
       ms[i][4*j]=un*mo[i-1][4*j-2]+np*mo[i][4*j-2]+up*mo[i][4*j+2]+(nn+no)*mo[i+1][4*j+2]
   return DQ.__Tailo(m,np,nn,no,up,un,k+1,o,mo,ms)
 def tailo(self,m=100): return DQ.__Tailo(m,self.__Nn,self.__Np,self.__No,self.__Un,self.__Up)
 #tailsP combines all decay rates and errors from Iterative Power, tails combines all decay rates from all methods
 def tailsP(self,m=100): return [self.tailn(2*m),self.tailo(m),self.tailp(2*m)]
 def tails(self,m=100): return [self.tailsK()]+self.tailsP(m)+[self.tailsD()]